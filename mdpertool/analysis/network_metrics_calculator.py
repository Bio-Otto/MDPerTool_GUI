
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import re

def generate_pc_bar_plot(df_metrics: pd.DataFrame, output_image: str):
    """
    Hesaplanan Propagation Coefficient (PC) skorlarını görsel (Bar Plot) haline getirir
    ve .png olarak kaydeder.
    """
    try:
        import os
        norm_output_image = os.path.normpath(output_image)
        #print(f"[DEBUG][PC Plot] PNG path: {norm_output_image}")
        #print(f"[DEBUG][PC Plot] DataFrame head:")
        #print(df_metrics.head())
        # Sadece pozitif PC'si olan en yüksek 30 rezidüyü gösterelim ki grafik okunabilir olsun
        df_plot = df_metrics[df_metrics['Propagation_Coefficient (PC)'] > 0].copy()
        df_plot = df_plot.sort_values('Propagation_Coefficient (PC)', ascending=False).head(30)

        #print(f"[DEBUG][PC Plot] Plot DataFrame head:")
        #print(df_plot.head())

        # Bar plot verilerini print et
        #print(f"[DEBUG][PC Plot] Bar X: {df_plot['Residue_ID'].tolist()}")
        #print(f"[DEBUG][PC Plot] Bar Y: {df_plot['Propagation_Coefficient (PC)'].tolist()}")

        plt.figure(figsize=(10, 6))
        if df_plot.empty:
            plt.text(0.5, 0.5, 'No data to plot', fontsize=16, ha='center', va='center')
            plt.axis('off')
        else:
            plt.bar(df_plot['Residue_ID'], df_plot['Propagation_Coefficient (PC)'], color='salmon', edgecolor='black')
            plt.xticks(rotation=45, ha='right', fontsize=9)
            plt.ylabel('Propagation Coefficient (PC)', fontsize=12, fontweight='bold')
            plt.xlabel('Residue', fontsize=12, fontweight='bold')
            plt.title('Top Allosteric Super-Hubs (Signal Propagation Capacity)', fontsize=14, fontweight='bold')
            plt.tight_layout()
        try:
            plt.savefig(norm_output_image, dpi=300)
            #print(f"[DEBUG][PC Plot] PNG saved: {norm_output_image}")
            # Kısa ve sade bir test yoluna da kaydet
            test_path = os.path.normpath('C:/PC_plot_test.png')
            plt.savefig(test_path, dpi=300)
            #print(f"[DEBUG][PC Plot] Test PNG saved: {test_path}")
        except Exception as e:
            #print(f"[PC Plot Hatası] PNG kaydedilemedi: {e}")
            pass
        plt.close()
    except Exception as e:
        pass
        #print(f"[PC Plot Hatası] Bar plot çizilemedi: {e}")

def generate_pymol_macro(df_metrics: pd.DataFrame, output_pml: str, source_node: str):
    """
    Makalelerdeki ısı haritası ve allosterik yolak renklendirme uygulamasını yapmak için
    PyMOL makro betiği (.pml) üretir. PyMOL'de 'run script' denilince PC skorları renklendirilir.
    """
    try:
        lines = []
        lines.append("# MDPerTool - Propagation Coefficient (PC) Coloring Macro")
        lines.append("hide everything, all")
        lines.append("show cartoon, all")
        lines.append("color white, all")
        lines.append("alter all, b=0.0") # Tüm B-factor'leri sıfırla

        # Her rezidü için numarasını bul (Örn: 'GLU45' -> 45)
        for _, row in df_metrics.iterrows():
            res_id = row['Residue_ID']
            pc_score = row['Propagation_Coefficient (PC)']
            
            match = re.search(r'\d+', str(res_id))
            if match and pc_score > 0:
                res_num = match.group()
                # PyMOL b-factor ataması
                lines.append(f"alter resi {res_num}, b={pc_score}")

        lines.append(f"spectrum b, white_red, all, minimum=0.001")
        
        # Source (Kaynak) düğümü ekstra vurgula
        source_match = re.search(r'\d+', str(source_node))
        if source_match:
            source_num = source_match.group()
            lines.append(f"show spheres, resi {source_num}")
            lines.append(f"color blue, resi {source_num}")
            lines.append(f"set sphere_scale, 0.8, resi {source_num}")

        with open(output_pml, "w") as f:
            f.write("\n".join(lines))
    except Exception as e:
        pass
        #print(f"[PC PyMOL Error] PyMOL makrosu oluşturulamadı: {e}")


def calculate_propagation_coefficient(G: nx.DiGraph, source_node: str, output_csv: str = "propagation_coefficients.csv"):
    """
    Kullanıcıya Allosterik Sinyal Ağı üzerindeki düğümlerin fonksiyonunu
    (Super-Hub, Kaynak, Hedef vb.) anlatan Propagation Coefficient (PC) analizi modülü.
    
    PC(i) = (m_i * n_i) / (sum(m_j * n_j) + l)
    
    Deep-preference hierarchical layout (BFS üzerinden Layer tespiti) kullanılarak
    her bir yönlü ağ düğümü (rezidü) için makaleye sadık hesaplama ve yapısal yorumlama yapar.
    
    :param G: networkx.DiGraph (Önceden oluşturulmuş sinyal ağı)
    :param source_node: str (Ağdaki sinyali başlatan ana kaynak düğüm/rezidü)
    :param output_csv: str (Sonuçların kaydedileceği dosya yolu)
    :return: Kaydedilen verileri içeren pandas.DataFrame veya hata durumunda None
    """
    
    if source_node not in G.nodes():
        #print(f"[PC HATA] {source_node} ağı içinde bulunamadı.")
        return None

    # 1. Hiyerarşik Katmanları Belirle (Shortest Path / BFS uzaklığı)
    try:
        layers = nx.single_source_shortest_path_length(G, source_node)
    except Exception as e:
        #print("[PC HATA] Kaynak düğümden ağa ulaşılamadı veya ağ bağlantısız!", e)
        return None

    # Düğümleri katmanlarına göre grupla
    layer_dict = {}
    for node, dist in layers.items():
        if dist not in layer_dict:
            layer_dict[dist] = []
        layer_dict[dist].append(node)
        G.nodes[node]['layer'] = dist

    metrics = []

    # 2. Her bir katman (N) için PC hesapla
    for N, nodes_in_layer_N in layer_dict.items():
        # L'yi hesapla (Katman N'yi transit geçen kenarlar / edges that bypass layer N)
        # Başlangıç düğümünün katmanı < N ve bitiş düğümünün katmanı > N olan kenarlar.
        l_edges = 0
        for u, v in G.edges():
            u_layer = G.nodes[u].get('layer', -1)
            v_layer = G.nodes[v].get('layer', -1)
            
            if u_layer != -1 and v_layer != -1:
                if u_layer < N and v_layer > N:
                    l_edges += 1
        
        # Paydadaki Düğüm Üzerinden Geçen Yolların Toplamı: sum(m_j * n_j)
        layer_pathways_sum = sum([G.in_degree(j) * G.out_degree(j) for j in nodes_in_layer_N])
        
        # Payda (Denominator): 
        # (Katmanda gerçekleşen tüm düğüm yolları + Katmanı pas geçen transit kollar)
        denominator = layer_pathways_sum + l_edges
        
        # Her düğüm için PC(i) Puanlarını Hesapla ve Yorumla
        for i in nodes_in_layer_N:
            m_i = G.in_degree(i)
            n_i = G.out_degree(i)
            numerator = m_i * n_i
            
            # Payda sıfır ise PC 0 alınır
            pc_i = (numerator / denominator) if denominator > 0 else 0.0
            
            # --- Kullanıcı için Yorumlama Mantığı (Interpretation logic) ---
            if m_i == 0 and n_i > 0:
                role = "Initiator / Source"
                explanation = "Starts the energy/signal perturbation."
            elif n_i == 0 and m_i > 0:
                role = "Sink / Target"
                explanation = "Final destination of the signal."
            elif pc_i > 0:
                if n_i > m_i:
                    role = "Amplifier / Super-Hub"
                    explanation = "Distributes incoming signal to a wider area. Crucial allosteric distribution hub."
                elif m_i > n_i:
                    role = "Funnel / Collector"
                    explanation = "Merges signals from different branches into a single channel. Prominent in signal isolation."
                else: 
                    role = "Relay / Bridge"
                    explanation = "Transmits the signal to the next layer without modification."
            else:
                role = "Isolated"
                explanation = "Does not play an active role in propagation."

            metrics.append({
                "Residue_ID": i,
                "Network_Layer": N,
                "In_Degree (m_i)": m_i,
                "Out_Degree (n_i)": n_i,
                "Bypassing_Edges (l)": l_edges,
                "Propagation_Coefficient (PC)": round(pc_i, 4),
                "Assigned_Role": role,
                "Biological_Interpretation": explanation
            })

    df_metrics = pd.DataFrame(metrics)
    # Katman ve PC değerine göre sırala. Kullanıcı sinyal yolunu adım adım PC skorlarına göre izleyebilir.
    df_metrics = df_metrics.sort_values(by=["Network_Layer", "Propagation_Coefficient (PC)"], ascending=[True, False])
    
    # İlgili dizin yoksa oluştur
    os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
    df_metrics.to_csv(output_csv, index=False)
    
    return df_metrics

def execute_network_metrics_workflow(G: nx.DiGraph, source_node: str, output_directory: str, prefix: str = ""):
    """
    Ana iş akışından kolayca çağırılabilmesi için yardımcı fonksiyon.
    Hesaplama sonrası matplotlib grafikleri ve PyMOL makrosu oluşturur.
    """
    csv_filename = os.path.join(output_directory, f"{prefix}propagation_metrics.csv")
    plot_filename = os.path.join(output_directory, f"{prefix}propagation_plot.png")
    pml_filename = os.path.join(output_directory, f"{prefix}pymol_heatmap.pml")
    
    #print(f"[*] Calculating Propagation Coefficients for source '{source_node}'. Output will be assigned to {csv_filename}")
    
    df = calculate_propagation_coefficient(G, source_node, csv_filename)
    if df is not None:
        #print("[PC] İlk 5 satır:\n", df.head())
        # Plot ve PyMOL senaryolarını tetikle
        generate_pc_bar_plot(df, plot_filename)
        generate_pymol_macro(df, pml_filename, source_node)
        #print(f"[+] PC analysis successfully completed and visual outputs saved.")
    #else:
        #print("[PC] DataFrame None döndü!")
    return df
