import os
import numpy as np
import pandas as pd
import networkx as nx
from scipy.optimize import curve_fit

"""
Bu script, MDPerTool ağ ve enerji analizini matematiksel eğri uydurma ve 
merkeziyet (hub) tespiti ile zenginleştirmek için deneysel bir test modülüdür.
"""

def boltzmann_two_state(x, A1, A2, x0, dx):
    """
    Two-state Boltzmann Dağılım Denklemi
    Çoğu Enerji Disipasyonu (Energy Dissipation) makalesinde kullanılan standart fonksiyon.
    :param x: Zaman (veya frame)
    :param A1: Başlangıç plato değeri (Initial state)
    :param A2: Bitiş plato değeri (Final state)
    :param x0: Eğrinin orta noktası, yani Yarı-Yanıt Süresi (t0 - half response time)
    :param dx: Zaman sabiti (Curve steepness / dissipation rate indicator)
    """
    return A2 + (A1 - A2) / (1 + np.exp((x - x0) / dx))

def fit_boltzmann_to_node_data(time_array, energy_array):
    """
    Tek bir amino asidin zaman-enerji verisine Boltzmann fonksiyonunu uydurur.
    t0 (yarı-yanıt süresi) değerini çeker.
    """
    try:
        # Başlangıç tahminleri (A1, A2, x0, dx) - Veriye göre dinamik ayarlanır
        p0 = [np.max(energy_array), np.min(energy_array), np.median(time_array), np.std(time_array)/2]
        
        popt, _ = curve_fit(boltzmann_two_state, time_array, energy_array, p0=p0, maxfev=10000)
        
        A1, A2, t0, dx = popt
        return {
            "A1": A1, 
            "A2": A2, 
            "t0_half_response": t0, 
            "dx_dissipation_rate": dx,
            "fitted": True
        }
    except Exception as e:
        return {"A1": np.nan, "A2": np.nan, "t0_half_response": np.nan, "dx_dissipation_rate": np.nan, "fitted": False}

def calculate_propagation_coefficient(G: nx.DiGraph, source_node, output_csv="propagation_coefficients.csv"):
    """
    Makaleye sadık kalınarak hazırlanmış Propagation Coefficient (PC) hesaplaması.
    Öncelikle ağ üzerinde 'Deep-preference hierarchical layout' (Katmanlı yapı) 
    oluşturulur (Source'dan uzaklığa göre BFS katmanları).
    Daha sonra her N katmanı için PC(i) hesaplanır.
    """
    # 1. Hiyerarşik Katmanları Belirle (Shortest Path / BFS uzaklığı)
    try:
        layers = nx.single_source_shortest_path_length(G, source_node)
    except Exception as e:
        print("[Hata] Source node ağda bulunamadı veya ağ bağlantısız!", e)
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
        # L'yi hesapla (Katman N'yi transit geçen kenarlar / edges that pass through layer N)
        # Yani başlangıç düğümünün katmanı < N ve bitiş düğümünün katmanı > N olan kenarlar.
        l_edges = 0
        for u, v in G.edges():
            u_layer = G.nodes[u].get('layer', -1)
            v_layer = G.nodes[v].get('layer', -1)
            # Eğer u daha önceki bir katmanda ve v daha sonraki bir katmandaysa bu 'l' sayılır
            if u_layer != -1 and v_layer != -1:
                if u_layer < N and v_layer > N:
                    l_edges += 1
        
        # Paydadaki Düğüm Üzerinden Geçen Yolların Toplamı: sum(m_j * n_j)
        layer_pathways_sum = sum([G.in_degree(j) * G.out_degree(j) for j in nodes_in_layer_N])
        
        # Payda (Denominator): (Katmanda gerçekleşen tüm yollar + Transit geçenler)
        denominator = layer_pathways_sum + l_edges
        
        # Her düğüm için PC(i) Puanlarını Hesapla ve Yorumla
        for i in nodes_in_layer_N:
            m_i = G.in_degree(i)
            n_i = G.out_degree(i)
            numerator = m_i * n_i
            
            # Payda sıfır ise (örn. son katmanda giden yol yoksa), PC 0 alınır
            pc_i = (numerator / denominator) if denominator > 0 else 0.0
            
            # Kullanıcı için Yorumlama Mantığı (Interpretation logic)
            # PC skoru ve oranlar incelenerek rezidünün ağdaki biyolojik/yapısal rolü belirlenir.
            if m_i == 0 and n_i > 0:
                role = "Initiator / Source (Sinyal Başlatıcı)"
                explanation = "Sisteme enerji/sinyal zerk eden başlangıç noktası."
            elif n_i == 0 and m_i > 0:
                role = "Sink / Target (Sinyal Sönümleyici)"
                explanation = "Sinyalin biriktiği ve sonlandığı hedef nokta."
            elif pc_i > 0:
                if n_i > m_i:
                    role = "Amplifier / Hub (Sinyal Çoğaltıcı Darboğaz)"
                    explanation = f"Gelen sinyali daha geniş bir alana dağıtıyor. PC Skoru: {pc_i:.3f}. Allosterik yolakta kritik bir dağıtım üssü (Super-Hub)."
                elif m_i > n_i:
                    role = "Funnel / Collector (Sinyal Toplayıcı Huni)"
                    explanation = f"Farklı kollardan gelen sinyalleri tek kanala indirgiyor. PC Skoru: {pc_i:.3f}. Sinyal izolasyonunda önemli."
                else: # m_i == n_i
                    role = "Relay / Transmission Bridge (Sinyal Aktarıcı Köprü)"
                    explanation = f"Sinyali değiştirmeden bir sonraki katmana iletiyor. PC Skoru: {pc_i:.3f}."
            else:
                role = "Isolated / Inactive (İzole/İnaktif)"
                explanation = "Ağda aktif bir iletim rolü üstlenmiyor."

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
    # Katman ve PC değerine göre sırala. Böylece kullanıcı sinyal yolunu adım adım PC skorlarına göre izleyebilir.
    df_metrics = df_metrics.sort_values(by=["Network_Layer", "Propagation_Coefficient (PC)"], ascending=[True, False])
    
    df_metrics.to_csv(output_csv, index=False)
    print(f"\n[BAŞARILI] Propagation Coefficient Analizi Hesaplandı -> {output_csv}")
    return df_metrics

if __name__ == "__main__":
    # --- 1) NETWORK HUB/MERKEZİYET TEST SENARYOSU ---
    print("--- Directed Graph Hub Açısından Propagation Coefficient (PC) Testi ---")
    G_test = nx.DiGraph()
    # Örnek bir Sinyal Transfer Ağı
    G_test.add_edges_from([
        ("Source", "A"), ("Source", "B"), ("Source", "C"), # Source = Layer 0, A/B/C = Layer 1
        ("A", "D"), ("B", "E"), ("C", "E"), ("B", "D"),    # D/E = Layer 2
        ("A", "Hub_Node"), ("E", "Hub_Node"),              # Hub_Node = Layer 3
        ("Hub_Node", "Target1"), ("Hub_Node", "Target2"),  # Target1/2 = Layer 4
        ("Source", "E") # Layer 0'dan Layer 2'ye giden (Layer 1'i baypas eden / l hücresini arttıran) kenar
    ])
    
    test_csv_path = "test_pc_metrics.csv"
    df = calculate_propagation_coefficient(G_test, source_node="Source", output_csv=test_csv_path)
    
    # Sonuçların Özeti
    print(df.to_string(index=False))