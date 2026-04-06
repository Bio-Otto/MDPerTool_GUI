
# --- Fonksiyonun bağımsız kopyası ---
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def generate_pc_bar_plot(df_metrics: pd.DataFrame, output_image: str):
    try:
        df_plot = df_metrics[df_metrics['Propagation_Coefficient (PC)'] > 0].copy()
        df_plot = df_plot.sort_values('Propagation_Coefficient (PC)', ascending=False).head(30)

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
        plt.savefig(output_image, dpi=300)
        plt.close()
    except Exception as e:
        print(f"[PC Plot Hatası] Bar plot çizilemedi: {e}")

def test_generate_pc_bar_plot():
    data = {
        'Residue_ID': ['ALA10', 'GLY20', 'SER30', 'THR40', 'LEU50'],
        'Propagation_Coefficient (PC)': [0.25, 0.0, 0.15, 0.4, 0.05]
    }
    df = pd.DataFrame(data)
    output_image = 'test_pc_plot.png'
    print('DataFrame to plot:')
    print(df)
    generate_pc_bar_plot(df, output_image)
    print(f'Plot saved to {output_image}')

if __name__ == '__main__':
    test_generate_pc_bar_plot()
