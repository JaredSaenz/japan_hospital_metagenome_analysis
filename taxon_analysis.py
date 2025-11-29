import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ===================== CONFIGURACIÓN =====================
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Ruta automática: misma carpeta que el CSV
csv_path = "DRR199648_focusTax/output_All_levels.csv"
df = pd.read_csv(csv_path)

# Directorio donde se guardarán las figuras (misma carpeta del CSV)
output_dir = os.path.dirname(csv_path)
if not output_dir:
    output_dir = "."  # si estás en la carpeta actual

# Renombrar columna de abundancia
df.rename(columns={"final.contigs.fasta": "Abundance"}, inplace=True)
df["Abundance"] = pd.to_numeric(df["Abundance"], errors='coerce')
df = df.dropna(subset=["Abundance"])

print(f"Total de genomas detectados: {len(df)}")
print(f"Abundancia total: {df['Abundance'].sum():.2f}%\n")

# ===================== FUNCIÓN PARA GRAFICAR Y GUARDAR =====================
def plot_top_n(level, n=15, title=None, filename=None):
    grouped = df.groupby(level)["Abundance"].sum().sort_values(ascending=False)
    top = grouped.head(n)
    others = pd.Series({"Others": grouped[n:].sum()}) if len(grouped) > n else pd.Series()
    data = pd.concat([top, others])
    
    plt.figure(figsize=(12, 8))
    ax = sns.barplot(x=data.values, y=data.index, palette="viridis")
    plt.title(title or f"Top {n} - {level}", fontsize=16, pad=20)
    plt.xlabel("Abundancia relativa (%)")
    plt.ylabel(level)
    
    # Añadir valores en las barras
    for i, v in enumerate(data.values):
        if v > 0.1:
            ax.text(v + 0.1, i, f"{v:.2f}%", va='center', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    
    # Guardar en PNG y PDF de alta calidad
    base_name = filename or f"Top_{n}_{level.replace(' ', '_')}"
    plt.savefig(os.path.join(output_dir, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f"{base_name}.pdf"), bbox_inches='tight')
    print(f"Guardado: {base_name}.png y {base_name}.pdf")
    
    plt.show()

# ===================== GRÁFICOS + GUARDADO AUTOMÁTICO =====================
print("=== PHYLA más abundantes ===")
plot_top_n("Phylum", n=10, title="Abundancia por Phylum", filename="01_Phylum_Top10")

print("\n=== CLASES más abundantes ===")
plot_top_n("Class", n=12, title="Abundancia por Clase", filename="02_Class_Top12")

print("\n=== FAMILIAS más abundantes ===")
plot_top_n("Family", n=15, title="Abundancia por Familia", filename="03_Family_Top15")

print("\n=== GÉNEROS más abundantes ===")
plot_top_n("Genus", n=15, title="Abundancia por Género", filename="04_Genus_Top15")

# Top 20 especies/cepas
print("\n=== TOP 20 ESPECIES/STRIANS ===")
top_species = df.sort_values("Abundance", ascending=False)[["Species", "Strain", "Abundance"]].head(20)
plt.figure(figsize=(12, 9))
ax = top_species.plot(kind='barh', x='Strain', y='Abundance', color='skyblue', legend=False)
plt.title("Top 20 genomas detectados por abundancia", fontsize=16, pad=20)
plt.xlabel("Abundancia relativa (%)")
plt.ylabel("Cepa")
for i, v in enumerate(top_species["Abundance"]):
    ax.text(v + 0.1, i, f"{v:.2f}%", va='center', fontsize=11, fontweight='bold')
plt.tight_layout()

# Guardar
plt.savefig(os.path.join(output_dir, "05_Top20_Species.png"), dpi=600, bbox_inches='tight')
plt.savefig(os.path.join(output_dir, "05_Top20_Species.pdf"), bbox_inches='tight')
print("Guardado: 05_Top20_Species.png y .pdf")
plt.show()

# ===================== RESUMEN FINAL =====================
print("\n" + "="*60)
print("RESUMEN TAXONÓMICO - FIGURAS GUARDADAS EN:")
print(output_dir)
print("="*60)
print(f"Phylum dominante: Proteobacteria ({df[df['Phylum']=='Proteobacteria']['Abundance'].sum():.2f}%)")
print(f" ├ Gammaproteobacteria: {df[df['Class']=='Gammaproteobacteria']['Abundance'].sum():.2f}%")
print(f" └ Familia dominante: Enterobacteriaceae ({df[df['Family']=='Enterobacteriaceae']['Abundance'].sum():.2f}%)")
enterob = df[df['Family']=='Enterobacteriaceae'].groupby("Genus")["Abundance"].sum().sort_values(ascending=False)
print(f" → Géneros principales en Enterobacteriaceae:")
for genus, abund in enterob.head(8).items():
    print(f"   • {genus}: {abund:.2f}%")
print("="*60)
