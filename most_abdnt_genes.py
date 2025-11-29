import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

# Nombre del archivo de entrada
archivo_entrada = 'DRR199648_diamond/card_proteins_dmnd.out.tsv'

# Verificar si el archivo existe
if not os.path.exists(archivo_entrada):
    print(f"ERROR: No encuentro el archivo '{archivo_entrada}'")
    # Si estás en Jupyter, puedes comentar la línea de sys.exit()
    # sys.exit(1) 
else:
    print("Archivo encontrado. Iniciando análisis...")

    # 1. Cargar los datos
    # DIAMOND genera un archivo sin encabezado, así que definimos los nombres de columna
    print("Cargando datos...")
    columnas = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    try:
        df = pd.read_csv(archivo_entrada, sep='\t', names=columnas)
        
        # 2. Limpiar los nombres de los genes
        # El nombre en CARD viene largo (ej: "gnl|CARD|ARO:3000123|NombreGen"). 
        # Queremos solo el nombre final.
        print("Procesando nombres de genes...")
        def limpiar_nombre(texto):
            if isinstance(texto, str) and '|' in texto:
                texto = texto.split('|')[-1]
                if isinstance(texto, str) and ':' in texto:
                	return texto.split(':')[-1] # Toma lo que está después del último ':'
                return texto.split('|')[-1] # Toma lo que está después del último '|'
            return texto

        df['Gen_Resistencia'] = df['sseqid'].apply(limpiar_nombre)

        # 3. Contar los 15 genes más frecuentes
        top_genes = df['Gen_Resistencia'].value_counts().head(15).reset_index()
        top_genes.columns = ['Gen', 'Conteo']

        # 4. Crear la gráfica
        print("Generando gráfica...")
        plt.figure(figsize=(12, 8)) # Tamaño de la imagen (ancho, alto)
        # Usamos la paleta 'viridis' que se ve muy profesional
        sns.barplot(data=top_genes, x='Conteo', y='Gen', palette='viridis')

        # Etiquetas y Título
        plt.title('Genes de Resistencia Antimicrobiana (ARGs) más abundantes', fontsize=16)
        plt.xlabel('Número de secuencias detectadas', fontsize=12)
        plt.ylabel('Gen de Resistencia', fontsize=12)
        plt.grid(axis='x', linestyle='--', alpha=0.7) # Cuadrícula suave para leer mejor
        plt.tight_layout()

        # 5. Guardar la imagen
        nombre_imagen = 'Figura_20genesResistencia.png'
        plt.savefig(nombre_imagen, dpi=300) # dpi 300 es alta calidad para tu reporte
        plt.show() # Muestra la gráfica en pantalla si usas Jupyter
        
        print(f"¡Listo! Gráfica guardada como: {nombre_imagen}")

    except Exception as e:
        print(f"Ocurrió un error: {e}")
