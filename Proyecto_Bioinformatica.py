import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import requests
import tempfile
from biopandas.pdb import PandasPdb
import py3Dmol  # Asegúrate de que py3Dmol esté instalado
import os

# Configuración de la página
st.set_page_config(page_title="Bioinformática en el area  de la Biomedica", 
                   page_icon=":cat:", layout="wide")

# Título de la app en Streamlit
st.title("Bioinformática en el area de la Biomédica")
st.markdown("""
Bienvenido a esta aplicacion donde gracias a la bioinformatica puedes visualizar de manera dinamica las propiedades principales de la proteina de su eleccion



""")

# Barra lateral para navegación
st.sidebar.title("Buscar")

# Caja de texto para ingresar el código de la proteína (PDB ID)
pdb_code = st.sidebar.text_input("Ingresa el código de la proteína (PDB ID)", value="1A8M")  # Valor por defecto "1A8M"

# Función para obtener el archivo PDB desde la URL de RCSB
def obtener_pdb(pdb_code):
    url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        st.sidebar.error(f"No se pudo obtener el archivo PDB para el código: {pdb_code}")
        return None

# Intentar obtener el archivo PDB basado en el código ingresado
pdb_data = obtener_pdb(pdb_code)

# Función para guardar el archivo PDB en un archivo temporal
def guardar_pdb_temporal(pdb_data):
    # Crear un archivo temporal
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
    with open(temp_file.name, 'w') as f:
        f.write(pdb_data)
    return temp_file.name

# Función para calcular el porcentaje de GC (guanina y citosina) en fragmentos
def calcular_gc_por_fragmento(pdb_data, fragment_size=10):
    gc_percentages = []
    secuencia = ""
    
    # PDB puede contener cadenas de ácidos nucleicos (ARN/ADN) en las líneas SEQRES
    for line in pdb_data.splitlines():
        if line.startswith("SEQRES"):  # Las líneas SEQRES contienen las secuencias de cadenas
            secuencia += line[19:].strip()  # Extraemos la secuencia después de la palabra SEQRES
    
    if len(secuencia) == 0:  # Si no se encontró secuencia de ADN/ARN
        return None

    # Dividir la secuencia en fragmentos y calcular el GC para cada fragmento
    for i in range(0, len(secuencia), fragment_size):
        fragment = secuencia[i:i+fragment_size]
        gc_count = fragment.count('G') + fragment.count('C')
        total_count = len(fragment)
        
        if total_count > 0:
            gc_percentage = (gc_count / total_count) * 100
            gc_percentages.append(gc_percentage)

    return gc_percentages

# Intentar obtener el archivo PDB y procesarlo
if pdb_data:
    # Guardar el archivo PDB en un archivo temporal
    pdb_file_path = guardar_pdb_temporal(pdb_data)
    
    # Leer el archivo PDB usando PandasPdb
    ppdb = PandasPdb().read_pdb(pdb_file_path)
    st.sidebar.success(f"Archivo PDB cargado correctamente: {pdb_code}")
else:
    ppdb = None

# Función para mostrar las secciones del archivo PDB
def datos_pdb():
    if ppdb is not None:
        # Mejoras visuales en la sección de información general
        st.markdown("---")
        st.subheader("Informacion del PDB")
        st.markdown("Explora las diferentes propieades del archivo PDB de tu eleccion:")

        for section, df in ppdb.df.items():
            st.write(f"### {section}")
            st.dataframe(df.head(), width=800)  # Estilo de tabla interactiva con un ancho adecuado

        # Información sobre los diferentes tipos de datos en el PDB
        st.markdown("#### Información sobre las dimensiones de las secciones del archivo PDB:")
        df_hetatm = ppdb.df["HETATM"].shape
        df_anistropic = ppdb.df["ANISOU"].shape
        df_others = ppdb.df["OTHERS"].shape
        st.write(f"Tamaño de HETATM: {df_hetatm}")
        st.write(f"Tamaño de ANISOU: {df_anistropic}")
        st.write(f"Tamaño de OTHERS: {df_others}")

        # Mostrar información sobre la sección "ATOM"
        st.markdown("---")
        df_atom = ppdb.df["ATOM"]
        st.subheader("Datos ATOM")
        st.write(f"Tamaño del DataFrame ATOM: {df_atom.shape}")
        st.write("Primeras filas de la sección ATOM:")
        st.dataframe(df_atom.head(), width=800)
    else:
        st.warning("No se pudo cargar el archivo PDB.")

# Función para la sección "Visualización 3D con Plotly"
def visualizacion_3d_plotly():
    if ppdb is not None:
        # Visualización 3D de las coordenadas atómicas con Plotly
        st.markdown("---")
        st.subheader("Visualización 3D de las coordenadas atómicas (Plotly)")
        df_atom = ppdb.df["ATOM"]
        fig = px.scatter_3d(df_atom, x="x_coord", y="y_coord", z="z_coord", color="element_symbol", 
                            template="plotly_dark", color_discrete_sequence=["white", "gray", "red", "yellow"])
        fig.update_coloraxes(showscale=True)
        fig.update_traces(marker=dict(size=3))
        st.plotly_chart(fig, use_container_width=True)

        # Visualización 3D de las coordenadas atómicas coloreadas por el nombre del residuo
        st.markdown("---")
        st.subheader("Visualización 3D por Nombre de Residuo (Plotly)")
        fig = px.scatter_3d(df_atom, x="x_coord", y="y_coord", z="z_coord", color="residue_name", template="plotly_dark")
        fig.update_traces(marker=dict(size=3))
        fig.update_coloraxes(showscale=True)
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning("No se pudo cargar la visualización 3D.")

# Función para la sección "Visualización 3D de Proteína (py3Dmol)"
def visualizacion_3d_proteina():
    if ppdb is not None:
        # Crear una visualización 3D con py3Dmol
        viewer = py3Dmol.view(width=800, height=600)
        viewer.addModel(pdb_data, "pdb")  # Cargar el modelo PDB
        viewer.setStyle({'cartoon': {}})  # Estilo de visualización en formato cartoon
        viewer.setBackgroundColor("white")  # Fondo blanco
        viewer.zoomTo()  # Ajuste para ver toda la proteína

        # Convertir la visualización a HTML
        viewer_html = viewer._make_html()

        # Mostrar el HTML en Streamlit usando st.components.v1.html
        st.components.v1.html(viewer_html, height=600)
    else:
        st.warning("No se pudo cargar la visualización 3D.")

# Llamar a las funciones para mostrar los datos y visualizaciones
datos_pdb()
visualizacion_3d_plotly()

# Mostrar el contenido en la misma pestaña
st.markdown("---")
st.markdown("### Análisis y Visualización de la proteína")

visualizacion_3d_proteina()


# Mostrar el porcentaje de GC en un gráfico de línea
if pdb_data:
    gc_percentages = calcular_gc_por_fragmento(pdb_data, fragment_size=10)  # Fragmento de 10 residuos

    if gc_percentages is not None:
        st.markdown("---")
        st.subheader(f"Patrón de contenido GC a lo largo de la secuencia")
        
        # Crear un DataFrame con los porcentajes de GC para graficar
        df_gc = pd.DataFrame({
            "Fragmento": range(len(gc_percentages)),
            "Porcentaje GC": gc_percentages
        })
        
        # Mostrar el gráfico de línea con el porcentaje de GC
        st.line_chart(df_gc.set_index("Fragmento"))
    else:
        st.warning("No se encontró secuencia de ARN/ADN en el archivo PDB para calcular el porcentaje de GC.")

# Pie de página con información adicional
st.markdown("---")
st.markdown("""
### Información adicional:
Si tienes alguna duda sobre el uso de esta aplicación, no dudes en hacer contacto con nosotros. 

#### Créditos de este trabajo:
* Datos obtenidos de la base de datos PDB.
* Visualización realizada con Plotly, py3Dmol y Streamlit.
* Este trabajo fue realizado en Visual Studio Code.
""")
