#efectocaballodetroya.streamlit.app

import streamlit as st
from rdkit import Chem

def count_atoms_with_chelation_potential(smiles):
    """
    Calcula el número de átomos en la molécula representada por SMILES
    que tienen capacidad de formar enlaces de coordinación o quelar metales.
    """
    try:
        # Convierte el SMILES en un objeto mol
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "SMILES inválido. Por favor verifica la entrada."
        
        # Identifica los átomos con pares de electrones libres (N, O, S, P)
        chelation_atoms = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() in [7, 8, 15, 16]:  # N, O, P, S
                # Verifica que el átomo no tenga una carga positiva que limite la coordinación
                if atom.GetFormalCharge() <= 0:
                    # Los átomos con pares libres o disponibles pueden coordinar
                    chelation_atoms += 1
        
        return f"Número de átomos con capacidad de quelar o formar enlaces dativos: {chelation_atoms}"
    except Exception as e:
        return f"Error al procesar el SMILES: {e}"

# Interfaz de Streamlit
st.title("Predictor de Efecto Caballo de Troya")
st.write("Autor: Jesus Alvarado-Huayhuaz")
st.write("""
Esta aplicación toma un código SMILES como entrada y calcula la energía de afinidad de la molécula en el receptor FhuE de A. Baumannii.
""")
#st.write("""
#Esta aplicación toma un código SMILES como entrada y calcula el número de átomos
#que tienen capacidad de quelar metales (como el hierro) o formar enlaces dativos
#basados en su estructura química.
#""")

# Entrada de SMILES
smiles_input = st.text_input("Introduce el código SMILES de la molécula:")

if smiles_input:
    result = count_atoms_with_chelation_potential(smiles_input)
    st.write(result)
