#efectocaballodetroya.streamlit.app
import streamlit as st
from rdkit import Chem

def count_negative_atoms(smiles):
    """
    A.
    """
    try:
        # Convierte el SMILES en un objeto mol
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "SMILES inválido. Por favor verifica la entrada."
        
        # Obtiene la lista de átomos y calcula los números de oxidación
        negative_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
        return f"Score: {negative_atoms}"
    except Exception as e:
        return f"Error al procesar el SMILES: {e}"

# Interfaz de Streamlit
st.title("Predictor de Efecto Caballo de Troya")
st.write("Autor: Jesus Alvarado-Huayhuaz")

st.write("""
Esta aplicación toma un código SMILES como entrada y calcula la energía de afinidad de la molécula en el receptor FhuE de A. Baumannii.
""")

# Entrada de SMILES
smiles_input = st.text_input("Introduce el código SMILES de la molécula:")

if smiles_input:
    result = count_negative_atoms(smiles_input)
    st.write(result)
