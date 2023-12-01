import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen




ch3_pattern = Chem.MolFromSmarts('[CH3]')
ch2_pattern = Chem.MolFromSmarts('[CH2]')
oh_pattern = Chem.MolFromSmarts('[OH]')
ch_pattern = Chem.MolFromSmarts('[CH]')
ch_pattern_ar = Chem.MolFromSmarts('[cH]')
bh_pattern = Chem.MolFromSmarts('[BH]')
bh_pattern_ar = Chem.MolFromSmarts('[bH]')
bh2_pattern = Chem.MolFromSmarts('[BH2]')
bh3_pattern = Chem.MolFromSmarts('[BH3]')
nh_pattern = Chem.MolFromSmarts('[NH]')
nh_pattern_ar = Chem.MolFromSmarts('[nH]')
nh2_pattern = Chem.MolFromSmarts('[NH2]')
nh3_pattern = Chem.MolFromSmarts('[NH3]')
sih_pattern = Chem.MolFromSmarts('[SiH]')
sih2_pattern = Chem.MolFromSmarts('[SiH2]')
sih3_pattern = Chem.MolFromSmarts('[SiH3]')
sih4_pattern = Chem.MolFromSmarts('[SiH4]')
ph_pattern = Chem.MolFromSmarts('[PH]')
ph_pattern_ar = Chem.MolFromSmarts('[pH]')
ph2_pattern = Chem.MolFromSmarts('[PH2]')
ph3_pattern = Chem.MolFromSmarts('[PH3]')
sh_pattern = Chem.MolFromSmarts('[SH]')
sh_pattern_ar = Chem.MolFromSmarts('[sH]')


# Function for counting bonds by type
def count_bond_types(smiles):
    mol = Chem.MolFromSmiles(smiles)
    bond_counts = {}

    if mol:
        for bond in mol.GetBonds():
            begin_atom = bond.GetBeginAtom().GetSymbol()
            end_atom = bond.GetEndAtom().GetSymbol()
            bond_type = bond.GetBondTypeAsDouble()

            bond_key = f"{begin_atom}-{end_atom} ({bond_type})"
            if bond_key in bond_counts:
                bond_counts[bond_key] += 1
            else:
                bond_counts[bond_key] = 1

    return bond_counts


def count_h_containing_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        group_counts = {
            'CH3': len(mol.GetSubstructMatches(ch3_pattern)),
            'CH2': len(mol.GetSubstructMatches(ch2_pattern)),
            'OH': len(mol.GetSubstructMatches(oh_pattern)),
            'CH': len(mol.GetSubstructMatches(ch_pattern)) + len(mol.GetSubstructMatches(ch_pattern_ar)),
            'BH': len(mol.GetSubstructMatches(bh_pattern)) + len(mol.GetSubstructMatches(bh_pattern_ar)),
            'Bh2': len(mol.GetSubstructMatches(bh2_pattern)),
            'Bh3': len(mol.GetSubstructMatches(bh3_pattern)),
            'NH': len(mol.GetSubstructMatches(nh_pattern)) + len(mol.GetSubstructMatches(nh_pattern_ar)),
            'NH2': len(mol.GetSubstructMatches(nh2_pattern)),
            'NH3': len(mol.GetSubstructMatches(nh3_pattern)),
            'SiH': len(mol.GetSubstructMatches(sih_pattern)),
            'SiH2': len(mol.GetSubstructMatches(sih2_pattern)),
            'SiH3': len(mol.GetSubstructMatches(sih3_pattern)),
            'SiH4': len(mol.GetSubstructMatches(sih4_pattern)),
            'PH': len(mol.GetSubstructMatches(ph_pattern)) + len(mol.GetSubstructMatches(ph_pattern_ar)),
            'PH2': len(mol.GetSubstructMatches(ph2_pattern)),
            'PH3': len(mol.GetSubstructMatches(ph3_pattern)),
            'SH': len(mol.GetSubstructMatches(sh_pattern))+ len(mol.GetSubstructMatches(sh_pattern_ar))
        }
        #I know that mixing the aliphatic and aromatic is not a very good idea, but let's keep it for now
        return group_counts

    else:
        return {}

def calculate_molecular_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Descriptors.MolWt(mol)
    else:
        return None


def count_aromatic_rings(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return len(Chem.GetSSSR(mol))
    return 0


def calculate_main_chain_length(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        # Find carbons
        carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
        if carbon_atoms:
            start_atom = max(carbon_atoms, key=lambda x: len(mol.GetAtomWithIdx(x).GetNeighbors()))
            visited_atoms = set()
            max_length = 0
            stack = [(start_atom, 0)]
            while stack:
                atom_idx, length = stack.pop()
                visited_atoms.add(atom_idx)
                max_length = max(max_length, length)

                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited_atoms:
                        stack.append((neighbor_idx, length + 1))

            return max_length
    return 0


def count_unique_elements(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return len(set(atom.GetSymbol() for atom in mol.GetAtoms()))
    return 0



def create_dataframe(SMILES):
    #Import df
    df = pd.DataFrame({"smiles":SMILES})

    #Create smiles to delete
    df["smiles_to_delete"] = df["smiles"]

    #Create one column per element
    df["Number of B"] = 0
    df["Number of C"] = 0
    df["Number of N"] = 0
    df["Number of O"] = 0
    df["Number of F"] = 0
    df["Number of Si"] = 0
    df["Number of P"] = 0
    df["Number of S"] = 0
    df["Number of Cl"] = 0
    df["Number of Br"] = 0
    df["Number of I"] = 0

    #Remove elements that appears less than 10 times
    element_to_remove = ["Na", "Mg", "Al", "Ca", "Ti", "Cr", "[Co]", "Zn", "Ge", "As", "Se", "Ag", "[Hg]"]

    indices_to_remove = []

    for index, row in df.iterrows():
        for element in element_to_remove:
            if element in row["smiles_to_delete"]:
                indices_to_remove.append(index)

    # Delete the corresponding rows
    df.drop(indices_to_remove, inplace=True)

    # Reset index
    df.reset_index(drop=True, inplace=True)

    #Remove special characters
    char_to_remove = ["(", ")", "=", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0", "@", "[", "]", "-", "+", "H",
                  "/", "#", ".", "\\", "%", " "]
                  #We remove H car it is not suppose to be in smile

    for index, row in df.iterrows():
        cleaned_smiles = ''.join([char for char in row["smiles_to_delete"] if char not in char_to_remove])
        df.loc[index, "smiles_to_delete"] = cleaned_smiles

    #Count the number of element
    element_to_count = ["Br", "B", "Cl", "C", "N", "O", "F", "Si", "P", "S", "I"]  # The order matters

    for index, row in df.iterrows():
        for element in element_to_count:
            count = row["smiles_to_delete"].count(element.lower()) + row["smiles_to_delete"].count(element)
            df.loc[index, "Number of " + element] = count
            row["smiles_to_delete"] = row["smiles_to_delete"].replace(element.lower(), "").replace(element, "")

        df.loc[index, "smiles_to_delete"] = row["smiles_to_delete"]

    #Count the type and number of bonds
    df['Bond Counts'] = df['smiles'].apply(count_bond_types)

    df['Bond Counts'] = df['Bond Counts'].apply(lambda x: {k: x.get(k, 0) for k in set().union(*df['Bond Counts'])})
    df_train_bondtomerge = pd.DataFrame(df['Bond Counts'].to_list())

    df = df.merge(df_train_bondtomerge, left_index=True, right_index=True)

    # Apply the function count_h_containing_groups to every smiles
    df['Functional Groups with H'] = df['smiles'].apply(count_h_containing_groups)

    # Replace missing values by 0
    df['Functional Groups with H'] = df['Functional Groups with H'].apply(lambda x: {k: x.get(k, 0) for k in set().union(*df['Functional Groups with H'])})

    # Create column for every key in the dictionnary
    df_train_bondtomerge = pd.DataFrame(df['Functional Groups with H'].to_list())
    df = df.merge(df_train_bondtomerge, left_index=True, right_index=True)

    df['Molecular_weight'] = df['smiles'].apply(calculate_molecular_weight)

    df['Aromatic Rings Count'] = df['smiles'].apply(count_aromatic_rings)

    df['Main Chain Length'] = df['smiles'].apply(calculate_main_chain_length)

    df['Nombre d\'éléments différents'] = df['smiles'].apply(count_unique_elements)

    df['Nombre de doubles liaisons'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x).GetNumBonds(Chem.BondType.DOUBLE) if Chem.MolFromSmiles(x) is not None else 0)

    df['XLogP'] = df['smiles'].apply(lambda x: Crippen.MolLogP(Chem.MolFromSmiles(x)) if Chem.MolFromSmiles(x) is not None else 0)

    df = df.drop(columns=["smiles_to_delete", "Bond Counts", "Functional Groups with H"
                                  ])

    return df
