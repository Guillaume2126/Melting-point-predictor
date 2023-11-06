import pandas as pd
import numpy as np
from rdkit import Chem


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


def create_detaframe(df):
    #Import df
    df = pd.read_csv("'"+df+"'")

    #TODO : creer df vide et faire feature engineering

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
    df_train['Bond Counts'] = df_train['smiles'].apply(count_bond_types)

    df_train['Bond Counts'] = df_train['Bond Counts'].apply(lambda x: {k: x.get(k, 0) for k in set().union(*df_train['Bond Counts'])})
    df_train_bondtomerge = pd.DataFrame(df_train['Bond Counts'].to_list())

    df_train = df_train.merge(df_train_bondtomerge, left_index=True, right_index=True)
