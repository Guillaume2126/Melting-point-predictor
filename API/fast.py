from fastapi import FastAPI
import pandas as pd
import joblib
from Preprocessing.Create_dataframe import create_dataframe

app = FastAPI()

model_path = "Model/best_full_data_model.joblib"
app.state.model = joblib.load(model_path)

@app.get("/")
def root():
    return {
        'message': "API is ON"
    }


@app.get("/predict")
def predict(SMILES: str):
    """
    Make a single melting point prediction.
    Assumes `SMILES` is provided as a string by the user.
    """

    df = create_dataframe(SMILES)

    #dataset cols
    to_keep_cols = [
        'name', 'smiles', 'mpC', 'Number of B', 'Number of C', 'Number of N', 'Number of O',
        'Number of F', 'Number of Si', 'Number of P', 'Number of S', 'Number of Cl', 'Number of Br',
        'Number of I', 'S-P (2.0)', 'N-S (1.0)', 'Si-O (1.0)', 'Si-N (1.0)', 'S-Br (1.0)', 'O-Si (1.0)',
        'P-S (1.0)', 'H-C (1.0)', 'F-P (1.0)', 'C-B (1.0)', 'C-S (1.0)', 'S-Cl (1.0)', 'F-B (1.0)',
        'N-Si (1.0)', 'P-C (1.0)', 'C-O (2.0)', 'I-O (1.0)', 'H-O (1.0)', 'P-Br (1.0)', 'S-O (1.0)',
        'F-C (1.0)', 'C-Si (1.0)', 'N-N (1.5)', 'N-H (1.0)', 'S-N (1.5)', 'N-C (3.0)', 'Si-Si (1.0)',
        'F-Si (1.0)', 'P-O (2.0)', 'P-S (2.0)', 'N-B (1.0)', 'Cl-I (1.0)', 'N-C (2.0)', 'P-N (1.0)',
        'F-S (1.0)', 'Si-P (1.0)', 'C-C (3.0)', 'Cl-P (1.0)', 'Br-Br (1.0)', 'Br-P (1.0)', 'N-N (2.0)',
        'C-C (1.5)', 'C-F (1.0)', 'Cl-Si (1.0)', 'N-I (1.0)', 'S-N (1.0)', 'B-Cl (1.0)', 'Br-B (1.0)',
        'C-H (1.0)', 'P-N (2.0)', 'B-C (1.0)', 'B-Br (1.0)', 'O-S (1.0)', 'B-S (1.0)', 'C-N (2.0)',
        'O-N (1.0)', 'N-P (2.0)', 'H-N (1.0)', 'S-S (1.0)', 'Cl-C (1.0)', 'N-N (3.0)', 'C-P (2.0)',
        'B-F (1.0)', 'C-I (1.0)', 'S-I (1.0)', 'C-N (1.5)', 'B-O (1.0)', 'O-N (1.5)', 'O-P (1.0)',
        'B-B (1.0)', 'O-C (1.5)', 'O-B (1.0)', 'C-Cl (1.0)', 'N-O (1.5)', 'O-C (3.0)', 'Si-Br (1.0)',
        'N-P (1.0)', 'Br-Si (1.0)', 'Br-C (1.0)', 'O-I (1.0)', 'Si-F (1.0)', 'C-O (1.5)', 'P-Cl (1.0)',
        'I-I (1.0)', 'S-C (1.5)', 'O-O (2.0)', 'Si-Cl (1.0)', 'O-H (1.0)', 'C-C (1.0)', 'S-S (1.5)',
        'N-S (2.0)', 'P-Si (1.0)', 'N-C (1.0)', 'N-N (1.0)', 'S-P (1.0)', 'Cl-S (1.0)', 'I-C (1.0)',
        'N-O (1.0)', 'O-C (2.0)', 'P-F (1.0)', 'O-N (2.0)', 'S-O (2.0)', 'Cl-O (1.0)', 'Cl-N (1.0)',
        'C-C (2.0)', 'S-F (1.0)', 'C-N (1.0)', 'O-C (1.0)', 'N-C (1.5)', 'C-Br (1.0)', 'N-Cl (1.0)',
        'O-P (2.0)', 'C-P (1.0)', 'C-S (1.5)', 'S-C (1.0)', 'O-S (2.0)', 'N-O (2.0)', 'Cl-B (1.0)',
        'O-O (1.0)', 'Si-C (1.0)', 'N-S (1.5)', 'N-Br (1.0)', 'O-Cl (1.0)', 'C-O (1.0)', 'C-N (3.0)',
        'C-O (3.0)', 'C-S (2.0)', 'Cl-Cl (1.0)', 'S-C (2.0)', 'B-N (1.0)', 'P-O (1.0)', 'SiH4', 'OH',
        'SiH2', 'PH3', 'CH', 'SiH', 'SiH3', 'Bh2', 'PH', 'SH', 'PH2', 'CH2', 'NH', 'Bh3', 'NH2', 'NH3',
        'CH3', 'BH', 'Molecular_weight', 'Aromatic Rings Count', 'Main Chain Length',
        'Nombre d\'éléments différents', 'Nombre de doubles liaisons', 'XLogP']

    # del generated col that we didnt use to train the model
    df = df[[col for col in df.columns if col in to_keep_cols]]

    # add missing cols that have not been generated creating the pred_df from smiles
    for col in to_keep_cols:
        if col not in df.columns:
            df[col] = 0

    df = df[to_keep_cols]


    # List of columns to remove (model pre-preprocessing)
    columns_to_remove = ['O-Si (1.0)', 'N-Si (1.0)', 'Br-P (1.0)', 'C-F (1.0)', 'H-N (1.0)',
                         'C-I (1.0)', 'C-Cl (1.0)', 'O-I (1.0)', 'C-O (1.5)', 'P-Si (1.0)',
                         'N-C (1.5)', 'C-Br (1.0)', 'C-S (1.5)', 'Si-C (1.0)', 'Bh2', 'BH',
                         'Nombre de doubles liaisons']

    #drop not used cols
    df.drop(columns=columns_to_remove, errors="ignore", inplace=True)

    X = df.drop(columns=["smiles"])

    y_pred = app.state.model.predict(X)

    return dict(melting_point=round(float(y_pred), 2))
