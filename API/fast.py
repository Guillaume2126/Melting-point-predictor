from fastapi import FastAPI
import pandas as pd
import joblib
from Preprocessing import Create_dataframe

app = FastAPI()

model_path = "Model/best_full_data_model.joblib"
app.state.model = joblib.load(model_path)

@app.get("/predict")
def predict(SMILES: str):
    """
    Make a single melting point prediction.
    Assumes `SMILES` is provided as a string by the user.
    """

    df = Create_dataframe.create_dataframe(SMILES)

    # COLONNES A SUPPRIMER CAR DANS LE PROCESSING DU NOTEBOOK DU MODELE
    # (VOIR NOTEBOOK LCCOPY_2) J'ENLEVE CES COLONNES MANUELLEMENT
    # CAR AU DESSUS DU SEUIL DE CORRELATION DE 90%.
    # LE RESTE DU PREPROCESSING EST DANS LA PIPELINE

    columns_to_remove = ['O-Si (1.0)', 'N-Si (1.0)', 'Br-P (1.0)', 'C-F (1.0)', 'H-N (1.0)',
                         'C-I (1.0)', 'C-Cl (1.0)', 'O-I (1.0)', 'C-O (1.5)', 'P-Si (1.0)',
                         'N-C (1.5)', 'C-Br (1.0)', 'C-S (1.5)', 'Si-C (1.0)', 'Bh2', 'BH',
                         'Nombre de doubles liaisons']

    df = df.drop(columns=columns_to_remove, errors='ignore')

    X = df.drop(columns=["smiles"])

    y_pred = app.state.model.predict(X)

    return dict(melting_point=float(y_pred))
