from fastapi import FastAPI
import pandas as pd
from Preprocessing.Create_dataframe import create_dataframe
from joblib import load

app = FastAPI()

model_path_joblib = "Model/best_full_data_model.joblib"
app.state.model = load(model_path_joblib)

@app.get("/predict")
def predict(SMILES):
    """
    Make a single melting point prediction.
    Assumes `SMILES` is provided as a string by the user
    """

    df = create_dataframe(SMILES)

    X = df.drop(columns=["smiles"])

    model = app.state.model
    assert model is not None
    print(type(X))
    print(X)
    y_pred = model.predict(X)

    return dict(melting_point=float(y_pred))
