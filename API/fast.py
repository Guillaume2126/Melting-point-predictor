from fastapi import FastAPI
import pandas as pd
from Preprocessing.Create_dataframe import create_dataframe

app = FastAPI()

model = #Charger le modele

@app.get("/predict")
def predict(SMILES):
    """
    Make a single melting point prediction.
    Assumes `SMILES` is provided as a string by the user
    """

    #Finish the preprocessing function before
    The code below is an example
    """X_pred = pd.DataFrame(locals(), index=[0])

    # Convert to US/Eastern TZ-aware!
    X_pred['pickup_datetime'] = pd.Timestamp(pickup_datetime, tz='US/Eastern')

    model = app.state.model
    assert model is not None

    X_processed = preprocess_features(X_pred)
    y_pred = model.predict(X_processed)"""

    # ⚠️ fastapi only accepts simple Python data types as a return value
    # among them dict, list, str, int, float, bool
    # in order to be able to convert the api response to JSON
    return dict(fare_amount=float(y_pred))
    # $CHA_END
