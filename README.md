# Melting-point-predictor

⏳ WORK IN PROGRESS ⏳

Name of collaborators here ?

The project aims to develop an advanced solution for estimating the melting point of organic compounds, a crucial characteristic in the field of organic chemistry. Accurate melting point determination is essential for compound characterization, with significant implications in various sectors of scientific research and industry, including the pharmaceutical and materials industries.

The estimator inputs a SMILES notation and outputs an estimated melting point.

![Image](Image_presentation.png)

To have more information about how this model was construct, please read the text above

## 1️⃣ Context

Interet du projet, bibliographie, peut-être une petite image ?



## 2️⃣ Data collection

## 3️⃣ Data visualization and data cleaning

Est-ce que j'ai des beaux graphiques a mettre ? si non, alors enlever "data vizualisation" du titre

## 4️⃣ Feature engineering and scaling

## 5️⃣ Find a model using machine learning and deep learning

### Baseline

Mettre belle image de la baseline

## 6️⃣ API

## 7️⃣ Create the website - Front end
https://github.com/lccopy/Melting-point-predictor-UI
___

More to come...

RMSE was used as metric.

List of models used:
Model 1: Linear Regression -> High RMSE (Around 8800000)

Model 2: RandomForest -> At the beginning: 46.2
                         After GridSearch: 45.9
        Hyperparameters:
    {'n_estimators': 300, 'min_samples_split': 5, 'min_samples_leaf': 1, 'max_depth': 30}

Model 3: MLPRegressor -> At the beginning: 45.5
                         After GridSearch:

## 8️⃣ Conclusion and future improvements
