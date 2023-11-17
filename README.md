# Melting-point-predictor

⏳ WORK IN PROGRESS ⏳ 

The aim of this project is to estimate the melting point of molecules
(explain more deeply (what is the input, output, how it works+
add a video about how it works).

To have more information about how this model was construct, please read the text above

## Data collection

## Data visualization and data cleaning

## Feature engineering and scaling

## Find a model using machine learning and deep learning

### Baseline


## API

## Create the website - Front end

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
