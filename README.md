# Melting-point-predictor

⏳ WORK IN PROGRESS ⏳

The project aims to develop an advanced solution for estimating the melting point of organic compounds, a crucial characteristic in the field of organic chemistry. Accurate melting point determination is essential for compound characterization, with significant implications in various sectors of scientific research and industry, including the pharmaceutical and materials industries.

The estimator inputs a SMILES notation and outputs an estimated melting point.

![Image_presentation](https://github.com/Guillaume2126/Melting-point-predictor/assets/111251905/fc922485-83ad-40c0-9c30-2098057bb488)


To have more information about how this model was construct, please read the text above

## 1️⃣ Context

The melting point of a molecule is a fundamental physical property with implications to its chemical stability and application. Traditional experimental methods for determining melting points are often time-consuming and costly. Machine learning models offer a fast, cost-effective alternative for predicting physical properties from molecular structures.

![Melting_point](https://github.com/Guillaume2126/Melting-point-predictor/assets/111251905/ea39b184-c869-44c6-9b15-60a306890f73)


Recent advances in machine learning and deep learning technologies have led to the development of new models for predicting the properties of molecules. One challenge, which is still of great importance today, is the prediction of the melting temperature of chemical compounds.

Many teams have tried to overcome this challenge, and the papers below are just a few examples:
1. Melting temperature prediction via first principles and deep learning, Qi-Jun Hong, Computational Materials Science, 2022, 214, (111684).
2. Building Machine Learning Small Molecule Melting Points and Solubility Models Using CCDC Melting Points Dataset Xiangwei Zhu, Valery R. Polyakov, Krishna Bajjuri, Huiyong Hu, Andreas Maderna, Clare A. Tovee, and Suzanna C. Ward, Journal of Chemical Information and Modeling, 2023, 63 (10), 2948-2959

In our case, the project was carried out independently, and we are not affiliated with any laboratory or organization.


## 2️⃣ Data collection

A dataset of molecular structures and associated melting points was gathered ([Jean-Claude Bradley Melting Point Dataset](https://figshare.com/articles/dataset/Jean_Claude_Bradley_Open_Melting_Point_Datset/1031637)).

## 3️⃣ Data visualization and data cleaning

We utilized Pearson correlation coefficients to visualize relationships between features, excluding those with coefficients greater than 0.90 to prevent multicollinearity from skewing the model's results.

A quick step of data cleaning was performed to remove duplicates or NaN values.

## 4️⃣ Feature engineering and scaling

Various stages of features engineering were carried out, mainly using the RDKit library where possible. As a result, the following columns were added to the initial dataframe:

***-> Type and number of bonds in each molecules (example: C-C or C=0)***

***-> Functional groups***

***-> Molecular weight***

***-> Number of different elements***

***-> Number of aromatic rings***

***-> Main chain lenght***

***-> Number of double bonds***

***-> XlogP***

Numerical features were scaled using an ensemble of scalers to normalize data and reduce the influence of outlier values. The scaled dataset was partitioned into a 90% training set and a 10% test set, balancing the need for model training and validation accuracy.

## 5️⃣ Find a model by using machine learning and deep learning

The LightGBM model was selected for its advantages in handling large datasets with high efficiency. A pipeline incorporating robust scaling, standard scaling and minmax normalizer was employed, ensuring that scaling would fit on the training data and could be applied consistently to the test data. The model's hyperparameters were fine-tuned using GridSearchCV, optimizing for minimal RMSE.

Post hyperparameter optimization, the LightGBM model achieved RMSE scores of 40.9°C and 39.5°C on two different test runs, indicating consistent predictive capabilities. Our linear regression baseline started at 45°C. The model's performance was further visualized in a scatter plot of actual versus predicted melting points, with deviations from the y=x line analyzed to understand prediction discrepancies. The best neural network model scored at RMSE = 43. We then decided to go on with classical machine learning.

The errors were further examined by creating bins based on actual melting points and calculating statistics within these bins, revealing patterns in the model's performance across different temperature ranges.

## 6️⃣ API

--- In progress ---

## 7️⃣ User Interface

Use the [User Interface](https://github.com/lccopy/Melting-point-predictor-UI) of the project.


<img width="700" alt="image" src="https://github.com/lccopy/Melting-point-predictor-UI/assets/111251905/003a0dfd-6932-4f07-bae8-420ee3c262a1">


The User Interface was created using Flask.

Using this estimator is relatively straightforward. First, log on to the site and enter the SMILES rating of any compound, then click on "Estimate". You can use an existing compound, but it's also possible to test the application with a compound that hasn't yet been synthesized or listed.

Note that the melting point is an estimate, not an exact prediction.

The authors would like to point out that this application has been created independently and should not be used for commercial purposes.

## 8️⃣ Conclusion and future improvements

The study demonstrates the potential of machine learning in predicting molecular melting points, with LightGBM providing a solid framework for development. While the current model performs adequately, further research into feature enhancement and model refinement is recommended to address areas of high predictive uncertainty. Our approach sets the stage for more sophisticated predictive tools in computational chemistry, aiding in the rapid assessment of molecular properties.

To improve these results, other avenues of improvement are possible, such as adding new columns through feature engineering, increasing the number of data items, using new models and modifying hyperparameters. In our case, the results obtained (RMSE of 39.5) are satisfactory, given that RMSE in the literature are currently around 35-40.

## Collaborators

PhD. Guillaume Bretel (Doctor of Chemistry, Data Scientist)

Luca Morel (Data Scientist, AI developer)
