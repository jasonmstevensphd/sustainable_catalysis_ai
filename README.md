# sustainable_catalysis_ai
This repository contains all accompanying code to support the Organic Process Research and Development publication titled "Accelerating the Development of Sustainable Catalytic Processes through Data Science"

# File descriptions
Borylation_EDBO.r - Loads the final EDBO predictions and visualizes the predictions  
2024-08-02_CatalysisML.xlsx - Dataset of ML predicted yields for substrates 1a, 1b, 1c  
CatalysisMLPaper.r - Contains all tables, data, and plots directly reported in the paper  
Borylation_Updated (6).csv - Contains all values from EDBO+ study  
ML_Mordred.ipynb - Notebook for training the gradient boost model  
Mordred_Aligned.csv - Dataset for machine learning training  
NiBorylationData.xlsx - Contains data from high-throughput screening  
pred_Borylation_Updated (1).csv Contains EDBO+ predicted values for yield and cost after Round 2  
pred_Borylation_Updated (2).csv Contains EDBO+ predicted values for yield and cost after Round 3  
pred_Borylation_Updated.csv - Contains EDBO+ predicted values for yield and cost after Round 1  
Prediction_RXN_NiBorylation.py - python script for executing predictions  
predictions.ipynb - Notebook that executes predictions and validates the model published in the paper  
Prediction_RXN_NiBorylationColnames.csv - Dependent file for Prediction_RXN_NiBorylation.py & predictions.ipynb  
Prediction_RXN_NiBorylationPredictionSet.csv - Dependent file for Prediction_RXN_NiBorylation.py & predictions.ipynb  

# Environment considerations
All work was conducted on a a virtual machine with Ubuntu 18 and Python 3.6 and R 3.5
RDKit version 2019.03.4
scikit-learn version 0.23.2

# Model deployment
Train the model by executing the ML_Mordred.ipynb notebook. It is recommended to use a minimum of 12 cores and 60 GB RAM. Training time with this configuration is approximately 150 minutes. The code will produce a .sav file that will be needed to query the model.

# Model query
To obtain predictions for the model load the predictions.ipynb file and enter the smiles string of your aryl halide in the second code chunk.  
Example; df = Pred_RXN_NiBor(['CCOC(=O)C(C)(C)OC1=CC=C(C=C1)Cl'])  

The output will be a dataframe containing all 48 possible conditions with predicted values
