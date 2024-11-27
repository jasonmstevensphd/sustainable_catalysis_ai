import os

import pandas as pd
import numpy as np
import plotly.express as px

import rdkit.Chem as Chem
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem import rdMolDescriptors
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.Draw import IPythonConsole

from mordred import Calculator, descriptors

import pickle

# Mordred Descriptor Calculation Function
 
def mordred(smiles_list, name='', dropna=False):
    """
    Compute all mordred descriptors for a list of smiles strings.
    """
    
    smiles_list = list(smiles_list)
    
    # Initialize descriptor calculator with all descriptors
 
    calc = Calculator(descriptors)
    
    output = []
    for entry in smiles_list:
        try:
            data_i = calc(Chem.MolFromSmiles(entry)).fill_missing()
        except:
            data_i = np.full(len(calc.descriptors),np.NaN)
            
        output.append(list(data_i))
        
    descriptor_names = list(calc.descriptors)
    columns = []
    for entry in descriptor_names:
        columns.append(name + str(entry))
        
    df = pd.DataFrame(data=output, columns=columns)
    df.insert(0, name + 'SMILES', smiles_list)
    
    if dropna == True:
        df = df.dropna(axis=1)
    
    return df

def Pred_RXN_NiBor(smiles_list):
    Pred_NiBor_MorVal = mordred(smiles_list, dropna = True)
    Pred_NiBor_MorVal.columns = Pred_NiBor_MorVal.columns.astype(str)
    tag = 'Electrophile_'
    Pred_NiBor_MorVal.columns = [tag + x for x in Pred_NiBor_MorVal.columns]
    Pred_NiBor_MorVal = pd.concat([Pred_NiBor_MorVal]*48, ignore_index=True)
    Pred_NiBor_MorVal_retain = Pred_NiBor_MorVal
    Pred_NiBor_MorVal.drop(['Electrophile_SMILES'], axis = 1, inplace = True)
    Pred_NiBor_PredSpace = pd.read_csv('Prediction/Prediction_RXN_NiBorylationPredictionSet.csv')
    Pred_NiBor_PredDF = Pred_NiBor_PredSpace[['Ligand_inchi', 'MeOH', 'Ligand']]
    Pred_Solvent = Pred_NiBor_PredDF.loc[:,'MeOH']
    Pred_NiBor_PredDF['Solvent'] = Pred_Solvent
    Pred_NiBor_PredDF.loc[Pred_NiBor_PredDF.Solvent == 1, 'Solvent'] = 'Methanol'
    Pred_NiBor_PredDF.loc[Pred_NiBor_PredDF.Solvent == 0, 'Solvent'] = 'Ethanol'
    Pred_NiBor_PredSpace = Pred_NiBor_PredSpace[Pred_NiBor_PredSpace.columns.drop(list(Pred_NiBor_PredSpace.filter(regex = 'Electrophile')))]
    Pred_NiBor_PredSpace = Pred_NiBor_PredSpace.drop(['Ligand_inchi', 'Ligand'], axis = 1)
    Pred_NiBor_PredSpace[[]] = Pred_NiBor_PredSpace[[]].apply(pd.to_numeric)
    ML = Pred_NiBor_MorVal.join(Pred_NiBor_PredSpace)
    Colnames = pd.read_csv('Prediction/Prediction_RXN_NiBorylationColnames.csv')
    Colnames = Colnames.drop(['Electrophile_inchi', 'Ligand_inchi'], axis = 1)
    Colnames = list(Colnames.columns.values)
    ML = ML[Colnames]
    ML[[]] = ML[[]].apply(pd.to_numeric)
    with open("Prediction/Prediction_RXN_NiBorylation_GBM.sav", 'rb') as file:  
        model = pickle.load(file)
    Pred_NiBor_PredDF['GBM_Prediction'] = model.predict(ML)
    return Pred_NiBor_PredDF
