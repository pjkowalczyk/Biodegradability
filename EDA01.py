"""
Cheminformatics | Biodegradability | EDA.py
2020-10-05
"""
import os
os.getcwd()

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
import pandas as pd
import numpy as np

from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops
from rdkit.Chem import inchi
from rdkit.Chem import MolStandardize
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import PandasTools
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 
import numpy as np
import math
from sklearn.ensemble import RandomForestRegressor

from molvs import Standardizer
s = Standardizer()
from molvs import fragment

df01 = pd.read_excel(open('JChemInfModel_52_655/Table S1.xls','rb'),
                       sheet_name='Training set and Test Set')
list(df01)
columns = ['Data Set', 'Appliation Domain (ID: In Domain, OD: Out of Domain)']
df01.drop(columns, inplace=True, axis=1)
df01 = df01.rename(columns={'Experimental Labels': 'EndPt'})

df02 = pd.read_excel(open('JChemInfModel_52_655/Table S1.xls','rb'),
                       sheet_name='817 compounds with BOD% Value')
list(df02)
def EndPt (row):
   if row['BOD%'] >= 60.0 :
      return 'RB'
   return 'NRB'
df02['EndPt'] = df02.apply (lambda row: EndPt(row), axis=1)
df02beta = df02[['CASRN', 'SMILES', 'EndPt']]

df03 = pd.read_excel(open('JCIM_53_867/MansouriData.xlsx','rb'),
                       sheet_name='Train+Test')
df03 = df03[['CAS-RN', 'Smiles', 'Class']]
df03 = df03.rename(columns={'CAS-RN': 'CASRN', 'Smiles': 'SMILES', 'Class': 'EndPt'})

df04 = pd.read_excel(open('JCIM_53_867/MansouriData.xlsx','rb'),
                       sheet_name='External validation')
df04 = df04[['CAS-RN', 'Smiles', 'class']]
df04 = df04.rename(columns={'CAS-RN': 'CASRN', 'Smiles': 'SMILES', 'class': 'EndPt'})

train_df = PandasTools.LoadSDF("OPERA/TR_RBioDeg_1197.sdf")
test_df = PandasTools.LoadSDF("OPERA/TST_RBioDeg_411.sdf")
df05 = pd.concat([train_df[['CAS', 'Canonical_QSARr', 'Ready_Biodeg']],
                 test_df[['CAS', 'Canonical_QSARr', 'Ready_Biodeg']]], ignore_index = True)
df05['Ready_Biodeg'] = pd.to_numeric(df05['Ready_Biodeg'])
def EndPtBeta (row):
   if row['Ready_Biodeg'] >= 0.5 :
      return 'RB'
   return 'NRB'
df05['EndPt'] = df05.apply (lambda row: EndPtBeta(row), axis=1)
df05 = df05.rename(columns={'CAS': 'CASRN', 'Canonical_QSARr': 'SMILES'})
df05beta = df05[['CASRN', 'SMILES', 'EndPt']]

# join data

df01['Source'] = 'Cheng'
df02beta['Source'] = 'Cheng'
df03['Source'] = 'Mansouri'
df04['Source'] = 'Mansouri'
df05beta['Source'] = 'OPERA'

frames = [df01, df02beta, df03, df04, df05beta]
df = pd.concat(frames, ignore_index = True)

for i in range(len(df)):
    mol = Chem.MolFromSmiles(df.loc[i, 'SMILES'])
    mol = s.standardize(mol)
    mols = rdmolops.GetMolFrags(mol, asMols = True)
    mol = max(mols, default = mol, key=lambda m: m.GetNumAtoms())
    df.at[i, 'mol'] = mol
    df.at[i, 'InChI'] = inchi.MolToInchi(mol)
    
# duplicate InChI recorsd
duplicateRowsDF = df[df.duplicated(['InChI'])]
DF = duplicateRowsDF.sort_values(by=['InChI'])
DF.to_csv('DF.csv')
