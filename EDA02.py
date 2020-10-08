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

train_df = PandasTools.LoadSDF("OPERA/TR_RBioDeg_1197.sdf")
test_df = PandasTools.LoadSDF("OPERA/TST_RBioDeg_411.sdf")
df = pd.concat([train_df[['CAS', 'Canonical_QSARr', 'Ready_Biodeg']],
                test_df[['CAS', 'Canonical_QSARr', 'Ready_Biodeg']]], ignore_index = True)
df['Ready_Biodeg'] = pd.to_numeric(df['Ready_Biodeg'])
def EndPtBeta (row):
   if row['Ready_Biodeg'] >= 0.5 :
      return 'RB'
   return 'NRB'
df['EndPt'] = df.apply (lambda row: EndPtBeta(row), axis=1)
df = df.rename(columns={'CAS': 'CASRN', 'Canonical_QSARr': 'SMILES'})
df = df[['CASRN', 'SMILES', 'EndPt']]

for i in range(len(df)):
    mol = Chem.MolFromSmiles(df.loc[i, 'SMILES'])
    mol = s.standardize(mol)
    mols = rdmolops.GetMolFrags(mol, asMols = True)
    mol = max(mols, default = mol, key=lambda m: m.GetNumAtoms())
    df.at[i, 'mol'] = mol
    
