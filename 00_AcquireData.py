import platform
print(platform.python_version())

import warnings
warnings.filterwarnings(action = 'once')

import csv
import pandas as pd
import numpy as np
import os

current_path = os.getcwd()
print(current_path)

# JCIM 2012, 52(3), 655-669

df01 = pd.read_excel('data\\JCIM_52_655.xls')

df01.shape

# JCIM 2013, 53(4), 867-878

df02 = pd.read_excel('data\\JCIM_53_867.xlsx')

df02.shape
