# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 20:40:11 2019

@author: US16120
"""

import numpy as np
import pandas as pd

import pyupset as pyu
from pickle import load
from upsetplot import from_memberships
from upsetplot import plot

example = from_memberships(
        [[],                             # 0
         ['Mansouri'],                   # 445
         ['Cheng'],                      # 388
         ['Cheng', 'Mansouri'],          # 578
         ['OPERA'],                      # 850
         ['OPERA', 'Mansouri'],          # 114
         ['OPERA', 'Cheng'],             # 52
         ['OPERA', 'Cheng', 'Mansouri'], # 563
         ],
         data=[0, 445, 388, 578, 850, 114, 52, 563]
         )
example 

plot(example)
