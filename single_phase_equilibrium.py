__author__ = 'Santiago Salas'

import pandas as pd
import sys, os
import csv
import numpy as np

# Import data
csv_file = open('./DATA/EX-1-parameters.csv')
dt = pd.read_csv(csv_file, quotechar='"').fillna(0)
csv_file.close()

# Initialize variables
NASA_coeffs_names = map(lambda x: 'a'+str(x), range(1,8))
NASA_coeffs_names.append('b1')
NASA_coeffs_names.append('b2')

for i in NASA_coeffs_names:
    dt[i] = 0

# Extract data
components_i = dt.filter(regex = 'Comp')
z_charges_i = dt.filter(regex = 'z.*i')
r_coeffs_ij = dt.filter(regex = 'nu_i[0-9]')
n_components = len(components_i)
n_reactions = len(r_coeffs_ij.columns)

r_coeffs_ij.sort_index(axis='columns', ascending=True, inplace=True)

# Calculate Cp, Delta_H_r and K

# Delta_Cp_T0 =
# Delta_H_T0 =
#

# Calculate equilibrium at T, P

def nasa_coeffs(component_name):
    coeffs = dict(map(lambda x: (x, 0.0), NASA_coeffs_names))
    return coeffs

def convert_thermo_inp_to_csv(file_path_name_ext):
    file_path = os.path.split(file_path_name_ext)[0]
    file_name, file_ext = os.path.splitext(os.path.split(file_path_name_ext)[1])
    from pyparsing import Word, Literal, alphas, Optional, OneOrMore
    caps       = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    lowers     = caps.lower()
    digits     = "0123456789"
    signs      = '+-'
    element    = Word( caps, lowers )
    elementRef = element + Optional( Word( digits ) )
    species    = OneOrMore( elementRef ) + Optional( OneOrMore( Word( signs ) ) )

    thermo_inp_file = open(file_path_name_ext)
    output_file = os.path.join(file_path, file_name + '.csv')
    with open(output_file, 'w') as outfile:
        writer = csv.writer(outfile, lineterminator = '\n')
        writer.writerow(['i', 'comp_i', 'z_i', ])
    thermo_inp_file.close()