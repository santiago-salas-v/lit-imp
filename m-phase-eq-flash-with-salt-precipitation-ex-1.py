__author__ = 'Santiago Salas'
__ref__ = 'Lucia, A., Henley, H., Thomas, E.,Multiphase Equilibrium' + \
          'Flash with Salt Precipitation in Systems with Multiple Salts, Chemical Engineering' + \
          'Research and Design (2014), http://dx.doi.org/10.1016/j.cherd.2014.04.034'
__ex__ = 'Example 1'

import pandas as pd

csv_file = open('./DATA/EX-1-parameters.csv')
dt = pd.read_csv(csv_file, quotechar='"')
csv_file.close()