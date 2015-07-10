import cantera as ct
import csv
import numpy as np
import sys
from PyQt4 import QtGui, QtCore
# define gui app and widget
app = QtGui.QApplication(sys.argv)
# Widget class with changed icon


class NSortableTableWidgetItem(QtGui.QTableWidgetItem):
    # Implement less than (<) for numeric table widget items.
    def __init__(self, text):
        QtGui.QTableWidgetItem.__init__(self, text)

    def __lt__(self, y):
        return float(self.text()) < float(y.text())


class MainForm(QtGui.QWidget):

    def __init__(self, parent=None):
        ## CALCULATIONS
        # Instantiate
        self.gas = ct.Solution('gri30.cti')
        # Properties T,P,X
        nitrogen_over_oxygen_mole_ratio = 78.85/21.15
        initial_nitrogen = 2*nitrogen_over_oxygen_mole_ratio
        self.gas.TPY = 300, 100000, 'CH4:1, O2:2, N2:' + \
                                    str(initial_nitrogen)
        # Equilibrate @ constant TP
        print '\n\n Equilibrating... \n\n'
        self.gas.equilibrate('TP')
        print '\n\n ...Done... \n\n'

        QtGui.QWidget.__init__(self, parent)

        self.setGeometry(300, 100, 550, 500)
        self.setFixedWidth(400)
        self.setMinimumHeight(500)
        self.setWindowTitle('Icon')
        self.setWindowIcon(QtGui.QIcon('./utils/ch_pot.png'))

        self.table_properties = QtGui.QTableWidget(4, 3, self)
        self.table_state_funcs = QtGui.QTableWidget(6, 4, self)
        self.table_composition = QtGui.QTableWidget(1, 4, self)

        self.table_properties.connect(self.table_properties,
                                      QtCore.SIGNAL('cellChanged(int,int)'),
                                      self.new_properties)

        self.clear_composition = QtGui.QPushButton('Clear Composition')
        self.clear_composition.connect(self.clear_composition,
                                       QtCore.SIGNAL('clicked()'),
                                       self.clear_composition_proc)

        # Widths
        for table in [self.table_state_funcs, self.table_composition, self.table_properties]:
            table.horizontalHeader().setStretchLastSection(True)
            table.verticalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
            if table != self.table_composition:
                table.verticalHeader().setStretchLastSection(True)
            else:
                table.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)

        v_box = QtGui.QVBoxLayout(self)
        v_box.addStrut(5)
        v_box.addWidget(self.table_properties)
        v_box.addWidget(self.table_state_funcs)
        v_box.addWidget(self.clear_composition)
        v_box.addWidget(self.table_composition)

        self.setLayout(v_box)

    def print_eq(self):
        ## PRINTING
        # Block signals to avoid recursion
        self.table_properties.blockSignals(True)
        self.table_state_funcs.blockSignals(True)
        self.table_composition.blockSignals(True)
        # Print to console
        print self.gas()
        # Save csv file in ./DATA/
        print '\n\n Saving csv... \n\n'
        csv_file = './DATA/equilibrium_simple_ex_01_result_ignorethis.csv'
        # Quickest output.
        with open(csv_file, 'w') as outfile:
            writer = csv.writer(outfile, lineterminator='\n')
            writer.writerow(['T (K)']+self.gas.species_names)
            writer.writerow([self.gas.T]+list(self.gas.X))
        # Redundant but more complete output.
        with open(csv_file, 'w') as outfile:
            sort_indexes = np.flipud(np.argsort(self.gas.Y))
            sorted_y = self.gas.Y[sort_indexes]
            sorted_x = self.gas.X[sort_indexes]
            sorted_names = np.array(self.gas.species_names)[sort_indexes]
            properties_list = [['T', 'K'],
                               ['P', 'Pa'],
                               ['density', 'kg/m^3'],
                               ['mean_molecular_weight', 'amu']]
            state_funcs_list = \
                [['enthalpy', 'J'], ['int_energy', 'J'],
                 ['entropy', 'J/K'], ['gibbs', 'J'],
                 ['cp', 'J/K'], ['cv', 'J/K']]
            writer = csv.writer(outfile, lineterminator='\n')
            for i in range(len(properties_list)):
                writer.writerow(
                    [properties_list[i][0],
                     getattr(self.gas, properties_list[i][0]),
                     properties_list[i][1]])
                self.table_properties.setItem(i, 0, QtGui.QTableWidgetItem(
                    properties_list[i][0]))
                self.table_properties.setItem(i, 1, QtGui.QTableWidgetItem(
                    str(getattr(self.gas, properties_list[i][0]))))
                self.table_properties.setItem(i, 2, QtGui.QTableWidgetItem(
                    properties_list[i][1]))
            writer.writerow([''])
            writer.writerow(['', '1 kg', '1 kmol', ''])
            for i in range(len(state_funcs_list)):
                writer.writerow(
                    [state_funcs_list[i][0],
                     getattr(self.gas, state_funcs_list[i][0]+'_mass'),
                     getattr(self.gas, state_funcs_list[i][0]+'_mole'),
                     state_funcs_list[i][1]])
                self.table_state_funcs.setItem(i, 0, QtGui.QTableWidgetItem(
                    state_funcs_list[i][0]))
                self.table_state_funcs.setItem(i, 1, QtGui.QTableWidgetItem(
                    str(getattr(self.gas, state_funcs_list[i][0]+'_mass'))))
                self.table_state_funcs.setItem(i, 2, QtGui.QTableWidgetItem(
                    str(getattr(self.gas, state_funcs_list[i][0]+'_mole'))))
                self.table_state_funcs.setItem(i, 3, QtGui.QTableWidgetItem(
                    state_funcs_list[i][1]))
            writer.writerow([''])
            writer.writerow(['', 'X', 'Y', 'Chem. Pot. / RT'])

            self.table_composition.setRowCount(sort_indexes.size)
            for i in range(sort_indexes.size):
                writer.writerow([sorted_names[i]] + [sorted_x[i]] + [sorted_y[i]] +
                                [self.gas.chemical_potentials[i]/self.gas.T/ct.gas_constant])
                self.table_composition.setItem(i, 0, QtGui.QTableWidgetItem(sorted_names[i]))
                self.table_composition.setItem(i, 1, NSortableTableWidgetItem(str(sorted_x[i])))
                self.table_composition.setItem(i, 2, NSortableTableWidgetItem(str(sorted_y[i])))
                self.table_composition.setItem(i, 3, NSortableTableWidgetItem(
                    str(self.gas.chemical_potentials[i]/self.gas.T/ct.gas_constant)))
                for j in range(2+1):
                    self.table_composition.item(i, j).setFlags(QtCore.Qt.ItemIsEnabled)

            for i in range(self.table_properties.rowCount()):
                for j in range(self.table_properties.columnCount()):
                    if j != 1 or i > 1:
                        self.table_properties.item(i, j).setFlags(QtCore.Qt.ItemIsEnabled)
            for i in range(self.table_state_funcs.rowCount()):
                for j in range(self.table_state_funcs.columnCount()):
                    self.table_state_funcs.item(i, j).setFlags(QtCore.Qt.ItemIsEnabled)

        self.table_state_funcs.setHorizontalHeaderLabels(['', '1 kg', '1 kgmol', ''])
        self.table_composition.setHorizontalHeaderLabels(['el. form.', 'X', 'Y', 'Chem. Pot. / RT'])
        self.table_properties.horizontalHeader().setVisible(False)
        self.table_properties.verticalHeader().setVisible(False)
        self.table_state_funcs.verticalHeader().setVisible(False)
        self.table_composition.verticalHeader().setVisible(False)
        self.table_composition.resizeColumnsToContents()
        self.table_state_funcs.resizeColumnsToContents()
        self.table_properties.resizeColumnsToContents()
        self.table_composition.setSortingEnabled(True)
        self.table_composition.sortByColumn(2)
        self.table_composition.horizontalHeader().setSortIndicator(2, QtCore.Qt.DescendingOrder)
        # Stop blocking signals
        self.table_properties.blockSignals(False)
        self.table_state_funcs.blockSignals(False)
        self.table_composition.blockSignals(False)

        print '\n\n ...Done... \n\n'

    def new_properties(self, row, column):
        selected_property = str(self.table_properties.item(row, 0).text())
        print 'new ' + selected_property + ' = ' + \
              str(self.table_properties.item(row, column).text())
        T = float(self.table_properties.item(0, 1).text())
        P = float(self.table_properties.item(1, 1).text())
        names_Y_dict = dict()

        for i in range(self.table_composition.rowCount()):
            names_Y_dict[str(self.table_composition.item(i, 0).text())] = \
                float(self.table_composition.item(i, 2).text())

        self.gas.TPY = T, P, names_Y_dict
        self.gas.equilibrate('TP')
        self.print_eq()

    def clear_composition_proc(self):
        self.table_composition.setSortingEnabled(False)
        for i in range(1, self.table_composition.rowCount()):
            self.table_composition.setItem(i, 2,
                                           QtGui.QTableWidgetItem(str(0.0)))
        self.table_composition.setSortingEnabled(True)
        self.new_properties(0, 1)
icon = MainForm()
icon.print_eq()
icon.show()
icon.setWindowTitle('CH4 / O2 Equilibrium composition')
# show the main widget, enter 'main loop'

# exit 'main loop'
sys.exit(app.exec_())