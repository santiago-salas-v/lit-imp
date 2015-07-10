import cantera as ct
import os
import numpy as np
import sys
import matplotlib
# Use PySide instead of Qt4
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
import matplotlib.pylab as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from PySide import QtGui, QtCore
from PySide.QtGui import QFontMetrics, QFont
# define gui app and widget
app = QtGui.QApplication(sys.argv)
# Widget class with changed icon

def reactor_solution(ct_mix, air, fuel_species,
                     stoich_ratio_phi, temp_init,
                     press_init, cond, nt, dt, t):
    # Arguments

    # Instantiate
    m = ct_mix.n_species
    stoich_oxygen = ct_mix.n_atoms(fuel_species, 'C') + \
                    1 / 4. * ct_mix.n_atoms(fuel_species, 'H')
    nitrogen_over_oxygen_mole_ratio = 78.85 / 21.15
    initial_nitrogen = stoich_oxygen * nitrogen_over_oxygen_mole_ratio

    ct_mix.TPX = temp_init, press_init, dict([
        (fuel_species, stoich_ratio_phi),
        ('O2', stoich_oxygen),
        ('N2', initial_nitrogen)])

    # ct.Reactor(ct_mix) does not converge easily (stiff solver?)
    r = ct.IdealGasReactor(ct_mix)

    # Default is 'UV'
    if cond == 'HP':
        # Dfine a flexible wall between the reaction mixture and
        # the environment, for pressure_inside = pressure_outside
        env = ct.Reservoir(air)
        wall = ct.Wall(left=r, right=env)
        wall.expansion_rate_coeff = 1.0e6  # set expansion parameter. dV/dt = KA(P_1 - P_2)
        wall.area = 1.0  # set wall area

    sim = ct.ReactorNet([r])
    time_progress = np.zeros(nt, dtype='float')
    temp_progress = np.zeros(nt, dtype='float')
    press_progress = np.zeros(nt, dtype='float')
    v_progress = np.zeros(nt, dtype='float')
    x_progress = np.zeros([nt, m], dtype='float')
    x_labels = ct_mix.species_names

    # loop nt number of dt seconds
    for n in range(nt):
        time_progress[n] = t
        temp_progress[n] = r.thermo.T
        press_progress[n] = r.thermo.P
        v_progress[n] = r.volume
        for i in range(m):
            x_progress[n, i] = r.thermo.X[i]
        t += dt
        sim.advance(t)

    outstruct = dict()
    outstruct['t'] =  time_progress
    outstruct['T'] =  temp_progress
    outstruct['P'] = press_progress
    outstruct['X'] = x_progress
    outstruct['x_labels'] = x_labels
    outstruct['V'] = v_progress
    return outstruct


class MainForm(QtGui.QWidget):
    def __init__(self, parent=None):
        # Superclass init call
        QtGui.QWidget.__init__(self, parent)
        #default_font = QFont('arial', 10)
        #cell_width = QFontMetrics(default_font).\
        #   width('+9.' + '9'*11 + 'e-999')
        outstruct = reactor_solution(
            ct_mix=ct.Solution('gri30.cti'),
            air=ct.Solution('air.cti'),
            fuel_species='CH4',  # select which is the fuel
            stoich_ratio_phi=1.0,  # stoichiometric ratio gas/O2
            temp_init=1000.0,  # K, initial temp
            press_init=101325.0,  # Pa , initial pressure
            cond='HP',  # cconstant H, P or U, V
            nt=100,  # number of time steps
            dt=(1.4-4e-1)/100, #*1/100 * 10e6 * 10e-6,  # us, time step size
            t=4e-1  # init. time
        )
        t = outstruct['t']
        temp = outstruct['T']
        p = outstruct['P']
        x = outstruct['X']
        x_labels = outstruct['x_labels']
        comp_labels_for_plots = \
            dict([(y, '$X_{' + y + '}$') for y in x_labels])
        comp_labels_for_table = \
            ['X[' + y + ']' for y in x_labels]
        v = outstruct['V']
        nColumns = len(p)
        nRows = x.shape[1] + 3
        self.table = QtGui.QTableWidget()
        self.table.setHorizontalHeaderLabels([])
        self.table.setRowCount(nRows)
        self.table.setColumnCount(nColumns)
        self.table.horizontalHeader().setVisible(False)
        self.table.setVerticalHeaderLabels(
            ['t, us', 'T, K', 'P, Pa'] + comp_labels_for_table)
        cell_width = self.table.fontMetrics(). \
            width('+9.' + '9' * 11 + 'e-999')
        del y

        for i in range(nColumns):
            self.table.setColumnWidth(i, cell_width)
            self.table.setItem(0, i, QtGui.QTableWidgetItem(str(t[i])))
            self.table.setItem(1, i, QtGui.QTableWidgetItem(str(temp[i])))
            self.table.setItem(2, i, QtGui.QTableWidgetItem(str(p[i])))
            for j in range(nRows - 3):
                self.table.setItem(3 + j, i,
                                   QtGui.QTableWidgetItem(
                                       str(x[i, j])))

        self.fig = plt.Figure((5.0, 4.0), dpi=100)
        self.canvas = FigureCanvas(self.fig)

        #super(self.canvas).__init__(FigureCanvas)
        #self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        #with plt.style.context('dark_background'):
        self.axes = self.fig.add_subplot(2, 2, 1)
        self.axes.plot(t, temp)
        self.axes.set_xlabel('t')
        self.axes.legend(['$T, K$']).draggable(True)
        self.axes = self.fig.add_subplot(2, 2, 2)
        self.axes.plot(t, v)
        self.axes.set_xlabel('t')
        self.axes.legend(['$V, m^3$']).draggable(True)
        self.axes = self.fig.add_subplot(2, 2, 3)
        self.axes.plot(t, x[:, x_labels.index('CH4')])
        self.axes.set_xlabel('t')
        self.axes.legend([comp_labels_for_plots['CH4']]).draggable(True)
        self.axes = self.fig.add_subplot(2, 2, 4)
        self.axes.plot(t, x[:, x_labels.index('O2')])
        self.axes.set_xlabel('t')
        self.axes.legend(['$X_{O2}$']).draggable(True)
        self.axes.hold(False)

        vertical_layout = QtGui.QVBoxLayout(self)
        vertical_layout.addWidget(self.table)
        vertical_layout.addWidget(self.canvas)
        self.setLayout(vertical_layout)


icon = MainForm()
icon.show()
icon.setWindowTitle('Costant P simple BR')
# show the main widget, enter 'main loop'

# exit 'main loop'
sys.exit(app.exec_())