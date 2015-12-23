# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 20:48:42 2015

@author: Santiago Salas
"""
import os, sys, numpy as np, scipy as sp, csv
from PySide import QtGui, QtCore
from functools import partial
from sympy import solve, nsolve, symbols
from numpy import log10, matlib

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8


    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)


class Ui_GroupBox(object):
    def setupUi(self, GroupBox):
        GroupBox.setObjectName(_fromUtf8("GroupBox"))
        # Default size
        # GroupBox.resize(394, 357)
        GroupBox.resize(500, 357)
        self.verticalLayout_2 = QtGui.QVBoxLayout(GroupBox)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.open_button = QtGui.QPushButton(GroupBox)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(
            _fromUtf8("utils/glyphicons-145-folder-open.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(
            _fromUtf8("utils/glyphicons-415-disk-save.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(
            _fromUtf8("utils/glyphicons-41-stats.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.open_button.setIcon(icon)
        self.open_button.setObjectName(_fromUtf8("open_button"))
        self.horizontalLayout_2.addWidget(self.open_button)
        self.save_button = QtGui.QPushButton(GroupBox)
        self.save_button.setIcon(icon1)
        self.save_button.setObjectName(_fromUtf8("save_button"))
        self.horizontalLayout_2.addWidget(self.save_button)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.label = QtGui.QLabel(GroupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.label.setAlignment(QtCore.Qt.AlignRight)
        self.horizontalLayout_5.addWidget(self.label)
        self.spinBox = QtGui.QSpinBox(GroupBox)
        self.spinBox.setProperty("value", 0)
        self.spinBox.setObjectName(_fromUtf8("spinBox"))
        self.horizontalLayout_5.addWidget(self.spinBox)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.Components = QtGui.QTableWidget(GroupBox)
        self.Components.setMinimumSize(QtCore.QSize(0, 210))
        self.Components.setObjectName(_fromUtf8("Components"))
        item = QtGui.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.Components.horizontalHeader().setCascadingSectionResizes(False)
        self.Components.horizontalHeader().setDefaultSectionSize(100)
        self.Components.horizontalHeader().setMinimumSectionSize(27)
        self.Components.horizontalHeader().setSortIndicatorShown(True)
        self.Components.verticalHeader().setVisible(False)
        self.verticalLayout_2.addWidget(self.Components)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.label_2 = QtGui.QLabel(GroupBox)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_2.setAlignment(QtCore.Qt.AlignRight)
        self.horizontalLayout_6.addWidget(self.label_2)
        self.spinBox_2 = QtGui.QSpinBox(GroupBox)
        self.spinBox_2.setProperty("value", 0)
        self.spinBox_2.setObjectName(_fromUtf8("spinBox_2"))
        self.horizontalLayout_6.addWidget(self.spinBox_2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_6)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.tableReacs = QtGui.QTableWidget(GroupBox)
        self.tableReacs.setObjectName(_fromUtf8("tableReacs"))
        self.tableReacs.horizontalHeader().setVisible(True)
        self.tableReacs.verticalHeader().setVisible(False)
        self.horizontalLayout.addWidget(self.tableReacs)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.comboBox = QtGui.QComboBox(GroupBox)
        self.comboBox.setObjectName(_fromUtf8("comboBox"))
        self.verticalLayout.addWidget(self.comboBox)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.doubleSpinBox = QtGui.QDoubleSpinBox(GroupBox)
        self.doubleSpinBox.setObjectName(_fromUtf8("doubleSpinBox"))
        self.horizontalLayout_3.addWidget(self.doubleSpinBox)
        self.doubleSpinBox_2 = QtGui.QDoubleSpinBox(GroupBox)
        self.doubleSpinBox_2.setObjectName(_fromUtf8("doubleSpinBox_2"))
        self.horizontalLayout_3.addWidget(self.doubleSpinBox_2)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.comboBox_2 = QtGui.QComboBox(GroupBox)
        self.comboBox_2.setObjectName(_fromUtf8("comboBox_2"))
        self.verticalLayout.addWidget(self.comboBox_2)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.doubleSpinBox_4 = QtGui.QDoubleSpinBox(GroupBox)
        self.doubleSpinBox_4.setObjectName(_fromUtf8("doubleSpinBox_4"))
        self.horizontalLayout_4.addWidget(self.doubleSpinBox_4)
        self.doubleSpinBox_3 = QtGui.QDoubleSpinBox(GroupBox)
        self.doubleSpinBox_3.setObjectName(_fromUtf8("doubleSpinBox_3"))
        self.horizontalLayout_4.addWidget(self.doubleSpinBox_3)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.pushButton = QtGui.QPushButton(GroupBox)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.pushButton.setIcon(icon2)
        self.verticalLayout.addWidget(self.pushButton)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.open_button.clicked.connect(partial(open_file, self))
        self.save_button.clicked.connect(partial(save_file, self))
        self.pushButton.clicked.connect(partial(plot_intervals, self))
        self.retranslateUi(GroupBox)
        QtCore.QMetaObject.connectSlotsByName(GroupBox)

    def retranslateUi(self, GroupBox):
        GroupBox.setWindowTitle(_translate("GroupBox", "Simplistic EC.", None))
        # GroupBox.setTitle(_translate("GroupBox", "Simplistic Eq.", None))
        __sortingEnabled = self.Components.isSortingEnabled()
        self.open_button.setText(_translate("GroupBox", "Open", None))
        self.save_button.setText(_translate("GroupBox", "Save", None))
        self.Components.setSortingEnabled(__sortingEnabled)
        self.pushButton.setText(_translate("GroupBox", "Plot", None))
        self.label_2.setText(_translate("GroupBox", "Nr (Reac.)", None))
        self.label.setText(_translate("GroupBox", "n (Comp.)", None))


def open_file(form):
    (filename, _) = \
        QtGui.QFileDialog.getOpenFileName(None,
                                          caption='Open file',
                                          dir=os.path.join(sys.path[0], 'DATA'),
                                          filter='*.csv')
    if os.path.isfile(filename):
        load_csv(filename, form)


def load_csv(filename, form):
    with open(filename) as csv_file:
        n = 0
        Nr = 0
        header_comps = []
        header_reacs = []
        reader = csv.reader(csv_file, dialect='excel')
        for row in reader:
            if 'COMP' in row:
                header_comps = next(reader)[0:4]
                next(reader)
            elif 'REAC' in row:
                header_reacs = next(reader)
                next(reader)
            if len(header_reacs) == 0 and len(header_comps) > 0:
                n = n + 1
            elif len(header_reacs) > 0 and len(header_comps) > 0:
                Nr = Nr + 1
        csv_file.close()
        form.spinBox.setProperty("value", n)
        form.spinBox_2.setProperty("value", Nr)
    with open(filename) as csv_file:
        reader = csv.reader(csv_file, dialect='excel')
        comps = np.empty([n, 4], dtype='S50')
        reacs = np.empty([Nr, n + 2], dtype='S50')
        i = 0
        j = 0
        for row in reader:
            if reader.line_num > 2 and reader.line_num <= n + 2:
                comps[i] = np.array(row[0:4])
                comps[i] = map(lambda x: '0' if x == '' else x, comps[i])
                i = i + 1
            elif reader.line_num > n + 2 + 2:
                reacs[j] = np.array(row)
                reacs[j] = map(lambda x: '0' if x == '' else x, reacs[j])
                j = j + 1
        csv_file.close()

        form.Components.setRowCount(n)
        form.Components.setColumnCount(len(header_comps) + 3)
        form.Components.setHorizontalHeaderLabels(
            header_comps + ['Cieq, mol/L', '-log10(Ci0)', '-log10(Cieq)'])

        form.tableReacs.setRowCount(Nr)
        form.tableReacs.setColumnCount(len(header_reacs) + 1)
        form.tableReacs.setHorizontalHeaderLabels(
            header_reacs + ['Xi_j'])

        i = range(0, n)
        j = range(0, len(header_comps))

        for column in j:
            for row in i:
                newItem = QtGui.QTableWidgetItem(str(comps[row][column]))
                if column != 1:  # Comp. i <Str>
                    newItem = NSortableTableWidgetItem(newItem)
                    form.Components.setItem(row, column, newItem)
                else:
                    form.Components.setItem(row, column, newItem)
                if not column in range(1, 3 + 1):
                    newItem.setFlags(QtCore.Qt.ItemIsEnabled)

        i = range(0, Nr)
        j = range(0, len(header_reacs))

        for column in j:
            for row in i:
                form.tableReacs.setItem(row, column, NSortableTableWidgetItem(str(reacs[row][column])))

        # Widths and heights
        form.Components.setSortingEnabled(True)
        form.Components.sortByColumn(0, QtCore.Qt.AscendingOrder)
        form.Components.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        form.Components.verticalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)

        form.tableReacs.setSortingEnabled(True)
        form.tableReacs.sortByColumn(0, QtCore.Qt.AscendingOrder)
        form.tableReacs.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        form.tableReacs.verticalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)

        C0_i = np.matrix([row[3] for row in comps], dtype=float).T
        z_i = np.matrix([row[2] for row in comps], dtype=float).T
        nu_ij = np.matrix([row[2:2 + n] for row in reacs], dtype=int).T
        pKa_j = np.matrix([row[1] for row in reacs], dtype=float).T

        calc_Xieq(C0_i, z_i, nu_ij, pKa_j)


def save_file(form):
    pass


def plot_intervals(form):
    pass


def calc_Xieq(C0_i, z_i, nu_ij, pKa_j, Xieq_0=0, Ceq_0=0, tol=np.finfo(float).eps, max_it=1000):
    """
    :return: tuple with Ceq_i, Xieq_j, f_0
    :param C0_i: np.matrix (n X 1) - Conc(i, alimentación)
    :param z_i: np.matrix (n X 1) - Carga(i, alimentación)
    :param nu_ij: np.matrix (n X Nr) - Coefs. esteq. componente i en reacción j
    :param pKa_j: np.matrix (n X 1) - (-1)*log10("Cte." de equilibrio en reacción j) = -log10 Kc_j(T)
    :param Xieq_0: np.matrix (n X 1) - avance de reacción j - estimado inicial
    :param Ceq_0: np.matrix (n X 1) - Conc(i, equilibrio)
    """
    n = nu_ij.shape[0]
    Nr = nu_ij.shape[1]
    Kc_j = np.power(10, -pKa_j) / (997 / (2 * 1.00794 + 15.9994))
    Ceq_0 = np.matrix(np.ones([n, 1])) * tol
    Xieq_j = np.matrix(np.ones([Nr, 1])) * tol
    X0 = np.concatenate([Ceq_0, Xieq_j])
    J0 = J(Ceq_0, C0_i, nu_ij, Xieq_j, Kc_j)
    k = 1
    while k < max_it:
        print 'hi'
        k += 1
    # Ceq_i = np.matrix(symbols('Ce0:' + str(n))).transpose()
    # Ceq_0 = Ceq_i
    # Ceq_0 = C0_i + np.finfo(float).eps
    # Xieq_j = np.matrix(symbols('xi0:' + str(Nr))).transpose()
    # Xieq_0 = np.matrix(np.zeros([Nr, 1]))
    # # f_0 = f_gl_0(Ceq_i, C0_i, nu_ij, Xieq_j, Kc_j)
    # f_0 = np.empty([n + Nr, 1], dtype=object)
    # f_0[0:n] = -Ceq_i + C0_i + nu_ij * Xieq_j
    # f_0[n:n + Nr] = -Kc_j + np.prod(np.power(Ceq_i, nu_ij), 0).T
    # soln = nsolve(f_0, np.concatenate([Ceq_i, Xieq_j]),
    #               np.concatenate([Ceq_0, Xieq_0]))
    # Ceq_i = soln[0:n]
    # Xieq_j = soln[n:len(soln)]
    # a = []
    # a.sort()
    return 0


def f_gl_0(Ceq_i, C0_i, nu_ij, Xieq_j, Kc_j):
    n = nu_ij.shape[0]
    Nr = nu_ij.shape[1]
    f_0 = np.matrix(np.empty([n + Nr, 1], dtype=object))
    f_0[0:n] = -Ceq_i + C0_i + nu_ij * Xieq_j
    f_0[n:n + Nr] = -Kc_j + np.prod(np.power(Ceq_i, nu_ij), 0).T


def J(Ceq_i, C0_i, nu_ij, Xieq_j, Kc_j):
    n = nu_ij.shape[0]
    Nr = nu_ij.shape[1]
    Eins_durch_C = np.diag(np.power(Ceq_i, -1).A1, 0)
    Quotient = np.diag(np.prod(np.power(Ceq_i, nu_ij), 0).A1)
    J = np.matrix(np.zeros([n + Nr, n + Nr], dtype=float))
    J[0:n, 0:n] = -1 * np.eye(n).astype(float)
    J[0:n, n:n + Nr] = nu_ij
    J[n:n + Nr, 0:n] = Quotient * nu_ij.T * Eins_durch_C
    return J


class NSortableTableWidgetItem(QtGui.QTableWidgetItem):
    # Implement less than (<) for numeric table widget items.
    def __init__(self, text):
        QtGui.QTableWidgetItem.__init__(self, text)

    def __lt__(self, y):
        return float(self.text()) < float(y.text())


class MainForm(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.setWindowIcon(QtGui.QIcon(
            os.path.join(sys.path[0], *['utils', 'icon_batch.png'])))
        self.ui = Ui_GroupBox()
        self.ui.setupUi(self)


app = QtGui.QApplication.instance()  # checks if QApplication already exists
if not app:  # create QApplication if it doesnt exist
    app = QtGui.QApplication(sys.argv)

main_form = MainForm()
main_form.show()
load_csv('./DATA/COMPONENTS_REACTIONS_EX_001.csv', main_form.ui)

sys.exit(app.exec_())