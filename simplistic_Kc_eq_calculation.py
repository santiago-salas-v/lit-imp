# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 20:48:42 2015

@author: Santiago Salas
"""
import os, sys, logging, pandas as pd, numpy as np, scipy as sp, csv
from PySide import QtGui, QtCore
from functools import partial
from mat_Zerlegungen import gausselimination
from datetime import datetime

# from sympy import solve, nsolve, symbols
# from numpy import log10, matlib

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
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
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
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(
            _fromUtf8("utils/glyphicons-82-refresh.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(
            _fromUtf8("utils/glyphicons-196-circle-info.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(
            _fromUtf8("utils/glyphicons-88-log-book.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.open_button.setIcon(icon)
        self.open_button.setObjectName(_fromUtf8("open_button"))
        self.horizontalLayout_2.addWidget(self.open_button)
        self.save_button = QtGui.QPushButton(GroupBox)
        self.save_button.setIcon(icon1)
        self.save_button.setObjectName(_fromUtf8("save_button"))
        self.info_button = QtGui.QPushButton(GroupBox)
        self.info_button.setIcon(icon4)
        self.info_button.setObjectName(_fromUtf8("info_button"))
        self.log_button = QtGui.QPushButton(GroupBox)
        self.log_button.setIcon(icon5)
        self.log_button.setObjectName(_fromUtf8("log_button"))
        self.horizontalLayout_2.addWidget(self.save_button)
        self.horizontalLayout_2.addWidget(self.log_button)
        self.horizontalLayout_2.addWidget(self.info_button)
        self.equilibrate_button = QtGui.QPushButton(GroupBox)
        self.equilibrate_button.setIcon(icon3)
        self.horizontalLayout_3.addWidget(self.equilibrate_button)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        # TODO: Add tol and max_it spinboxes functionality
        # TODO: Add log button
        self.label_3 = QtGui.QLabel(GroupBox)
        self.label_3.setObjectName(_fromUtf8("max_it_label"))
        self.label_3.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label_3)
        self.spinBox_3 = QtGui.QSpinBox(GroupBox)
        self.spinBox_3.setMaximum(2000)
        self.spinBox_3.setMinimum(2)
        self.spinBox_3.setProperty("value", 1000)
        self.spinBox_3.setObjectName(_fromUtf8("max_it_spinbox"))
        self.horizontalLayout_5.addWidget(self.spinBox_3)
        self.label_4 = QtGui.QLabel(GroupBox)
        self.label_4.setObjectName(_fromUtf8("tol_label"))
        self.label_4.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label_4)
        self.doubleSpinBox_5 = QtGui.QDoubleSpinBox(GroupBox)
        self.doubleSpinBox_5.setDecimals(int(-np.log10(np.finfo(float).eps) + 1))
        self.doubleSpinBox_5.setMaximum(float(1))
        self.doubleSpinBox_5.setMinimum(np.finfo(float).eps * 2)
        self.doubleSpinBox_5.setSingleStep(0.1 / 100.0 * (1.0 - np.finfo(float).eps))
        self.doubleSpinBox_5.setProperty("value", float(1.0e-08))
        self.doubleSpinBox_5.setObjectName(_fromUtf8("tol_spinbox"))
        self.horizontalLayout_5.addWidget(self.doubleSpinBox_5)
        self.label = QtGui.QLabel(GroupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label)
        self.spinBox = QtGui.QSpinBox(GroupBox)
        self.spinBox.setProperty("value", 0)
        self.spinBox.setObjectName(_fromUtf8("spinBox"))
        self.horizontalLayout_5.addWidget(self.spinBox)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.tableComps = QtGui.QTableWidget(GroupBox)
        self.tableComps.setMinimumSize(QtCore.QSize(0, 210))
        self.tableComps.setObjectName(_fromUtf8("tableComps"))
        item = QtGui.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.tableComps.horizontalHeader().setCascadingSectionResizes(False)
        self.tableComps.horizontalHeader().setDefaultSectionSize(100)
        self.tableComps.horizontalHeader().setMinimumSectionSize(27)
        self.tableComps.horizontalHeader().setSortIndicatorShown(True)
        self.tableComps.verticalHeader().setVisible(False)
        self.verticalLayout_2.addWidget(self.tableComps)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.label_5 = QtGui.QLabel(GroupBox)
        self.label_5.setObjectName(_fromUtf8("solvent_label"))
        self.label_5.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_6.addWidget(self.label_5)
        self.comboBox_3 = QtGui.QComboBox(GroupBox)
        self.horizontalLayout_6.addWidget(self.comboBox_3)
        self.label_6 = QtGui.QLabel(GroupBox)
        self.label_6.setObjectName(_fromUtf8("C_solvent_Tref"))
        self.label_6.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_6.addWidget(self.label_6)
        self.doubleSpinBox_6 = QtGui.QDoubleSpinBox(GroupBox)
        self.doubleSpinBox_6.setDecimals(int(-np.log10(np.finfo(float).eps) + 1))
        self.doubleSpinBox_6.setMaximum(float(1000))
        self.doubleSpinBox_6.setMinimum(np.finfo(float).eps * 1.1)
        self.doubleSpinBox_6.setSingleStep(1.0e-2)
        self.doubleSpinBox_6.setObjectName(_fromUtf8("C_solvent_Tref_doublespinbox"))
        self.horizontalLayout_6.addWidget(self.doubleSpinBox_6)
        self.label_2 = QtGui.QLabel(GroupBox)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
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
        self.label_7 = QtGui.QLabel(GroupBox)
        self.label_7.setObjectName(_fromUtf8("horizontalAxisLabel"))
        self.label_7.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignCenter)
        self.label_8 = QtGui.QLabel(GroupBox)
        self.label_8.setObjectName(_fromUtf8("verticalAxisLabel"))
        self.label_8.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignCenter)
        self.verticalLayout.addWidget(self.label_7)
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
        self.verticalLayout.addWidget(self.label_8)
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
        self.equilibrate_button.clicked.connect(partial(gui_equilibrate, self))
        self.tableComps.cellChanged.connect(partial(recalculate_after_cell_edit, self))
        self.info_button.clicked.connect(partial(display_about_info, GroupBox))
        self.log_button.clicked.connect(partial(show_log, GroupBox))
        self.retranslateUi(GroupBox)
        QtCore.QMetaObject.connectSlotsByName(GroupBox)

    def retranslateUi(self, GroupBox):
        GroupBox.setWindowTitle(_translate("GroupBox", "Simplistic EC.", None))
        # GroupBox.setTitle(_translate("GroupBox", "Simplistic Eq.", None))
        __sortingEnabled = self.tableComps.isSortingEnabled()
        self.open_button.setText(_translate("GroupBox", "Open", None))
        self.save_button.setText(_translate("GroupBox", "Save", None))
        self.log_button.setText(_translate("GroupBox", "Log", None))
        self.info_button.setText(_translate("GroupBox", "About", None))
        self.equilibrate_button.setText(_translate("GroupBox", "Equilibrate", None))
        self.tableComps.setSortingEnabled(__sortingEnabled)
        self.pushButton.setText(_translate("GroupBox", "Plot", None))
        self.label_2.setText(_translate("GroupBox", "Nr (Reac.)", None))
        self.label.setText(_translate("GroupBox", "n (Comp.)", None))
        self.label_3.setText(_translate("GroupBox", "max. it", None))
        self.label_4.setText(_translate("GroupBox", "tol", None))
        self.label_5.setText(_translate("GroupBox", "solvent", None))
        self.label_6.setText(_translate("GroupBox", "C_solvent (25C)", None))
        self.label_7.setText(_translate("GroupBox", "Horizontal 'X' axis", None))
        self.label_8.setText(_translate("GroupBox", "Vertical 'Y' axis", None))


def open_file(form):
    (filename, _) = \
        QtGui.QFileDialog.getOpenFileName(None,
                                          caption='Open file',
                                          dir=os.path.join(sys.path[0], 'DATA'),
                                          filter='*.csv')
    if os.path.isfile(filename):
        header_comps, comps, header_reacs, reacs = \
            load_csv(form, filename)
        equilibrate(form, header_comps, comps, header_reacs, reacs)
        form.tableComps.sortByColumn(0, QtCore.Qt.AscendingOrder)
        form.tableReacs.sortByColumn(0, QtCore.Qt.AscendingOrder)


def load_csv(form, filename):
    global n, Nr
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
    form.tableComps.setRowCount(n)
    form.tableComps.setColumnCount(len(header_comps) + 3)
    form.tableComps.setHorizontalHeaderLabels(
        header_comps + ['Ceq_i, mol/L', '-log10(C0_i)', '-log10(Ceq_i)'])

    form.tableReacs.setRowCount(Nr)
    form.tableReacs.setColumnCount(len(header_reacs) + 1)
    form.tableReacs.setHorizontalHeaderLabels(
        header_reacs + ['Xieq_j'])
    return header_comps, comps, header_reacs, reacs


def load_QTableWidget(form):
    global component_order_in_table
    n = form.tableComps.rowCount()
    Nr = form.tableReacs.rowCount()
    comps = np.empty([n, 4], dtype='S50')
    reacs = np.empty([Nr, n + 2], dtype='S50')
    header_comps = []
    header_reacs = []
    component_order_in_table = []
    for i in range(form.tableComps.rowCount()):
        component_order_in_table.append(int(form.tableComps.item(i, 0).text()) - 1)
    for i in range(comps.shape[0]):
        header_comps.append(form.tableComps.horizontalHeaderItem(i).data(0))
    for i in range(reacs.shape[1]):
        header_reacs.append(form.tableReacs.horizontalHeaderItem(i).data(0))
    for j in range(comps.shape[1]):
        for i in range(comps.shape[0]):
            comps[component_order_in_table[i], j] = form.tableComps.item(i, j).text()
    for j in range(reacs.shape[1]):
        for i in range(reacs.shape[0]):
            reacs[i, j] = form.tableReacs.item(i, j).text()
    return header_comps, comps, header_reacs, reacs


def recalculate_after_cell_edit(form, row, column):
    gui_equilibrate(form)


def gui_equilibrate(form):
    header_comps, comps, header_reacs, reacs = load_QTableWidget(form)
    equilibrate(form, header_comps, comps, header_reacs, reacs)


def equilibrate(form, header_comps, comps, header_reacs, reacs):
    # Solve
    global C0_i, z_i, nu_ij, pKa_j
    global max_it, tol, index_of_solvent, C_solvent_Tref
    global component_order_in_table
    # Setup logging
    if not os.path.exists('./logs'):
        os.mkdir('./logs')
    logging.basicConfig(filename='./logs/calculation_results.log', level=logging.DEBUG,
                        format='%(asctime)s;%(message)s')

    # Collect variables
    n = len(comps)
    Nr = len(reacs)
    C0_i = np.matrix([row[3] for row in comps], dtype=float).T
    z_i = np.matrix([row[2] for row in comps], dtype=float).T
    nu_ij = np.matrix([row[2:2 + n] for row in reacs], dtype=int).T
    pKa_j = np.matrix([row[1] for row in reacs], dtype=float).T
    max_it = int(form.spinBox_3.value())
    tol = float(form.doubleSpinBox_5.value())

    # Gui setup with calculated values
    form.comboBox.clear()
    form.comboBox_2.clear()
    form.comboBox_3.clear()
    for item in comps[:, 1].T:
        form.comboBox.addItem('C_{' + item + '} | 0')
        form.comboBox_2.addItem('C_{' + item + '} | eq')
        form.comboBox_3.addItem(item)

    highestC0Indexes = np.argpartition(C0_i.A1, (-1, -2))
    index_of_solvent = highestC0Indexes[-1]
    if len(C0_i) > 1:
        index_of_second_highest_C0 = highestC0Indexes[-2]
    else:
        index_of_second_highest_C0 = highestC0Indexes[-1]
    C_solvent_Tref = C0_i[index_of_solvent].item()
    form.comboBox_3.setCurrentIndex(index_of_solvent)
    form.doubleSpinBox_6.setValue(C_solvent_Tref)
    form.doubleSpinBox_6.setPrefix('(mol/L)')

    C_second_highest_C0_Tref = C0_i[index_of_second_highest_C0].item()
    form.comboBox.setCurrentIndex(index_of_second_highest_C0)
    form.comboBox_2.setCurrentIndex(0)

    # Calculate equilibrium composition
    Ceq_i, Xieq_j = calc_Xieq()

    if 'component_order_in_table' in globals():
        i = component_order_in_table
    else:
        i = range(0, n)
    j = range(0, 4 + 3)

    # As usual, problems occurr when sorting is combined with setting QTableWidgetItems
    form.tableComps.setSortingEnabled(False)
    form.tableReacs.setSortingEnabled(False)

    form.tableComps.blockSignals(True)
    form.tableReacs.blockSignals(True)

    for column in j:
        for row in i:
            if column < 4:
                newItem = QtGui.QTableWidgetItem(str(comps[row, column]))
            elif column == 4:
                newItem = QtGui.QTableWidgetItem(str(Ceq_i[row].item()))
            elif column == 5:
                if C0_i[row] <= 0:
                    newItem = QtGui.QTableWidgetItem(str(np.nan))
                else:
                    newItem = QtGui.QTableWidgetItem(str(-np.log10(C0_i[row].item())))
            elif column == 6:
                if Ceq_i[row].item() <= 0:
                    newItem = QtGui.QTableWidgetItem(str(np.nan))
                else:
                    newItem = QtGui.QTableWidgetItem(str(-np.log10(Ceq_i[row].item())))
            # sortierbar machen
            if column != 1:  # Comp. i <Str>
                newItem = NSortableTableWidgetItem(newItem)
                form.tableComps.setItem(row, column, newItem)
            else:
                form.tableComps.setItem(row, column, newItem)
            if not column in range(1, 3 + 1):
                newItem.setFlags(QtCore.Qt.ItemIsEnabled)

    i = range(0, Nr)
    j = range(0, len(header_reacs) + 1)

    for column in j:
        for row in i:
            if column != len(header_reacs):
                form.tableReacs.setItem(row, column, NSortableTableWidgetItem(str(reacs[row][column])))
            elif column == len(header_reacs):
                form.tableReacs.setItem(row, column, NSortableTableWidgetItem(str(Xieq_j[row].item())))

    # Widths and heights, re-enable sorting
    form.tableComps.setSortingEnabled(True)
    form.tableComps.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
    form.tableComps.verticalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)

    form.tableReacs.setSortingEnabled(True)
    form.tableReacs.horizontalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
    form.tableReacs.verticalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)

    form.tableComps.blockSignals(False)
    form.tableReacs.blockSignals(False)


def save_file(form):
    pass


def plot_intervals(form):
    pass


def calc_Xieq():
    """Steepest descent for good initial estimate, then Newton method for non-linear algebraic system
    :return: tuple with Ceq_i, Xieq_j, f_0
    :param C0_i: np.matrix (n X 1) - Conc(i, alimentación)
    :param z_i: np.matrix (n X 1) - Carga(i, alimentación)
    :param nu_ij: np.matrix (n X Nr) - Coefs. esteq. componente i en reacción j
    :param pKa_j: np.matrix (n X 1) - (-1)*log10("Cte." de equilibrio en reacción j) = -log10 Kc_j(T)
    :param Xieq_0: np.matrix (n X 1) - avance de reacción j - estimado inicial
    :param Ceq_0: np.matrix (n X 1) - Conc(i, equilibrio)
    """
    global C0_i, z_i, nu_ij, pKa_j, max_it, tol, n, Nr, Kc_j, index_of_solvent, C_solvent_Tref
    n = nu_ij.shape[0]
    Nr = nu_ij.shape[1]
    Kc_j = np.multiply(np.power(10, -pKa_j), np.power(C_solvent_Tref, nu_ij[index_of_solvent, :]).T)
    Xieq_j = np.matrix(np.zeros([Nr, 1]))
    X0 = np.concatenate([C0_i + abs(nu_ij * np.matrix(np.ones([Nr, 1])) * tol), Xieq_j])
    X = X0
    # Steepest descent: min(g(X))=min(f(X).T*f(X))
    X, F_val = steepest_descent(X, f_gl_0, J, g, 1.0e-3)
    """
    if any(map(lambda x: np.sign(x)==-1, X[0:n])):
        X = X0
        minus_f = lambda x: -1*f_gl_0(x)
        minus_J = lambda x: -1*J(x)
        X, F_val = steepest_descent(X0, minus_f, minus_J, g, 1.0e-3)
    """
    # Newton method: G(X) = J(X)^-1 * F(X)
    k = 0
    J_val = J(X)
    Y = np.matrix(np.ones(len(X))).T * tol / (np.sqrt(len(X)) * tol)
    diff = np.matrix(np.empty([len(X), 1]))
    diff.fill(np.nan)
    stop = False
    while k <= max_it and not stop:
        new_log_Entry('Newton', k, X, diff, F_val, Y, np.nan, np.nan, stop)
        X_k_m_1 = X
        Y = gausselimination(J_val, -F_val)
        X = X + Y
        diff = X - X_k_m_1
        J_val = J(X)
        F_val = f_gl_0(X)
        if np.sqrt((Y.T * Y).item()) < tol:
            stop = True  # Procedure successful
        k += 1
    new_log_Entry('Newton', k, X, diff, F_val, Y, np.nan, np.nan, stop)
    Ceq_i = X[0:n]
    Xieq_j = X[n:n + Nr]
    return Ceq_i, Xieq_j


def f_gl_0(X):
    global C0_i, nu_ij, n, Nr, Kc_j
    Ceq_i = X[0:n, 0]
    Xieq_j = X[n:n + Nr, 0]
    f_gl_0 = np.matrix(np.empty([n + Nr, 1], dtype=float))
    f_gl_0[0:n] = -Ceq_i + C0_i + nu_ij * Xieq_j
    f_gl_0[n:n + Nr] = -Kc_j + np.prod(np.power(Ceq_i, nu_ij), 0).T
    return f_gl_0


def J(X):
    global C0_i, nu_ij, n, Nr, Kc_j
    Ceq_i = X[0:n, 0]
    Eins_durch_C = np.diag(np.power(Ceq_i, -1).A1, 0)
    Quotient = np.diag(np.prod(np.power(Ceq_i, nu_ij), 0).A1)
    J = np.matrix(np.zeros([n + Nr, n + Nr], dtype=float))
    J[0:n, 0:n] = -1 * np.eye(n).astype(float)
    J[0:n, n:n + Nr] = nu_ij
    J[n:n + Nr, 0:n] = Quotient * nu_ij.T * Eins_durch_C
    return J


def g(X):
    f_0 = f_gl_0(X)
    return (f_0.T * f_0).item()


def steepest_descent(X0, f, J, g, tol):
    X = X0
    f_val = f(X)
    k = 0
    stop = False
    diff = np.matrix(np.empty([len(X), 1]))
    diff.fill(np.nan)
    Y = np.matrix(np.empty([len(X), 1]))
    Y.fill(np.nan)
    while k < max_it and not stop:
        z = 2 * J(X).T * f(X)  # z(X) = nabla(g(X)) = 2*J(X).T*F(X)
        z0 = np.sqrt((z.T * z).item())
        if z0 == 0:
            # Zero gradient
            stop = True
            break
        z = z / z0
        alpha1 = 0
        alpha3 = 1
        g1 = g(X - alpha1 * z)
        g3 = g(X - alpha3 * z)
        while g3 >= g1 and alpha3 > tol / 2.0:
            alpha3 = alpha3 / 2.0
            g3 = g(X - alpha3 * z)
        alpha2 = alpha3 / 2.0
        g2 = g(X - alpha2 * z)
        """
        (Note: Newton’s forward divided-difference formula is used to find
        the quadratic P(α) = g1 + h1α + h3α(α − α2) that interpolates
        h(α) at α = 0, α = α2, α = α3.)
        """
        h1 = (g2 - g1) / alpha2
        h2 = (g3 - g2) / (alpha3 - alpha2)
        h3 = (h2 - h1) / alpha3
        alpha0 = 0.5 * (alpha2 - h1 / h3)  # (The critical point of P occurs at α0.)
        g0 = g(X - alpha0 * z)
        if g0 < g3:
            alpha = alpha0
            g_min = g0
        else:
            alpha = alpha3
            g_min = g3
        if abs(g_min - g1) < tol:
            stop = True  # Procedure successful
        new_log_Entry('Steepest descent', k, X, diff, f_val, Y, g_min, g1, stop)
        X_k_m_1 = X
        X = X - alpha * z
        diff = X_k_m_1 - X
        f_val = f(X)
        k += 1
    new_log_Entry('Steepest descent', k, X, diff, f_val, Y, g_min, g1, stop)
    return X, f_val


def new_log_Entry(method, k, X, diff, f_val, Y, g_min, g1, stop):
    logging.debug(method + ' method loop;' +
                  'k=' + str(k) +
                  ';X=' + '[' + ','.join(map(str, X.T.A1)) + ']' +
                  ';||X(k)-X(k-1)||=' + str((diff.T * diff).item()) +
                  ';f(X)=' + '[' + ','.join(map(str, f_val.T.A1)) + ']' +
                  ';||f(X)||=' + str(np.sqrt((f_val.T * f_val).item())) +
                  ';Y=' + '[' + ','.join(map(str, Y.T.A1)) + ']' +
                  ';||Y||=' + str(np.sqrt((Y.T * Y).item())) +
                  ';g=' + str(g_min) +
                  ';|g-g1|=' + str(abs(g_min - g1)) +
                  ';stop=' + str(stop))


def f_test(X):
    X1, X2, X3 = X[0].item(), X[1].item(), X[2].item()
    f = np.matrix([ \
        3 * X1 - np.cos(X2 * X3) - 1.0 / 2.0, \
        X1 ** 2 - 81.0 * (X2 + 0.1) ** 2 + np.sin(X3) + 1.06, \
        np.exp(-X1 * X2) + 20.0 * X3 + (10 * np.pi - 3) / 3.0
    ]).T
    return f


def J_test(X):
    X1, X2, X3 = X[0].item(), X[1].item(), X[2].item()
    J = np.matrix([ \
        [3, X2 * np.sin(X3), X3 * np.sin(X2)], \
        [2 * X1, -81.0 * 2 * (X2 + 0.1), np.cos(X3)], \
        [-X2 * np.exp(-X1 * X2), -X1 * np.exp(-X1 * X2), 20.0]
    ])
    return J


def g_test(X):
    f_0 = f_test(X)
    return (f_0.T * f_0).item()


def display_about_info(form):
    textStreamOutput = ''
    with open('README.md') as readme_file:
        for row in readme_file:
            textStreamOutput += row
    readme_file.close()
    aboutBox_1 = QtGui.QMessageBox.about(form, 'simplistic Kc eq calculation', textStreamOutput)


def show_log(form):
    headers_and_types = np.array(
        (('date', str),
         ('method', str),
         ('k', int),
         ('X', list),
         ('||X(k)-X(k-1)||', float),
         ('f(X)', list),
         ('||f(X)||', float),
         ('Y', list),
         ('||Y||', float),
         ('g', float),
         ('|g-g1|', float),
         ('stop', bool)))

    headers_and_types_dict = dict(headers_and_types)
    col_numbers_with_float = \
        np.argwhere(map(lambda x: x == float, headers_and_types[:, 1]))
    col_numbers_with_list = \
        np.argwhere(map(lambda x: x == list, headers_and_types[:, 1]))
    col_numbers_with_int = \
        np.argwhere(map(lambda x: x == int, headers_and_types[:, 1]))
    col_numbers_with_str = \
        np.argwhere(map(lambda x: x == str, headers_and_types[:, 1]))
    col_numbers_with_bool = \
        np.argwhere(map(lambda x: x == bool, headers_and_types[:, 1]))

    take_float = lambda x: float(x.rpartition('=')[-1])
    take_list = lambda x: \
        np.fromstring(x.rpartition('=')[-1]
                      .replace('[', '').replace(']', ''),
                      sep=',')
    take_int = lambda x: int(x.rpartition('=')[-1])
    take_bool = lambda x: x.rpartition('=')[-1] == 'True'
    take_date = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S,%f')

    cell_conversions = dict.fromkeys(headers_and_types_dict.keys())

    for i in range(len(headers_and_types)):
        if i == 0:
            cell_conversions[headers_and_types[i, 0]] = take_date
        elif i in col_numbers_with_float:
            cell_conversions[headers_and_types[i, 0]] = take_float
        elif i in col_numbers_with_list:
            cell_conversions[headers_and_types[i, 0]] = take_list
        elif i in col_numbers_with_int:
            cell_conversions[headers_and_types[i, 0]] = take_int
        elif i in col_numbers_with_bool:
            cell_conversions[headers_and_types[i, 0]] = take_bool
        elif i in col_numbers_with_str:
            pass
        i += 1

    log = pd.read_csv(
        filepath_or_buffer='./logs/calculation_results.log',
        delimiter=';',
        names=headers_and_types[:, 0],
        index_col=False,
        converters=cell_conversions)
    form.logWidget = LogWidget(log)



class LogWidget(QtGui.QWidget):
    def __init__(self, _log, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.log = _log
        self.setupUi()
        self.group_2.show()

    def setupUi(self):
        minLogHeight = 400
        minLogWidth = minLogHeight*16/9
        displayItemsByPage = 50
        self.group_2 = QtGui.QGroupBox()
        self.group_2.resize(minLogWidth,minLogHeight)
        self.verticalLayout = QtGui.QVBoxLayout(self.group_2)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.group_2.setLayout(self.verticalLayout)
        self.pandasView = QtGui.QTableView(self.group_2)
        self.firstButton = QtGui.QPushButton(self.group_2)
        self.lastButton = QtGui.QPushButton(self.group_2)
        self.nextButton = QtGui.QPushButton(self.group_2)
        self.previousButton = QtGui.QPushButton(self.group_2)
        self.pageLabel = QtGui.QLabel(self.group_2)
        self.totPagesLabel = QtGui.QLabel(self.group_2)
        self.pageBox = QtGui.QLineEdit(self.group_2)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout.addWidget(self.pandasView)
        self.horizontalLayout.addWidget(self.firstButton)
        self.horizontalLayout.addWidget(self.previousButton)
        self.horizontalLayout.addWidget(self.pageLabel)
        self.horizontalLayout.addWidget(self.pageBox)
        self.horizontalLayout.addWidget(self.totPagesLabel)
        self.horizontalLayout.addWidget(self.nextButton)
        self.horizontalLayout.addWidget(self.lastButton)
        self.pandasModel = PandasModel(self.log.tail(displayItemsByPage))
        self.pandasView.setModel(self.pandasModel)

        # To ensure full display, first set resize modes, then resize columns to contents
        self.pandasView.horizontalHeader().setResizeMode(QtGui.QHeaderView.Interactive)
        self.pandasView.verticalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.pandasView.resizeColumnsToContents()
        self.pageBox.setMaximumWidth(int(round(minLogHeight/float(5))))
        self.pageBox.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignCenter)

        lastPage = int(round(len(self.log.values)/float(displayItemsByPage)))
        currentPageFirstEntry = len(self.log.values)-50
        currentPageLastEntry = len(self.log.values)
        self.totPagesLabel.setText(' / ' + str(lastPage))
        self.pageLabel.setText('Entries ' + str(currentPageFirstEntry) +
                               ' to ' + str(currentPageLastEntry)+ '; Page: ')
        self.pageBox.setText(str(lastPage))
        self.firstButton.setText('<< First')
        self.lastButton.setText('Last >>')
        self.previousButton.setText('< Previous')
        self.nextButton.setText('Next >')


class PandasModel(QtCore.QAbstractTableModel):
    """
    Used to populate a QTableView with the pandas model
    """

    def __init__(self, data, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return len(self._data.values)

    def columnCount(self, parent=None):
        return self._data.columns.size

    def data(self, index, role=QtCore.Qt.DisplayPropertyRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                return str(self._data.values[index.row()][index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self._data.columns[col]
        return None


class aboutBox(QtGui.QMessageBox):
    def __init__(self, parent=None):
        QtGui.QMessageBox.__init__(self, parent)

    def about(self, title_text, contained_text):
        QtGui.QMessageBox.__init__(self, None)
        self.setText(title_text)
        self.setDetailedText(contained_text)
        self.show()


class NSortableTableWidgetItem(QtGui.QTableWidgetItem):
    # Implement less than (<) for numeric table widget items.
    def __init__(self, text):
        QtGui.QTableWidgetItem.__init__(self, text)

    def __lt__(self, y):
        float_self = float(self.text())
        if np.isnan(float_self):
            return True
        else:
            return float(self.text()) < float(y.text())


class MainForm(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.setWindowIcon(QtGui.QIcon(
            os.path.join(sys.path[0], *['utils', 'icon_batch.png'])))
        self.ui = Ui_GroupBox()
        self.ui.setupUi(self)


if __name__ == '__main__':
    app = QtGui.QApplication.instance()  # checks if QApplication already exists
    if not app:  # create QApplication if it doesnt exist
        app = QtGui.QApplication(sys.argv)

    main_form = MainForm()
    main_form.show()
    header_comps, comps, header_reacs, reacs = \
        load_csv(main_form.ui, './DATA/COMPONENTS_REACTIONS_EX_001.csv')
    equilibrate(main_form.ui, header_comps, comps, header_reacs, reacs)
    main_form.ui.tableComps.sortByColumn(0, QtCore.Qt.AscendingOrder)
    main_form.ui.tableReacs.sortByColumn(0, QtCore.Qt.AscendingOrder)
    sys.exit(app.exec_())
