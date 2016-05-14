# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 20:48:42 2015

@author: Santiago Salas
@ref: Denbigh, p. 298
"""
import os
import sys
import logging
import re
import pandas as pd
import numpy as np
import csv
import bisect
import uuid
import urllib
import matplotlib
import colormaps
import ctypes  # Needed to set the app icon correctly
from functools import partial
from mat_Zerlegungen import gausselimination
from datetime import datetime

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PySide import QtGui, QtCore, QtWebKit
from mpldatacursor import datacursor

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromutf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8

    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)


def take_float(x):
    return float(x.rpartition('=')[-1])


def take_list(x):
    return np.fromstring(x.rpartition('=')[-1]
                         .replace('[', '').replace(']', ''),
                         sep=',')


def take_int(x):
    return int(x.rpartition('=')[-1])


def take_bool(x):
    return x.rpartition('=')[-1] == 'True'


def take_date(x):
    return datetime.strptime(x, '%Y-%m-%d %H:%M:%S,%f')

# Variables used
colormap_colors = colormaps.viridis.colors + colormaps.inferno.colors
markers = matplotlib.markers.MarkerStyle.filled_markers
fillstyles = matplotlib.markers.MarkerStyle.fillstyles
# Structure of input files
# Header of components will match this expression, only need to find
# indexes in file.
header_comps_input_model = [
    'i', 'Comp.', 'z', 'M/(g/mol)', 'w0/g', 'xw0',
    'n0/mol', 'x0', 'c0/(mol/L)',
    'm0/(mol/kg_{solvent})'
]
comp_variable_input_names = [
    'index', 'comp_id', 'z', 'molar_mass',
    'w0', 'xw0', 'n0', 'x0', 'c0', 'm0'
]
header_comps_output_model = [
    'weq/g', 'xweq', 'neq/mol', 'xeq',
    'ceq/(mol/L)', 'meq/(mol/kg_{solvent})',
    'rhoeq/(g/mL)',
    '\gamma_{eq}^{II}', '\gamma_{eq}^{III}',
    'aeq',
    '-log10(xw0)', '-log10(x0)', '-log10(c0)',
    '-log10(m0)', '-log10(a0)',
    '-log10(\gamma_{eq}^{II})',
    '-log10(\gamma_{eq}^{III})',
    '-log10(xweq)', '-log10(xeq)', '-log10(ceq)',
    '-log10(meq)', '-log10(aeq)'
]
comp_variable_output_names = [
    'weq', 'xweq', 'neq', 'xeq',
    'ceq', 'meq',
    'rhoeq',
    'gammaeq_ii', 'gammaeq_iii',
    'aeq',
    'mlog10xw0', 'mlog10x0', 'mlog10c0',
    'mlog10m0', 'mlog10a0',
    'mlog10gammaeq_ii',
    'mlog10gammaeq_iii',
    'mlog10xweq', 'mlog10xeq', 'mlog10ceq',
    'mlog10meq', 'mlog10aeq',
]
# Store programatic name vs. header name in a dict.
comp_names_headers = \
    dict(zip(header_comps_input_model + header_comps_output_model,
             comp_variable_input_names + comp_variable_output_names))
# Regular expression compilations
float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')
matchingHLine = re.compile('=+')
# Default component and reaction headers according to structure
reac_headers_re = re.compile(
    r'(\bj)' +
    r'|(^pKaj$)' +
    r'|nu_?i?([0-9]+)?j(\(i=([0-9]+)\))?')  # For now, no more than 999 reactions
comp_headers_re = re.compile(
    r'(\bi)|(Comp\.?i?)|([z|Z]_?i?)|(M_?i?)|(w_?0_?i?)|(x_?w_?0?)|' +
    r'(n_?0?_?i?)|(x_?0)|([c|C]_?0_?i?)|(m_?0)')
doc_hline_re = re.compile(
    r'(\s*-{3,})')
html_title = re.compile('<title>(.*?)</title>',
                        re.IGNORECASE | re.DOTALL)


class UiGroupBox(QtGui.QWidget):
    _was_canceled = False

    def __init__(self, parent):
        QtGui.QWidget.__init__(self, parent)
        # Default size
        parent.resize(500, 357)
        # Assignments
        self.verticalLayout_2 = QtGui.QVBoxLayout(parent)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.open_button = QtGui.QPushButton()
        self.save_button = QtGui.QPushButton()
        self.info_button = QtGui.QPushButton()
        self.log_button = QtGui.QPushButton()
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.label_3 = QtGui.QLabel()
        self.equilibrate_button = QtGui.QPushButton()
        self.spinBox_3 = QtGui.QSpinBox()
        self.label_4 = QtGui.QLabel()
        self.label = QtGui.QLabel()
        self.spinBox = QtGui.QSpinBox()
        self.tableComps = QtGui.QTableWidget()
        self.label_9 = QtGui.QLabel()
        self.progress_var = QtGui.QProgressBar(parent)
        self.cancelButton = QtGui.QPushButton(parent)
        self.doubleSpinBox_5 = ScientificDoubleSpinBox()
        item = QtGui.QTableWidgetItem()
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.label_5 = QtGui.QLabel()
        self.comboBox_3 = QtGui.QComboBox()
        self.label_6 = QtGui.QLabel()
        self.doubleSpinBox_6 = QtGui.QDoubleSpinBox()
        self.label_2 = QtGui.QLabel()
        self.spinBox_2 = QtGui.QSpinBox()
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.tableReacs = QtGui.QTableWidget()
        self.verticalLayout = QtGui.QVBoxLayout()
        self.label_7 = QtGui.QLabel()
        self.comboBox = QtGui.QComboBox()
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.doubleSpinBox = ScientificDoubleSpinBox()
        self.doubleSpinBox_2 = ScientificDoubleSpinBox()
        self.plotButton = QtGui.QPushButton()
        self.radio_group = QtGui.QHBoxLayout()
        self.radio_b_1 = QtGui.QRadioButton()
        self.radio_b_2 = QtGui.QRadioButton()
        self.radio_b_3 = QtGui.QRadioButton()
        self.groupBox = None  # self.groupBox will contain either the plotBox or logBox
        # Object names
        parent.setObjectName(_fromutf8("GroupBox"))
        self.verticalLayout_2.setObjectName(_fromutf8("verticalLayout_2"))
        self.horizontalLayout_2.setObjectName(_fromutf8("horizontalLayout_2"))
        self.horizontalLayout_3.setObjectName(_fromutf8("horizontalLayout_3"))
        self.open_button.setObjectName(_fromutf8("open_button"))
        self.save_button.setObjectName(_fromutf8("save_button"))
        self.info_button.setObjectName(_fromutf8("info_button"))
        self.log_button.setObjectName(_fromutf8("log_button"))
        self.horizontalLayout_5.setObjectName(_fromutf8("horizontalLayout_5"))
        self.label_3.setObjectName(_fromutf8("max_it_label"))
        self.spinBox_3.setObjectName(_fromutf8("max_it_spinbox"))
        self.label_4.setObjectName(_fromutf8("tol_label"))
        self.doubleSpinBox_5.setObjectName(_fromutf8("tol_spinbox"))
        self.label.setObjectName(_fromutf8("label"))
        self.spinBox.setObjectName(_fromutf8("spinBox"))
        self.tableComps.setObjectName(_fromutf8("tableComps"))
        self.horizontalLayout_6.setObjectName(_fromutf8("horizontalLayout_6"))
        self.label_5.setObjectName(_fromutf8("solvent_label"))
        self.label_6.setObjectName(_fromutf8("c_solvent_tref"))
        self.doubleSpinBox_6.setObjectName(
            _fromutf8("c_solvent_tref_doublespinbox"))
        self.label_2.setObjectName(_fromutf8("label_2"))
        self.spinBox_2.setObjectName(_fromutf8("spinBox_2"))
        self.horizontalLayout.setObjectName(_fromutf8("horizontalLayout"))
        self.tableReacs.setObjectName(_fromutf8("tableReacs"))
        self.verticalLayout.setObjectName(_fromutf8("verticalLayout"))
        self.label_7.setObjectName(_fromutf8("horizontalAxisLabel"))
        self.comboBox.setObjectName(_fromutf8("comboBox"))
        self.horizontalLayout_3.setObjectName(_fromutf8("horizontalLayout_3"))
        self.doubleSpinBox.setObjectName(_fromutf8("doubleSpinBox"))
        self.doubleSpinBox_2.setObjectName(_fromutf8("doubleSpinBox_2"))
        self.plotButton.setObjectName(_fromutf8("plotButton"))
        self.radio_group.setObjectName(_fromutf8("radio_group"))
        self.radio_b_1.setObjectName(_fromutf8("radio_b_1"))
        self.radio_b_2.setObjectName(_fromutf8("radio_b_2"))
        self.radio_b_3.setObjectName(_fromutf8("radio_b_3"))
        # Operations
        self.horizontalLayout_2.addWidget(self.open_button)
        self.horizontalLayout_2.addWidget(self.save_button)
        self.horizontalLayout_2.addWidget(self.log_button)
        self.horizontalLayout_2.addWidget(self.info_button)
        self.verticalLayout_2.addWidget(self.equilibrate_button)
        self.radio_group.addWidget(self.radio_b_1)
        self.radio_group.addWidget(self.radio_b_2)
        self.radio_group.addWidget(self.radio_b_3)
        self.radio_b_1.setChecked(True)
        self.radio_b_2.setEnabled(False)
        self.radio_b_3.setEnabled(False)
        self.radio_b_1.setToolTip('<b>%s</b><br><img src="%s">' %
                                  ('Ideal solution',
                                   'utils/Ideal_solution.png'))
        self.radio_b_2.setToolTip('<b>%s</b><br><img src="%s">' %
                                  ('Debye-Hueckel',
                                   'utils/Debye-Hueckel.png'))
        self.radio_b_3.setToolTip('<b>%s</b><br><img src="%s">' %
                                  ('Davies (DH ext.)',
                                   'utils/Davies.png'))
        self.verticalLayout_2.addLayout(self.radio_group)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.label_3.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label_3)
        self.spinBox_3.setMaximum(2000)
        self.spinBox_3.setMinimum(2)
        self.spinBox_3.setProperty("value", 1000)
        self.horizontalLayout_5.addWidget(self.spinBox_3)
        self.label_4.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label_4)
        self.doubleSpinBox_5.setDecimals(
            int(-np.log10(np.finfo(float).eps) + 1))
        self.doubleSpinBox_5.setMaximum(float(1))
        self.doubleSpinBox_5.setMinimum(np.finfo(float).eps * 2)
        self.doubleSpinBox_5.setSingleStep(
            0.1 / 100.0 * (1.0 - np.finfo(float).eps))
        self.doubleSpinBox_5.setProperty("value", float(1.0e-08))
        self.horizontalLayout_5.addWidget(self.doubleSpinBox_5)
        self.label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label)
        self.spinBox.setProperty("value", 0)
        self.horizontalLayout_5.addWidget(self.spinBox)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.tableComps.setMinimumSize(QtCore.QSize(0, 210))
        item.setTextAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.tableComps.horizontalHeader().setCascadingSectionResizes(False)
        self.tableComps.horizontalHeader().setDefaultSectionSize(100)
        self.tableComps.horizontalHeader().setMinimumSectionSize(27)
        self.tableComps.horizontalHeader().setSortIndicatorShown(True)
        self.tableComps.verticalHeader().setVisible(False)
        self.verticalLayout_2.addWidget(self.tableComps)
        self.label_9.setAlignment(QtCore.Qt.AlignTop)
        self.label_9.setFrameStyle(QtGui.QFrame.Box | QtGui.QFrame.Raised)
        self.verticalLayout_2.addWidget(self.label_9)
        self.horizontalLayout_7.setAlignment(QtCore.Qt.AlignRight)
        self.horizontalLayout_7.setSizeConstraint(QtGui.QLayout.SetFixedSize)
        self.horizontalLayout_7.addStrut(max(
            [self.progress_var.frameSize().height(),
             self.cancelButton.frameSize().height()]))
        self.horizontalLayout_7.setAlignment(QtCore.Qt.AlignLeft)
        self.horizontalLayout_7.addWidget(self.cancelButton)
        self.horizontalLayout_7.addWidget(self.progress_var)
        self.cancelButton.setEnabled(False)
        self.progress_var.setEnabled(False)
        self.verticalLayout_2.addLayout(self.horizontalLayout_7)
        self.label_5.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_6.addWidget(self.label_5)
        self.horizontalLayout_6.addWidget(self.comboBox_3)
        self.label_6.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_6.addWidget(self.label_6)
        self.doubleSpinBox_6.setDecimals(
            int(-np.log10(np.finfo(float).eps) + 1))
        self.doubleSpinBox_6.setMaximum(float(1000))
        self.doubleSpinBox_6.setMinimum(np.finfo(float).eps * 1.1)
        self.doubleSpinBox_6.setSingleStep(1.0e-2)
        self.horizontalLayout_6.addWidget(self.doubleSpinBox_6)
        self.label_2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
        self.horizontalLayout_6.addWidget(self.label_2)
        self.spinBox_2.setProperty("value", 0)
        self.horizontalLayout_6.addWidget(self.spinBox_2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_6)
        self.tableReacs.horizontalHeader().setVisible(True)
        self.tableReacs.verticalHeader().setVisible(False)
        self.horizontalLayout.addWidget(self.tableReacs)
        self.label_7.setAlignment(
            QtCore.Qt.AlignHCenter | QtCore.Qt.AlignCenter)
        self.verticalLayout.addWidget(self.label_7)
        self.verticalLayout.addWidget(self.comboBox)
        self.doubleSpinBox.setMinimum(0.0)
        self.horizontalLayout_3.addWidget(self.doubleSpinBox)
        self.doubleSpinBox_2.setMinimum(0.0)
        self.horizontalLayout_3.addWidget(self.doubleSpinBox_2)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.verticalLayout.addWidget(self.plotButton)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        # Events
        self.open_button.clicked.connect(partial(self.open_file))
        self.save_button.clicked.connect(partial(self.save_file))
        self.plotButton.clicked.connect(partial(self.solve_intervals))
        self.equilibrate_button.clicked.connect(
            partial(self.recalculate_after_cell_edit, 0, 0))
        self.tableComps.cellChanged.connect(
            partial(self.recalculate_after_cell_edit))
        self.info_button.clicked.connect(partial(self.display_about_info))
        self.log_button.clicked.connect(partial(self.show_log))
        self.cancelButton.clicked.connect(partial(self.cancel_loop))
        self.comboBox.currentIndexChanged.connect(
            partial(self.populate_input_spinboxes))

        # Icons
        icon0 = QtGui.QIcon()
        icon0.addPixmap(QtGui.QPixmap(
            _fromutf8("utils/glyphicons-145-folder-open.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(
            _fromutf8("utils/glyphicons-415-disk-save.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(
            _fromutf8("utils/glyphicons-41-stats.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(
            _fromutf8("utils/glyphicons-82-refresh.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(
            _fromutf8("utils/glyphicons-196-circle-info.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        icon5 = QtGui.QIcon()
        icon5.addPixmap(QtGui.QPixmap(
            _fromutf8("utils/glyphicons-88-log-book.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.open_button.setIcon(icon0)
        self.save_button.setIcon(icon1)
        self.plotButton.setIcon(icon2)
        self.equilibrate_button.setIcon(icon3)
        self.info_button.setIcon(icon4)
        self.log_button.setIcon(icon5)
        # Retranslate, connect
        self.retranslate_ui(parent)
        QtCore.QMetaObject.connectSlotsByName(parent)

    def retranslate_ui(self, parent):
        parent.setWindowTitle(_translate("parent", "Homogeneous EC.", None))
        parent.setTitle(QtGui.QApplication.translate("parent", "EC", None))
        __sortingEnabled = self.tableComps.isSortingEnabled()
        self.open_button.setText(_translate("parent", "Open", None))
        self.save_button.setText(_translate("parent", "Save", None))
        self.log_button.setText(_translate("parent", "Log", None))
        self.info_button.setText(_translate("parent", "About", None))
        self.equilibrate_button.setText(
            _translate("parent", "Equilibrate", None))
        self.radio_b_1.setText(
            _translate("parent", "Ideal solution", None))
        self.radio_b_2.setText(
            _translate("parent", "Debye-Hückel", None))
        self.radio_b_3.setText(
            _translate("parent", "Davies (DH ext.)", None))
        self.tableComps.setSortingEnabled(__sortingEnabled)
        self.plotButton.setText(_translate("parent", "Plot", None))
        self.label_2.setText(_translate("parent", "nr (Reac.)", None))
        self.label.setText(_translate("parent", "n (Comp.)", None))
        self.label_3.setText(_translate("parent", "max. it", None))
        self.label_4.setText(_translate("parent", "tol", None))
        self.label_5.setText(_translate("parent", "solvent", None))
        self.label_6.setText(_translate("parent", "C_solvent (25C)", None))
        self.label_7.setText(_translate("parent", "Horizontal 'X' axis", None))
        self.label_9.setText(
            _translate(
                "parent",
                'Currently unequilibrated',
                None))
        self.cancelButton.setText('cancel')

    def cancel_loop(self):
        self._was_canceled = True
        self.progress_var.setValue(0)

    def populate_input_spinboxes(self, index):
        comps = self.comps
        c0_component = self.c0[index]
        self.doubleSpinBox.setValue(c0_component / 10.0 ** 7)
        self.doubleSpinBox_2.setValue(c0_component * (1 + 20 / 100.0))

    def remove_canceled_status(self):
        self._was_canceled = False

    def was_canceled(self):
        return self._was_canceled

    def open_file(self):
        (filename, _) = \
            QtGui.QFileDialog.getOpenFileName(None,
                                              caption='Open file',
                                              dir=os.path.join(
                                                  sys.path[0], 'DATA'),
                                              filter='*.csv')
        if os.path.isfile(filename):
            # Reset solution state and order of items
            if hasattr(self, 'acceptable_solution'):
                delattr(self, 'acceptable_solution')
            if hasattr(self, 'component_order_in_table'):
                delattr(self, 'component_order_in_table')
            # Load csv data into form variables
            self.load_csv(filename)

            # Continue with typical solution and table population procedure
            self.gui_equilibrate()
            self.tableComps.sortByColumn(0, QtCore.Qt.AscendingOrder)
            self.tableReacs.sortByColumn(0, QtCore.Qt.AscendingOrder)

    def load_csv(self, filename):
        with open(filename) as csv_file:
            n = 0
            nr = 0
            header_comps = []
            header_reacs = []
            reader = csv.reader(csv_file, dialect='excel')
            reading_comps = False
            reading_reacs = False
            comps = []
            reacs = []
            valid_columns_reacs = []
            valid_columns_comps = []
            # First two columns of the header of reacions will match this
            # expression, need to find indexes in file and append coefficient
            # matrix.
            header_reacs_model = ['j', 'pKa']
            max_row_length = len([])
            for row in reader:
                row_without_whitespace = [x.replace(' ', '') for x in row]
                row_without_blanks = [
                    x for x in row_without_whitespace if len(x) > 0]
                if len(row_without_blanks) == 0:
                    pass  # skip empty line
                elif 'COMP' in row:
                    reading_comps = True
                    reading_reacs = False
                    header_comps = next(reader)
                    # Get column mappings of form
                    # [
                    # [output_index_1, input_index_1],
                    # [output_index_2, input_index_2],
                    # ...]
                    io_column_index_map = []
                    for (col_no, column) in enumerate(header_comps):
                        old_index = col_no
                        matches = comp_headers_re.match(column)
                        if matches is not None:
                            new_index = [j for j, match \
                                         in enumerate(matches.groups()) \
                                         if  match is not None]
                            io_column_index_map.append(
                                [new_index[0], old_index]
                            )
                elif 'REAC' in row:
                    reading_reacs = True
                    reading_comps = False
                    header_reacs = next(reader)
                    # Enumerate (index, column)
                    # Remove spaces.
                    header_reacs_enum = [(i, j.replace(' ', '')) for i, j in
                                         enumerate(header_reacs) if
                                         j.replace(' ', '') != '']
                    valid_columns_reacs = [x[0] for x in header_reacs_enum]
                elif reading_comps:
                    n += 1
                    # put 0 instead of blank and keep all columns to add in
                    # model
                    row_to_add = ['']*len(header_comps_input_model)
                    for new_index, old_index  in io_column_index_map:
                        row_to_add[new_index] = \
                            row_without_whitespace[old_index]
                    comps.append(row_to_add)
                elif reading_reacs:
                    nr += 1
                    # put 0 instead of blank and keep only columns to add
                    row_to_add = ['0' if row_without_whitespace[x] == '' else
                                  row_without_whitespace[x]
                                  for x in valid_columns_reacs]
                    reacs.append(row_to_add)
        csv_file.close()
        # Rearrange reacs to pKaj nu_ij (i=1) nu_ij(i=2) ... nu_ij(i=n)
        # Rearrange comps rearrange to Comp. i zi c0
        header_reacs_groups_enum = \
            map(lambda y: [y[0], reac_headers_re.match(y[1]).groups()],
                header_reacs_enum)
        priorities_reacs = []
        # Get reaction set structure based on same order as regex groups.
        # Group 2 corresponds to nu_2j or nu_2j(i=2). First is priority.
        # Group 4 corresponds to nu_ij(i=2). If redundant with 2, 2 used.
        for x in header_reacs_groups_enum:
            if x is not None:
                index_origin = x[0]
                if x[1][2] is None \
                        and x[1][4] is None:
                    variable = [y for y in x[1] if y is not None][0]
                    index_destination = [i for i, j in enumerate(x[1]) if j is not None][
                        0]
                elif x[1][2] is not None:
                    variable = 'nu_ij(i=' + x[1][2] + ')'
                    # add space for [j, pKaj]
                    index_destination = int(x[1][2]) + 0 + 1
                elif x[1][2] is None \
                        and x[1][4] is not None:
                    variable = 'nu_ij(i=' + x[1][4] + ')'
                    # add space for [j, pKaj]
                    index_destination = int(x[1][4]) + 0 + 1
                priorities_reacs.append(
                    (variable, index_origin, index_destination))
        # reacs header, form
        # [None, ..., None, 'pKaj', i=1', 'i=2', ... 'i=n']
        # comps header, form
        # [None, ..., None, 'Comp. i', 'z', 'M', 'n0', 'w0']
        # Sort priorities arrays by destination
        sorted_priorities_reacs = sorted(priorities_reacs, key=lambda y: y[2])
        header_reacs = \
            [str(x[0]) for x in sorted_priorities_reacs]
        # Add components index form. All have been worked on in
        # scaffold.
        header_comps = comp_variable_input_names
        comps = np.array(comps)
        reacs = np.array(reacs)
        self.spinBox.setProperty("value", n)
        self.spinBox_2.setProperty("value", nr)
        self.tableComps.setRowCount(n)
        self.tableComps.setColumnCount(len(header_comps_input_model) +
                                       len(header_comps_output_model))
        self.tableComps.setHorizontalHeaderLabels(
            header_comps_input_model + header_comps_output_model)
        self.tableReacs.setRowCount(nr)
        self.tableReacs.setColumnCount(n + 2 + 1)
        self.tableReacs.setHorizontalHeaderLabels(
            header_reacs + ['xieq'])

        # Pass variables to self before loop start
        variables_to_pass = ['header_comps',
                             'comps', 'header_reacs', 'reacs',
                             'n', 'nr'
                             ]
        for var in variables_to_pass:
            setattr(self, var, locals()[var])

    def initialize_variables(self):
        header_comps = self.header_comps
        header_reacs = self.header_reacs
        comps = self.comps
        reacs = self.reacs

    def load_variables_from_form(self):
        n = self.n
        nr = self.nr
        comps = np.empty([n, 4], dtype='S50')
        reacs = np.empty([nr, n + 2], dtype='S50')
        header_comps = []
        header_reacs = []
        component_order_in_table = []
        for i in range(self.tableComps.rowCount()):
            component_order_in_table.append(
                int(self.tableComps.item(i, 0).text()) - 1)
        for i in range(comps.shape[1]):
            header_comps.append(
                self.tableComps.horizontalHeaderItem(i).data(0))
        for i in range(reacs.shape[1]):
            header_reacs.append(
                self.tableReacs.horizontalHeaderItem(i).data(0))
        for j in range(comps.shape[1]):
            for i in range(comps.shape[0]):
                comps[
                    component_order_in_table[i],
                    j] = self.tableComps.item(
                    i,
                    j).text()
        for j in range(reacs.shape[1]):
            for i in range(reacs.shape[0]):
                reacs[i, j] = self.tableReacs.item(i, j).text()
        # Pass variables to self before loop start
        variables_to_pass = ['header_comps', 'comps', 'header_reacs',
                             'reacs', 'component_order_in_table']
        for var in variables_to_pass:
            setattr(self, var, locals()[var])

    def gui_setup_and_variables(self):
        # Collect variables
        n = self.n
        nr = self.nr
        comps = self.comps
        reacs = self.reacs
        header_comps = self.header_comps
        header_reacs = self.header_reacs
        molar_masses_valid = False
        unset_variables = []
        for col, name in enumerate(header_comps):
            if name in ['comp_id']:
                data_type = str
            elif name in ['index']:
                data_type = int
            else:
                data_type = float
            try:
                column_vector = \
                    np.matrix(
                        [row[col] for row in comps],
                        dtype=data_type).T
            except ValueError as detail:
                unset_variables.append(name)
                column_vector = np.empty([n, 1])
                column_vector[:] = np.nan
                if name in ['index', 'comp_id', 'z']:
                    print detail
                    raise Exception('Input field missing: '
                                    + name)
            # Put values of each column name into self by name
            setattr(self, name, column_vector)
        index = self.index
        comp_id = self.comp_id
        z = self.z
        molar_mass = self.molar_mass
        n0 = self.n0
        w0 = self.w0
        xw0 = self.xw0
        x0 = self.x0
        c0 = self.c0
        m0 = self.m0
        rho_solvent = 0.997
        positive_molar_masses = \
            all(map(lambda x: x > 0, molar_mass))
        can_calculate_n_from_w = \
            all(map(lambda x: x not in unset_variables,
                    ['w0', 'molar_mass']))
        mol_variables_empty = \
            all(map(lambda x: x in unset_variables,
                    ['n0', 'c0']))
        # Gui setup with calculated values
        # First determine concentration variables from available data, and
        # determine the index of the solvent.
        highest_n0_indexes = []
        index_of_solvent = []
        c_solvent_tref = []
        if mol_variables_empty:
            # If there moles cannot be determined, reaction quotients
            # cannot be determined either, request molar mass inputs.
            if not can_calculate_n_from_w or \
                    not positive_molar_masses:
                raise Exception('Need positive molar masses defined')
        elif 'n0' in unset_variables and \
                'c0' not in unset_variables:
            # moles not given, but molarity given:
            # overwrite mole number based on 1L of molarity
            n0 = c0
        elif 'n0' not in unset_variables:
            # moles and molar masses given:
            # use as main generators
            pass
        highest_n0_indexes = np.argpartition(n0.A1, (-1, -2))
        index_of_solvent = highest_n0_indexes[-1]
        # Calculate concentration types based on n0, M
        x0 = n0 / sum(n0)
        if positive_molar_masses:
            mm0 = molar_mass[index_of_solvent]
            n0_mm0 = n0[index_of_solvent] * mm0
            # always calculate weight values & fract.
            w0 = np.multiply(n0, molar_mass)
            xw0 = w0 / sum(w0)
        else:
            # Default molar mass of solvent
            molar_mass[index_of_solvent] = 18.01528
            mm0 = molar_mass[index_of_solvent]
            n0_mm0 = n0[index_of_solvent] * mm0
        # overwrite molality and molarity based on available n0.
        # molarity has units mol/g, leave conversion to mol/kg
        # for the end.
        m0 = n0 / (n0_mm0)
        c0 = m0 * rho_solvent * 1000
        # \frac{x_i}{m_i} = \frac{M_0}{1 + sum_{j\neq0}{m_j M_0}}
        mi_mm0_over_xi = \
            1 + sum([m_j for j, m_j in enumerate(m0 * mm0) if
                     j != index_of_solvent])
        rho0 = (mi_mm0_over_xi) * rho_solvent
        rho0_i = np.multiply(c0, mm0)
        c_solvent_tref = c0[index_of_solvent].item()
        z = np.matrix([row[2] for row in comps], dtype=float).T
        nu_ij = np.matrix([row[2:2 + n] for row in reacs], dtype=int).T
        pka = np.matrix([row[1] for row in reacs], dtype=float).T
        max_it = int(self.spinBox_3.value())
        tol = float(self.doubleSpinBox_5.value())
        # Determine second highest component for plotting possibilities
        if len(n0) > 1:
            index_of_second_highest_n0 = highest_n0_indexes[-2]
        else:
            index_of_second_highest_n0 = highest_n0_indexes[-1]
        n_second_highest_n0_tref = n0[index_of_second_highest_n0].item()
        # Calculated variables to output.
        variables_to_pass = [
            'index_of_solvent', 'molar_mass',
            'n0', 'x0', 'w0', 'xw0', 'c0', 'z', 'm0',
            'rho_solvent', 'rho0', 'c_solvent_tref',
            'nu_ij', 'pka', 'max_it', 'tol',
            'n_second_highest_n0_tref',
            'positive_molar_masses',
            'can_calculate_n_from_w'
        ]
        for var in variables_to_pass:
            setattr(self, var, locals()[var])
        # Setup plotting tools
        self.comboBox.clear()
        self.comboBox_3.clear()
        for item in comps[:, 0:2]:
            self.comboBox.addItem('c0_' + item[0] + ' {' + item[1] + '}')
            self.comboBox_3.addItem(item[1])
        self.comboBox.setCurrentIndex(index_of_second_highest_n0)
        self.doubleSpinBox.setValue(n_second_highest_n0_tref / 10.0 ** 7)
        self.doubleSpinBox_2.setValue(
            n_second_highest_n0_tref * (1 + 20 / 100.0))
        self.comboBox_3.setCurrentIndex(index_of_solvent)
        self.doubleSpinBox_6.setValue(c_solvent_tref)
        self.doubleSpinBox_6.setPrefix('(mol/L)')

    def retabulate(self):
        n = self.n
        nr = self.nr
        header_comps = self.header_comps
        reacs = self.reacs
        xieq = self.xieq

        if hasattr(self, 'component_order_in_table'):
            i = getattr(self, 'component_order_in_table')
        else:
            i = range(0, n)
        j = range(0, len(comp_names_headers))

        self.tableComps.blockSignals(True)
        self.tableReacs.blockSignals(True)
        self.comboBox.blockSignals(True)
        # As usual, problems occurr when sorting is combined with setting QTableWidgetItems.
        # Therefore disable sorting, then set QTableWidgetItems and finally
        # reenable sorting.
        self.tableComps.setSortingEnabled(False)
        self.tableReacs.setSortingEnabled(False)

        for column in j:
            column_name = self.tableComps.model().headerData(
                column, QtCore.Qt.Horizontal)
            var_name = comp_names_headers[column_name]
            column_data = getattr(self, var_name)
            if var_name in ['m0', 'meq']:
                # Present units converted: mol/gsolvent to mol/kgsolvent
                column_data = 1000 * column_data
            elif var_name in ['mlog10m0', 'mlog10meq']:
                # Present units converted: mol/gsolvent to mol/kgsolvent.
                # Add -log10(10^3)
                column_data = column_data - 3
            for row in i:
                if column <= len(comp_names_headers):
                    new_item = \
                        QtGui.QTableWidgetItem(str(column_data[row].item()))
                # sortierbar machen
                if column != 1:  # Comp. i <Str>
                    new_item = NSortableTableWidgetItem(new_item)
                    self.tableComps.setItem(row, column, new_item)
                else:
                    self.tableComps.setItem(row, column, new_item)
                if column not in range(1, len(header_comps)):
                    new_item.setFlags(QtCore.Qt.ItemIsEnabled)

        i = range(0, nr)
        j = range(0, n + 2 + 1)

        for column in j:
            for row in i:
                if column != n + 2:
                    self.tableReacs.setItem(
                        row, column, NSortableTableWidgetItem(str(reacs[row][column])))
                elif column == n + 2:
                    self.tableReacs.setItem(
                        row, column, NSortableTableWidgetItem(str(xieq[row].item())))

        # Widths and heights, re-enable sorting
        self.tableComps.setSortingEnabled(True)
        self.tableComps.horizontalHeader().setResizeMode(
            QtGui.QHeaderView.ResizeToContents)
        self.tableComps.verticalHeader().setResizeMode(
            QtGui.QHeaderView.ResizeToContents)

        self.tableReacs.setSortingEnabled(True)
        self.tableReacs.horizontalHeader().setResizeMode(
            QtGui.QHeaderView.ResizeToContents)
        self.tableReacs.verticalHeader().setResizeMode(
            QtGui.QHeaderView.ResizeToContents)

        self.tableComps.blockSignals(False)
        self.tableReacs.blockSignals(False)
        self.comboBox.blockSignals(False)

    def save_file(self):
        pass

    def solve_intervals(self):
        variables_to_check = [
            'ceq_series',
            'xieq_series',
            'indep_var_series',
            'dep_var_series',
            'index_of_variable']
        for var in variables_to_check:
            if hasattr(self, var):
                delattr(self, var)
        comps = self.comps
        reacs = self.reacs
        c0 = self.c0
        n = self.n
        nr = self.nr
        index_of_variable = self.comboBox.currentIndex()
        indep_var_label = 'c0_{' + comps[index_of_variable, 0] + ', ' + \
                          comps[index_of_variable, 1] + '}/(mol/L)'
        dep_var_labels = \
            ['Ceq_' + '{' + item[0] + ', ' + item[1] + '}/(mol/L)' for item in comps[:, 0:2]] + \
            ['\\xi eq_' +
                '{' + str(item) + '}/(mol/L)' for item in range(1, nr + 1, 1)]
        min_value = self.doubleSpinBox.value()
        max_value = self.doubleSpinBox_2.value()
        n_points = 20
        indep_var_series_single = \
            [min_value + x
             for x in np.arange(n_points + 1) * (max_value - min_value) / n_points]
        c0_variable_comp = c0[index_of_variable]
        mid_index = bisect.bisect(
            indep_var_series_single,
            c0_variable_comp) - 1
        xieq = self.xieq
        ceq = self.ceq
        ceq_series = np.matrix(np.zeros([n_points + 1, len(ceq)]))
        xieq_series = np.matrix(np.zeros([n_points + 1, len(xieq)]))
        dep_var_series = dict(
            zip(dep_var_labels, np.empty(n + nr, dtype=np.ndarray)))
        indep_var_series = dict.fromkeys(
            dep_var_labels, indep_var_series_single)
        # Keep current solution intact for after plotting range
        self.stored_solution_ceq = self.ceq
        self.stored_solution_xieq = self.xieq
        for j in range(mid_index, -1, -1):
            self.c0[index_of_variable] = indep_var_series_single[j]
            self.equilibrate()
            ceq_series[j, :] = self.ceq.T
            xieq_series[j, :] = self.xieq.T
        self.ceq = self.stored_solution_ceq
        self.xieq = self.stored_solution_xieq
        for j in range(mid_index + 1, n_points + 1, +1):
            self.c0[index_of_variable] = indep_var_series_single[j]
            self.equilibrate()
            ceq_series[j, :] = self.ceq.T
            xieq_series[j, :] = self.xieq.T
        for j in range(n):
            dep_var_series[dep_var_labels[j]] = ceq_series[:, j]
        for j in range(nr):
            dep_var_series[dep_var_labels[n + j]] = xieq_series[:, j]
        self.ceq = self.stored_solution_ceq
        self.xieq = self.stored_solution_xieq
        self.Ceq_series = ceq_series
        self.Xieq_series = xieq_series
        self.indep_var_series = indep_var_series
        self.dep_var_series = dep_var_series
        self.dep_var_labels = dep_var_labels
        self.indep_var_label = indep_var_label
        self.index_of_variable = index_of_variable
        self.initiate_plot()

    def initiate_plot(self):
        n = self.n
        nr = self.nr
        dep_var_labels = self.dep_var_labels
        labels_to_plot = [x for x in self.dep_var_labels if x.find('Ceq') >= 0]
        # dict, keys:ceq_labels; bindings: plottedseries
        plotted_series = dict(
            zip(dep_var_labels, np.empty(n + nr, dtype=object)))
        dep_var_labels = self.dep_var_labels
        dep_var_series = self.dep_var_series
        indep_var_series = self.indep_var_series
        indep_var_label = self.indep_var_label
        self.groupBox = QtGui.QGroupBox()
        self.groupBox.plotBox = UiGroupBoxPlot(self.groupBox)
        self.groupBox.plotBox.set_variables(
            plotted_series=plotted_series,
            dep_var_series=dep_var_series,
            dep_var_labels=dep_var_labels,
            dep_var_labels_to_plot=labels_to_plot,
            indep_var_label=indep_var_label,
            indep_var_series=indep_var_series,
            add_path_arrows=False)
        self.groupBox.show()
        self.groupBox.plotBox.plot_intervals(labels_to_plot)

    def recalculate_after_cell_edit(self, row, column):
        self.load_variables_from_form()
        self.gui_equilibrate()

    def gui_equilibrate(self):
        if self.groupBox is not None:
            self.groupBox.hide()  # would need to refresh plotBox or logBox
        self.tableComps.blockSignals(True)
        self.tableReacs.blockSignals(True)
        self.comboBox.blockSignals(True)
        self.gui_setup_and_variables()
        self.equilibrate()
        self.retabulate()

        self.tableComps.blockSignals(False)
        self.tableReacs.blockSignals(False)
        self.comboBox.blockSignals(False)

    def equilibrate(self):
        # Collect variables
        n = self.n
        nr = self.nr
        n0 = self.n0
        w0 = self.w0
        xw0 = self.xw0
        x0 = self.x0
        c0 = self.c0
        m0 = self.m0
        z = self.z
        mm = self.molar_mass
        rho_solvent = self.rho_solvent
        comps = self.comps
        reacs = self.reacs
        nu_ij = self.nu_ij
        pka = self.pka
        max_it = self.max_it
        tol = self.tol
        c_solvent_tref = self.c_solvent_tref
        index_of_solvent = self.index_of_solvent

        # Init. calculations
        kc = np.multiply(
            np.power(10, -pka), np.power(c_solvent_tref, nu_ij[index_of_solvent, :]).T)
        self.kc = kc

        # Setup logging
        if not os.path.exists('./logs'):
            os.mkdir('./logs')
        logging.basicConfig(
            filename='./logs/calculation_results.log',
            level=logging.DEBUG,
            format='%(asctime)s;%(message)s')

        # First estimates for eq. Composition Ceq and Reaction extent Xieq
        if not hasattr(self, 'acceptable_solution'):
            neq_0 = n0
            # replace 0 by 10^-6*smallest value: Smith, Missen 1988 DOI:
            # 10.1002/cjce.5450660409
            neq_0[n0 == 0] = min(n0[n0 != 0].A1) * np.finfo(float).eps
            xieq_0 = np.matrix(np.zeros([nr, 1]))
        else:
            # Use previous solution as initial estimate, if it was valid.
            neq_0 = self.neq_0
            xieq_0 = self.xieq_0

        # Pass variables to self before loop start
        variables_to_pass = ['n0', 'z', 'nu_ij', 'pka',
                             'max_it', 'tol',
                             'neq_0', 'xieq_0']
        for var in variables_to_pass:
            setattr(self, var, locals()[var])

        k = 1
        stop = False
        self.cancelButton.setEnabled(True)
        self.progress_var.setEnabled(True)
        self.remove_canceled_status()
        self.acceptable_solution = False
        self.initialEstimateAttempts = 1
        self.methodLoops = [0, 0]  # loop numbers: [Line search, Newton]
        # Calculate equilibrium composition: Newton method
        # TODO: Implement global homotopy-continuation method
        while not self.acceptable_solution \
                and k < max_it and stop is False \
                and not self.was_canceled():
            neq, xieq = calc_xieq(self)
            k += 1
            # TODO: if progress_var.wasCanceled() == True then stop
            if all(neq >= 0) and not any(np.isnan(neq)):
                self.acceptable_solution = True
            else:
                # Set reactions to random extent and recalculate
                # TODO: scale to concentration sizes
                self.xieq_0 = np.matrix(
                    np.random.normal(0.0, 1.0 / 3.0, nr)).T
                # Set aequilibrium composition to initial value + estimated
                # conversion
                self.neq_0 = n0  # + nu_ij * self.xieq_0
                # replace 0 by 10^-6*smallest value: Smith, Missen 1988 DOI:
                # 10.1002/cjce.5450660409
                self.neq_0[self.neq_0 == 0] = min(
                    n0[n0 != 0].A1) * np.finfo(float).eps
                self.initialEstimateAttempts += 1
                self.methodLoops = [0, 0]

        if not self.acceptable_solution:
            delattr(self, 'acceptable_solution')
            self.label_9.setText(self.label_9.text() + '\n')
        else:
            self.neq_0 = neq
            self.xieq_0 = xieq
            self.label_9.setText(self.label_9.text() +
                                 '\nsum(n0*z) = ' + str((z.T * n0).item()) +
                                 ' \t\t\t sum(Ceq*z) = ' + str((z.T * neq).item()) +
                                 '\nI_0 = ' + str((1 / 2.0 * np.power(z, 2).T * n0).item()) +
                                 '\t\t\t\t I_eq = ' + str((1 / 2.0 * np.power(z, 2).T * neq).item()))

        # Once solved, calculate conc. variables at equilibrium
        mm0 = mm[index_of_solvent]
        n0_mm0 = neq[index_of_solvent] * mm0
        meq = neq / (n0_mm0)
        ceq = meq * rho_solvent * 1000
        # \frac{x_i}{m_i} = \frac{M_0}{1 + sum_{j\neq0}{m_j M_0}}
        mi_mm0_over_xi = \
            1 + sum([m_j for j, m_j in enumerate(meq * mm0) if
                     j != index_of_solvent])
        rho = (mi_mm0_over_xi) * rho_solvent
        rhoeq = np.multiply(ceq, mm)
        weq = np.multiply(neq, mm)
        xweq = weq / sum(weq)
        xeq = neq / sum(neq)
        # TODO: Implement activity coefficients
        gammaeq_ii = np.ones_like(meq)
        gammaeq_iii = np.ones_like(meq)
        aeq = np.multiply(gammaeq_iii, meq)*np.nan
        a0 = np.multiply(gammaeq_iii, m0)
        mlog10xw0 = -np.log10(xw0)
        mlog10x0 = -np.log10(x0)
        mlog10c0 = -np.log10(c0)
        mlog10m0 = -np.log10(m0)
        mlog10a0 = -np.log10(a0)
        mlog10gammaeq_ii = -np.log10(np.ones_like(meq))
        mlog10gammaeq_iii = -np.log10(np.ones_like(meq))
        mlog10xweq = -np.log10(xweq)
        mlog10xeq = -np.log10(xeq)
        mlog10ceq = -np.log10(ceq)
        mlog10meq = -np.log10(meq)
        mlog10aeq = -np.log10(aeq)
        variables_to_pass = [
            'xieq',
            'weq', 'xweq', 'neq', 'xeq',
            'ceq', 'meq',
            'rhoeq',
            'gammaeq_ii', 'gammaeq_iii',
            'aeq',
            'mlog10xw0', 'mlog10x0', 'mlog10c0',
            'mlog10m0', 'mlog10a0',
            'mlog10gammaeq_ii',
            'mlog10gammaeq_iii',
            'mlog10xweq', 'mlog10xeq', 'mlog10ceq',
            'mlog10meq', 'mlog10aeq',
        ]
        for var in variables_to_pass:
            setattr(self, var, locals()[var])
        self.cancelButton.setEnabled(False)
        self.progress_var.setEnabled(False)

    def display_about_info(self):
        row_string = unicode('', 'utf_8')
        self.browser_window = QtGui.QWidget()
        self.browser_window.setWindowIcon(QtGui.QIcon(
            os.path.join(sys.path[0], *['utils', 'icon_batch_16X16.png'])))
        self.browser_window.setWindowTitle('About')
        verticalLayout = QtGui.QVBoxLayout(self.browser_window)
        aboutBox_1 = QtWebKit.QWebView()
        toolbar_1 = QtGui.QToolBar()
        back_icon = QtGui.QIcon(
            os.path.join(
                sys.path[0],
                *['utils', 'glyphicons-217-circle-arrow-left.png']))
        forward_icon = QtGui.QIcon(
            os.path.join(
                sys.path[0],
                *['utils', 'glyphicons-218-circle-arrow-right.png']))
        comboBox_1 = QtGui.QComboBox()
        action_select = toolbar_1.addWidget(comboBox_1)
        action_back = toolbar_1.addAction(back_icon, 'Back')
        action_forward = toolbar_1.addAction(forward_icon, 'Forward')

        verticalLayout.addWidget(toolbar_1)
        verticalLayout.addWidget(aboutBox_1)
        aboutBox_1.setWindowTitle('About')

        adding_table = False
        html_stream = unicode('<!DOCTYPE html>', 'utf_8')
        html_stream += unicode('<html>', 'utf_8')
        html_stream += unicode(
            '<head><title>00 - README</title>' +
            '<meta name="qrichtext" content="1">' +
            '<meta charset="utf-8">' +
            '<style type="text/css">' +
            '\\np,li {white-space: pre-wrap;}\n' +
            '\\np,br {line-height: 10%;}\n' +
            '</style></head>',
            'utf_8')

        html_stream += unicode('<body style=' +
                               '"' +
                               ' font-family:' +
                               "'" +
                               'MS Shell Dlg 2' +
                               "'" +
                               '; font-size:8.25pt; font-weight:400; font-style:normal;' +
                               '"' +
                               '>', 'utf_8')
        string_to_add = unicode('', 'utf_8')
        starting_p = unicode(
            "<p style=' margin-top:0px; margin-bottom:0px; margin-left:0px;" +
            "margin-right:0px; -qt-block-indent:0; text-indent:0px;'>",
            'utf_8')
        ending_p = unicode('</p>', 'utf_8')
        with open('README.md') as readme_file:
            for row in readme_file:
                row_string = unicode(row, 'utf_8')
                if not adding_table and row_string.find('|') != -1:
                    string_to_add = ''.join(
                        ['<table>', '<tr><td>',
                         row_string.replace(
                             '|', '</td><td>').replace('\n', ''),
                         '</td></tr>'])
                    adding_table = True
                elif adding_table and row_string.find('|') == -1:
                    string_to_add = '</table>' + starting_p + row_string + ending_p
                    adding_table = False
                elif adding_table and row_string.find('|') != -1:
                    if doc_hline_re.match(row_string) is not None:
                        string_to_add = ''
                    else:
                        string_to_add = ''.join(['<tr><td>', row_string.replace(
                            '|', '</td><td>').replace('\n', ''), '</td></tr>'])
                elif not adding_table and row_string.find('|') == -1:
                    if doc_hline_re.match(row_string) is not None:
                        string_to_add = '<hr>'
                    else:
                        string_to_add = starting_p + row_string + ending_p
                if len(row_string.replace('\n', '')) == 0:
                    html_stream += string_to_add
                elif matchingHLine.match(row_string):
                    html_stream += '<hr>'
                else:
                    html_stream += string_to_add
            if adding_table:
                html_stream += unicode('</table>', 'utf_8')
        html_stream += unicode('<hr>', 'utf_8')
        html_stream += unicode('<ul>', 'utf_8')
        file_name_list = []
        for file in os.listdir(os.path.abspath('docs')):
            filename_ext = os.path.splitext(os.path.basename(file))
            ext = filename_ext[-1]
            file_name = filename_ext[0]
            if ext == '.html':
                file_path = os.path.join('docs', file)
                with open(file_path) as opened_file:
                    read_file = '\n'.join(opened_file.readlines())
                    file_title = \
                        html_title.search(read_file).groups()[0]
                    opened_file.close()
                    comboBox_1.addItem(file_title, file_path)
                if file_name.lower() != 'readme':
                    file_name_list.append((file_title, file_path))
        comboBox_1.model().sort(0)
        comboBox_1.setCurrentIndex(0)
        for file_name_path in sorted(file_name_list, key=lambda x: x[0]):
            html_stream += unicode(
                '<li><a href=' +
                'file://' +
                urllib.pathname2url(os.path.abspath(file_name_path[1])) +
                '>' +
                file_name_path[0] +
                '</a></li>', 'utf-8')
        html_stream += unicode('</ul>', 'utf_8')
        html_stream += unicode('<hr>', 'utf_8')
        html_stream += unicode(
            "<footer><p>" +
            '<a href=' +
            'https://github.com/santiago-salas-v/lit-impl-py' + '>' +
            'https://github.com/santiago-salas-v/lit-impl-py</a></p></footer>',
            'utf_8')
        html_stream += unicode('</body></html>', 'utf_8')
        readme_file.close()
        html_file_path = os.path.join('docs', 'README.html')
        html_file = open(html_file_path, 'w')
        html_file.write(html_stream.encode('utf-8'))
        html_file.close()
        aboutBox_1.load(html_file_path)
        aboutBox_1.setMinimumWidth(500)
        aboutBox_1.setMinimumHeight(400)
        aboutBox_1.show()
        self.browser_window.show()
        # Make connections at the end
        action_back.triggered.connect(partial(aboutBox_1.back))
        action_forward.triggered.connect(partial(aboutBox_1.forward))
        comboBox_1.currentIndexChanged.connect(
            partial(
                lambda x: aboutBox_1.load(
                    comboBox_1.itemData(x))
            ))

    def show_log(self):
        headers_and_types = np.array(
            (('date', str),
             ('method', str),
             ('k', int),
             ('backtrack', int),
             ('lambda_ls', float),
             ('accum_step', float),
             ('X', list),
             ('||X(k)-X(k-1)||', float),
             ('f(X)', list),
             ('||f(X)||', float),
             ('Y', list),
             ('||Y||', float),
             ('g', float),
             ('|g-g1|', float),
             ('stop', bool),
             ('series_id', str)))

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
        self.groupBox = QtGui.QGroupBox()
        self.groupBox.log_widget = LogWidget(log, parent=self.groupBox)
        self.groupBox.show()


class UiGroupBoxPlot(QtGui.QWidget):

    def __init__(self, parent):
        QtGui.QWidget.__init__(self, parent)
        # Default size
        parent.resize(741, 421)
        # Assignments
        log_scale_func_list = [
            ('-log10', lambda x: -1.0 * np.log10(x)), ('', lambda x: 10.0 ** (-1.0 * x))]
        self.icon_down = QtGui.QIcon(os.path.join(
            sys.path[0], *['utils', 'glyphicons-602-chevron-down.png']))
        self.icon_up = QtGui.QIcon(os.path.join(
            sys.path[0], *['utils', 'glyphicons-601-chevron-up.png']))
        self.verticalLayout_1 = QtGui.QVBoxLayout(parent)
        self.horizontalTools = QtGui.QHBoxLayout()
        self.toolsFrame = QtGui.QFrame()
        self.toggleLogButtonX = QtGui.QPushButton()
        self.toggleLogButtonY = QtGui.QPushButton()
        self.eraseAnnotationsB = QtGui.QPushButton()
        self.navigation_frame = QtGui.QFrame()
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.figure = Figure(dpi=72, facecolor=(1, 1, 1),
                             edgecolor=(0, 0, 0))
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.label = QtGui.QLabel(parent)
        self.listWidget_2 = QtGui.QListWidget(parent)
        self.label_2 = QtGui.QLabel(parent)
        self.listWidget = QtGui.QListWidget(parent)
        self.plotButton = QtGui.QPushButton(parent)
        self.toolbar = NavigationToolbar(self.canvas, self.navigation_frame)
        self.horizontalLayout_1 = QtGui.QHBoxLayout()
        self.all_to_displayed = QtGui.QPushButton()
        self.all_to_available = QtGui.QPushButton()
        # Default log_log_scale_func_list
        self.set_log_scale_func_list(log_scale_func_list)
        # Object names
        self.verticalLayout_1.setObjectName("verticalLayout_1")
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.canvas.setObjectName("canvas")
        self.verticalLayout.setObjectName("verticalLayout")
        self.label.setObjectName("label")
        self.listWidget_2.setObjectName("listWidget_2")
        self.label_2.setObjectName("label_2")
        self.listWidget.setObjectName("listWidget")
        self.plotButton.setObjectName("plotButton")
        # Operations
        self.eraseAnnotationsB.setIcon(QtGui.QIcon(os.path.join(
            sys.path[0], *['utils', 'glyphicons-551-erase.png'])))
        self.figure.subplots_adjust(bottom=0.15)
        self.ax.grid('on')
        self.horizontalLayout.addWidget(self.canvas)
        self.verticalLayout_1.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalTools.setContentsMargins(0, 0, 0, 0)
        self.horizontalTools.setAlignment(QtCore.Qt.AlignVCenter)
        self.toolsFrame.setLayout(self.horizontalTools)
        self.toggleLogButtonX.setCheckable(True)
        self.toggleLogButtonY.setCheckable(True)
        self.verticalLayout_1.addWidget(self.navigation_frame)
        self.horizontalTools.addWidget(self.toggleLogButtonY)
        self.horizontalTools.addWidget(self.toggleLogButtonX)
        self.horizontalTools.addWidget(self.eraseAnnotationsB)
        self.verticalLayout_1.addWidget(self.toolsFrame)
        self.verticalLayout_1.addLayout(self.horizontalLayout)
        self.listWidget_2.setSelectionMode(
            QtGui.QAbstractItemView.MultiSelection)
        self.listWidget.setSelectionMode(
            QtGui.QAbstractItemView.MultiSelection)
        self.listWidget_2.setViewMode(QtGui.QListView.ListMode)
        self.listWidget_2.setSortingEnabled(False)
        self.listWidget.setSortingEnabled(False)
        self.verticalLayout.addWidget(self.label)
        self.verticalLayout.addWidget(self.listWidget_2)
        self.verticalLayout.addLayout(self.horizontalLayout_1)
        self.all_to_displayed.setIcon(self.icon_up)
        self.all_to_available.setIcon(self.icon_down)
        self.horizontalLayout_1.addWidget(self.all_to_available)
        self.horizontalLayout_1.addWidget(self.all_to_displayed)
        self.verticalLayout.addWidget(self.label_2)
        self.verticalLayout.addWidget(self.listWidget)
        self.verticalLayout.addWidget(self.plotButton)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.navigation_frame.setMinimumHeight(self.toolbar.height())
        # Events
        self.listWidget_2.itemDoubleClicked.connect(
            partial(self.move_to_available))
        self.listWidget.itemDoubleClicked.connect(
            partial(self.move_to_displayed))
        # self.all_to_available.clicked.connect(
        #    partial(self.move_to_available(self.listWidget_2.items)))
        self.toggleLogButtonX.toggled.connect(
            partial(self.toggled_toggle_log_button_x))
        self.toggleLogButtonY.toggled.connect(
            partial(self.toggled_toggle_log_button_y))
        self.eraseAnnotationsB.clicked.connect(partial(self.erase_annotations))
        self.plotButton.clicked.connect(partial(self.force_update_plot))
        # Retranslate, connect
        self.retranslate_ui(parent)
        QtCore.QMetaObject.connectSlotsByName(parent)

    def set_variables(
            self,
            plotted_series,
            dep_var_series,
            dep_var_labels,
            indep_var_label,
            indep_var_series,
            dep_var_labels_to_plot=None,
            log_x_checked=True,
            log_y_checked=True,
            log_scale_func_list=None,
            add_path_arrows=False):
        self.toggleLogButtonX.blockSignals(True)
        self.toggleLogButtonY.blockSignals(True)
        self.toggleLogButtonX.setChecked(log_x_checked)
        self.toggleLogButtonY.setChecked(log_y_checked)
        self.toggleLogButtonX.blockSignals(False)
        self.toggleLogButtonY.blockSignals(False)

        # Variables to be set
        self.dc = dict()  # space to store data cursors for each series
        self.plotted_series = plotted_series
        self.dep_var_series = dep_var_series
        self.dep_var_labels = dep_var_labels
        self.indep_var_label = indep_var_label
        self.indep_var_series = indep_var_series
        self.log_x_checked = log_x_checked
        self.log_y_checked = log_y_checked
        self.log_scale_func_list = log_scale_func_list
        self.log_dep_var_series = dict.fromkeys(dep_var_series.keys())
        self.log_indep_var_series = dict.fromkeys(indep_var_series.keys())
        self.add_path_arrows = add_path_arrows
        # Default log scale to -log10, but enable use of other log scales depending on setting of this array.
        # Form:
        # [(log_scale_string,log_scale_func),(invlog_scale_string,invlog_scale_func)]
        if log_scale_func_list is None:
            pass
        else:
            self.set_log_scale_func_list(log_scale_func_list)
        # Variabeln in Log-Skala
        for k in self.log_dep_var_series.iterkeys():
            self.log_dep_var_series[k] = \
                self.log_scale_func(self.dep_var_series[k])
            self.log_indep_var_series[k] = \
                self.log_scale_func(self.indep_var_series[k])
        # Populate lists with displayed / available functions
        for label in dep_var_labels:
            if label in dep_var_labels_to_plot:
                new_item = QtGui.QListWidgetItem(label, self.listWidget_2)
                new_item.setIcon(self.icon_down)
            else:
                new_item = QtGui.QListWidgetItem(label, self.listWidget)
                new_item.setIcon(self.icon_up)
        # Sort them
        self.listWidget_2.sortItems(QtCore.Qt.AscendingOrder)
        self.listWidget.sortItems(QtCore.Qt.DescendingOrder)

    def force_update_plot(self):
        ylabel = self.ax.get_ylabel()
        self.ax.clear()
        self.ax.grid('on')
        self.plot_intervals()
        self.ax.set_ylabel(ylabel)

    def set_log_scale_func_list(self, log_scale_func_list):
        self.log_scale_string = log_scale_func_list[0][0]
        self.log_scale_func = log_scale_func_list[0][1]
        self.invlog_scale_string = log_scale_func_list[1][0]
        self.invlog_scale_func = log_scale_func_list[1][1]
        self.find_log_variable = re.compile(
            '\$?(?P<log>' +
            self.log_scale_string +
            '\()(?P<id>[^\$]*)(\))\$?|\$?(?P<id2>[^\$]*)\$?')
        self.toggleLogButtonX.setText(
            self.log_scale_string + "(x) - horizontal")
        self.toggleLogButtonY.setText(self.log_scale_string + "(y) - vertical")

    def retranslate_ui(self, parent):
        parent.setWindowTitle(
            QtGui.QApplication.translate(
                "parent",
                "Plot",
                None,
                QtGui.QApplication.UnicodeUTF8))
        parent.setTitle(
            QtGui.QApplication.translate(
                "parent",
                "Plot",
                None,
                QtGui.QApplication.UnicodeUTF8))
        self.label.setText(
            QtGui.QApplication.translate(
                "parent",
                "Displayed",
                None,
                QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(
            QtGui.QApplication.translate(
                "parent",
                "Available",
                None,
                QtGui.QApplication.UnicodeUTF8))
        self.plotButton.setText(
            QtGui.QApplication.translate(
                "parent",
                "Plot",
                None,
                QtGui.QApplication.UnicodeUTF8))
        self.toggleLogButtonX.setText(
            QtGui.QApplication.translate(
                "parent",
                self.log_scale_string +
                "(x) - horizontal",
                None,
                QtGui.QApplication.UnicodeUTF8))
        self.toggleLogButtonY.setText(
            QtGui.QApplication.translate(
                "parent",
                self.log_scale_string +
                "(y) - vertical",
                None,
                QtGui.QApplication.UnicodeUTF8))
        self.eraseAnnotationsB.setText(
            QtGui.QApplication.translate(
                "parent",
                "Erase annotations",
                None,
                QtGui.QApplication.UnicodeUTF8))
        self.all_to_available.setText(
            QtGui.QApplication.translate(
                "parent",
                "All",
                None,
                QtGui.QApplication.UnicodeUTF8))
        self.all_to_displayed.setText(
            QtGui.QApplication.translate(
                "parent",
                "All",
                None,
                QtGui.QApplication.UnicodeUTF8))

    def toggled_toggle_log_button_x(self, checked):
        self.delete_arrows()
        match = self.find_log_variable.search(self.ax.get_xlabel())
        for line in self.ax.lines:
            line_xdata = line.get_xdata()
            line_label_match = self.find_log_variable.search(line.get_label())
            line_label = line_label_match.group('id')
            if line_label is None:
                line_label = line_label_match.group('id2')
            if not checked and (
                    match.group('log') is not None):  # just unchecked
                if line_label in self.indep_var_series.keys():
                    line.set_xdata(self.indep_var_series[line_label])
                else:
                    line.set_xdata(self.invlog_scale_func(line_xdata))
                self.ax.set_xlabel(u'$' + match.group('id') + u'$')
            else:  # just checked
                if line_label in self.indep_var_series.keys():
                    line.set_xdata(self.log_indep_var_series[line_label])
                else:
                    # TODO: Check first if any data are nan due to conversion.
                    line.set_xdata(self.log_scale_func(line_xdata))
                self.ax.set_xlabel(
                    '$' + self.log_scale_string + '(' + match.group('id2') + ')' + '$')
        if self.add_path_arrows:
            add_arrow_to_line2d(
                self.ax,
                self.ax.get_lines(),
                arrow_locs=np.linspace(
                    0.,
                    1.,
                    15),
                arrow_style='->',
                arrow_size=2)
        self.ax.legend(
            loc='best',
            fancybox=True,
            borderaxespad=0.,
            framealpha=0.5).draggable(True)
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()

    def toggled_toggle_log_button_y(self, checked):
        self.delete_arrows()
        for line in self.ax.lines:
            line_label = line.get_label()
            line_ydata = line.get_ydata()
            match = self.find_log_variable.search(line_label)
            if not checked and (match.group('log') is not None):
                if match.group('id') in self.dep_var_series.keys():
                    line.set_ydata(
                        self.dep_var_series[
                            match.group('id')].A1.tolist())
                else:
                    line.set_ydata(self.invlog_scale_func(line_ydata))
                line.set_label(u'$' + match.group('id') + u'$')
            else:
                if match.group('id2') in self.dep_var_series.keys():
                    line.set_ydata(
                        self.log_dep_var_series[
                            match.group('id2')].A1.tolist())
                else:
                    # TODO: Check first if any data are nan due to conversion.
                    line.set_ydata(self.invlog_scale_func(line_ydata))
                line.set_label('$' + self.log_scale_string +
                               '(' + match.group('id2') + ')' + '$')
        if self.add_path_arrows:
            add_arrow_to_line2d(
                self.ax,
                self.ax.get_lines(),
                arrow_locs=np.linspace(
                    0.,
                    1.,
                    15),
                arrow_style='->',
                arrow_size=2)
        self.ax.legend(
            loc='best',
            fancybox=True,
            borderaxespad=0.,
            framealpha=0.5).draggable(True)
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()

    def plot_intervals(self, item_texts=None):
        self.erase_annotations()
        self.delete_arrows()
        plotted_series = self.plotted_series
        dc = dict(zip(self.dep_var_labels, np.empty(len(self.dep_var_labels))))
        if item_texts is None:
            item_texts = []
            for i in range(self.listWidget_2.count()):
                item_texts.append(self.listWidget_2.item(i).text())
        # Nur die Linien hinzufügen, die noch nicht aufgezeichnet wurden.
        done_series = []
        series_to_plot = []
        for x in self.ax.get_lines():
            match = self.find_log_variable.search(x.get_label())
            if match is not None and match.group('id') is not None:
                done_series.append(match.group('id'))
            elif match is not None and match.group('id2') is not None:
                done_series.append(match.group('id2'))
        for x in item_texts:
            if x not in done_series:
                series_to_plot.append(x)
        if not self.toggleLogButtonX.isChecked():
            indep_var_label = '$' + self.indep_var_label + '$'
        else:
            indep_var_label = '$' + self.log_scale_string + \
                              '(' + self.indep_var_label + ')$'
        for label in series_to_plot:
            if not self.toggleLogButtonX.isChecked():
                indep_var_values = self.indep_var_series[label]
            else:
                indep_var_values = self.log_indep_var_series[label]
            if not self.toggleLogButtonY.isChecked():
                dep_var_values = self.dep_var_series[label]
                series_label = '$' + label + '$'
            else:
                dep_var_values = self.log_dep_var_series[label]
                series_label = '$' + self.log_scale_string + '(' + label + ')$'
            plotted_series[label] = self.ax.plot(
                indep_var_values, dep_var_values.A1.tolist(), 'go-', label=series_label,
                color=colormap_colors[
                    np.random.randint(
                        0, len(colormap_colors), 1)],
                markerfacecolor=colormap_colors[
                    np.random.randint(0, len(colormap_colors), 1)],
                marker=markers[np.random.randint(0, len(markers) - 1)],
                fillstyle=fillstyles[np.random.randint(0, len(fillstyles) - 1)])
            dc[label] = datacursor(
                plotted_series[label], draggable=True, display='multiple', arrowprops=dict(
                    arrowstyle='simple', fc='white', alpha=0.5), bbox=dict(
                    fc='white', alpha=0.5), formatter='x: {x:0.3g},y: {y:0.3g}\n{label}'.format)
        if self.add_path_arrows:
            add_arrow_to_line2d(
                self.ax,
                self.ax.get_lines(),
                arrow_locs=np.linspace(
                    0.,
                    1.,
                    15),
                arrow_style='->',
                arrow_size=2)
        self.ax.legend(
            loc='best',
            fancybox=True,
            borderaxespad=0.,
            framealpha=0.5).draggable(True)
        self.plotted_series = plotted_series
        self.dc = dc
        self.ax.set_xlabel(indep_var_label, fontsize=14)
        self.listWidget_2.setMinimumWidth(
            self.listWidget_2.sizeHintForColumn(0))
        self.canvas.draw()

    def erase_annotations(self, text_list=None):
        if text_list is None:
            text_list = [x.get_text() for x in self.figure.texts]
        for text in text_list:
            indexes_of_text = np.where(
                [x.get_text().find(text) >= 0 for x in self.figure.texts])[0]
            if len(indexes_of_text) > 1:
                l = self.figure.texts.pop(indexes_of_text[0].item())
                del l
            else:
                l = self.figure.texts.pop(indexes_of_text.item())
                del l
        self.canvas.draw()

    def move_to_available(self, item):
        self.delete_arrows()
        if self.listWidget_2.count() <= 1:
            return
        selected_items = self.listWidget_2.selectedItems()
        for selected_item in selected_items:
            name = selected_item.text()
            new_item = self.listWidget_2.takeItem(
                self.listWidget_2.indexFromItem(selected_item).row())
            new_item.setIcon(QtGui.QIcon(os.path.join(
                sys.path[0], *['utils', 'glyphicons-601-chevron-up.png'])))
            self.listWidget.insertItem(self.listWidget.count(), new_item)
            associated_annotations = [
                x.get_text() for x in self.figure.texts if x.get_text().find(name) >= 0]
            self.erase_annotations(associated_annotations)
            l = self.ax.lines.pop(np.where(
                [x.properties()['label'].find(name) >= 0 for x in self.ax.lines])[0].item())
            del l
        if len(self.ax.lines) > 0:
            self.ax.legend(
                loc='best',
                fancybox=True,
                borderaxespad=0.,
                framealpha=0.5).draggable(True)
        else:
            del self.ax.legend
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()
        self.listWidget_2.sortItems(QtCore.Qt.AscendingOrder)
        self.listWidget.sortItems(QtCore.Qt.DescendingOrder)

    def move_to_displayed(self, item):
        self.delete_arrows()
        selected_items = self.listWidget.selectedItems()
        selected_items_names = [x.text() for x in selected_items]
        for selected_item in selected_items:
            name = item.text()
            new_item = self.listWidget.takeItem(
                self.listWidget.indexFromItem(selected_item).row())
            new_item.setIcon(QtGui.QIcon(os.path.join(
                sys.path[0], *['utils', 'glyphicons-602-chevron-down.png'])))
            self.listWidget_2.insertItem(self.listWidget_2.count(), new_item)
        self.plot_intervals(selected_items_names)
        if len(self.ax.lines) > 0:
            self.ax.legend(
                loc='best',
                fancybox=True,
                borderaxespad=0.,
                framealpha=0.5).draggable(True)
        else:
            del self.ax.legend
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()
        self.listWidget_2.sortItems(QtCore.Qt.AscendingOrder)
        self.listWidget.sortItems(QtCore.Qt.DescendingOrder)

    def clear_all(self):
        self.ax.clear()
        self.listWidget.clear()
        self.listWidget_2.clear()

    def delete_arrows(self):
        while self.ax.patches:
            l = self.ax.patches.pop(0)
            del l


class LogWidget(QtGui.QWidget):

    def __init__(self, _log, parent):
        QtGui.QWidget.__init__(self, parent)
        self.log = _log
        self.icon = QtGui.QIcon()
        self.icon.addPixmap(QtGui.QPixmap(
            _fromutf8("utils/glyphicons-88-log-book.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        parent.setWindowIcon(self.icon)
        self.setup_ui(parent)

    def setup_ui(self, parent):
        # Default size
        self.minLogHeight = 400
        self.minLogWidth = self.minLogHeight * 16 / 9
        parent.resize(self.minLogWidth, self.minLogHeight)
        # Assignments
        self.verticalLayout = QtGui.QVBoxLayout(parent)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.pandasView = QtGui.QTableView()
        self.firstButton = QtGui.QPushButton()
        self.lastButton = QtGui.QPushButton()
        self.nextButton = QtGui.QPushButton()
        self.previousButton = QtGui.QPushButton()
        self.pageLabel = QtGui.QLabel()
        self.totPagesLabel = QtGui.QLabel()
        self.page_box = QtGui.QLineEdit()
        self.exportButton = QtGui.QPushButton()
        self.plotButton = QtGui.QPushButton()
        self.display_items_by_page = 50
        self.current_page_first_entry = len(self.log.values) \
            - self.display_items_by_page
        # Operations
        parent.setWindowTitle('calculation log')
        self.firstButton.setText('<< First')
        self.lastButton.setText('Last >>')
        self.previousButton.setText('< Previous')
        self.nextButton.setText('Next >')
        self.page_box.setMaximumWidth(int(round(self.minLogHeight / float(5))))
        self.page_box.setAlignment(
            QtCore.Qt.AlignHCenter | QtCore.Qt.AlignCenter)
        self.exportButton.setText('Export (csv)')
        self.plotButton.setText('Plot solution paths')
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout.addWidget(self.pandasView)
        self.verticalLayout.addWidget(self.exportButton)
        self.verticalLayout.addWidget(self.plotButton)
        self.horizontalLayout.addWidget(self.firstButton)
        self.horizontalLayout.addWidget(self.previousButton)
        self.horizontalLayout.addWidget(self.pageLabel)
        self.horizontalLayout.addWidget(self.page_box)
        self.horizontalLayout.addWidget(self.totPagesLabel)
        self.horizontalLayout.addWidget(self.nextButton)
        self.horizontalLayout.addWidget(self.lastButton)

        # To ensure full display, first set resize modes, then resize columns
        # to contents
        self.pandasView.horizontalHeader().setResizeMode(QtGui.QHeaderView.Interactive)
        self.pandasView.verticalHeader().setResizeMode(
            QtGui.QHeaderView.ResizeToContents)
        self.pandasView.setWordWrap(True)
        self.pandasView.resizeColumnsToContents()

        # Events
        self.firstButton.clicked.connect(partial(self.first_page))
        self.lastButton.clicked.connect(partial(self.last_page))
        self.nextButton.clicked.connect(partial(self.next_page))
        self.previousButton.clicked.connect(partial(self.previous_page))
        self.exportButton.clicked.connect(partial(self.export_data))
        self.plotButton.clicked.connect(partial(self.plot_data))
        self.page_box.editingFinished.connect(partial(self.go_to_page_no))
        # Display
        self.display()

    def first_page(self):
        self.current_page_first_entry = 0
        self.display_items_by_page = self.display_items_by_page
        self.display()

    def last_page(self):
        self.current_page_first_entry = len(self.log.values) \
            - self.display_items_by_page
        self.display()

    def next_page(self):
        if self.currentPageLastEntry + \
                self.display_items_by_page > len(self.log.values):
            self.current_page_first_entry = len(self.log.values) \
                - self.display_items_by_page
        else:
            self.current_page_first_entry = self.current_page_first_entry + \
                self.display_items_by_page
        self.display()

    def previous_page(self):
        if self.current_page_first_entry - \
                self.display_items_by_page < 0:
            self.current_page_first_entry = 0
        else:
            self.current_page_first_entry = self.current_page_first_entry - \
                self.display_items_by_page
        self.display()

    def go_to_page_no(self):
        page_text = self.page_box.text()
        try:
            page_no = int(round(float(page_text)))
            self.last_page = int(round(len(self.log.values) /
                                       float(self.display_items_by_page)))
            if page_no < self.last_page:
                self.current_page_first_entry = \
                    (page_no - 1) * self.display_items_by_page + 1
            elif page_no >= self.last_page:
                self.current_page_first_entry = len(self.log.values) \
                    - self.display_items_by_page
            elif page_no <= 1:
                self.current_page_first_entry = 1
            self.display()
        except ValueError:
            pass

    def display(self):
        # Page number
        self.currentPageLastEntry = self.current_page_first_entry + \
            self.display_items_by_page
        self.currentPage = int(round(self.currentPageLastEntry /
                                     float(self.display_items_by_page)))
        self.last_page = int(round(len(self.log.values) /
                                   float(self.display_items_by_page)))

        self.pandas_model = PandasModel(
            self.log[self.current_page_first_entry: self.currentPageLastEntry])
        self.pandasView.setModel(self.pandas_model)
        self.pandasView.resizeColumnsToContents()
        self.pandasView.verticalHeaders = \
            map(str, range(self.current_page_first_entry,
                           self.currentPageLastEntry + 1, 1))
        self.totPagesLabel.setText(' / ' + str(self.last_page))
        self.pageLabel.setText('Entries ' +
                               str(self.current_page_first_entry) +
                               ' to ' +
                               str(self.currentPageLastEntry) +
                               '; Page: ')
        self.page_box.blockSignals(True)
        self.page_box.setText(str(self.currentPage))
        self.page_box.blockSignals(False)

    def export_data(self):
        supported_filters = ['CSV file (*.csv)']
        # TODO: supported_filters = ['CSV file (*.csv)', 'XLSX (2010) (*.xlsx)',
        # 'XLS (2007) (*.xls)']
        (file_name, selected_filter) = QtGui.QFileDialog.getSaveFileName(
            parent=self.group_2,
            caption='enter file name to save...',
            filter=';;'.join(supported_filters))
        if selected_filter == supported_filters[0]:
            self.log.to_csv(file_name)
            # elif selected_filter == supported_filters[1] or \
            #                selected_filter == supported_filters[2]:
            #    self.log.to_excel(file_name)

    def plot_data(self):
        grouped = self.log.groupby('series_id')
        # dict, keys:ceq_labels; bindings: plottedseries
        dep_var_labels = grouped.head(1)['date'].apply(lambda x: str(x)).values
        indep_var_label = 'accum step'
        dep_var_series = dict(
            zip(dep_var_labels, np.empty(len(dep_var_labels))))
        indep_var_series = dict(
            zip(dep_var_labels, np.empty(len(dep_var_labels))))
        plotted_series = dict(
            zip(dep_var_labels, np.empty(len(dep_var_labels))))
        log_scale_func_list = [
            ('log10', lambda x: +1.0 * np.log10(x)), ('', lambda x: 10.0 ** (+1.0 * x))]

        for name, group in grouped:
            index = group['date'].head(1).apply(lambda x: str(x)).values.item()
            indep_var_series[index] = group['accum_step'].values
            dep_var_series[index] = np.matrix(group['||f(X)||'].values)
        # Generate the plot
        self.group_3 = QtGui.QGroupBox()
        self.plotBox = UiGroupBoxPlot(self.group_3)
        self.plotBox.set_variables(plotted_series=plotted_series,
                                   dep_var_series=dep_var_series,
                                   dep_var_labels=dep_var_labels,
                                   dep_var_labels_to_plot=[dep_var_labels[-1]],
                                   indep_var_label=indep_var_label,
                                   indep_var_series=indep_var_series,
                                   log_x_checked=False, log_y_checked=False,
                                   log_scale_func_list=log_scale_func_list,
                                   add_path_arrows=True)
        self.plotBox.plot_intervals([dep_var_labels[-1]])
        # FIXME: 2 Punkte, Logx ein- und ausschalten
        self.group_3.show()
        self.plotBox.ax.set_ylabel('||f(X)||')


def calc_xieq(form):
    """Newton method for non-linear algebraic system, with line-search
    :return: tuple with neq, xieq, f_0
    :param n0: np.matrix (n x 1) - Mol(i, alimentación)
    :param z: np.matrix (n x 1) - Carga(i, alimentación)
    :param nu_ij: np.matrix (n x nr) - Coefs. esteq. componente i en reacción j
    :param pka: np.matrix (n x 1) - (-1)*log10("Cte." de equilibrio en reacción j) = -log10 kc(T)
    :param xieq_0: np.matrix (n x 1) - avance de reacción j - estimado inicial
    :param neq_0: np.matrix (n x 1) - Mol(i, equilibrio)
    """
    n = form.n
    nr = form.nr
    n0 = form.n0
    mm = form.molar_mass
    s_index = form.index_of_solvent
    pka = form.pka
    nu_ij = form.nu_ij
    neq_0 = form.neq_0
    xieq_0 = form.xieq_0
    max_it = form.max_it
    tol = form.tol
    z = form.z
    kc = form.kc

    f = lambda x: f_gl_0(x, n0, nu_ij, n, nr, kc, mm, s_index)
    j = lambda x: jac(x, n0, nu_ij, n, nr, kc, mm, s_index)

    x0 = np.concatenate([neq_0, xieq_0])
    x = x0
    # Newton method: G(x) = J(x)^-1 * F(x)
    k = 0
    j_val = j(x)
    f_val = f(x)
    y = np.matrix(np.ones(len(x))).T * tol / (np.sqrt(len(x)) * tol)
    magnitude_f = np.sqrt((f_val.T * f_val).item())
    # Line search variable lambda
    lambda_ls = 0.0
    accum_step = 0.0
    # For progress bar, use log scale to compensate for quadratic convergence
    log10_to_o_max_magnitude_f = np.log10(tol / magnitude_f)
    progress_k = (1.0 - np.log10(tol / magnitude_f) /
                  log10_to_o_max_magnitude_f) * 100.0
    diff = np.matrix(np.empty([len(x), 1]))
    diff.fill(np.nan)
    stop = False
    divergent = False
    # Add progress bar & variable
    form.progress_var.setValue(0)
    update_status_label(form, k, 'solving...' if not stop else 'solved.')
    # Unique identifier for plotting logged solutions
    series_id = str(uuid.uuid1())
    new_log_entry(
        'Newton',
        k,
        0,
        0,
        accum_step,
        x,
        diff,
        f_val,
        0 * y,
        np.nan,
        np.nan,
        stop,
        series_id)
    while k <= max_it and not stop:
        k += 1
        form.methodLoops[1] += 1
        j_it = 0
        lambda_ls = 1.0
        accum_step += lambda_ls
        x_k_m_1 = x
        progress_k_m_1 = progress_k
        y = gausselimination(j_val, -f_val)
        # First attempt without backtracking
        x = x + lambda_ls * y  # FIXME: Autoassignment throws 1001 iterations
        diff = x - x_k_m_1
        j_val = j(x)
        f_val = f(x)
        magnitude_f = np.sqrt((f_val.T * f_val).item())
        new_log_entry(
            'Newton',
            k,
            j_it,
            lambda_ls,
            accum_step,
            x,
            diff,
            f_val,
            lambda_ls * y,
            np.nan,
            np.nan,
            stop,
            series_id)
        if magnitude_f < tol and all(x[0:n] >= 0):
            stop = True  # Procedure successful
            form.progress_var.setValue(100.0)
        else:
            # For progress bar, use log scale to compensate for quadratic
            # convergence
            update_status_label(
                form, k, 'solving...' if not stop else 'solved.')
            progress_k = (1.0 - np.log10(tol / magnitude_f) /
                          log10_to_o_max_magnitude_f) * 100.0
            if np.isnan(magnitude_f) or np.isinf(magnitude_f):
                stop = True  # Divergent method
                divergent = True
                form.progress_var.setValue(
                    (1.0 -
                     np.log10(
                         np.finfo(float).eps) /
                     log10_to_o_max_magnitude_f) *
                    100.0)
            else:
                form.progress_var.setValue(
                    (1.0 -
                     np.log10(
                         tol /
                         magnitude_f) /
                     log10_to_o_max_magnitude_f) *
                    100.0)
            if round(progress_k) == round(progress_k_m_1):
                QtGui.QApplication.processEvents()
                # if form.progress_var.wasCanceled():
                # stop = True
        while j_it <= max_it and not all(x[0:n] >= 0):
            # Backtrack if any conc < 0. Line search method.
            # Ref. http://dx.doi.org/10.1016/j.compchemeng.2013.06.013
            j_it += 1
            lambda_ls = lambda_ls / 2.0
            accum_step += -lambda_ls
            x = x_k_m_1
            progress_k = progress_k_m_1
            x = x + lambda_ls * y
            diff = x - x_k_m_1
            j_val = j(x)
            f_val = f(x)
            new_log_entry(
                'Newton',
                k,
                j_it,
                lambda_ls,
                accum_step,
                x,
                diff,
                f_val,
                lambda_ls * y,
                np.nan,
                np.nan,
                stop,
                series_id)
            form.methodLoops[0] += 1
            update_status_label(
                form, k, 'solving...' if not stop else 'solved.')
    update_status_label(
        form,
        k,
        'solved.' if stop and not divergent else 'solution not found.')
    neq = x[0:n]
    xieq = x[n:n + nr]
    return neq, xieq


def f_gl_0(x, n0, nu_ij, n, nr, kc, mm, s_index):
    neq = x[0:n, 0]
    n0_mm0 = n0[s_index] * mm[s_index]
    meq = neq / n0_mm0
    xieq = x[n:n + nr, 0]
    result = np.matrix(np.empty([n + nr, 1], dtype=float))
    result[0:n] = -neq + n0 + nu_ij * xieq
    result[n:n + nr] = -kc + np.prod(np.power(meq, nu_ij), 0).T
    return result


def jac(x, n0, nu_ij, n, nr, kc, mm, s_index):
    neq = x[0:n, 0]
    n0_mm0 = n0[s_index] * mm[s_index]
    meq = neq / n0_mm0
    eins_durch_m = np.diag(np.power(meq, -1).A1, 0)
    quotient = np.diag(np.prod(np.power(meq, nu_ij), 0).A1)
    result = np.matrix(np.zeros([n + nr, n + nr], dtype=float))
    result[0:n, 0:n] = -1 * np.eye(n).astype(float)
    result[0:n, n:n + nr] = nu_ij
    result[n:n + nr, 0:n] = quotient * nu_ij.T * eins_durch_m * 1 / n0_mm0
    return result


def update_status_label(form, k, solved):
    form.label_9.setText('Loops: Newton \t' +
                         str(form.methodLoops[1]) +
                         ' \t Line search (backtrack) \t' +
                         str(form.methodLoops[0]) +
                         ' \t Initial estimate attempts \t' +
                         str(form.initialEstimateAttempts) +
                         '\n' +
                         'Iteration (k) \t' +
                         str(k) +
                         '\n' +
                         str(solved))


def new_log_entry(
        method,
        k,
        backtrack,
        lambda_ls,
        accum_step,
        x,
        diff,
        f_val,
        y,
        g_min,
        g1,
        stop,
        series_id):
    logging.debug(method + ' ' +
                  ';k=' + str(k) +
                  ';backtrack=' + str(backtrack) +
                  ';lambda_ls=' + str(lambda_ls) +
                  ';accum_step=' + str(accum_step) +
                  ';X=' + '[' + ','.join(map(str, x.T.A1)) + ']' +
                  ';||X(k)-X(k-1)||=' + str((diff.T * diff).item()) +
                  ';f(X)=' + '[' + ','.join(map(str, f_val.T.A1)) + ']' +
                  ';||f(X)||=' + str(np.sqrt((f_val.T * f_val).item())) +
                  ';Y=' + '[' + ','.join(map(str, y.T.A1)) + ']' +
                  ';||Y||=' + str(np.sqrt((y.T * y).item())) +
                  ';g=' + str(g_min) +
                  ';|g-g1|=' + str(abs(g_min - g1)) +
                  ';stop=' + str(stop) +
                  ';' + series_id)


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
        # TODO: Implement vertical header indicating row numbers
        return None


class AboutBox(QtGui.QMessageBox):

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


class ScientificDoubleSpinBox(QtGui.QDoubleSpinBox):

    def __init__(self, parent=None):
        QtGui.QDoubleSpinBox.__init__(self, parent)
        self.setMinimum(-np.inf)
        self.setMaximum(np.inf)
        self.validator = FloatValidator()
        self.setDecimals(1000)

    def validate(self, text, position):
        return self.validator.validate(text, position)

    def fixup(self, text):
        return self.validator.fixup(text)

    def valueFromText(self, text):
        return float(text)

    def textFromValue(self, value):
        return format_float(value)

    def stepBy(self, steps):
        text = self.cleanText()
        groups = float_re.search(text).groups()
        decimal = float(groups[1])
        decimal += steps
        new_string = "{:g}".format(decimal) + (groups[3] if groups[3] else "")
        self.lineEdit().setText(new_string)


class FloatValidator(QtGui.QValidator):

    def validate(self, string, position):
        if valid_float_string(string):
            return self.State.Acceptable
        if string == "" or string[position - 1] in 'e.-+':
            return self.State.Intermediate
        return self.State.Invalid

    def fixup(self, text):
        match = float_re.search(text)
        return match.groups()[0] if match else ""


def add_arrow_to_line2d(
        axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8],
        arrow_style='-|>', arrow_size=1, transform=None):
    """
    Add arrows to a matplotlib.lines.Line2D at selected locations.

    Parameters:
    -----------
    axes:
    line: list of 1 Line2D obbject as returned by plot command
    arrow_locs: list of locations where to insert arrows, % of total length
    arrowstyle: style of the arrow
    arrowsize: size of the arrow
    transform: a matplotlib transform instance, default to data coordinates

    Returns:
    --------
    arrows: list of arrows
    """
    if (not (isinstance(line, list)) or not (
            isinstance(line[0], matplotlib.lines.Line2D))):
        raise ValueError("expected a matplotlib.lines.Line2D object")
    x, y = line[0].get_xdata(), line[0].get_ydata()
    finite_indexes = np.bitwise_and(np.isfinite(x), np.isfinite(y))
    x = np.array(x)[finite_indexes]
    y = np.array(y)[finite_indexes]

    if sum(finite_indexes) < 2:
        return

    arrow_kw = dict(arrowstyle=arrow_style, mutation_scale=10 * arrow_size)

    color = line[0].get_color()
    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        raise NotImplementedError("multicolor lines not supported")
    else:
        arrow_kw['color'] = color

    linewidth = line[0].get_linewidth()
    if isinstance(linewidth, np.ndarray):
        raise NotImplementedError("multiwidth lines not supported")
    else:
        arrow_kw['linewidth'] = linewidth

    if transform is None:
        axes.relim()
        axes.autoscale_view()
        transform = axes.transData

    arrows = []
    for loc in arrow_locs:
        s = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
        n = np.searchsorted(s, s[-1] * loc)
        arrow_tail = (x[n], y[n])
        arrow_head = (np.mean(x[n:n + 2]), np.mean(y[n:n + 2]))
        p = matplotlib.patches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)
    return arrows


def format_float(value):
    """Modified form of the 'g' format specifier."""
    string = "{:g}".format(value).replace("e+", "e")
    string = re.sub("e(-?)0*(\d+)", r"e\1\2", string)
    return string


def valid_float_string(string):
    match = float_re.search(string)
    return match.groups()[0] == string if match else False


def main():
    app = QtGui.QApplication.instance()  # checks if QApplication already exists
    if not app:  # create QApplication if it doesnt exist
        app = QtGui.QApplication(sys.argv)
    # following 2 lines for setting app icon correctly
    myappid = u'mycompany.myproduct.subproduct.version'  # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    # ended lines for setting app icon correctly
    main_form = QtGui.QGroupBox()
    icon = QtGui.QIcon(
        os.path.join(sys.path[0], *['utils', 'icon_batch_32X32.png']))
    main_form.setWindowIcon(icon)
    main_form.ui = UiGroupBox(main_form)
    main_form.show()
    main_form.ui.load_csv('./DATA/COMPONENTS_REACTIONS_EX_001.csv')
    main_form.ui.initialize_variables()
    main_form.ui.gui_equilibrate()
    main_form.ui.tableComps.sortByColumn(0, QtCore.Qt.AscendingOrder)
    main_form.ui.tableReacs.sortByColumn(0, QtCore.Qt.AscendingOrder)
    app.exec_()


if __name__ == '__main__':
    main()
