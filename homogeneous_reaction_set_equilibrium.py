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
from reaction_equilibrium import calc_xieq
from datetime import datetime

#matplotlib.use('Qt4Agg')
#matplotlib.rcParams['backend.qt5'] = 'PySide'
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from qtpy import QtGui, QtCore, QtWebEngineWidgets, QtWidgets
from mplcursors import cursor as datacursor

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromutf8(s):
        return s

try:
    _encoding = QtWidgets.QApplication.UnicodeUTF8

    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig)


def take_float(x):
    return float(x.rpartition('=')[-1])


def take_list(x):
    separator = ','
    raw_list = x.rpartition('=')[-1].replace('[', '').replace(']', '')
    output_string = np.fromstring(raw_list, dtype=float, sep=separator)
    if x.find('j') > -1:
        square_dim = int(round(np.sqrt(len(output_string))))
        output_string = output_string.reshape(square_dim, square_dim)
    return output_string


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
    'm0/(mol/kg_{solvent})',
    '-log10(xw0)', '-log10(x0)', '-log10(c0)',
    '-log10(m0)', '-log10(a0)'
]
comp_variable_input_names = [
    'index', 'comp_id', 'z', 'molar_mass',
    'w0', 'xw0', 'n0', 'x0', 'c0', 'm0',
    'mlog10xw0', 'mlog10x0', 'mlog10c0',
    'mlog10m0', 'mlog10a0'
]
header_comps_output_model = [
    'weq/g', 'xweq', 'neq/mol', 'xeq',
    'ceq/(mol/L)', 'meq/(mol/kg_{solvent})',
    'rhoeq/(g/mL)',
    '\gamma_{eq}^{II}', '\gamma_{eq}^{III}',
    'aeq',
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
    'mlog10gammaeq_ii',
    'mlog10gammaeq_iii',
    'mlog10xweq', 'mlog10xeq', 'mlog10ceq',
    'mlog10meq', 'mlog10aeq',
]
# First two columns of the header of reacions will match this
# expression, need to find indexes in file and append coefficient
# matrix.
header_reacs_model = ['j', 'pKa']
# Store programatic name vs. header name in a dict.
comp_names_headers = \
    dict(zip(header_comps_input_model + header_comps_output_model,
             comp_variable_input_names + comp_variable_output_names))
# Regular expression compilations
float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')
matchingHLine = re.compile('=+')
# Default component and reaction headers according to structure
# [j, pKa, nu_ij(i=1), nu_ij(i=2),...]
# Sample resulting regex groups when applied to text values:
# 'j'           ['j', None,  None, None, None]
# 'pKaj'        [None, 'pKaj', None, None, None]
# 'nu2000j'     [None, None, 2000, None, None]
# 'nu2j(i=99)'  [None, None, 2, (i=99), 99]
reac_headers_re = re.compile(
    r'(\bj$)' +
    r'|(^pKaj$)' +
    r'|nu_?i?([0-9]+)?j(\(i=([0-9]+)\))?')
comp_headers_re = re.compile(
    r'(\bi$)|(Comp\.?i?)|([z|Z]_?i?)|(M_?i?)|(w_?0_?i?)|(x_?w_?0?)|' +
    r'(n_?0?_?i?)|(x_?0)|([c|C]_?0_?i?)|(m_?0)')
doc_hline_re = re.compile(
    r'(\s*-{3,})')
html_title = re.compile('<title>(.*?)</title>',
                        re.IGNORECASE | re.DOTALL)
main_window_title = 'Homogeneous EC.'


class UiGroupBox(QtWidgets.QWidget):
    _was_canceled = False

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        # Assignments
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(parent)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.open_button = QtWidgets.QPushButton()
        self.save_button = QtWidgets.QPushButton()
        self.info_button = QtWidgets.QPushButton()
        self.log_button = QtWidgets.QPushButton()
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.label_3 = QtWidgets.QLabel()
        self.equilibrate_button = QtWidgets.QPushButton()
        self.spinBox_3 = QtWidgets.QSpinBox()
        self.label_4 = QtWidgets.QLabel()
        self.label = QtWidgets.QLabel()
        self.spinBox = QtWidgets.QSpinBox()
        self.tableComps = QtWidgets.QTableView()
        self.label_9 = QtWidgets.QLabel()
        self.progress_var = QtWidgets.QProgressBar(parent)
        self.cancelButton = QtWidgets.QPushButton(parent)
        self.doubleSpinBox_5 = ScientificDoubleSpinBox()
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.label_5 = QtWidgets.QLabel()
        self.comboBox_3 = QtWidgets.QComboBox()
        self.label_6 = QtWidgets.QLabel()
        self.doubleSpinBox_6 = QtWidgets.QDoubleSpinBox()
        self.label_2 = QtWidgets.QLabel()
        self.spinBox_2 = QtWidgets.QSpinBox()
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.tableReacs = QtWidgets.QTableView()
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.label_7 = QtWidgets.QLabel()
        self.comboBox = QtWidgets.QComboBox()
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.doubleSpinBox = ScientificDoubleSpinBox()
        self.doubleSpinBox_2 = ScientificDoubleSpinBox()
        self.plotButton = QtWidgets.QPushButton()
        self.radio_group = QtWidgets.QHBoxLayout()
        self.radio_b_1 = QtWidgets.QRadioButton()
        self.radio_b_2 = QtWidgets.QRadioButton()
        self.radio_b_3 = QtWidgets.QRadioButton()
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
        #item.setTextAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.tableComps.horizontalHeader().setCascadingSectionResizes(False)
        self.tableComps.horizontalHeader().setDefaultSectionSize(100)
        self.tableComps.horizontalHeader().setMinimumSectionSize(27)
        self.tableComps.horizontalHeader().setSortIndicatorShown(True)
        self.tableComps.verticalHeader().setVisible(False)
        self.verticalLayout_2.addWidget(self.tableComps)
        self.label_9.setAlignment(QtCore.Qt.AlignTop)
        self.label_9.setFrameStyle(QtWidgets.QFrame.Box | QtWidgets.QFrame.Raised)
        self.verticalLayout_2.addWidget(self.label_9)
        self.horizontalLayout_7.setAlignment(QtCore.Qt.AlignRight)
        self.horizontalLayout_7.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
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
        parent.setWindowTitle(_translate("parent", main_window_title, None))
        parent.setTitle(QtWidgets.QApplication.translate(
            "parent", main_window_title[-3:], None))
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
            QtWidgets.QFileDialog.getOpenFileName(None,
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
            if hasattr(self, 'reaction_order_in_table'):
                delattr(self, 'reaction_order_in_table')
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
                        matches = \
                            comp_headers_re.match(column.replace(' ', ''))
                        if matches is not None:
                            new_index = [j for j, match
                                         in enumerate(matches.groups())
                                         if match is not None][0]
                            io_column_index_map.append(
                                [new_index, old_index]
                            )
                elif 'REAC' in row:
                    reading_reacs = True
                    reading_comps = False
                    header_reacs = next(reader)
                    # Get column mappings of form
                    # [
                    # [output_index_1, input_index_1],
                    # [output_index_2, input_index_2],
                    # ...]
                    io_column_index_map = []
                    for (col_no, column) in enumerate(header_reacs):
                        old_index = col_no
                        matches = \
                            reac_headers_re.match(column.replace(' ', ''))
                        if matches is not None:
                            re_index = [(j, match) for j, match
                                        in enumerate(matches.groups())
                                        if match is not None]
                            # When index is
                            # 0: 'j', 1: 'pKa'
                            # 2: X , 3: '(i=Y)',
                            # 4: Y
                            if len(re_index) == 1:
                                first_index = re_index[0][0]
                                if first_index in [0, 1]:
                                    # 0: 'j', 1: 'pKa'
                                    new_index = first_index
                                else:
                                    # In 'nu_Xj'
                                    # 2: X
                                    # add 0 + 1 for j, pKa
                                    new_index = \
                                        int(re_index[0][1]) + 0 + 1
                            elif len(re_index) > 2:  # In 'nu_iX(i=Y)'
                                                    # 2: X , 3: '(i=Y)',
                                                    # 4: Y
                                # Prioritize 4
                                # add 0 + 1 for j, pKa
                                new_index = \
                                    int(re_index[-1][1]) + 0 + 1
                            elif len(re_index) == 2:  # In 'nu_iX(i=Y)'
                                                    # 3: '(i=Y)',
                                                    # 4: Y
                                # Prioritize 4
                                # add 0 + 1 for j, pKa
                                new_index = \
                                    int(re_index[-1][1]) + 0 + 1
                            io_column_index_map.append(
                                [new_index, old_index]
                            )
                elif reading_comps:
                    n += 1
                    # put 0 instead of blank and keep all columns to add in
                    # model
                    row_to_add = [''] * len(header_comps_input_model)
                    for new_index, old_index in io_column_index_map:
                        text_with_number = row_without_whitespace[old_index]
                        if new_index!=1 and len(text_with_number) == 0:
                            row_to_add[new_index] = float(0)
                        elif new_index!=1:
                            row_to_add[new_index] = float(text_with_number)
                        else:
                            row_to_add[new_index] = row_without_whitespace[old_index]
                    comps.append(row_to_add)
                elif reading_reacs:
                    nr += 1
                    sorted_io_column_index_map =\
                        sorted(io_column_index_map, key=lambda x: x[1])
                    max_column_no = max([item[0]
                                         for item in io_column_index_map]) + 1
                    row_to_add = [''] * max_column_no
                    # put 0 instead of blank and keep only columns to add
                    for new_index, old_index in sorted_io_column_index_map:
                        text_with_number = row_without_whitespace[old_index]
                        if len(text_with_number) == 0:
                            row_to_add[new_index] = float(0)
                        else:
                            row_to_add[new_index] = float(text_with_number)
                    reacs.append(row_to_add)
        csv_file.close()
        self.parentWidget().setWindowTitle(
            main_window_title[:-1] + ' - ' + os.path.basename(filename))
        column_of_index_comps = comp_variable_input_names.index('index')
        column_of_index_reacs = header_reacs_model.index('j')
        # First, sort by existing order, if available
        sorted_comps = sorted(comps, key=lambda x: x[column_of_index_comps])
        sorted_reacs = sorted(reacs, key=lambda x: x[column_of_index_reacs])
        # Add indexes or replace existing indexes with simple ones.
        for k, row in enumerate(sorted_comps):
            row[column_of_index_comps] = k + 1
        for k, row in enumerate(sorted_reacs):
            row[column_of_index_reacs] = k + 1
        header_comps = comp_variable_input_names \
            + comp_variable_output_names
        header_reacs = header_reacs_model \
            + ['nu_' + str(x + 1) + 'j' for x in range(n)]
        comps = np.array(sorted_comps, dtype=object)  # do not convert yet
        reacs = np.array(sorted_reacs, dtype=object)  # do not convert yet
        # Grow comps and reacs to their final widths from the start.
        # Width of reactions matrix could not be known before this step.
        header_comps_complete = \
            header_comps_input_model + header_comps_output_model
        header_reacs_complete = \
            header_reacs + ['xieq']
        comps_column_width = len(header_comps_complete)
        reacs_column_width = len(header_reacs_complete)
        comps_completed_matrix = np.empty(
            [comps.shape[0], comps_column_width], dtype=object
        )
        comps_completed_matrix[:, 0:comps.shape[1]] = comps
        reacs_completed_matrix = np.empty(
            [reacs.shape[0], reacs_column_width], dtype=object
        )
        reacs_completed_matrix[:, 0:reacs.shape[1]] = reacs
        self.spinBox.setProperty("value", n)
        self.spinBox_2.setProperty("value", nr)

        self.comps_model = MatrixModel(
            comps_completed_matrix,
            header_comps_complete,
            editable_columns=range(len(header_comps_input_model))
        )
        self.reacs_model = MatrixModel(
            reacs_completed_matrix,
            header_reacs_complete,
            editable_columns=range(len(header_reacs_complete) - 1)
        )
        self.tableComps.setSortingEnabled(False)
        self.tableReacs.setSortingEnabled(False)
        # Set modelsfor tables and connect datachanged signals
        self.tableComps.setModel(self.comps_model)
        self.tableReacs.setModel(self.reacs_model)
        self.tableComps.model().change_data.connect(
            self.recalculate_after_cell_edit)
        self.tableReacs.model().change_data.connect(
            self.recalculate_after_cell_edit)
        # Pass variables to self before loop start
        variables_to_pass = ['header_comps', 'comps',
                             'header_comps_complete',
                             'comps_completed_matrix',
                             'header_reacs', 'reacs',
                             'header_reacs_complete',
                             'reacs_completed_matrix',
                             'n', 'nr'
                             ]
        for var in variables_to_pass:
            setattr(self, var, locals()[var])

    def load_variables_from_form(self):
        comps_completed_matrix = \
            self.tableComps.model().return_data()
        reacs_completed_matrix = \
            self.tableReacs.model().return_data()
        header_comps = self.tableComps.model().return_headers()
        header_reacs = self.tableReacs.model().return_headers()
        n = len(comps_completed_matrix)
        nr = len(reacs_completed_matrix)
        index_of_component_order_in_table = \
            header_comps.index('i')
        index_of_reaction_order_in_table = \
            header_reacs.index('j')
        component_order_in_table = \
            (comps_completed_matrix[
                :, index_of_component_order_in_table
            ] - 1).tolist()
        reaction_order_in_table = \
            (reacs_completed_matrix[
                :, index_of_reaction_order_in_table
            ] - 1).tolist()
        # Pass variables to self before loop start
        variables_to_pass = ['header_comps',
                             'comps_completed_matrix',
                             'header_reacs',
                             'reacs_completed_matrix',
                             'component_order_in_table',
                             'reaction_order_in_table']
        for var in variables_to_pass:
            setattr(self, var, locals()[var])
        self.gui_setup_and_variables()

    def gui_setup_and_variables(self):
        # Collect variables
        n = self.n
        nr = self.nr
        comps = np.array(sorted(
            self.comps_completed_matrix, key=lambda x: x[0]))
        reacs = np.array(sorted(
            self.reacs_completed_matrix, key=lambda x: x[0]))
        header_comps = self.header_comps
        header_reacs = self.header_reacs
        molar_masses_valid = False
        unset_variables = []
        for col, name in enumerate(comp_variable_input_names):
            if name in ['comp_id']:
                data_type = str
            elif name in ['index']:
                data_type = int
            else:
                data_type = float
                # take floats, replace empty strings for 0.0
                column_vector = np.array([float(x) if  not isinstance(x,str) else 0.0 for x in comps[:,col] ])
                if all(column_vector == 0):
                    unset_variables.append(name)
            if data_type != float:
                try:
                    column_vector = np.array(comps[:, col], dtype=data_type)
                except ValueError as detail:
                    # print detail
                    unset_variables.append(name)
                    column_vector = np.empty([n, 1])
                    column_vector[:] = np.nan
                    if name in ['index', 'comp_id', 'z']:
                        raise Exception('Input field missing: '
                                        + name)
            # Put / Reset values of each column name into self by name
            if hasattr(self, name):
                delattr(self, name)
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
        rho_solvent = 0.997  # TODO: Adjust density with temperature
        positive_molar_masses = \
            all(molar_mass > 0)
        can_calculate_n_from_w = \
            'molar_mass' not in unset_variables and \
            ('xw0' not in unset_variables or
             'w0' not in unset_variables)
        mol_variables_empty = \
            all(map(lambda x: x in unset_variables,
                    ['n0', 'x0', 'c0']))
        # Gui setup with calculated values
        # First determine concentration variables from available data, and
        # determine the index of the solvent.
        highest_n0_indexes = []
        index_of_solvent = []
        c_solvent_tref = []
        if mol_variables_empty:
            # If number of moles cannot be determined, reaction quotients
            # cannot be determined either, request molar mass inputs.
            if not can_calculate_n_from_w or \
                    not positive_molar_masses:
                raise Exception('Need positive molar masses defined')
            else:
                # calculate n0 from w0 or xw0
                if 'w0' in unset_variables:
                    w0 = xw0
                elif 'xw0' in unset_variables:
                    xw0 = w0 / sum(w0)
                n0 = np.divide(w0, molar_mass)
        elif 'n0' in unset_variables and \
                'c0' not in unset_variables:
            # moles not given, but molarity given:
            # overwrite mole number based on 1L of molarity
            n0 = c0
        elif 'w0' in unset_variables and \
                'xw0' not in unset_variables and \
                'molar_mass' not in unset_variables:
            # only w% given, use weight in base 1.
            w0 = xw0
            n0 = np.divide(w0, molar_mass)
        elif 'n0' not in unset_variables:
            # moles and molar masses given:
            # use as main generators
            pass
        highest_n0_indexes = np.argpartition(n0.flatten(), (-1, -2))
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
        gamma0_ii = np.ones_like(m0)
        gamma0_iii = np.ones_like(m0)
        a0 = np.multiply(gamma0_iii, m0) * np.nan
        rho0_i = np.multiply(c0, mm0)
        c_solvent_tref = c0[index_of_solvent].item()
        # Overwrite logarithmic calculations, no need to keep
        # supplied values
        mlog10xw0 = -np.log10(xw0)
        mlog10x0 = -np.log10(x0)
        mlog10c0 = -np.log10(c0)
        mlog10m0 = -np.log10(m0)
        mlog10a0 = -np.log10(a0)
        nu_ij = np.array([row[2:2 + n] for row in reacs], dtype=float).T
        pka = np.array([row[1] for row in reacs], dtype=float).T
        max_it = int(self.spinBox_3.value())
        tol = float(self.doubleSpinBox_5.value())
        # Determine second highest component for plotting possibilities
        if len(n0) > 1:
            index_of_second_highest_n0 = highest_n0_indexes[-2]
        else:
            index_of_second_highest_n0 = highest_n0_indexes[-1]
        n_second_highest_n0_tref = n0[index_of_second_highest_n0].item()
        # Calculated variables to output for other functions.
        variables_to_pass = [
            'index_of_solvent', 'molar_mass',
            'rho0', 'rho_solvent', 'c_solvent_tref',
            'max_it', 'tol',
            'n_second_highest_n0_tref',
            'positive_molar_masses',
            'can_calculate_n_from_w',
            'comps', 'reacs'
        ]   + comp_variable_input_names \
            + ['nu_ij', 'pka']
        for var in variables_to_pass:
            setattr(self, var, locals()[var])
        # Setup plotting tools
        self.comboBox.clear()
        self.comboBox_3.clear()
        for item in comps[:, 0:2]:
            self.comboBox.addItem('c0_' + str(item[0]) +
                                  ' {' + item[1] + '}')
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
        header_comps_complete = self.header_comps_complete
        header_reacs_complete = self.header_reacs_complete
        reacs = self.reacs
        xieq = self.xieq

        if hasattr(self, 'component_order_in_table'):
            i = getattr(self, 'component_order_in_table')
        else:
            i = range(0, n)

        if hasattr(self, 'reaction_order_in_table'):
            j = getattr(self, 'reaction_order_in_table')
        else:
            j = range(0, nr)

        self.tableComps.blockSignals(True)
        self.tableReacs.blockSignals(True)
        self.comboBox.blockSignals(True)
        # As usual, problems occurr when sorting is combined with setting QTableView.
        # Therefore disable sorting, then set model and finally
        # reenable sorting.
        self.tableComps.setSortingEnabled(False)
        self.tableReacs.setSortingEnabled(False)

        for column, column_name in enumerate(header_comps_complete):
            var_name = comp_names_headers[column_name]
            # Keep original order from table
            column_data = getattr(self, var_name)[i]
            if var_name in ['m0', 'meq']:
                # Present units converted: mol/gsolvent to mol/kgsolvent
                column_data = 1000 * column_data
            elif var_name in ['mlog10m0', 'mlog10meq']:
                # Present units converted: mol/gsolvent to mol/kgsolvent.
                # Add -log10(10^3)
                column_data = column_data - 3
            self.comps_model.set_column(column, column_data)

        for column, column_name in enumerate(header_reacs_complete):
            if column != n + 0 + 1 + 1:
                # Keep original order from table in sorted_reacs
                column_data = reacs[j, column]
            elif column == n + 0 + 1 + 1:
                # Keep original order from table
                column_data = getattr(self, column_name)[j]
            self.reacs_model.set_column(column, column_data)

        # Widths and heights, re-enable sorting
        self.tableComps.setSortingEnabled(True)
        self.tableComps.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeToContents)
        self.tableComps.verticalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeToContents)

        self.tableReacs.setSortingEnabled(True)
        self.tableReacs.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeToContents)
        self.tableReacs.verticalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeToContents)

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
        n_points = 20
        index_of_variable = self.comboBox.currentIndex()
        comps = self.comps
        c0 = self.c0
        n = self.n
        nr = self.nr
        c0_variable_comp = c0[index_of_variable]
        rho_solvent = self.rho_solvent
        xieq = self.xieq
        ceq = self.ceq
        ceq_series = np.array(np.zeros([n_points + 1, n]))
        xieq_series = np.array(np.zeros([n_points + 1, nr]))
        # Keep current solution intact for after plotting range
        self.stored_solution_ceq = self.ceq
        self.stored_solution_xieq = self.xieq
        # TODO: Get plotting to work with format QTableView model (comps
        # matrix)
        indep_var_label = 'c0_{' + str(comps[index_of_variable, 0]) + ', ' + \
                          comps[index_of_variable, 1] + '}/(mol/L)'
        dep_var_labels = \
            ['ceq_' + '{' + str(item[0]) + ', ' + item[1] + '}/(mol/L)' for item in comps[:, 0:2]] + \
            ['\\xi eq_' +
                '{' + str(item) + '}/(mol/L)' for item in range(1, nr + 1, 1)]
        min_value = self.doubleSpinBox.value()
        max_value = self.doubleSpinBox_2.value()
        indep_var_series_single = \
            [min_value + x
             for x in np.arange(n_points + 1) * (max_value - min_value) / n_points]
        mid_index = bisect.bisect(
            indep_var_series_single,
            c0_variable_comp) - 1
        dep_var_series = dict(
            zip(dep_var_labels, np.empty(n + nr, dtype=np.ndarray)))
        indep_var_series = dict.fromkeys(
            dep_var_labels, indep_var_series_single)
        for j in range(mid_index, -1, -1):
            # input is in molar conc. Get molal, n, x
            self.c0[index_of_variable] = indep_var_series_single[j]
            self.m0[index_of_variable] = \
                self.c0[index_of_variable] / rho_solvent
            self.x0 = self.c0 / sum(self.c0)
            tot_n0_const = sum(self.n0)
            self.n0[index_of_variable] = \
                tot_n0_const * self.x0[index_of_variable]
            self.equilibrate()
            ceq_series[j, :] = self.ceq.T
            xieq_series[j, :] = self.xieq.T
        self.ceq = self.stored_solution_ceq
        self.xieq = self.stored_solution_xieq
        for j in range(mid_index + 1, n_points + 1, +1):
            # input is in molar conc. Get molal, n, x
            self.c0[index_of_variable] = indep_var_series_single[j]
            self.m0[index_of_variable] = \
                self.c0[index_of_variable] / rho_solvent
            self.x0 = self.c0 / sum(self.c0)
            tot_n0_const = sum(self.n0)
            self.n0[index_of_variable] = \
                tot_n0_const * self.x0[index_of_variable]
            self.equilibrate()
            ceq_series[j, :] = self.ceq.T
            xieq_series[j, :] = self.xieq.T
        for j in range(n):
            dep_var_series[dep_var_labels[j]] = ceq_series[:, j]
        for j in range(nr):
            dep_var_series[dep_var_labels[n + j]] = xieq_series[:, j]
        self.ceq = self.stored_solution_ceq
        self.xieq = self.stored_solution_xieq
        self.ceq_series = ceq_series
        self.xieq_series = xieq_series
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
        labels_to_plot = [x for x in self.dep_var_labels if x.find('ceq') >= 0]
        # dict, keys:ceq_labels; bindings: plottedseries
        plotted_series = dict(
            zip(dep_var_labels, np.empty(n + nr, dtype=object)))
        dep_var_labels = self.dep_var_labels
        dep_var_series = self.dep_var_series
        indep_var_series = self.indep_var_series
        indep_var_label = self.indep_var_label
        self.groupBox = QtWidgets.QGroupBox()
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
        s_index = self.index_of_solvent
        m0_ref = 1 / 1000.0  # ref. 1mol/kgsolvent conv. to mol/gsolvent
        method = 'ideal_solution'

        if self.radio_b_2.isChecked():
            method = 'debye-hueckel'
        elif self.radio_b_3.isChecked():
            method = 'davies'

        # Init. calculations
        kc = np.multiply(np.power(10, -pka),
                         np.power(c_solvent_tref, nu_ij[s_index, :]).T)
        self.kc = kc

        # Setup logging
        if not os.path.exists('./logs'):
            os.mkdir('./logs')
        logging.basicConfig(
            filename='./logs/calculation_results.log',
            level=logging.DEBUG,
            format='%(asctime)s;%(message)s')

        # First estimates for eq. Composition ceq and Reaction extent xieq
        if not hasattr(self, 'acceptable_solution'):
            neq_0 = n0
            # replace 0 by 10^-6*smallest value: Smith, Missen 1988 DOI:
            # 10.1002/cjce.5450660409
            neq_0[n0 == 0] = min(n0[n0 != 0].flatten()) * np.finfo(float).eps
            xieq_0 = np.zeros(nr)
        else:
            # Use previous solution as initial estimate, if it was valid.
            neq_0 = self.neq_0
            xieq_0 = self.xieq_0

        # Pass variables to self before loop start
        variables_to_pass = ['nu_ij', 'pka',
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
        self.method_loops = [0, 0]  # loop numbers: [Line search, Newton]
        # Unique identifier for plotting logged solutions
        series_id = str(uuid.uuid1())
        # Calculate equilibrium composition: Newton method
        # TODO: Implement global homotopy-continuation method
        while not self.acceptable_solution \
                and k < max_it and stop is False \
                and not self.was_canceled():
            neq, meq, xieq, gammaeq_iii, ionic_str_eq, method_loops = \
                calc_xieq(
                    n0, mm, z, s_index, kc, nu_ij,
                    neq_0, xieq_0, 298.15, method, max_it, tol,
                    self.method_loops,
                    partial(update_status_label,
                            series_id,
                            self,
                            'Newton'
                            ),
                    lambda: QtWidgets.QApplication.processEvents()
                )
            self.method_loops = method_loops
            k += 1
            # TODO: if progress_var.wasCanceled() == True then stop
            if all(neq >= 0) and not any(np.isnan(neq)):
                self.acceptable_solution = True
            else:
                # Set reactions to random extent and recalculate
                self.xieq_0 = np.array(
                    np.random.normal(0.0, 1.0 / 3.0, nr)).T
                # Set aequilibrium composition to initial value + estimated
                # conversion
                self.neq_0 = n0  # + nu_ij * self.xieq_0
                # replace 0 by 10^-6*smallest value: Smith, Missen 1988 DOI:
                # 10.1002/cjce.5450660409
                self.neq_0[self.neq_0 == 0] = min(
                    n0[n0 != 0].flatten()) * np.finfo(float).eps
                self.initialEstimateAttempts += 1
                self.method_loops = [0, 0]
                stop = True

        if not self.acceptable_solution:
            delattr(self, 'acceptable_solution')
            self.label_9.setText(self.label_9.text() + '\n')
        else:
            self.neq_0 = neq
            self.xieq_0 = xieq
            self.label_9.setText(self.label_9.text() +
                                 '\nsum(n0*z) = ' +
                                 str((z.T.dot(n0)).item()) +
                                 ' \t\t\t sum(neq*z) = ' +
                                 str((z.T.dot(neq)).item()) +
                                 '\nI_0 = ' +
                                 str((1 / 2.0 * np.power(z, 2).T.dot(m0) * 1000).item()) +
                                 '(mol/kg_solv)' +
                                 '\t\t\t\t I_eq = ' +
                                 str(ionic_str_eq.item() / m0_ref) +
                                 '(mol/kg_solv)')

        # Once solved, calculate conc. variables at equilibrium
        xeq = neq / sum(neq)
        ceq = meq * rho_solvent * 1000
        rhoeq = np.multiply(ceq, mm)
        weq = np.multiply(neq, mm)
        xweq = weq / sum(weq)
        aeq = np.multiply(gammaeq_iii, meq * 1000)
        gammaeq_ii = mm[s_index].item() * np.multiply(
            gammaeq_iii, np.multiply(
                meq, np.power(xeq, -1)))
        mlog10gammaeq_ii = -np.log10(gammaeq_ii)
        mlog10gammaeq_iii = -np.log10(gammaeq_iii)
        mlog10xweq = -np.log10(xweq)
        mlog10xeq = -np.log10(xeq)
        mlog10ceq = -np.log10(ceq)
        mlog10meq = -np.log10(meq)
        mlog10aeq = -np.log10(aeq)
        # Join passing variables to self with comp matrix
        for col, name in enumerate(comp_variable_output_names):
            column_in_comps_matrix = \
                len(comp_variable_input_names) + col
            comps[:, column_in_comps_matrix] = \
                locals()[name].reshape(1, -1)
            setattr(self, name, locals()[name])
        for col, name in enumerate(['xieq']):
            column_in_reacs_matrix = \
                reacs.shape[1] - 1 + col
            reacs[:, column_in_reacs_matrix] = \
                locals()[name].reshape(1, -1)
            setattr(self, name, locals()[name])
        self.comps = comps
        self.cancelButton.setEnabled(False)
        self.progress_var.setEnabled(False)

    def display_about_info(self):
        row_string = unicode('', 'utf_8')
        self.browser_window = QtWidgets.QWidget()
        self.browser_window.setWindowIcon(QtGui.QIcon(
            os.path.join(sys.path[0], *['utils', 'icon_batch_16X16.png'])))
        self.browser_window.setMinimumWidth(1000)
        self.browser_window.setWindowTitle('About')
        verticalLayout = QtWidgets.QVBoxLayout(self.browser_window)
        # Handle fixed 96 dpi for Webkit bugreport
        # https://bugreports.qt-project.org/browse/QTBUG-29571
        window = QtWidgets.QApplication.desktop().screen()
        horizontalDpi = window.logicalDpiX()
        aboutBox_1 = QtWebKit.QWebView()
        aboutBox_1.setZoomFactor(horizontalDpi / 96.0)
        toolbar_1 = QtWidgets.QToolBar()
        back_icon = QtGui.QIcon(
            os.path.join(
                sys.path[0],
                *['utils', 'glyphicons-217-circle-arrow-left.png']))
        forward_icon = QtGui.QIcon(
            os.path.join(
                sys.path[0],
                *['utils', 'glyphicons-218-circle-arrow-right.png']))
        comboBox_1 = QtWidgets.QComboBox()
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
             ('stop', bool),
             ('X', list),
             ('||X(k)-X(k-1)||', float),
             ('f(X)', list),
             ('||f(X)||', float),
             ('j(X)', list),
             ('Y', list),
             ('||Y||', float),
             ('g', float),
             ('|g-g1|', float),
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
        self.groupBox = QtWidgets.QGroupBox()
        self.groupBox.log_widget = LogWidget(log, parent=self.groupBox)
        self.groupBox.show()


class UiGroupBoxPlot(QtWidgets.QWidget):

    def __init__(self, parent):
        QtWidgets.QWidget.__init__(self, parent)
        # Assignments
        log_scale_func_list = [
            ('-log10', lambda x: -1.0 * np.log10(x)), ('', lambda x: 10.0 ** (-1.0 * x))]
        self.icon_down = QtGui.QIcon(os.path.join(
            sys.path[0], *['utils', 'glyphicons-602-chevron-down.png']))
        self.icon_up = QtGui.QIcon(os.path.join(
            sys.path[0], *['utils', 'glyphicons-601-chevron-up.png']))
        self.verticalLayout_1 = QtWidgets.QVBoxLayout(parent)
        self.horizontalTools = QtWidgets.QHBoxLayout()
        self.toolsFrame = QtWidgets.QFrame()
        self.toggleLogButtonX = QtWidgets.QPushButton()
        self.toggleLogButtonY = QtWidgets.QPushButton()
        self.eraseAnnotationsB = QtWidgets.QPushButton()
        self.navigation_frame = QtWidgets.QFrame()
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        # Handle fixed 96 dpi for Webkit bugreport
        # https://bugreports.qt-project.org/browse/QTBUG-29571
        window = QtWidgets.QApplication.desktop().screen()
        horizontalDpi = window.logicalDpiX()
        canvas_factor = 72.0 / 96.0
        self.figure = Figure(dpi=horizontalDpi * canvas_factor,
                             facecolor=(1, 1, 1),
                             edgecolor=(0, 0, 0))
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.label = QtWidgets.QLabel(parent)
        self.listWidget_2 = QtWidgets.QListWidget(parent)
        self.label_2 = QtWidgets.QLabel(parent)
        self.listWidget = QtWidgets.QListWidget(parent)
        self.plotButton = QtWidgets.QPushButton(parent)
        self.toolbar = NavigationToolbar(self.canvas, self.navigation_frame)
        self.horizontalLayout_1 = QtWidgets.QHBoxLayout()
        self.move_10_to_displayed = QtWidgets.QPushButton()
        self.move_10_to_available = QtWidgets.QPushButton()
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
            QtWidgets.QAbstractItemView.MultiSelection)
        self.listWidget.setSelectionMode(
            QtWidgets.QAbstractItemView.MultiSelection)
        self.listWidget_2.setViewMode(QtWidgets.QListView.ListMode)
        self.listWidget_2.setSortingEnabled(False)
        self.listWidget.setSortingEnabled(False)
        self.verticalLayout.addWidget(self.label)
        self.verticalLayout.addWidget(self.listWidget_2)
        self.verticalLayout.addLayout(self.horizontalLayout_1)
        self.move_10_to_displayed.setIcon(self.icon_up)
        self.move_10_to_available.setIcon(self.icon_down)
        self.horizontalLayout_1.addWidget(self.move_10_to_available)
        self.horizontalLayout_1.addWidget(self.move_10_to_displayed)
        self.verticalLayout.addWidget(self.label_2)
        self.verticalLayout.addWidget(self.listWidget)
        self.verticalLayout.addWidget(self.plotButton)
        self.horizontalLayout.addLayout(self.verticalLayout)
        # Toolbar height adjusted once added.
        self.toolbar.adjustSize()
        self.navigation_frame.setMinimumHeight(self.toolbar.height())
        # Events
        self.listWidget_2.itemDoubleClicked.connect(
            lambda: self.move_to_available(item_no=-1))
        self.listWidget.itemDoubleClicked.connect(
            lambda: self.move_to_displayed(item_no=-1))
        self.toggleLogButtonX.toggled.connect(
            partial(self.toggled_toggle_log_button_x))
        self.toggleLogButtonY.toggled.connect(
            partial(self.toggled_toggle_log_button_y))
        self.eraseAnnotationsB.clicked.connect(partial(self.erase_annotations))
        self.plotButton.clicked.connect(partial(self.force_update_plot))
        self.move_10_to_available.clicked.connect(
            lambda: self.move_to_available(item_no=10))
        self.move_10_to_displayed.clicked.connect(
            lambda: self.move_to_displayed(item_no=10))
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
                new_item = QtWidgets.QListWidgetItem(label, self.listWidget_2)
                new_item.setIcon(self.icon_down)
            else:
                new_item = QtWidgets.QListWidgetItem(label, self.listWidget)
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
            QtWidgets.QApplication.translate(
                "parent",
                "Plot",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        parent.setTitle(
            QtWidgets.QApplication.translate(
                "parent",
                "Plot",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        self.label.setText(
            QtWidgets.QApplication.translate(
                "parent",
                "Displayed",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        self.label_2.setText(
            QtWidgets.QApplication.translate(
                "parent",
                "Available",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        self.plotButton.setText(
            QtWidgets.QApplication.translate(
                "parent",
                "Plot",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        self.toggleLogButtonX.setText(
            QtWidgets.QApplication.translate(
                "parent",
                self.log_scale_string +
                "(x) - horizontal",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        self.toggleLogButtonY.setText(
            QtWidgets.QApplication.translate(
                "parent",
                self.log_scale_string +
                "(y) - vertical",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        self.eraseAnnotationsB.setText(
            QtWidgets.QApplication.translate(
                "parent",
                "Erase annotations",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        self.move_10_to_available.setText(
            QtWidgets.QApplication.translate(
                "parent",
                "(10)",
                None,
                QtWidgets.QApplication.UnicodeUTF8))
        self.move_10_to_displayed.setText(
            QtWidgets.QApplication.translate(
                "parent",
                "(10)",
                None,
                QtWidgets.QApplication.UnicodeUTF8))

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
                            match.group('id')].flatten().tolist())
                else:
                    line.set_ydata(self.invlog_scale_func(line_ydata))
                line.set_label(u'$' + match.group('id') + u'$')
            else:
                if match.group('id2') in self.dep_var_series.keys():
                    line.set_ydata(
                        self.log_dep_var_series[
                            match.group('id2')].flatten().tolist())
                else:
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
                indep_var_values, dep_var_values.flatten().tolist(), 'go-', label=series_label,
                color=colormap_colors[
                    np.random.randint(
                        0, len(colormap_colors), 1).item()],
                markerfacecolor=colormap_colors[
                    np.random.randint(0, len(colormap_colors), 1).item()],
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

    def move_to_available(self, item_no=-1):
        self.delete_arrows()
        selected_items = self.listWidget_2.selectedItems()
        if self.listWidget_2.count() <= 1:
            return  # Stop if already showing only one item
        if item_no < 1:  # Move selected
            pass  # Default -1 for selected itemd
        else:  # Move item_no items
            no_items_to_move = self.listWidget_2.count()
            if item_no >= no_items_to_move:
                no_items_to_move = no_items_to_move - 1
            elif item_no < no_items_to_move:
                # enough available, to remove item_no
                no_items_to_move = item_no
            selected_items = \
                [self.listWidget_2.item(x)
                 for x in range(no_items_to_move)]
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

    def move_to_displayed(self, item_no=-1):
        self.delete_arrows()
        selected_items = self.listWidget.selectedItems()
        if item_no < 1:  # Move selected
            pass  # Default -1 for selected items
        else:  # Move item_no items
            no_items_to_move = self.listWidget.count()
            if item_no >= no_items_to_move:
                # Possible to move all to displayef
                no_items_to_move = no_items_to_move
            elif item_no < no_items_to_move:
                # enough available, to move item_no
                no_items_to_move = item_no
            selected_items = \
                [self.listWidget.item(x)
                 for x in range(no_items_to_move)]
        selected_items_names = [x.text() for x in selected_items]
        for selected_item in selected_items:
            name = selected_item.text()
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


class LogWidget(QtWidgets.QWidget):

    def __init__(self, _log, parent):
        QtWidgets.QWidget.__init__(self, parent)
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
        # Assignments
        self.verticalLayout = QtWidgets.QVBoxLayout(parent)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.pandasView = QtWidgets.QTableView()
        self.firstButton = QtWidgets.QPushButton()
        self.lastButton = QtWidgets.QPushButton()
        self.nextButton = QtWidgets.QPushButton()
        self.previousButton = QtWidgets.QPushButton()
        self.pageLabel = QtWidgets.QLabel()
        self.totPagesLabel = QtWidgets.QLabel()
        self.page_box = QtWidgets.QLineEdit()
        self.exportButton = QtWidgets.QPushButton()
        self.plotButton = QtWidgets.QPushButton()
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
        self.pandasView.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Interactive)
        self.pandasView.verticalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeToContents)
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
        (file_name, selected_filter) = QtWidgets.QFileDialog.getSaveFileName(
            parent=self,
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
            dep_var_series[index] = np.array(group['||f(X)||'].values)
        # Generate the plot
        self.group_3 = QtWidgets.QGroupBox()
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
        self.group_3.show()
        self.plotBox.ax.set_ylabel('||f(X)||')


def update_status_label(
        series_id,
        form,
        nle_method,
        progress,
        stop_value,
        k,
        j_it_backtrack,
        lambda_ls,
        accum_step,
        x,
        diff,
        f_val,
        j_val,
        lambda_ls_y,
        method_loops):
    status = 'solving...'
    if stop_value and progress == 100.0:
        status = 'solved.'
    else:
        status = 'divergent'
    # Steepest descent method deprecated.
    g_min = np.nan
    g1 = np.nan
    y = lambda_ls_y
    # Add progress bar & variable
    form.progress_var.setValue(int(progress))
    form.label_9.setText('Loops: Newton \t' +
                         str(method_loops[1]) +
                         ' \t Line search (backtrack) \t' +
                         str(method_loops[0]) +
                         ' \t Initial estimate attempts \t' +
                         str(form.initialEstimateAttempts) +
                         '\n' +
                         'Iteration (k) \t' +
                         str(k) +
                         '\n' +
                         str('solving...' if not stop_value
                             else status))
    logging.debug(nle_method + ' ' +
                  ';k=' + str(k) +
                  ';backtrack=' + str(j_it_backtrack) +
                  ';lambda_ls=' + str(lambda_ls) +
                  ';accum_step=' + str(accum_step) +
                  ';stop=' + str(stop_value) +
                  ';X=' + '[' + ','.join(map(str, x.T.flatten())) + ']' +
                  ';||X(k)-X(k-1)||=' + str((diff.T.dot(diff)).item()) +
                  ';f(X)=' + '[' + ','.join(map(str, f_val.T.flatten())) + ']' +
                  ';||f(X)||=' + str(np.sqrt((f_val.T.dot(f_val)).item())) +
                  ';j(X)=' + str(j_val.tolist()) +
                  ';Y=' + '[' + ','.join(map(str, y.T.flatten())) + ']' +
                  ';||Y||=' + str(np.sqrt((y.T.dot(y)).item())) +
                  ';g=' + str(g_min) +
                  ';|g-g1|=' + str(abs(g_min - g1)) +
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

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            return str(self._data.values[index.row()][index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self._data.columns[col]
        return None


class MatrixModel(QtCore.QAbstractTableModel):
    """
    Used to populate a QTableView with an np.Matrix
    """

    change_data=QtCore.Signal(QtCore.QModelIndex,QtCore.QModelIndex,name='dataChanged')

    def __init__(self, data, column_names,
                 editable_columns=None, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = data
        if editable_columns is None:
            self._editable_columns = [None]
        elif len(editable_columns) == 1:
            self._editable_columns = [editable_columns]
        elif len(editable_columns) > 1:
            self._editable_columns = editable_columns
        width = data.shape[1]
        if len(column_names) == width:
            self._column_names = column_names
        else:
            self._column_names = [''] * width

    def rowCount(self, *args, **kwargs):
        return self._data.shape[0]

    def columnCount(self, *args, **kwargs):
        return self._data.shape[1]

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            return str(self._data[index.row(), index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self._column_names[col]
        return None

    def set_column(self, index, column_array):
        height = len(self._data)
        # match whether it is column or row vector
        if column_array.shape[0] == height:
            self._data[:, index] = column_array.reshape(1, -1)
        elif column_array.shape[1] == height:
            self._data[:, index] = column_array

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        # if self._data.shape == value.shape:
        #     self._data = new_data
        if index.isValid() and 0 <= index.row() < len(self._data):
            row = index.row()
            column = index.column()
            self._data[row, column] = value
            self.dirty = True
            self.change_data.emit(
                index, index
            )
            return True
        return False

    def sort(self, column, order=QtCore.Qt.AscendingOrder):
        if order == QtCore.Qt.AscendingOrder:
            reverse = False
        else:
            reverse = True
        self._data = np.array(
            sorted(self._data,
                   key=lambda x: x[column],
                   reverse=reverse)
        )
        self.beginResetModel()
        # print 'sorted by col: ' + str(column)

    def flags(self, index):
        if index.column() in self._editable_columns:
            return QtCore.Qt.ItemIsEnabled | \
                QtCore.Qt.ItemIsSelectable | \
                QtCore.Qt.ItemIsEditable
        else:
            return QtCore.Qt.ItemIsEnabled | \
                QtCore.Qt.ItemIsSelectable

    def return_data(self):
        return self._data

    def return_headers(self):
        return self._column_names


class ScientificDoubleSpinBox(QtWidgets.QDoubleSpinBox):

    def __init__(self, parent=None):
        QtWidgets.QDoubleSpinBox.__init__(self, parent)
        self.setMinimum(-np.inf)
        self.setMaximum(np.inf)
        self.validator = QtGui.QDoubleValidator()
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
    app = QtWidgets.QApplication.instance()  # checks if QApplication already exists
    if not app:  # create QApplication if it doesnt exist
        app = QtWidgets.QApplication(sys.argv)
    # following 2 lines for setting app icon correctly
    myappid = u'mycompany.myproduct.subproduct.version'  # arbitrary string
    #ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)
    # ended lines for setting app icon correctly
    main_form = QtWidgets.QGroupBox()
    icon = QtGui.QIcon(
        os.path.join(sys.path[0], *['utils', 'icon_batch_chp_96X96.png']))
    main_form.setWindowIcon(icon)
    main_form.ui = UiGroupBox(main_form)
    main_form.show()
    main_form.ui.load_csv('./DATA/COMPONENTS_REACTIONS_EX_001.csv')
    main_form.ui.gui_equilibrate()
    main_form.ui.tableComps.sortByColumn(0, QtCore.Qt.AscendingOrder)
    main_form.ui.tableReacs.sortByColumn(0, QtCore.Qt.AscendingOrder)
    app.exec_()


if __name__ == '__main__':
    main()
