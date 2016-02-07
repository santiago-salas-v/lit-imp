# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 20:48:42 2015

@author: Santiago Salas
@ref: Denbigh, p. 298
"""
import os, sys, logging, re, pandas as pd, numpy as np, scipy as sp, csv, bisect
import matplotlib, colormaps

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PySide import QtGui, QtCore
from functools import partial
from mat_Zerlegungen import gausselimination
from datetime import datetime
from mpldatacursor import datacursor

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


class UiGroupBox(object):
    _was_canceled = False

    def __init__(self, parent):
        parent.setObjectName(_fromUtf8("GroupBox"))
        # Default size
        parent.resize(500, 357)
        self.verticalLayout_2 = QtGui.QVBoxLayout(parent)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.open_button = QtGui.QPushButton(parent)
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
        self.save_button = QtGui.QPushButton(parent)
        self.save_button.setIcon(icon1)
        self.save_button.setObjectName(_fromUtf8("save_button"))
        self.info_button = QtGui.QPushButton(parent)
        self.info_button.setIcon(icon4)
        self.info_button.setObjectName(_fromUtf8("info_button"))
        self.log_button = QtGui.QPushButton(parent)
        self.log_button.setIcon(icon5)
        self.log_button.setObjectName(_fromUtf8("log_button"))
        self.horizontalLayout_2.addWidget(self.save_button)
        self.horizontalLayout_2.addWidget(self.log_button)
        self.horizontalLayout_2.addWidget(self.info_button)
        self.equilibrate_button = QtGui.QPushButton(parent)
        self.equilibrate_button.setIcon(icon3)
        self.horizontalLayout_3.addWidget(self.equilibrate_button)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        # TODO: Add log button
        self.label_3 = QtGui.QLabel(parent)
        self.label_3.setObjectName(_fromUtf8("max_it_label"))
        self.label_3.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label_3)
        self.spinBox_3 = QtGui.QSpinBox(parent)
        self.spinBox_3.setMaximum(2000)
        self.spinBox_3.setMinimum(2)
        self.spinBox_3.setProperty("value", 1000)
        self.spinBox_3.setObjectName(_fromUtf8("max_it_spinbox"))
        self.horizontalLayout_5.addWidget(self.spinBox_3)
        self.label_4 = QtGui.QLabel(parent)
        self.label_4.setObjectName(_fromUtf8("tol_label"))
        self.label_4.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label_4)
        self.doubleSpinBox_5 = ScientificDoubleSpinBox(parent)
        self.doubleSpinBox_5.setDecimals(int(-np.log10(np.finfo(float).eps) + 1))
        self.doubleSpinBox_5.setMaximum(float(1))
        self.doubleSpinBox_5.setMinimum(np.finfo(float).eps * 2)
        self.doubleSpinBox_5.setSingleStep(0.1 / 100.0 * (1.0 - np.finfo(float).eps))
        self.doubleSpinBox_5.setProperty("value", float(1.0e-08))
        self.doubleSpinBox_5.setObjectName(_fromUtf8("tol_spinbox"))
        self.horizontalLayout_5.addWidget(self.doubleSpinBox_5)
        self.label = QtGui.QLabel(parent)
        self.label.setObjectName(_fromUtf8("label"))
        self.label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
        self.horizontalLayout_5.addWidget(self.label)
        self.spinBox = QtGui.QSpinBox(parent)
        self.spinBox.setProperty("value", 0)
        self.spinBox.setObjectName(_fromUtf8("spinBox"))
        self.horizontalLayout_5.addWidget(self.spinBox)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.tableComps = QtGui.QTableWidget(parent)
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
        self.label_9 = QtGui.QLabel(parent)
        self.label_9.setAlignment(QtCore.Qt.AlignTop)
        self.label_9.setFrameStyle(QtGui.QFrame.Box | QtGui.QFrame.Raised)
        self.verticalLayout_2.addWidget(self.label_9)
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setAlignment(QtCore.Qt.AlignRight)
        self.cancelButton = QtGui.QPushButton(parent)
        self.progressBar = QtGui.QProgressBar(parent, visible=True)
        self.cancelButton.setEnabled(False)
        self.progressBar.setEnabled(False)
        self.horizontalLayout_7.setSizeConstraint(QtGui.QLayout.SetFixedSize)
        self.horizontalLayout_7.addStrut(max( \
            [self.progressBar.frameSize().height(),
             self.cancelButton.frameSize().height()]))
        self.horizontalLayout_7.setAlignment(QtCore.Qt.AlignLeft)
        self.horizontalLayout_7.addWidget(self.cancelButton)
        self.horizontalLayout_7.addWidget(self.progressBar)
        self.verticalLayout_2.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.label_5 = QtGui.QLabel(parent)
        self.label_5.setObjectName(_fromUtf8("solvent_label"))
        self.label_5.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_6.addWidget(self.label_5)
        self.comboBox_3 = QtGui.QComboBox(parent)
        self.horizontalLayout_6.addWidget(self.comboBox_3)
        self.label_6 = QtGui.QLabel(parent)
        self.label_6.setObjectName(_fromUtf8("C_solvent_Tref"))
        self.label_6.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        self.horizontalLayout_6.addWidget(self.label_6)
        self.doubleSpinBox_6 = QtGui.QDoubleSpinBox(parent)
        self.doubleSpinBox_6.setDecimals(int(-np.log10(np.finfo(float).eps) + 1))
        self.doubleSpinBox_6.setMaximum(float(1000))
        self.doubleSpinBox_6.setMinimum(np.finfo(float).eps * 1.1)
        self.doubleSpinBox_6.setSingleStep(1.0e-2)
        self.doubleSpinBox_6.setObjectName(_fromUtf8("C_solvent_Tref_doublespinbox"))
        self.horizontalLayout_6.addWidget(self.doubleSpinBox_6)
        self.label_2 = QtGui.QLabel(parent)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.label_2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignCenter)
        self.horizontalLayout_6.addWidget(self.label_2)
        self.spinBox_2 = QtGui.QSpinBox(parent)
        self.spinBox_2.setProperty("value", 0)
        self.spinBox_2.setObjectName(_fromUtf8("spinBox_2"))
        self.horizontalLayout_6.addWidget(self.spinBox_2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_6)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.tableReacs = QtGui.QTableWidget(parent)
        self.tableReacs.setObjectName(_fromUtf8("tableReacs"))
        self.tableReacs.horizontalHeader().setVisible(True)
        self.tableReacs.verticalHeader().setVisible(False)
        self.horizontalLayout.addWidget(self.tableReacs)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_7 = QtGui.QLabel(parent)
        self.label_7.setObjectName(_fromUtf8("horizontalAxisLabel"))
        self.label_7.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignCenter)
        self.verticalLayout.addWidget(self.label_7)
        self.comboBox = QtGui.QComboBox(parent)
        self.comboBox.setObjectName(_fromUtf8("comboBox"))
        self.verticalLayout.addWidget(self.comboBox)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.doubleSpinBox = ScientificDoubleSpinBox(parent)
        self.doubleSpinBox.setObjectName(_fromUtf8("doubleSpinBox"))
        self.doubleSpinBox.setMinimum(0.0)
        self.horizontalLayout_3.addWidget(self.doubleSpinBox)
        self.doubleSpinBox_2 = ScientificDoubleSpinBox(parent)
        self.doubleSpinBox_2.setObjectName(_fromUtf8("doubleSpinBox_2"))
        self.doubleSpinBox_2.setMinimum(0.0)
        self.horizontalLayout_3.addWidget(self.doubleSpinBox_2)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        self.pushButton = QtGui.QPushButton(parent)
        self.pushButton.setObjectName(_fromUtf8("pushButton"))
        self.pushButton.setIcon(icon2)
        self.verticalLayout.addWidget(self.pushButton)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        # Events
        self.open_button.clicked.connect(partial(open_file, self))
        self.save_button.clicked.connect(partial(save_file, self))
        self.pushButton.clicked.connect(partial(solve_intervals, self))
        self.equilibrate_button.clicked.connect(partial(recalculate_after_cell_edit, self, 0, 0))
        self.tableComps.cellChanged.connect(partial(recalculate_after_cell_edit, self))
        self.info_button.clicked.connect(partial(display_about_info, self))
        self.log_button.clicked.connect(partial(show_log))
        self.cancelButton.clicked.connect(partial(self.cancel_loop))
        self.comboBox.currentIndexChanged.connect(partial(self.populate_input_spinboxes))
        self.retranslateUi(parent)
        QtCore.QMetaObject.connectSlotsByName(parent)

    def retranslateUi(self, parent):
        parent.setWindowTitle(_translate("parent", "Simplistic EC.", None))
        parent.setTitle(QtGui.QApplication.translate("parent", "EC", None))
        __sortingEnabled = self.tableComps.isSortingEnabled()
        self.open_button.setText(_translate("parent", "Open", None))
        self.save_button.setText(_translate("parent", "Save", None))
        self.log_button.setText(_translate("parent", "Log", None))
        self.info_button.setText(_translate("parent", "About", None))
        self.equilibrate_button.setText(_translate("parent", "Equilibrate", None))
        self.tableComps.setSortingEnabled(__sortingEnabled)
        self.pushButton.setText(_translate("parent", "Plot", None))
        self.label_2.setText(_translate("parent", "Nr (Reac.)", None))
        self.label.setText(_translate("parent", "n (Comp.)", None))
        self.label_3.setText(_translate("parent", "max. it", None))
        self.label_4.setText(_translate("parent", "tol", None))
        self.label_5.setText(_translate("parent", "solvent", None))
        self.label_6.setText(_translate("parent", "C_solvent (25C)", None))
        self.label_7.setText(_translate("parent", "Horizontal 'X' axis", None))
        self.label_9.setText(_translate("parent", 'Currently unequilibrated', None))
        self.cancelButton.setText('cancel')

    def cancel_loop(self):
        self._was_canceled = True
        self.progressBar.setValue(0)

    def populate_input_spinboxes(self, index):
        comps = self.comps
        C0_component = self.C0_i[index]
        self.doubleSpinBox.setValue(C0_component / 10.0 ** 7)
        self.doubleSpinBox_2.setValue(C0_component * (1 + 20 / 100.0))

    def remove_canceled_status(self):
        self._was_canceled = False

    def was_canceled(self):
        return self._was_canceled


class UiGroupBoxPlot(object):
    def __init__(self, parent, associated_form):
        parent.setObjectName("GroupBox")
        parent.resize(762, 450)
        # parent.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_1Widget = QtGui.QWidget(parent)
        self.verticalLayout_1Widget.setGeometry(QtCore.QRect(9, 19, 741, 421))
        self.verticalLayout_1Widget.setObjectName("verticalLayout_1Widget")
        self.verticalLayout_1Widget.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_1 = QtGui.QVBoxLayout(self.verticalLayout_1Widget)
        self.verticalLayout_1.setObjectName("verticalLayout_1")
        self.verticalLayout_1.setContentsMargins(0, 0, 0, 0)
        self.horizontalTools = QtGui.QHBoxLayout()
        self.horizontalTools.setAlignment(QtCore.Qt.AlignVCenter)
        self.toolsFrame = QtGui.QFrame()
        self.toolsFrame.setLayout(self.horizontalTools)
        self.toggleLogButtonX = QtGui.QPushButton()
        self.toggleLogButtonY = QtGui.QPushButton()
        self.eraseAnnotationsB = QtGui.QPushButton()
        self.toggleLogButtonX.setCheckable(True)
        self.toggleLogButtonY.setCheckable(True)
        self.toggleLogButtonX.setChecked(True)
        self.toggleLogButtonY.setChecked(True)
        self.eraseAnnotationsB.setIcon(QtGui.QIcon(os.path.join(sys.path[0],
                                                                *['utils', 'glyphicons-551-erase.png'])))
        self.navigation_frame = QtGui.QFrame()
        self.verticalLayout_1.addWidget(self.navigation_frame)
        self.horizontalTools.addWidget(self.toggleLogButtonY)
        self.horizontalTools.addWidget(self.toggleLogButtonX)
        self.horizontalTools.addWidget(self.eraseAnnotationsB)
        self.verticalLayout_1.addWidget(self.toolsFrame)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout_1.addLayout(self.horizontalLayout)
        # Generate the plot
        self.figure = Figure(figsize=(600, 450), dpi=72, facecolor=(1, 1, 1),
                             edgecolor=(0, 0, 0))
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(bottom=0.15)
        self.ax.grid('on')
        # Generate the canvas to display the plot
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setObjectName("canvas")
        self.horizontalLayout.addWidget(self.canvas)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtGui.QLabel(self.verticalLayout_1Widget)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.listWidget_2 = QtGui.QListWidget(self.verticalLayout_1Widget)
        self.listWidget_2.setObjectName("listWidget_2")
        self.listWidget_2.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.listWidget_2.setViewMode(QtGui.QListView.ListMode)
        self.listWidget_2.setSortingEnabled(True)
        self.verticalLayout.addWidget(self.listWidget_2)
        self.label_2 = QtGui.QLabel(self.verticalLayout_1Widget)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.listWidget = QtGui.QListWidget(self.verticalLayout_1Widget)
        self.listWidget.setObjectName("listWidget")
        self.listWidget.setSortingEnabled(True)
        self.verticalLayout.addWidget(self.listWidget)
        self.pushButton = QtGui.QPushButton(self.verticalLayout_1Widget)
        self.pushButton.setObjectName("pushButton")
        self.verticalLayout.addWidget(self.pushButton)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.listWidget_2.itemDoubleClicked.connect(partial(self.move_to_available))
        self.listWidget.itemDoubleClicked.connect(partial(self.move_to_displayed))
        self.toggleLogButtonX.toggled.connect(partial(self.plot_intervals, associated_form))
        self.toggleLogButtonY.toggled.connect(partial(self.plot_intervals, associated_form))
        self.eraseAnnotationsB.clicked.connect(partial(self.erase_annotations))
        self.pushButton.clicked.connect(partial(plot_intervals,associated_form,None))
        self.toolbar = NavigationToolbar(self.canvas, self.navigation_frame)
        self.navigation_frame.setMinimumHeight(self.toolbar.height())

        self.retranslateUi(parent)
        QtCore.QMetaObject.connectSlotsByName(parent)

    def retranslateUi(self, parent):
        parent.setWindowTitle(QtGui.QApplication.translate("parent", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        parent.setTitle(QtGui.QApplication.translate("parent", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("parent", "Displayed", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(
            QtGui.QApplication.translate("parent", "Available", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("parent", "Plot", None, QtGui.QApplication.UnicodeUTF8))
        self.toggleLogButtonX.setText(
            QtGui.QApplication.translate("parent", "-log10(x) - horizontal", None, QtGui.QApplication.UnicodeUTF8))
        self.toggleLogButtonY.setText(
            QtGui.QApplication.translate("parent", "-log10(y) - vertical", None, QtGui.QApplication.UnicodeUTF8))
        self.eraseAnnotationsB.setText(
            QtGui.QApplication.translate("parent", "Erase annotations", None, QtGui.QApplication.UnicodeUTF8))

    def plot_intervals(self, form, toggled):
        self.erase_annotations()
        plot_intervals(form)

    def erase_annotations(self, text_list=None):
        if text_list is None:
            text_list = [x.get_text() for x in self.figure.texts]
        for text in text_list:
            indexes_of_text = \
                np.where([x.get_text().find(text)>=0 for x in self.figure.texts])[0]
            if len(indexes_of_text)>1:
                l = self.figure.texts.pop(indexes_of_text[0].item())
                del l
            else:
                l = self.figure.texts.pop(indexes_of_text.item())
                del l
        self.canvas.draw()


    def move_to_available(self, item):
        if self.listWidget_2.count() <= 1:
            return
        name = item.text()
        new_item = QtGui.QListWidgetItem(name, self.listWidget)
        new_item.setIcon(QtGui.QIcon(os.path.join(sys.path[0],
                                                  *['utils', 'glyphicons-601-chevron-up.png'])))
        self.listWidget_2.takeItem(self.listWidget_2.currentRow())
        associated_annotations = \
            [self.figure.texts[x].get_text() for x in
             np.where([x.get_text().find(name)>=0 for x in self.figure.texts])[0]]
        self.erase_annotations(associated_annotations)
        l = self.ax.lines.pop(np.where(
            [x.properties()['label'].find(name) >= 0 for x in self.ax.lines])[0].item())
        del l
        self.ax.relim()
        y_lim = self.ax.get_ylim()
        y_min = y_lim[0]
        y_max = y_min
        for line in self.ax.findobj(lambda x: \
                                                                    x.properties()['visible'] == True and \
                                                                    type(x) == matplotlib.lines.Line2D and \
                                                            len(x.properties()['label']) > 0):
            y_max = max([y_max, max(line.get_ydata())])
        if y_min == y_max:
            pass
        else:
            self.ax.set_ylim(y_min, y_max * 1.05)
        if len(self.ax.lines) > 0:
            self.ax.legend(loc='best', fancybox=True, borderaxespad=0., framealpha=0.5).draggable(True)
        else:
            del self.ax.legend
        self.canvas.draw()

    def move_to_displayed(self, item):
        name = item.text()
        new_item = QtGui.QListWidgetItem(name, self.listWidget_2)
        new_item.setIcon(QtGui.QIcon(os.path.join(sys.path[0],
                                                  *['utils', 'glyphicons-602-chevron-down.png'])))
        self.listWidget.takeItem(self.listWidget.currentRow())
        plot_intervals(self.main_form, [name])
        self.ax.relim()
        y_lim = self.ax.get_ylim()
        y_min = y_lim[0]
        y_max = y_min
        for line in self.ax.findobj(lambda x: \
                                                                    x.properties()['visible'] == True and \
                                                                    type(x) == matplotlib.lines.Line2D and \
                                                            len(x.properties()['label']) > 0):
            y_max = max([y_max, max(line.get_ydata())])
        if y_min == y_max:
            pass
        else:
            self.ax.set_ylim(y_min, y_max * 1.05)
        if len(self.ax.lines) > 0:
            self.ax.legend(loc='best', fancybox=True, borderaxespad=0., framealpha=0.5).draggable(True)
        else:
            del self.ax.legend
        self.canvas.draw()


def open_file(form):
    (filename, _) = \
        QtGui.QFileDialog.getOpenFileName(None,
                                          caption='Open file',
                                          dir=os.path.join(sys.path[0], 'DATA'),
                                          filter='*.csv')
    if os.path.isfile(filename):
        # Reset solution state and order of items
        del form.acceptable_solution
        if hasattr(form, 'component_order_in_table'):
            delattr(form, 'component_order_in_table')
        # Load csv data into form variables
        load_csv(form, filename)

        # Continue with typical solution and table population procedure
        gui_equilibrate(form)
        form.tableComps.sortByColumn(0, QtCore.Qt.AscendingOrder)
        form.tableReacs.sortByColumn(0, QtCore.Qt.AscendingOrder)


def load_csv(form, filename):
    with open(filename) as csv_file:
        n = 0
        Nr = 0
        header_comps = []
        header_reacs = []
        reader = csv.reader(csv_file, dialect='excel')
        reading_comps = False
        reading_reacs = False
        comps = []
        reacs = []
        for row in reader:
            row_without_blanks = filter(lambda x: len(x.replace(' ', '')) > 0, row)
            if len(row_without_blanks) == 0:
                next(reader) #skip empty line
            elif 'COMP' in row:
                reading_comps = True
                reading_reacs = False
                header_comps = next(reader)[0:4]
            elif 'REAC' in row:
                reading_reacs = True
                reading_comps = False
                header_reacs = next(reader)
                header_reacs = filter(lambda x: len(x.replace(' ', '')) > 0, header_reacs)
            elif reading_comps:
                n += 1
                comps.append(map(lambda x: '0' if x=='' or x==' ' else x, row_without_blanks))
            elif reading_reacs:
                Nr += 1
                reacs.append(map(lambda x: '0' if x=='' or x==' ' else x, row[0 : n + 2]))
        csv_file.close()
    comps = np.array(comps)
    reacs = np.array(reacs)

    form.spinBox.setProperty("value", n)
    form.spinBox_2.setProperty("value", Nr)
    form.tableComps.setRowCount(n)
    form.tableComps.setColumnCount(len(header_comps) + 3)
    form.tableComps.setHorizontalHeaderLabels(
        header_comps + ['Ceq_i, mol/L', '-log10(C0_i)', '-log10(Ceq_i)'])

    form.tableReacs.setRowCount(Nr)
    form.tableReacs.setColumnCount(n + 2 + 1)
    form.tableReacs.setHorizontalHeaderLabels(
        header_reacs + ['Xieq_j'])

    # Pass variables to form before loop start
    variables_to_pass = ['header_comps', 'comps', 'header_reacs', 'reacs',
                         'n', 'Nr']
    for var in variables_to_pass:
        setattr(form, var, locals()[var])


def load_variables_from_form(form):
    n = form.n
    Nr = form.Nr
    comps = np.empty([n, 4], dtype='S50')
    reacs = np.empty([Nr, n + 2], dtype='S50')
    header_comps = []
    header_reacs = []
    component_order_in_table = []
    for i in range(form.tableComps.rowCount()):
        component_order_in_table.append(int(form.tableComps.item(i, 0).text()) - 1)
    for i in range(comps.shape[1]):
        header_comps.append(form.tableComps.horizontalHeaderItem(i).data(0))
    for i in range(reacs.shape[1]):
        header_reacs.append(form.tableReacs.horizontalHeaderItem(i).data(0))
    for j in range(comps.shape[1]):
        for i in range(comps.shape[0]):
            comps[component_order_in_table[i], j] = form.tableComps.item(i, j).text()
    for j in range(reacs.shape[1]):
        for i in range(reacs.shape[0]):
            reacs[i, j] = form.tableReacs.item(i, j).text()
    # Pass variables to form before loop start
    variables_to_pass = ['header_comps', 'comps', 'header_reacs', 'reacs',
                         'component_order_in_table']
    for var in variables_to_pass:
        setattr(form, var, locals()[var])


def gui_setup_and_variables(form):
    # Collect variables
    n = form.n
    Nr = form.Nr
    comps = form.comps
    reacs = form.reacs
    # Gui setup with calculated values

    C0_i = np.matrix([row[3] for row in comps], dtype=float).T
    highestC0Indexes = np.argpartition(C0_i.A1, (-1, -2))
    index_of_solvent = highestC0Indexes[-1]
    C_solvent_Tref = C0_i[index_of_solvent].item()
    if len(C0_i) > 1:
        index_of_second_highest_C0 = highestC0Indexes[-2]
    else:
        index_of_second_highest_C0 = highestC0Indexes[-1]
    C_second_highest_C0_Tref = C0_i[index_of_second_highest_C0].item()

    form.C0_i = C0_i
    form.index_of_solvent = index_of_solvent
    form.C_solvent_Tref = C_solvent_Tref
    form.z_i = np.matrix([row[2] for row in comps], dtype=float).T
    form.nu_ij = np.matrix([row[2:2 + n] for row in reacs], dtype=int).T
    form.pKa_j = np.matrix([row[1] for row in reacs], dtype=float).T
    form.max_it = int(form.spinBox_3.value())
    form.tol = float(form.doubleSpinBox_5.value())
    form.C_second_highest_C0_Tref = C_second_highest_C0_Tref

    form.comboBox.clear()
    form.comboBox_3.clear()
    for item in comps[:, 0:2]:
        form.comboBox.addItem('C0_' + item[0] + ' {' + item[1] + '}')
        form.comboBox_3.addItem(item[1])
    form.comboBox.setCurrentIndex(index_of_second_highest_C0)
    form.doubleSpinBox.setValue(C_second_highest_C0_Tref / 10.0 ** 7)
    form.doubleSpinBox_2.setValue(C_second_highest_C0_Tref * (1 + 20 / 100.0))
    form.comboBox_3.setCurrentIndex(index_of_solvent)
    form.doubleSpinBox_6.setValue(C_solvent_Tref)
    form.doubleSpinBox_6.setPrefix('(mol/L)')


def retabulate(form):
    # Collect variables
    n = form.n
    Nr = form.Nr
    C0_i = form.C0_i
    comps = form.comps
    reacs = form.reacs
    header_comps = form.header_comps
    header_reacs = form.header_reacs
    Ceq_i = form.Ceq_i
    Xieq_j = form.Xieq_j
    if hasattr(form, 'component_order_in_table'):
        i = getattr(form, 'component_order_in_table')
    else:
        i = range(0, n)
    j = range(0, 4 + 3)

    form.tableComps.blockSignals(True)
    form.tableReacs.blockSignals(True)
    form.comboBox.blockSignals(True)
    # As usual, problems occurr when sorting is combined with setting QTableWidgetItems.
    # Therefore disable sorting, then set QTableWidgetItems and finally reenable sorting.
    form.tableComps.setSortingEnabled(False)
    form.tableReacs.setSortingEnabled(False)

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
    j = range(0, n + 2 + 1)

    for column in j:
        for row in i:
            if column != n + 2:
                form.tableReacs.setItem(row, column, NSortableTableWidgetItem(str(reacs[row][column])))
            elif column == n + 2:
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
    form.comboBox.blockSignals(False)


def save_file(form):
    pass


def solve_intervals(form):
    variables_to_check = ['Ceq_series', 'Xieq_series', 'indep_var_series', 'dep_var_series',
                          'index_of_variable']
    for var in variables_to_check:
        if hasattr(form, var):
            delattr(form, var)
    comps = form.comps
    reacs = form.reacs
    C0_i = form.C0_i
    n = form.n
    Nr = form.Nr
    index_of_variable = form.comboBox.currentIndex()
    indep_var_label = 'C0_{' + comps[index_of_variable, 0] + ', ' + \
                      comps[index_of_variable, 1] + '}/(mol/L)'
    dep_var_labels = \
        ['Ceq_' + '{' + item[0] + ', ' + item[1] + '}/(mol/L)' for item in comps[:, 0:2]] + \
        ['\\xi eq_' + '{' + str(item) + '}/(mol/L)' for item in range(1, Nr + 1, 1)]
    min_value = form.doubleSpinBox.value()
    max_value = form.doubleSpinBox_2.value()
    n_points = 20
    indep_var_series = \
        [min_value + x
         for x in np.arange(n_points + 1) * (max_value - min_value) / (n_points)]
    C0_variable_comp = C0_i[index_of_variable]
    mid_index = bisect.bisect(indep_var_series, C0_variable_comp) - 1
    Xieq_j = form.Xieq_j
    Ceq_i = form.Ceq_i
    Ceq_series = np.matrix(np.zeros([n_points + 1, len(Ceq_i)]))
    Xieq_series = np.matrix(np.zeros([n_points + 1, len(Xieq_j)]))
    dep_var_series = dict(zip(dep_var_labels, np.empty(n + Nr, dtype=np.ndarray)))
    # Keep current solution intact for after plotting range
    form.stored_solution_Ceq_i = form.Ceq_i
    form.stored_solution_Xieq_j = form.Xieq_j
    for j in range(mid_index, -1, -1):
        form.C0_i[index_of_variable] = indep_var_series[j]
        equilibrate(form)
        Ceq_series[j, :] = form.Ceq_i.T
        Xieq_series[j, :] = form.Xieq_j.T
    form.Ceq_i = form.stored_solution_Ceq_i
    form.Xieq_j = form.stored_solution_Xieq_j
    for j in range(mid_index + 1, n_points + 1, +1):
        form.C0_i[index_of_variable] = indep_var_series[j]
        equilibrate(form)
        Ceq_series[j, :] = form.Ceq_i.T
        Xieq_series[j, :] = form.Xieq_j.T
    for j in range(n):
        dep_var_series[dep_var_labels[j]] = Ceq_series[:, j]
    for j in range(Nr):
        dep_var_series[dep_var_labels[n + j]] = Xieq_series[:, j]
    form.Ceq_i = form.stored_solution_Ceq_i
    form.Xieq_j = form.stored_solution_Xieq_j
    form.Ceq_series = Ceq_series
    form.Xieq_series = Xieq_series
    form.indep_var_series = indep_var_series
    form.dep_var_series = dep_var_series
    form.dep_var_labels = dep_var_labels
    form.indep_var_label = indep_var_label
    form.index_of_variable = index_of_variable
    form.groupBox = QtGui.QGroupBox()
    form.groupBox.plotBox = UiGroupBoxPlot(parent=form.groupBox, associated_form=form)
    form.groupBox.plotBox.main_form = form
    form.groupBox.show()
    initiate_plot(form)
    plot_intervals(form, None)


def initiate_plot(form):
    n = form.n
    Nr = form.Nr
    dep_var_labels = form.dep_var_labels
    form.groupBox.plotBox.ax.clear()
    form.groupBox.plotBox.listWidget.clear()
    form.groupBox.plotBox.listWidget_2.clear()
    for label in dep_var_labels:
        if label.find('Ceq') >= 0:
            new_item = QtGui.QListWidgetItem(label, form.groupBox.plotBox.listWidget_2)
            new_item.setIcon(QtGui.QIcon(os.path.join(sys.path[0],
                                                      *['utils', 'glyphicons-602-chevron-down.png'])))
        else:
            new_item = QtGui.QListWidgetItem(label, form.groupBox.plotBox.listWidget)
            new_item.setIcon(QtGui.QIcon(os.path.join(sys.path[0],
                                                      *['utils', 'glyphicons-601-chevron-up.png'])))
    # dict, keys:ceq_labels; bindings: plottedseries
    form.groupBox.plotBox.plotted_series = dict(zip(dep_var_labels, np.empty(n + Nr, dtype=object)))
    form.groupBox.plotBox.dc = dict(zip(dep_var_labels, np.empty(n + Nr, dtype=object)))


def plot_intervals(form, item_texts=None):
    form.groupBox.plotBox.erase_annotations()
    dep_var_series = form.dep_var_series
    dep_var_labels = form.dep_var_labels
    indep_var_label = form.indep_var_label
    indep_var_series = form.indep_var_series
    colormap_colors = colormaps.viridis.colors + colormaps.inferno.colors
    markers = matplotlib.markers.MarkerStyle.filled_markers
    fillstyles = matplotlib.markers.MarkerStyle.fillstyles
    dc = form.groupBox.plotBox.dc
    plotted_series = form.groupBox.plotBox.plotted_series
    if not form.groupBox.plotBox.toggleLogButtonX.isChecked():
        indep_var_series = indep_var_series
        indep_var_label = '$' + indep_var_label + '$'
    else:
        indep_var_series = -np.log10(indep_var_series)
        indep_var_label = '$-log10(' + indep_var_label + ')$'
    if item_texts is None:
        item_texts = []
        for i in range(form.groupBox.plotBox.listWidget_2.count()):
            item_texts.append(form.groupBox.plotBox.listWidget_2.item(i).text())
        form.groupBox.plotBox.ax.clear()
    for label in item_texts:
        if not form.groupBox.plotBox.toggleLogButtonY.isChecked():
            dep_var_values = dep_var_series[label]
            series_label = '$' + label + '$'
        else:
            dep_var_values = -np.log10(dep_var_series[label])
            series_label = '$-log10(' + label + ')$'
        plotted_series[label] = form.groupBox.plotBox.ax.plot(
            indep_var_series, dep_var_values.A1.tolist(), 'go-', label=series_label,
            color=colormap_colors[np.random.randint(0, len(colormap_colors), 1)],
            markerfacecolor=colormap_colors[np.random.randint(0, len(colormap_colors), 1)],
            marker=markers[np.random.randint(0, len(markers) - 1)],
            fillstyle=fillstyles[np.random.randint(0, len(fillstyles) - 1)])
        dc[label] = datacursor(plotted_series[label], draggable=True, display='multiple',
                               arrowprops=dict(arrowstyle='simple', fc='white', alpha=0.5),
                               bbox=dict(fc='white', alpha=0.5),
                               formatter='x: {x:0.3g},y: {y:0.3g}\n{label}'.format)
    form.groupBox.plotBox.ax.legend(
        loc='best', fancybox=True, borderaxespad=0., framealpha=0.5).draggable(True)
    form.groupBox.plotBox.plotted_series = plotted_series
    form.groupBox.plotBox.dc = dc
    form.groupBox.plotBox.ax.set_xlabel(indep_var_label, fontsize=14)
    form.groupBox.plotBox.listWidget_2.setMinimumWidth(
        form.groupBox.plotBox.listWidget_2.sizeHintForColumn(0))
    form.groupBox.plotBox.canvas.draw()


def recalculate_after_cell_edit(form, row, column):
    load_variables_from_form(form)
    gui_equilibrate(form)


def gui_equilibrate(form):
    form.tableComps.blockSignals(True)
    form.tableReacs.blockSignals(True)
    form.comboBox.blockSignals(True)
    gui_setup_and_variables(form)
    equilibrate(form)
    retabulate(form)

    form.tableComps.blockSignals(False)
    form.tableReacs.blockSignals(False)
    form.comboBox.blockSignals(False)


def equilibrate(form):
    # Collect variables
    n = form.n
    Nr = form.Nr
    C0_i = form.C0_i
    z_i = form.z_i
    comps = form.comps
    reacs = form.reacs
    nu_ij = form.nu_ij
    pKa_j = form.pKa_j
    max_it = form.max_it
    tol = form.tol
    C_solvent_Tref = form.C_solvent_Tref
    index_of_solvent = form.index_of_solvent

    # Init. calculations
    Kc_j = np.multiply(np.power(10, -pKa_j), np.power(C_solvent_Tref, nu_ij[index_of_solvent, :]).T)
    form.Kc_j = Kc_j

    # Setup logging
    if not os.path.exists('./logs'):
        os.mkdir('./logs')
    logging.basicConfig(filename='./logs/calculation_results.log', level=logging.DEBUG,
                        format='%(asctime)s;%(message)s')

    # First estimates for eq. Composition Ceq and Reaction extent Xieq
    if not hasattr(form, 'acceptable_solution'):
        Ceq_i_0 = C0_i
        # replace 0 by 10^-6*smallest value: Smith, Missen 1988 DOI: 10.1002/cjce.5450660409
        Ceq_i_0[C0_i == 0] = min(C0_i[C0_i != 0].A1) * 10 ** -6
        Xieq_j_0 = np.matrix(np.zeros([Nr, 1]))
    else:
        # Use previous solution as initial estimate, if it was valid.
        Ceq_i_0 = form.Ceq_i_0
        Xieq_j_0 = form.Xieq_j_0

    # Pass variables to form before loop start
    variables_to_pass = ['C0_i', 'z_i', 'nu_ij', 'pKa_j',
                         'max_it', 'tol',
                         'Ceq_i_0', 'Xieq_j_0']
    for var in variables_to_pass:
        setattr(form, var, locals()[var])

    k = 1
    stop = False
    form.cancelButton.setEnabled(True)
    form.progressBar.setEnabled(True)
    form.remove_canceled_status()
    form.acceptable_solution = False
    form.initialEstimateAttempts = 1
    form.methodLoops = [0, 0]  # loop numbers: [Line search, Newton]
    # Calculate equilibrium composition: Newton method
    # TODO: Implement global homotopy-continuation method
    while not form.acceptable_solution \
            and k < max_it and stop == False \
            and not form.was_canceled():
        Ceq_i, Xieq_j = calc_Xieq(form)
        k += 1
        # TODO: if progressBar.wasCanceled() == True then stop
        if all(Ceq_i >= 0) and not any(np.isnan(Ceq_i)):
            form.acceptable_solution = True
        else:
            # Set reactions to random extent and recalculate
            # TODO: scale to concentration sizes
            form.Xieq_j_0 = np.matrix(np.random.normal(0.0, 1.0 / 3.0, Nr)).T
            # Set aequilibrium composition to initial value + estimated conversion
            form.Ceq_i_0 = C0_i  # + nu_ij * form.Xieq_j_0
            # replace 0 by 10^-6*smallest value: Smith, Missen 1988 DOI: 10.1002/cjce.5450660409
            form.Ceq_i_0[form.Ceq_i_0 == 0] = min(C0_i[C0_i != 0].A1) * 10 ** -6
            form.initialEstimateAttempts += 1
            form.methodLoops = [0, 0]

    if not form.acceptable_solution:
        delattr(form, 'acceptable_solution')
        form.label_9.setText(form.label_9.text() + '\n')
    else:
        form.Ceq_i_0 = Ceq_i
        form.Xieq_j_0 = Xieq_j
        form.label_9.setText(form.label_9.text() +
                             '\nsum(C0*z_i) = ' + str((z_i.T * C0_i).item()) +
                             ' \t\t\t sum(Ceq_i*z_i) = ' + str((z_i.T * Ceq_i).item()) +
                             '\nI_0 = ' + str((1/2.0*np.power(z_i,2).T * C0_i).item()) +
                             '\t\t\t\t I_eq = ' + str((1/2.0*np.power(z_i,2).T * Ceq_i).item()))

    form.Ceq_i = Ceq_i
    form.Xieq_j = Xieq_j
    form.cancelButton.setEnabled(False)
    form.progressBar.setEnabled(False)


def calc_Xieq(form):
    """Steepest descent for good initial estimate, then Newton method for non-linear algebraic system
    :return: tuple with Ceq_i, Xieq_j, f_0
    :param C0_i: np.matrix (n X 1) - Conc(i, alimentación)
    :param z_i: np.matrix (n X 1) - Carga(i, alimentación)
    :param nu_ij: np.matrix (n X Nr) - Coefs. esteq. componente i en reacción j
    :param pKa_j: np.matrix (n X 1) - (-1)*log10("Cte." de equilibrio en reacción j) = -log10 Kc_j(T)
    :param Xieq_j_0: np.matrix (n X 1) - avance de reacción j - estimado inicial
    :param Ceq_i_0: np.matrix (n X 1) - Conc(i, equilibrio)
    """
    n = form.n
    Nr = form.Nr
    C0_i = form.C0_i
    pKa_j = form.pKa_j
    nu_ij = form.nu_ij
    Ceq_i_0 = form.Ceq_i_0
    Xieq_j_0 = form.Xieq_j_0
    max_it = form.max_it
    tol = form.tol
    z_i = form.z_i
    Kc_j = form.Kc_j

    f = lambda x: f_gl_0(x, C0_i, nu_ij, n, Nr, Kc_j)
    j = lambda x: jac(x, C0_i, nu_ij, n, Nr, Kc_j)
    g = lambda x: g_vec(x, C0_i, nu_ij, n, Nr, Kc_j)

    X0 = np.concatenate([Ceq_i_0, Xieq_j_0])
    X = X0
    # Newton method: G(X) = J(X)^-1 * F(X)
    k = 0
    J_val = j(X)
    F_val = f(X)
    Y = np.matrix(np.ones(len(X))).T * tol / (np.sqrt(len(X)) * tol)
    magnitude_Y = np.sqrt((Y.T * Y).item())
    # For progress bar, use log scale to compensate for quadratic convergence
    log10_to_o_max_magnityde_y = np.log10(tol / magnitude_Y)
    progress_k = \
        (1.0 - np.log10(tol / magnitude_Y) / log10_to_o_max_magnityde_y) * 100.0
    diff = np.matrix(np.empty([len(X), 1]))
    diff.fill(np.nan)
    stop = False
    divergent = False
    # Add progress bar & variable
    form.progressBar.setValue(0)
    update_status_label(form, k, stop)
    lambda_ls = 1.0
    j_it = 1
    while k <= max_it and not stop:
        new_log_entry('Newton', k, X, diff, F_val, Y, np.nan, np.nan, stop)
        X_k_m_1 = X
        progress_k_m_1 = progress_k
        Y = gausselimination(J_val, -F_val)
        magnitude_Y = np.sqrt((Y.T * Y).item())
        # First attempt without backtracking
        X = X + 1.0 * Y
        diff = X - X_k_m_1
        if magnitude_Y < tol and all(X[0:n] >= 0):
            stop = True  # Procedure successful
            form.progressBar.setValue(100.0)
        else:
            # For progress bar, use log scale to compensate for quadratic convergence
            update_status_label(form, k, stop)
            progress_k = \
                (1.0 - np.log10(tol / magnitude_Y) / log10_to_o_max_magnityde_y) * 100.0
            # TODO: Fix case in which magnitude_Y == inf (divergent)
            if np.isnan(magnitude_Y) or np.isinf(magnitude_Y):
                stop = True  # Divergent method
                divergent = True
                form.progressBar.setValue(
                    (1.0 - np.log10(np.finfo(float).eps) / log10_to_o_max_magnityde_y) * 100.0)
            else:
                form.progressBar.setValue(
                    (1.0 - np.log10(tol / magnitude_Y) / log10_to_o_max_magnityde_y) * 100.0)
            if round(progress_k) == round(progress_k_m_1):
                QtGui.QApplication.processEvents()
                # if form.progressBar.wasCanceled():
                # stop = True
        while j_it <= max_it and not all(X[0:n] >= 0):
            # Backtrack if any conc < 0. Line search method.
            # Ref. http://dx.doi.org/10.1016/j.compchemeng.2013.06.013
            lambda_ls = lambda_ls / 2.0
            X = X_k_m_1
            progress_k = progress_k_m_1
            X = X + lambda_ls * Y
            diff = X - X_k_m_1
            new_log_entry('Newton it. ' + str(k + 1) + ' - lambda: ' + str(lambda_ls) + ' backtrack: k',
                          j_it, X, diff, F_val, lambda_ls * Y, np.nan, np.nan, stop)
            j_it += 1
            form.methodLoops[0] += 1
            update_status_label(form, k, stop)
        j_it = 1
        lambda_ls = 1.0
        J_val = j(X)
        F_val = f(X)
        k += 1
        form.methodLoops[1] += 1
    new_log_entry('Newton', k, X, diff, F_val, Y, np.nan, np.nan, stop and not divergent)
    update_status_label(form, k, stop and not divergent)
    Ceq_i = X[0:n]
    Xieq_j = X[n:n + Nr]
    return Ceq_i, Xieq_j


def f_gl_0(X, C0_i, nu_ij, n, Nr, Kc_j):
    Ceq_i = X[0:n, 0]
    Xieq_j = X[n:n + Nr, 0]
    f_gl_0 = np.matrix(np.empty([n + Nr, 1], dtype=float))
    f_gl_0[0:n] = -Ceq_i + C0_i + nu_ij * Xieq_j
    f_gl_0[n:n + Nr] = -Kc_j + np.prod(np.power(Ceq_i, nu_ij), 0).T
    return f_gl_0


def jac(X, C0_i, nu_ij, n, Nr, Kc_j):
    Ceq_i = X[0:n, 0]
    eins_durch_c = np.diag(np.power(Ceq_i, -1).A1, 0)
    quotient = np.diag(np.prod(np.power(Ceq_i, nu_ij), 0).A1)
    jac = np.matrix(np.zeros([n + Nr, n + Nr], dtype=float))
    jac[0:n, 0:n] = -1 * np.eye(n).astype(float)
    jac[0:n, n:n + Nr] = nu_ij
    jac[n:n + Nr, 0:n] = quotient * nu_ij.T * eins_durch_c
    return jac


def g_vec(X, C0_i, nu_ij, n, Nr, Kc_j):
    f_0 = f_gl_0(X, C0_i, nu_ij, n, Nr, Kc_j)
    return (f_0.T * f_0).item()


def steepest_descent(X0, f, J, g, tol, form):
    # Skipping this method in favor of line search with Newton
    max_it = form.max_it
    X = X0
    f_val = f(X)
    k = 0
    stop = False
    diff = np.matrix(np.empty([len(X), 1]))
    diff.fill(np.nan)
    Y = np.matrix(np.empty([len(X), 1]))
    Y.fill(np.nan)
    abs_gmin_minus_g1 = np.nan
    progressBar = form.progressBar
    label_9 = form.label_9
    progressBar.setValue(0)
    update_status_label(form, k, False)
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
        abs_gmin_minus_g1 = abs(g_min - g1)
        if abs_gmin_minus_g1 < tol:
            stop = True  # Procedure successful
            progressBar.setValue(100.0)
        else:
            progressBar.setValue(tol / abs_gmin_minus_g1 * 100)
            update_status_label(form, k, False)
        new_log_entry('Steepest descent', k, X, diff, f_val, Y, g_min, g1, stop)
        update_status_label(form, k, False)
        X_k_m_1 = X
        X = X - alpha * z
        diff = X_k_m_1 - X
        f_val = f(X)
        k += 1
        form.methodLoops[0] += 1
    new_log_entry('Steepest descent', k, X, diff, f_val, Y, g_min, g1, stop)
    return X, f_val


def update_status_label(form, k, solved):
    if solved:
        solved = 'solved.'
    else:
        solved = 'solution not found.'

    form.label_9.setText('Loops: Newton \t' +
                         str(form.methodLoops[1]) + ' \t Line search (backtrack) \t' +
                         str(form.methodLoops[0]) +
                         ' \t Initial estimate attempts \t' +
                         str(form.initialEstimateAttempts) + '\n' +
                         'Iteration (k) \t' + str(k) +
                         '\n' + str(solved))


def new_log_entry(method, k, X, diff, f_val, Y, g_min, g1, stop):
    logging.debug(method + ' ' +
                  ';k=' + str(k) +
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
    rowString = unicode('', 'utf_8')
    form.aboutBox_1 = QtGui.QTextBrowser()
    form.aboutBox_1.setWindowTitle('About')
    form.aboutBox_1.setWindowIcon(QtGui.QIcon(
        os.path.join(sys.path[0], *['utils', 'icon_batch.png'])))
    form.aboutBox_1.setOpenExternalLinks(True)
    addingTable = False

    htmlStream = \
        unicode('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN" "http://www.w3.org/TR/REC-html40/strict.dtd">\n',
                'utf_8')
    htmlStream += unicode('<html>', 'utf_8')
    htmlStream += unicode('<head><meta name="qrichtext" content="1" /><style type="text/css">' +
                          '\\np,li {white-space: pre-wrap;}\n' +
                          '\\np,br {line-height: 10%;}\n' +
                          '</style></head>', 'utf_8')

    htmlStream += unicode('<body style=' + '"' + ' font-family:' + "'" + 'MS Shell Dlg 2' + "'" +
                          '; font-size:8.25pt; font-weight:400; font-style:normal;' + '"' + '>', 'utf_8')
    stringToAdd = unicode('', 'utf_8')
    startingP = unicode("<p style=' margin-top:0px; margin-bottom:0px; margin-left:0px;" +
                        "margin-right:0px; -qt-block-indent:0; text-indent:0px;'>", 'utf_8')
    endingP = unicode('</p>', 'utf_8')
    matchingHLine = re.compile('=+')
    with open('README.md') as readme_file:
        for row in readme_file:
            rowString = unicode(row, 'utf_8')
            if not addingTable and rowString.find('|') != -1:
                stringToAdd = ''.join(
                    ['<table>', '<tr><td>',
                     rowString.replace('|', '</td><td>').replace('\n', ''),
                     '</td></tr>'])
                addingTable = True
            elif addingTable and rowString.find('|') == -1:
                stringToAdd = '</table>' + startingP + rowString + endingP
                addingTable = False
            elif addingTable and rowString.find('|') != -1:
                stringToAdd = ''.join(
                    ['<tr><td>',
                     rowString.replace('|', '</td><td>').replace('\n', ''),
                     '</td></tr>'])
            elif not addingTable and rowString.find('|') == -1:
                stringToAdd = startingP + rowString + endingP
            if len(rowString.replace('\n', '')) == 0:
                htmlStream += stringToAdd + '<br>'
            elif matchingHLine.match(rowString):
                htmlStream += '<hr />'
            else:
                htmlStream += stringToAdd
    htmlStream += unicode('<hr />', 'utf_8')
    htmlStream += unicode("<footer><p>" +
                          '<a href=' + '"' + 'https://github.com/santiago-salas-v/lit-impl-py' + '"' + '>' +
                          'https://github.com/santiago-salas-v/lit-impl-py</a></p></footer>',
                          'utf_8')
    htmlStream += unicode('</body></html>', 'utf_8')
    readme_file.close()
    form.aboutBox_1.setHtml(htmlStream)
    form.aboutBox_1.setMinimumWidth(500)
    form.aboutBox_1.setMinimumHeight(400)
    form.aboutBox_1.show()


def show_log():
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
    return LogWidget(log)


class LogWidget(QtGui.QWidget):
    def __init__(self, _log, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.icon = QtGui.QIcon()
        self.icon.addPixmap(QtGui.QPixmap(
            _fromUtf8("utils/glyphicons-88-log-book.png")),
            QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.log = _log
        self.setupUi()
        self.group_2.setWindowIcon(self.icon)
        self.group_2.show()

    def setupUi(self):
        self.minLogHeight = 400
        self.minLogWidth = self.minLogHeight * 16 / 9
        self.group_2 = QtGui.QGroupBox()
        self.group_2.resize(self.minLogWidth, self.minLogHeight)
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
        self.exportButton = QtGui.QPushButton(self.group_2)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout.addWidget(self.pandasView)
        self.verticalLayout.addWidget(self.exportButton)
        self.horizontalLayout.addWidget(self.firstButton)
        self.horizontalLayout.addWidget(self.previousButton)
        self.horizontalLayout.addWidget(self.pageLabel)
        self.horizontalLayout.addWidget(self.pageBox)
        self.horizontalLayout.addWidget(self.totPagesLabel)
        self.horizontalLayout.addWidget(self.nextButton)
        self.horizontalLayout.addWidget(self.lastButton)

        self.firstButton.clicked.connect(partial(self.firstPage))
        self.lastButton.clicked.connect(partial(self.lastPage))
        self.nextButton.clicked.connect(partial(self.nextPage))
        self.previousButton.clicked.connect(partial(self.previousPage))
        self.exportButton.clicked.connect(partial(self.exportData))
        self.pageBox.editingFinished.connect(partial(self.goToPageNo))

        # To ensure full display, first set resize modes, then resize columns to contents
        self.pandasView.horizontalHeader().setResizeMode(QtGui.QHeaderView.Interactive)
        self.pandasView.verticalHeader().setResizeMode(QtGui.QHeaderView.ResizeToContents)
        self.pandasView.setWordWrap(True)
        self.pandasView.resizeColumnsToContents()

        self.group_2.setWindowTitle('calculation log')
        self.firstButton.setText('<< First')
        self.lastButton.setText('Last >>')
        self.previousButton.setText('< Previous')
        self.nextButton.setText('Next >')
        self.pageBox.setMaximumWidth(int(round(self.minLogHeight / float(5))))
        self.pageBox.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignCenter)
        self.exportButton.setText('Export (csv)')

        self.displayItemsByPage = 50
        self.currentPageFirstEntry = len(self.log.values) \
                                     - self.displayItemsByPage
        self.display()

    def firstPage(self):
        self.currentPageFirstEntry = 1
        self.displayItemsByPage = self.displayItemsByPage
        self.display()

    def lastPage(self):
        self.currentPageFirstEntry = len(self.log.values) \
                                     - self.displayItemsByPage
        self.display()

    def nextPage(self):
        if self.currentPageLastEntry + \
                self.displayItemsByPage > len(self.log.values):
            self.currentPageFirstEntry = len(self.log.values) \
                                         - self.displayItemsByPage
        else:
            self.currentPageFirstEntry = self.currentPageFirstEntry + \
                                         self.displayItemsByPage
        self.display()

    def previousPage(self):
        if self.currentPageFirstEntry - \
                self.displayItemsByPage < 0:
            self.currentPageFirstEntry = 0
        else:
            self.currentPageFirstEntry = self.currentPageFirstEntry - \
                                         self.displayItemsByPage
        self.display()

    def goToPageNo(self):
        pageText = self.pageBox.text()
        try:
            pageNo = int(round(float(pageText)))
            self.lastPage = int(round(len(self.log.values) / \
                                      float(self.displayItemsByPage)))
            if pageNo < self.lastPage:
                self.currentPageFirstEntry = \
                    (pageNo - 1) * self.displayItemsByPage + 1
            elif pageNo >= self.lastPage:
                self.currentPageFirstEntry = len(self.log.values) \
                                             - self.displayItemsByPage
            elif pageNo <= 1:
                self.currentPageFirstEntry = 1
            self.display()
        except ValueError:
            pass

    def display(self):
        # Page number
        self.currentPageLastEntry = self.currentPageFirstEntry + \
                                    self.displayItemsByPage
        self.currentPage = int(round(self.currentPageLastEntry / \
                                     float(self.displayItemsByPage)))
        self.lastPage = int(round(len(self.log.values) / \
                                  float(self.displayItemsByPage)))

        self.pandasModel = PandasModel(
            self.log[self.currentPageFirstEntry: self.currentPageLastEntry])
        self.pandasView.setModel(self.pandasModel)
        self.pandasView.resizeColumnsToContents()
        self.pandasView.verticalHeaders = \
            map(str, range(self.currentPageFirstEntry,
                           self.currentPageLastEntry + 1, 1))
        self.totPagesLabel.setText(' / ' + str(self.lastPage))
        self.pageLabel.setText('Entries ' + str(self.currentPageFirstEntry) +
                               ' to ' + str(self.currentPageLastEntry) + '; Page: ')
        self.pageBox.blockSignals(True)
        self.pageBox.setText(str(self.currentPage))
        self.pageBox.blockSignals(False)

    def exportData(self):
        supportedFilters = ['CSV file (*.csv)']
        # TODO: supportedFilters = ['CSV file (*.csv)', 'XLSX (2010) (*.xlsx)', 'XLS (2007) (*.xls)']
        (fileName, selectedFilter) = QtGui.QFileDialog.getSaveFileName(
            parent=self.group_2,
            caption='enter file name to save...',
            filter=';;'.join(supportedFilters))
        if selectedFilter == supportedFilters[0]:
            self.log.to_csv(fileName)
            # elif selectedFilter == supportedFilters[1] or \
            #                selectedFilter == supportedFilters[2]:
            #    self.log.to_excel(fileName)


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


# Following classes from git@gist.github.com:0be2e44981159d0854f5.git
# Regular expression to find floats. Match groups are the whole string, the
# whole coefficient, the decimal part of the coefficient, and the exponent
# part.
_float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')


def valid_float_string(string):
    match = _float_re.search(string)
    return match.groups()[0] == string if match else False


class FloatValidator(QtGui.QValidator):
    def validate(self, string, position):
        if valid_float_string(string):
            return self.State.Acceptable
        if string == "" or string[position - 1] in 'e.-+':
            return self.State.Intermediate
        return self.State.Invalid

    def fixup(self, text):
        match = _float_re.search(text)
        return match.groups()[0] if match else ""


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
        groups = _float_re.search(text).groups()
        decimal = float(groups[1])
        decimal += steps
        new_string = "{:g}".format(decimal) + (groups[3] if groups[3] else "")
        self.lineEdit().setText(new_string)


def format_float(value):
    """Modified form of the 'g' format specifier."""
    string = "{:g}".format(value).replace("e+", "e")
    string = re.sub("e(-?)0*(\d+)", r"e\1\2", string)
    return string


def main():
    app = QtGui.QApplication.instance()  # checks if QApplication already exists
    if not app:  # create QApplication if it doesnt exist
        app = QtGui.QApplication(sys.argv)

    main_form = QtGui.QGroupBox()
    icon = QtGui.QIcon(
        os.path.join(sys.path[0], *['utils', 'icon_batch.png']))
    main_form.setWindowIcon(icon)
    main_form.ui = UiGroupBox(main_form)
    main_form.show()
    load_csv(main_form.ui, './DATA/COMPONENTS_REACTIONS_EX_001.csv')
    gui_equilibrate(main_form.ui)
    main_form.ui.tableComps.sortByColumn(0, QtCore.Qt.AscendingOrder)
    main_form.ui.tableReacs.sortByColumn(0, QtCore.Qt.AscendingOrder)
    app.exec_()


if __name__ == '__main__':
    main()
