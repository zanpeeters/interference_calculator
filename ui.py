#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" GUI for interference calculator. """

try:
    from PyQt5 import QtCore, QtGui
    from PyQt5 import QtWidgets as widgets
except ImportError:
    try:
        from PyQt4 import QtCore, QtGui
        from PyQt4 import QtGui as widgets
        widgets.QStyleOptionViewItem = widgets.QStyleOptionViewItemV4
    except ImportError:
        raise ImportError('You need to have either PyQt4 or PyQt5 installed.')

try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
except ImportError:
    try:
        from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    except ImportError:
        raise ImportError('You need to have either the qt5agg or the qt4agg matplotlib backend installed.')

import matplotlib as mpl
import numpy as np
import sys, re, platform
from pyparsing import ParseException
from itertools import combinations
from interference_calculator.main import interference, standard_ratio
from interference_calculator.molecule import Molecule, periodic_table
from interference_calculator.ui_help import *

_isotope_rx = re.compile(r'(\d*[A-Z][a-z]{0,2})')
_charges_rx = re.compile(r'(\d+)')

# Qt uses 0-255 ints, Matplotlib uses 0-1 floats for RGB
_red = (193, 24, 78)   #c1184e, fuchsia
_blue = (31, 119, 180) #1f77b4, blue
_blueF = [c/255 for c in _blue]
_redF = [c/255 for c in _red]

mpl.rc('font', family='sans-serif', size=14)

class TableModel(QtCore.QAbstractTableModel):
    """ Take a pandas DataFrame and set data in a QTableModel (read-only). """
    def __init__(self, data, table='interference', parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent=parent)
        self._data = data
        self.table = table

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                if self.table == 'interference':
                    if index.column() == 0:
                        # formula
                        molec = self._data.iloc[index.row(), index.column()]
                        try:
                            m = Molecule(molec)
                        except ParseException:
                            return molec
                        else:
                            return m.formula(style='html', all_isotopes=True)
                    elif index.column() == 1:
                        # mass
                        return '{:.6f}'.format(self._data.iloc[index.row(), index.column()])
                    elif index.column() == 2:
                        # mass difference
                        return '{:.7f}'.format(self._data.iloc[index.row(), index.column()])
                    elif index.column() == 3:
                        # MRP
                        return '{:.2f}'.format(self._data.iloc[index.row(), index.column()])
                    elif index.column() == 4:
                        # probability
                        return '{:.5g}'.format(self._data.iloc[index.row(), index.column()])
                    else:
                        return '{}'.format(self._data.iloc[index.row(), index.column()])
                elif self.table == 'std_ratios':
                    if index.column() == 0:
                        # formula
                        m = Molecule(self._data.iloc[index.row(), index.column()])
                        return m.formula(style='html', all_isotopes=True)
                    elif index.column() == 1:
                        # mass
                        return '{:.6f}'.format(self._data.iloc[index.row(), index.column()])
                    elif index.column() in (2, 3):
                        # rel. abundance and ratio
                        return '{:.5g}'.format(self._data.iloc[index.row(), index.column()])
                    elif index.column() == 4:
                        # inverse ratio
                        return '{:.2f}'.format(self._data.iloc[index.row(), index.column()])
                    else:
                        return '{}'.format(self._data.iloc[index.row(), index.column()])
            elif role == QtCore.Qt.TextAlignmentRole:
                if index.column() == 0:
                    return QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
                else:
                    return QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter
            elif role == QtCore.Qt.EditRole:
                return self._data.iloc[index.row(), index.column()]
            elif role == QtCore.Qt.BackgroundRole:
                if self._data['target'].iloc[index.row()]:
                    return QtGui.QColor(*_red, alpha=32)

    def headerData(self, rowcol, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self._data.columns[rowcol]
        if orientation == QtCore.Qt.Vertical and role == QtCore.Qt.DisplayRole:
            return self._data.index[rowcol]

    def sort(self, column, order):
        # QtCore.Qt.AscendingOrder = 0
        # QtCore.Qt.DescendingOrder = 1
        ascending = not bool(order)
        colname = self._data.columns[column]
        self._data.sort_values(colname, ascending=ascending, inplace=True)
        self.beginResetModel()
        self.endResetModel()

    def copy(self, selection):
        mask = np.zeros(self._data.shape, dtype=bool)
        for s in selection:
            mask[s.row(), s.column()] = True
        output = self._data.where(mask)
        output = output.dropna(how='all', axis=(0,1))
        pasteboard = widgets.QApplication.clipboard()
        pasteboard.setText(output.to_csv(index=False))


class TableView(widgets.QTableView):
    """ Implement a QTableView which can display HTML in arbitrary columns """
    def __init__(self, html_cols=None, parent=None):
        widgets.QTableView.__init__(self, parent=parent)
        self.setSortingEnabled(True)
        if html_cols is not None:
            if isinstance(html_cols, int):
                html_cols = [html_cols]
            [self.setItemDelegateForColumn(c, HTMLDelegate(parent=parent)) for c in html_cols]

    def copy(self):
        self.model().copy(self.selectedIndexes())

    def contextMenuEvent(self, event):
        menu = widgets.QMenu(self)
        copy_action = menu.addAction('Copy')
        copy_action.setShortcut('Ctrl+C')
        select_all_action = menu.addAction('Select All')
        select_all_action.setShortcut('Ctrl+A')
        action = menu.exec_(self.mapToGlobal(event.pos()))

        if action == copy_action:
            self.copy()
        elif action == select_all_action:
            self.selectAll()


class HTMLDelegate(widgets.QStyledItemDelegate):
    """ Display HTML in a table cell. """
    def __init__(self, parent=None):
        widgets.QStyledItemDelegate.__init__(self, parent=parent)

    def createEditor(self, parent, option, index):
        """ disable editing """
        return None

    def paint(self, painter, option, index):
        """ paint QTextDocument """
        options = widgets.QStyleOptionViewItem(option)
        self.initStyleOption(options, index)

        style = widgets.QApplication.style()
        textbox = QtGui.QTextDocument()
        textbox.setHtml(options.text)
        textbox.setTextWidth(options.rect.width())
        options.text = ''
        style.drawControl(widgets.QStyle.CE_ItemViewItem, options, painter)
        context = QtGui.QAbstractTextDocumentLayout.PaintContext()
        textrect = style.subElementRect(widgets.QStyle.SE_ItemViewItemText, options)

        painter.save()
        painter.translate(textrect.topLeft())
        painter.setClipRect(textrect.translated(-textrect.topLeft()))
        textbox.documentLayout().draw(painter, context)
        painter.restore()

    def sizeHint(self, option, index):
        """ Set size hint for HTMLDelegate. """
        options = widgets.QStyleOptionViewItem(option)
        self.initStyleOption(options, index)

        textbox = QtGui.QTextDocument()
        textbox.setHtml(options.text)
        textbox.setTextWidth(options.rect.width())

        return QtCore.QSize(textbox.idealWidth(), textbox.size().height())


class Spectrum(widgets.QWidget):
    def __init__(self, data=None, parent=None):
        widgets.QWidget.__init__(self, parent=parent)
        self.setWindowTitle('Mass interference spectrum')
        self.setWindowFlags(QtCore.Qt.Window)
        self.setAttribute(QtCore.Qt.WA_ShowWithoutActivating)

        self.fig = mpl.figure.Figure(figsize=(720/72,560/72), dpi=72)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.layout = widgets.QVBoxLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.setSpacing(0)
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)

        self._data = None
        self.x = None
        self.y = None
        self.label_offset = (0, 24)
        self.minimum = 0

        self.MAX_LABELS = 15

        self.ax = self.fig.add_subplot(111)

        self.plot_spectrum(data)

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, newdata):
        self._data = newdata.copy().sort_values('probability', ascending=False)
        self.x = self._data['mass/charge'].values
        self.y = self._data['probability'].values

    def plot_spectrum(self, data=None):
        """ Plot the spectrum. """
        if data is not None:
            self.data = data

        if self._data is None:
            return

        self.ax.clear()
        self.ax.set_xlabel('mass/charge')
        self.ax.set_ylabel('probability')
        self.ax.set_title('Mass interference spectrum')
        self.ax.minorticks_on()
        self.ax.grid(True)

        self.colours = [_redF if c else _blueF for c in self.data['target']]

        self.data_points = self.ax.scatter(self.x, self.y, marker='D', c=self.colours)

        self.data_lines = self.ax.vlines(self.x, [self.minimum] * self.data.shape[0],
                                    self.y, colors=self.colours, linewidth=3)

        self.renderer = self.fig.canvas.get_renderer()

        self.labels = []
        for molec, x, y in zip(self.data['molecule'].iloc[:self.MAX_LABELS], self.x, self.y):
            m = Molecule(molec)
            l = m.formula(all_isotopes=True, style='latex')

            label = self.ax.text(x, y, l, horizontalalignment='center')

            # Set initial offset
            self.shift_label(label)
            self.labels.append(label)

        # Iterate over labels until none overlap
        done = False
        while not done:
            done = True
            for label1, label2 in combinations(self.labels, 2):
                box1 = label1.get_window_extent(self.renderer)
                box2 = label2.get_window_extent(self.renderer)
                if box1.intersection(box1, box2):
                    self.shift_label(label2)
                    # There was an overlap, check again.
                    done = False

        # Draw thin lines from datapoints to labels
        for label, x, y in zip(self.labels, self.x, self.y):
            lx, ly = label.get_position()
            self.ax.plot((x,lx), (y, ly), linewidth=0.5, color='black')

        # For log plot
        self.minimum = self.data['probability'].min()/100
        self.ax.autoscale()
        self.ax.set_ybound(lower=self.minimum)
        self.canvas.draw()

    def shift_label(self, label, offset=None):
        """ Shift text label 'label' up by 'offset' amount (in pixels).
            Label must be a matplotlib.text.Text object. Offset defaults
            to Spectrum.label_offset.
        """
        if not offset:
            offset = self.label_offset

        xpos, ypos = label.get_position()
        xpos_px, ypos_px = self.ax.transData.transform((xpos, ypos))
        xpos_px += offset[0]
        ypos_px += offset[1]
        xpos, ypos = self.ax.transData.inverted().transform((xpos_px, ypos_px))
        label.set_position((xpos, ypos))


class MainWindow(widgets.QMainWindow):
    """ Main window for interference calculator ui. """
    def __init__(self):
        widgets.QMainWindow.__init__(self)
        self.setWindowTitle('Mass interference calculator')
        self.resize(600, 720)
        self.setMinimumSize(600, 720)
        self.setCentralWidget(MainWidget(parent=self))


class MainWidget(widgets.QWidget):
    """ Central widget class for interference calculator ui. """
    def __init__(self, parent=None):
        widgets.QWidget.__init__(self, parent=parent)
        self.atoms = []
        self.charges = [1]
        self.mz = ''
        self.mzrange = 0.3
        self.maxsize = 5

        # Inputs
        self.atoms_input = widgets.QLineEdit(parent=self)
        self.atoms_input.setPlaceholderText('list of atoms in sample: C H Si Ca')

        self.maxsize_label = widgets.QLabel('max size', parent=self)
        self.maxsize_input = widgets.QSpinBox(parent=self)
        self.maxsize_input.setValue(self.maxsize)
        self.maxsize_input.setMinimumWidth(80)

        self.charges_label = widgets.QLabel('charge(s)', parent=self)
        self.charges_input = widgets.QLineEdit(' '.join([str(c) for c in self.charges]), parent=self)

        self.chargesign_input = widgets.QComboBox(parent=self)
        self.chargesign_input.addItem('-', '-')
        self.chargesign_input.addItem('o', 'o')
        self.chargesign_input.addItem('+', '+')

        self.mz_input = widgets.QLineEdit(parent=self)
        self.mz_input.setPlaceholderText('number or formula: 26.0 or 12C 14N')

        self.mzrange_label = widgets.QLabel('±', parent=self)
        self.mzrange_input = widgets.QDoubleSpinBox(parent=self)
        self.mzrange_input.setValue(self.mzrange)
        self.mzrange_input.setSingleStep(0.1)
        self.mzrange_input.setDecimals(3)

        # Action button
        self.interference_button = widgets.QPushButton('calculate interference', parent=self)
        self.standard_ratio_button = widgets.QPushButton('standard ratio', parent=self)
        self.help_button = widgets.QPushButton('?', parent=self)
        self.help_button.setFixedSize(20, 20)
        ss = 'background-color: rgb({},{},{}); color: white; border-radius: 10;'
        self.help_button.setStyleSheet(ss.format(*_blue))
        self.spectrum_button = widgets.QPushButton('▶︎', parent=self)
        self.spectrum_button.setFixedSize(20,20)
        ss = 'border: none;'
        self.spectrum_button.setStyleSheet(ss)

        # Table and spectrum output
        self.table_output = TableView(html_cols=0)
        self.spectrum_window = Spectrum(parent=self)

        # Qt4 on Mac Snow Leopard and older has a problem with the focus rectangle
        # around input boxes. (error: unlockFocus called too many times)
        if sys.platform == 'darwin':
            if int(platform.mac_ver()[0].split('.')[1]) < 7 and int(QtCore.QT_VERSION) < 0x50000:
                self.atoms_input.setAttribute(QtCore.Qt.WA_MacShowFocusRect, on=False)
                self.mz_input.setAttribute(QtCore.Qt.WA_MacShowFocusRect, on=False)
                self.mzrange_input.setAttribute(QtCore.Qt.WA_MacShowFocusRect, on=False)
                self.charges_input.setAttribute(QtCore.Qt.WA_MacShowFocusRect, on=False)
                self.maxsize_input.setAttribute(QtCore.Qt.WA_MacShowFocusRect, on=False)
                self.minabundance_input.setAttribute(QtCore.Qt.WA_MacShowFocusRect, on=False)
                self.table_output.setAttribute(QtCore.Qt.WA_MacShowFocusRect, on=False)

        # Show input errors on statusbar
        self.statusbar = self.parent().statusBar()
        self.statusbar.setStyleSheet('color: rgb({},{},{});'.format(*_red))

        # Layout objects
        pol = widgets.QSizePolicy()
        pol.setHorizontalPolicy(widgets.QSizePolicy.Ignored)
        pol.setVerticalPolicy(widgets.QSizePolicy.Maximum)

        self.atoms_group = widgets.QGroupBox('sample composition')
        self.atoms_group.setSizePolicy(pol)
        self.atoms_group_layout = widgets.QHBoxLayout()
        self.atoms_group_layout.addWidget(self.atoms_input)
        self.atoms_group.setLayout(self.atoms_group_layout)

        self.target_group = widgets.QGroupBox('target')
        self.target_group.setSizePolicy(pol)
        self.target_group_layout = widgets.QHBoxLayout()
        self.target_group_layout.addWidget(self.mz_input)
        self.target_group_layout.addWidget(self.mzrange_label)
        self.target_group_layout.addWidget(self.mzrange_input)
        self.target_group.setLayout(self.target_group_layout)

        self.parameters_group = widgets.QGroupBox('interference parameters')
        self.parameters_group.setSizePolicy(pol)
        self.parameters_group_layout = widgets.QHBoxLayout()
        self.parameters_group_layout.addWidget(self.charges_label)
        self.parameters_group_layout.addWidget(self.charges_input)
        self.parameters_group_layout.addWidget(self.chargesign_input)
        self.parameters_group_layout.addSpacing(25)
        self.parameters_group_layout.addWidget(self.maxsize_label)
        self.parameters_group_layout.addWidget(self.maxsize_input)
        self.parameters_group.setLayout(self.parameters_group_layout)

        self.button_layout = widgets.QHBoxLayout()
        self.button_layout.addWidget(self.interference_button)
        self.button_layout.addWidget(self.standard_ratio_button)
        self.button_layout.addWidget(self.help_button)
        self.button_layout.addWidget(self.spectrum_button)

        self.output_layout = widgets.QHBoxLayout()
        self.output_layout.addWidget(self.table_output)

        self.layout = widgets.QVBoxLayout()
        self.layout.addWidget(self.atoms_group)
        self.layout.addSpacing(10)
        self.layout.addWidget(self.target_group)
        self.layout.addSpacing(10)
        self.layout.addWidget(self.parameters_group)
        self.layout.addLayout(self.button_layout)
        self.layout.addLayout(self.output_layout)
        self.setLayout(self.layout)

        # Connect
        self.interference_button.clicked.connect(self.calculate_interference)
        self.standard_ratio_button.clicked.connect(self.show_standard_ratio)
        self.help_button.clicked.connect(self.show_help)
        self.atoms_input.editingFinished.connect(self.check_atoms_input)
        self.charges_input.editingFinished.connect(self.check_charges_input)
        self.mz_input.editingFinished.connect(self.check_mz_input)
        self.spectrum_button.clicked.connect(self.toggle_spectrum)

        # Set jump order for tab
        self.setTabOrder(self.atoms_input, self.mz_input)
        self.setTabOrder(self.mz_input, self.mzrange_input)
        self.setTabOrder(self.mzrange_input, self.charges_input)
        self.setTabOrder(self.charges_input, self.chargesign_input)
        self.setTabOrder(self.chargesign_input, self.maxsize_input)
        self.setTabOrder(self.maxsize_input, self.interference_button)
        self.setTabOrder(self.interference_button, self.standard_ratio_button)
        self.setTabOrder(self.standard_ratio_button, self.help_button)
        self.setTabOrder(self.help_button, self.spectrum_button)
        self.setTabOrder(self.spectrum_button, self.table_output)

        # Set tooltip help
        self.atoms_input.setToolTip(atoms_input_tooltip)
        self.charges_input.setToolTip(charges_input_tooltip)
        self.chargesign_input.setToolTip(chargesign_input_tooltip)
        self.mz_input.setToolTip(mz_input_tooltip)
        self.mzrange_input.setToolTip(mzrange_input_tooltip)
        self.maxsize_input.setToolTip(maxsize_input_tooltip)
        self.interference_button.setToolTip(interference_button_tooltip)
        self.standard_ratio_button.setToolTip(standard_ratio_button_tooltip)
        self.spectrum_button.setToolTip(spectrum_button_tooltip)
        self.help_button.setToolTip(help_button_tooltip)

    def warn(self, text, time=5000):
        """ Display a warning message in the status bar. """
        self.statusbar.showMessage(text, msecs=time)

    def check_atoms_input(self):
        """ Validate input for atoms_input.
            Returns True on proper validation, False on error.
        """
        atoms = re.findall(_isotope_rx, self.atoms_input.text())
        if not atoms:
            self.warn('Enter at least one element or isotope.')
            return False
        for a in atoms:
            if not (periodic_table['element'] == a).any(): 
                self.warn('{} is not an element or missing from the periodic table.'.format(a))
                return False
        self.atoms = atoms
        return True

    def check_charges_input(self):
        """ Validate input for charges_input.
            Returns True on correct input, False on error.
        """
        charges = re.findall(_charges_rx, self.charges_input.text())
        if not charges:
            self.warn('Enter at least one charge value.')
            return False
        icharges = []
        for c in charges:
            try:
                c = int(c)
            except ValueError:
                self.warn('{} is not a valid charge.'.format(c))
                return False
            icharges.append(c)
        self.charges = icharges
        return True

    def check_mz_input(self):
        """ Validate input for mz_input.
            Returns True on correct input, False on error.
        """
        if self.mz_input.text() == '':
            self.mz = None
        else:
            try:
                self.mz = float(self.mz_input.text())
            except ValueError:
                try:
                    m = Molecule(self.mz_input.text())
                except:
                    self.warn('Enter target as a number or as a molecular formula.')
                    return False
                else:
                    self.mz = self.mz_input.text()
        return True

    def keyPressEvent(self, event):
        """ Link enter/return to calculate button,
            cmd/ctrl-enter/return to standard ratio,
            cmd/ctrl-c to table data copy,
            cmd/ctrl-a to select all in table,
            cmd/ctrl-h to open help,
            cmd/ctrl-d to display spectrum.
        """
        key = event.key()
        mod = event.modifiers()
        if (key == QtCore.Qt.Key_Enter or key == QtCore.Qt.Key_Return):
            if mod == QtCore.Qt.ControlModifier:
                self.show_standard_ratio()
            else:
                self.calculate_interference()
        elif (key == QtCore.Qt.Key_C and mod == QtCore.Qt.ControlModifier):
            self.table_output.copy()
        elif (key == QtCore.Qt.Key_A and mod == QtCore.Qt.ControlModifier):
            self.table_output.selectAll()
        elif (key == QtCore.Qt.Key_H and mod == QtCore.Qt.ControlModifier):
            self.show_help()
        elif (key == QtCore.Qt.Key_D and mod == QtCore.Qt.ControlModifier):
            self.toggle_spectrum()
        else:
            super().keyPressEvent(event)

    @QtCore.pyqtSlot()
    def calculate_interference(self):
        """ Take input, calculate mass spectrum, display in table. """
        if not (self.check_atoms_input() and
                self.check_charges_input() and
                self.check_mz_input()):
            return

        if not self.mz:
            qmsg = widgets.QMessageBox(self)
            qmsg.setText('Long calculation warning')
            qmsg.setInformativeText(mz_warning)
            qmsg.setIcon(widgets.QMessageBox.Warning)
            qmsg.setStandardButtons(widgets.QMessageBox.Ok|widgets.QMessageBox.Cancel)
            if qmsg.exec_() == widgets.QMessageBox.Cancel:
                return

        self.maxsize = self.maxsize_input.value()
        self.chargesign = self.chargesign_input.currentText()
        self.mzrange = self.mzrange_input.value()

        data = interference(self.atoms, self.mz, targetrange=self.mzrange,
            maxsize=self.maxsize, charge=self.charges, chargesign=self.chargesign)
        data.pop('charge')
        data.columns = ['molecule', 'mass/charge', 'Δmass/charge', 'mz/Δmz (MRP)', 'probability', 'target']
        data.index = range(1, data.shape[0] + 1)

        model = TableModel(data, table='interference')
        self.table_output.setModel(model)
        self.table_output.setColumnHidden(5, True)
        try:
            # PyQt5
            self.table_output.horizontalHeader().setSectionResizeMode(widgets.QHeaderView.Stretch)
        except AttributeError:
            # PyQt4
            self.table_output.horizontalHeader().setResizeMode(widgets.QHeaderView.Stretch)

        self.spectrum_window.plot_spectrum(data)

    @QtCore.pyqtSlot()
    def show_standard_ratio(self):
        """ Show the standard ratios. """
        if not self.check_atoms_input():
            return

        data = standard_ratio(self.atoms)
        data['target'] = False
        if self.check_mz_input() and isinstance(self.mz, str):
            m = Molecule(self.mz)
            target_data = standard_ratio(m.elements)
            target_data['target'] = True
            data = data.append(target_data)
        data.index = range(1, data.shape[0] + 1)

        model = TableModel(data, table='std_ratios')
        self.table_output.setModel(model)
        self.table_output.setColumnHidden(5, False)
        self.table_output.setColumnHidden(6, True)
        try:
            # PyQt5
            self.table_output.horizontalHeader().setSectionResizeMode(widgets.QHeaderView.Stretch)
        except AttributeError:
            # PyQt4
            self.table_output.horizontalHeader().setResizeMode(widgets.QHeaderView.Stretch)

    @QtCore.pyqtSlot()
    def show_help(self):
        """ Display help window. """
        dialog = widgets.QDialog(parent=self)
        dialog.resize(600,620)
        text = widgets.QLabel(parent=dialog)
        text.setText(help_text)
        text.setWordWrap(True)
        text.setMargin(20)
        dialog.exec_()

    @QtCore.pyqtSlot()
    def toggle_spectrum(self):
        """ Show interference data in a spectrum. """
        if self.spectrum_window.isHidden():
            self.spectrum_window.show()
        else:
            self.spectrum_window.hide()


def run():
    """ Run the gui. """
    app = widgets.QApplication([])
    mainwindow = MainWindow()
    mainwindow.move(200,100)
    mainwindow.show()
    pos = mainwindow.pos()
    pos.setX(pos.x() + mainwindow.width())
    mainwindow.centralWidget().spectrum_window.move(pos)
    sys.exit(app.exec_())

if __name__ == '__main__':
    run()
