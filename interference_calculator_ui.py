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

import sys, re, platform
from interference_calculator import interference, standard_ratio
from molecule import Molecule, periodic_table

_float_qrx = QtCore.QRegExp(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$')
_float_validator = QtGui.QRegExpValidator(_float_qrx)
_isotope_rx = re.compile(r'(\d*[A-Z][a-z]{0,2})')
_charges_rx = re.compile(r'(\d+)')

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
                        m = Molecule(self._data.iloc[index.row(), index.column()])
                        return m.formula(style='html')
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
                        return m.formula(style='html')
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


class TableView(widgets.QTableView):
    """ Implement a QTableView which can display HTML in arbitrary columns """
    def __init__(self, html_cols=None, parent=None):
        widgets.QTableView.__init__(self, parent=parent)
        self.setSortingEnabled(True)
        if html_cols is not None:
            if isinstance(html_cols, int):
                html_cols = [html_cols]
            [self.setItemDelegateForColumn(c, HTMLDelegate(parent=parent)) for c in html_cols]

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

    def keyPressEvent(self, event):
        key = event.key()
        mod = event.modifiers()
        if (mod == QtCore.Qt.ControlModifier and key == QtCore.Qt.Key_C):
            self.copy()
        elif (mod == QtCore.Qt.ControlModifier and key == QtCore.Qt.Key_A):
            self.selectAll()

    def copy(self):
        pasteboard = widgets.QApplication.clipboard()
        selected = self.selectedIndexes()
        output = ''
        prev_col = 0
        for s in selected:
            if s.column() < prev_col:
                output += '\n'
            elif output:
                output += '    '
            output += str(s.data(role=QtCore.Qt.EditRole))
            prev_col = s.column()
        pasteboard.setText(output)


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


class MainWindow(widgets.QMainWindow):
    """ Main window for interference calculator ui. """
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Mass interference calculator')
        self.resize(600, 720)
        self.setMinimumSize(600, 720)
        self.setCentralWidget(MainWidget(parent=self))


class MainWidget(widgets.QWidget):
    """ Central widget class for interference calculator ui. """
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.atoms = []
        self.charges = [1]
        self.mz = 0.0
        self.mzrange = 0.3
        self.maxsize = 5

        # Inputs
        self.atoms_input = widgets.QLineEdit(parent=self)
        self.atoms_input.setPlaceholderText('list of atoms in sample: C H Si Ca')

        self.maxsize_label = widgets.QLabel('max size', parent=self)
        self.maxsize_input = widgets.QSpinBox(parent=self)
        self.maxsize_input.setValue(self.maxsize)

        self.charges_label = widgets.QLabel('target charge(s)', parent=self)
        self.charges_input = widgets.QLineEdit(' '.join([str(c) for c in self.charges]), parent=self)

        self.chargesign_input = widgets.QComboBox(parent=self)
        self.chargesign_input.addItem('-', '-')
        self.chargesign_input.addItem('o', 'o')
        self.chargesign_input.addItem('+', '+')

        self.mz_label = widgets.QLabel('target mass/charge or molecule', parent=self)
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

        # Table output
        self.table_output = TableView(html_cols=0)

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
        self.statusbar.setStyleSheet('color: #FF0000;')

        # Layout objects
        pol = widgets.QSizePolicy()
        pol.setHorizontalPolicy(widgets.QSizePolicy.Ignored)
        pol.setVerticalPolicy(widgets.QSizePolicy.Maximum)

        self.atoms_group = widgets.QGroupBox('sample composition')
        self.atoms_group.setSizePolicy(pol)
        self.atoms_group_layout = widgets.QHBoxLayout()
        self.atoms_group_layout.addWidget(self.atoms_input)
        self.atoms_group.setLayout(self.atoms_group_layout)

        self.parameters_group = widgets.QGroupBox('interference parameters')
        self.parameters_group.setSizePolicy(pol)
        self.parameters_group_layout = widgets.QVBoxLayout()
        self.parameters_group.setLayout(self.parameters_group_layout)

        self.param_top_layout = widgets.QGridLayout()
        self.param_top_layout.addWidget(self.mz_label, 0, 0)
        self.param_top_layout.addWidget(self.mz_input, 1, 0)
        self.param_top_layout.addWidget(self.mzrange_label, 1, 1)
        self.param_top_layout.addWidget(self.mzrange_input, 1, 2)

        self.param_bottom_layout = widgets.QGridLayout()
        self.param_bottom_layout.addWidget(self.charges_label, 0, 0)
        self.param_bottom_layout.addWidget(self.charges_input, 1, 0)
        self.param_bottom_layout.addWidget(self.chargesign_input, 1, 1)
        self.param_bottom_layout.addWidget(self.maxsize_label, 0, 2)
        self.param_bottom_layout.addWidget(self.maxsize_input, 1, 2)

        self.parameters_group_layout.addLayout(self.param_top_layout)
        self.parameters_group_layout.addLayout(self.param_bottom_layout)

        self.button_layout = widgets.QHBoxLayout()
        self.button_layout.addWidget(self.interference_button)
        self.button_layout.addWidget(self.standard_ratio_button)

        self.output_layout = widgets.QHBoxLayout()
        self.output_layout.addWidget(self.table_output)

        self.layout = widgets.QVBoxLayout()
        self.layout.addWidget(self.atoms_group)
        self.layout.addItem(widgets.QSpacerItem(10,10))
        self.layout.addWidget(self.parameters_group)
        self.layout.addLayout(self.button_layout)
        self.layout.addLayout(self.output_layout)
        self.setLayout(self.layout)

        # Connect
        self.interference_button.clicked.connect(self.calculate_interference)
        self.standard_ratio_button.clicked.connect(self.show_standard_ratio)
        self.atoms_input.editingFinished.connect(self.check_atoms_input)
        self.charges_input.editingFinished.connect(self.check_charges_input)
        self.mz_input.editingFinished.connect(self.check_mz_input)

        # Set jump order for tab
        self.setTabOrder(self.atoms_input, self.mz_input)
        self.setTabOrder(self.mz_input, self.mzrange_input)
        self.setTabOrder(self.mzrange_input, self.charges_input)
        self.setTabOrder(self.charges_input, self.chargesign_input)
        self.setTabOrder(self.chargesign_input, self.maxsize_input)
        self.setTabOrder(self.maxsize_input, self.interference_button)
        self.setTabOrder(self.interference_button, self.standard_ratio_button)
        self.setTabOrder(self.standard_ratio_button, self.table_output)

        # Set tooltip help
        self.atoms_input.setToolTip('''<html><head/><body>
            <p><b>List of atoms to include.</b></p>
            <p>Give the composition of the sample. List all atoms separated by space.
            All isotopes for each element are included automatically.</p>
            </body></html>''')
        self.charges_input.setToolTip('''<html><head/><body>
            <p><b>List of charges to apply.</b></p>
            <p>Give a list of ion charges, separated by space, to be considered for the
            target ion. Only use numbers, do not include the charge sign. For example,
            \'1 2 3\' will calculate mass-to-charge ratios for +, 2+, and 3+ (or -)
            charged ions. Default value is 1.</p>
            </body></html>''')
        self.chargesign_input.setToolTip('''<html><head/><body>
            <p><b>Select sign of ionic charge.</b></p>
            <p>The mass of the resulting ion depends on the charge sign and is corrected
            for the mass of the extra electron (&ndash;), missing electron (+), or not
            corrected (o).</p>
            </body></html>''')
        self.mz_input.setToolTip('''<html><head/><body>
            <p><b>Target mass-to-charge ratio.</b></p>
            <p>Give a mass-to-charge (m/z) value to filter the results. Only molecules
            with m/z ± range and up to max size atoms will be shown.</p>
            <p>If empty, all possible combinations of all isotopes of the selected
            atoms will be displayed. <i>This may be a very long list!</i></p>
            </body></html>''')
        self.mzrange_input.setToolTip('''<html><head/><body>
            <p><b>Target mass-to-charge range.</b></p>
            <p>Give a range of mass-to-charge ratios to filter the results.
            Default value is 0.3</p>
            </body></html>''')
        self.maxsize_input.setToolTip('''<html><head/><body>
            <p><b>Maximum molecule size.</b></p>
            <p>Give the maximum number of atoms in a molecule. The number of
            possible combinations of <i>n</i> atoms in a molecule grows exponentially
            with <i>n</i>, so keep this number low to avoid obscenely large lists.</p>
            </body></html>''')
        self.interference_button.setToolTip('''<html><head/><body>
            <p>(enter)</p>
            </body></html>''')
        if sys.platform == 'darwin':
            modifier = '&#8984;'
        else:
            modifier = 'ctrl'
        self.standard_ratio_button.setToolTip('''<html><head/><body>
            <p>({}-enter)</p>
            </body></html>'''.format(modifier))

    @QtCore.pyqtSlot()
    def check_atoms_input(self):
        """ Validate input for atoms_input.
            Returns True on proper validation, False on error.
        """
        atoms = re.findall(_isotope_rx, self.atoms_input.text())
        if not atoms:
            msg = 'Enter at least one element or isotope.'
            self.statusbar.showMessage(msg, msecs=5000)
            return False
        for a in atoms:
            if not (periodic_table['element'] == a).any(): 
                msg = '{} is not an element or missing from the internal periodic table.'.format(a)
                self.statusbar.showMessage(msg, msecs=5000)
                return False
        self.atoms = atoms
        return True

    @QtCore.pyqtSlot()
    def check_charges_input(self):
        """ Validate input for charges_input.
            Returns True on correct input, False on error.
        """
        charges = re.findall(_charges_rx, self.charges_input.text())
        if not charges:
            msg = 'Enter at least one charge value.'
            self.statusbar.showMessage(msg, msecs=5000)
            return False
        icharges = []
        for c in charges:
            try:
                c = int(c)
            except ValueError:
                msg = '{} is not a valid charge.'.format(c)
                self.statusbar.showMessage(msg, msecs=5000)
                return False
            icharges.append(c)
        self.charges = icharges
        return True

    @QtCore.pyqtSlot()
    def check_mz_input(self):
        """ Validate input for mz_input.
            Returns True on correct input, False on error.
        """
        if self.mz_input.text() == '':
            mz = None
        else:
            try:
                mz = float(self.mz_input.text())
            except ValueError:
                try:
                    m = Molecule(self.mz_input.text())
                    if m.charge == 0:
                        mz = m.mass
                    else:
                        mz = m.mass/m.charge
                except:
                    msg = 'Enter mz as a number or as a molecular formula.'
                    self.statusbar.showMessage(msg, msecs=5000)
                    return False
        self.mz = mz
        return True

    def keyPressEvent(self, event):
        """ Link enter/return to calculate button,
            cmd/ctrl-enter/return to standard ratio.
        """
        key = event.key()
        mod = event.modifiers()
        if (key == QtCore.Qt.Key_Enter or key == QtCore.Qt.Key_Return):
            if mod == QtCore.Qt.ControlModifier:
                self.show_standard_ratio()
            else:
                self.calculate_interference()
        else:
            super().keyPressEvent(event)

    @QtCore.pyqtSlot()
    def calculate_interference(self):
        """ Take input, calculate mass spectrum, display in table. """
        if not (self.check_atoms_input() and
                self.check_charges_input() and
                self.check_mz_input()):
            return

        self.maxsize = self.maxsize_input.value()
        self.chargesign = self.chargesign_input.currentText()
        self.mzrange = self.mzrange_input.value()

        data = interference(self.atoms, self.mz, mzrange=self.mzrange, 
            maxsize=self.maxsize, charge=self.charges, chargesign=self.chargesign, style='plain')
        data.pop('charge')
        data.columns = ['molecule', 'mass/charge', 'Δmass/charge', 'mz/Δmz (MRP)', 'probability']
        data.index = range(1, data.shape[0] + 1)

        model = TableModel(data, table='interference')
        self.table_output.setModel(model)
        try:
            # PyQt5
            self.table_output.horizontalHeader().setSectionResizeMode(widgets.QHeaderView.Stretch)
        except AttributeError:
            # PyQt4
            self.table_output.horizontalHeader().setResizeMode(widgets.QHeaderView.Stretch)

    @QtCore.pyqtSlot()
    def show_standard_ratio(self):
        """ Show the standard ratios. """
        if not self.check_atoms_input():
            return

        data = standard_ratio(self.atoms)
        data.index = range(1, data.shape[0] + 1)

        model = TableModel(data, table='std_ratios')
        self.table_output.setModel(model)
        try:
            # PyQt5
            self.table_output.horizontalHeader().setSectionResizeMode(widgets.QHeaderView.Stretch)
        except AttributeError:
            # PyQt4
            self.table_output.horizontalHeader().setResizeMode(widgets.QHeaderView.Stretch)


def run():
    """ Run the gui. """
    app = widgets.QApplication([])
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    run()
