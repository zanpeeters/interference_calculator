#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" GUI for mass spectrum program. """

try:
    from PyQt5 import QtCore, QtGui
    from PyQt5 import QtWidgets as widgets
except ImportError:
    try:
        from PyQt4 import QtCore, QtGui
        from PyQt4 import QtGui as widgets
    except ImportError:
        raise ImportError('You need to have either PyQt4 or PyQt5 installed.')

from pandas import DataFrame
import sys, re, platform
import masstable

_float_qrx = QtCore.QRegExp(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$')
_float_validator = QtGui.QRegExpValidator(_float_qrx)
_isotope_rx = re.compile(r'(\d*[A-Z][a-z]{0,2})')
_charges_rx = re.compile(r'(\d+)')

class TableModel(QtCore.QAbstractTableModel):
    """ Take a pandas DataFrame and set data in a QTableModel (read-only). """

    def __init__(self, data, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                if index.column() in (1, 2):
                    return '{:.5f}'.format(self._data.iloc[index.row(), index.column()])
                elif index.column() == 3:
                    return '{:.1f}'.format(self._data.iloc[index.row(), index.column()])
                elif index.column() == 4:
                    return '{:.5g}'.format(self._data.iloc[index.row(), index.column()])
                else:
                    return '{}'.format(self._data.iloc[index.row(), index.column()])
            elif role == QtCore.Qt.TextAlignmentRole:
                if index.column() == 0:
                    return QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
                else:
                    return QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter

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

        if options.widget:
            style = options.widget.style()
        else:
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


class LogSpinBox(widgets.QDoubleSpinBox):
    """ Like QDoubleSpinBox but steps by powers of 10 

        Default minimum and maximum are 10^-99 and 10^+99,
        although these limits are arbitrary. The real limits
        are sys.float_info.min and sys.float_info.max.

        The number of decimals is set to 99, so that the
        entered value is stored without rounding, while
        displaying at a convenient length.
    """
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setMaximum(1.0e99)
        self.setMinimum(1.0e-99)
        self.setDecimals(99)
        self.setValue(1.0)

    def stepBy(self, steps):
        """ Step by powers of 10. """
        self.setValue(self.value() * 10 ** steps)

    def validate(self, inputs, pos):
        """ Validate input based on float regex. """
        return _float_validator.validate(inputs, pos)

    def textFromValue(self, value):
        """ Override QDoubleSpinBox.textFromValue to display floats. """
        return '{:n}'.format(value)


class MainWindow(widgets.QMainWindow):
    """ Main window for masstable ui. """
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Mass interference calculator')
        self.resize(600, 720)
        self.setMinimumSize(600, 720)
        self.setCentralWidget(MainWidget(parent=self))


class MainWidget(widgets.QWidget):
    """ Central widget class for masstable ui. """
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.atoms = []
        self.charges = [1]
        self.mz = 0.0
        self.mzrange = 0.3
        self.maxsize = 5
        self.minabundance = 1.0e-8
        self.has_minabundance = False

        # Inputs
        self.atoms_input = widgets.QLineEdit(parent=self)
        self.atoms_input.setPlaceholderText('list of atoms in sample: C H Si Ca')

        self.maxsize_label = widgets.QLabel('max size', parent=self)
        self.maxsize_input = widgets.QSpinBox(parent=self)
        self.maxsize_input.setValue(self.maxsize)

        self.charges_label = widgets.QLabel('target charge(s)', parent=self)
        self.charges_input = widgets.QLineEdit(' '.join([str(c) for c in self.charges]), parent=self)

        self.chargesign_input = widgets.QComboBox(parent=self)
        self.chargesign_input.addItem('+', '+')
        self.chargesign_input.addItem('-', '-')

        self.mz_label = widgets.QLabel('target mass/charge', parent=self)
        self.mz_input = widgets.QLineEdit(parent=self)
        self.mz_input.setPlaceholderText('number or formula: 26.0 or 12C 14N')

        self.mzrange_label = widgets.QLabel('±', parent=self)
        self.mzrange_input = widgets.QDoubleSpinBox(parent=self)
        self.mzrange_input.setValue(self.mzrange)
        self.mzrange_input.setSingleStep(0.1)
        self.mzrange_input.setDecimals(3)

        self.minabundance_label = widgets.QLabel('missing abundance value', parent=self)
        self.minabundance_input = LogSpinBox(parent=self)
        self.minabundance_input.setValue(self.minabundance)
        self.minabundance_input.setEnabled(self.has_minabundance)
        self.minabundance_select = widgets.QCheckBox()
        self.minabundance_select.setChecked(self.has_minabundance)

        # Action button
        self.interference_button = widgets.QPushButton('calculate interference', parent=self)
        self.isotopes_button = widgets.QPushButton('natural abundance', parent=self)

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
        self.param_bottom_layout.addWidget(self.minabundance_label, 0, 3)
        self.param_bottom_layout.addWidget(self.minabundance_input, 1, 3)
        self.param_bottom_layout.addWidget(self.minabundance_select, 1, 4)

        self.parameters_group_layout.addLayout(self.param_top_layout)
        self.parameters_group_layout.addLayout(self.param_bottom_layout)

        self.button_layout = widgets.QHBoxLayout()
        self.button_layout.addWidget(self.interference_button)
        self.button_layout.addWidget(self.isotopes_button)

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
        self.interference_button.clicked.connect(self.calculate)
        self.isotopes_button.clicked.connect(self.show_isotopes)
        self.minabundance_select.stateChanged.connect(self.toggle_minabundance)
        self.atoms_input.editingFinished.connect(self.check_atoms_input)
        self.charges_input.editingFinished.connect(self.check_charges_input)
        self.mz_input.editingFinished.connect(self.check_mz_input)

        # Set jump order for tab
        self.setTabOrder(self.atoms_input, self.mz_input)
        self.setTabOrder(self.mz_input, self.mzrange_input)
        self.setTabOrder(self.mzrange_input, self.charges_input)
        self.setTabOrder(self.charges_input, self.chargesign_input)
        self.setTabOrder(self.chargesign_input, self.maxsize_input)
        self.setTabOrder(self.maxsize_input, self.minabundance_input)
        self.setTabOrder(self.minabundance_input, self.minabundance_select)
        self.setTabOrder(self.minabundance_select, self.interference_button)
        self.setTabOrder(self.interference_button, self.isotopes_button)
        self.setTabOrder(self.isotopes_button, self.table_output)

        # Set tooltip help
        self.atoms_input.setToolTip('''<html><head/><body>
            <p><b>List of atoms to include.</b></p>
            <p>Give the composition of the sample. List all atoms separated by space.
            All isotopes for each atom are included automatically.</p>
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
            <p>This sign is only used for display in the molecular formula, the value
            mass-to-charge ratio does not change with this selection.</p>
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
        self.minabundance_input.setToolTip('''<html><head/><body>
            <p><b>Missing abundance replacement value.</b></p>
            <p>The mass table downloaded from the NIST website has missing abundance
            values for some isotopes, possibly because the abundance is too low to measure
            accurately.</p>
            <p>Use this input field to give a (low) replacement value for the abundance.
            With a replacement value, a probability can be calculated even though the actual
            value may be unrealistic. This way, two molecules with unknown abundances can
            still be compared, relative to eachother.</p>
            <p>Use the tick box to disable this input field. When this input field is
            disabled, NaN (Not a Number) will be used and the probability of the molecule
            will not be calculated (NaN is respected).</p>
            </body></html>''')
        self.minabundance_select.setToolTip('''<html><head/><body>
            <p><b>Enable/disable missing abundance replacement.</b></p>
            <p>With this option selected, the given missing abundance replacement value
            will be used. Unselected, the missing abundance input field will be disabled
            and NaN (Not a Number) will be displayed. All molecules that include one or
            more isotopes with an unknown abundance will have NaN probability.</p>
            </body></html>''')
        self.interference_button.setToolTip('''<html><head/><body>
            <p>Calculate mass interference.</p>
            </body></html>''')
        self.isotopes_button.setToolTip('''<html><head/><body>
            <p>Show natural abundances for isotopes.</p>
            </body></html>''')

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
            if not (a in masstable.elements or a in masstable.isotopes):
                msg = '{} is not an element or isotope.'.format(a)
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
                    m = masstable.Molecule(self.mz_input.text())
                    ch = m.charge
                    if ch == 0:
                        ch = 1
                    mz = m.mass/ch
                except:
                    msg = 'Enter mz as a number or as a molecular formula.'
                    self.statusbar.showMessage(msg, msecs=5000)
                    return False
        self.mz = mz
        return True

    @QtCore.pyqtSlot()
    def toggle_minabundance(self):
        """ Dis/enable minabundance on toggle checkbox. """
        self.has_minabundance = bool(self.minabundance_select.checkState())
        self.minabundance_input.setEnabled(self.has_minabundance)

    def keyPressEvent(self, event):
        """ Link enter/return to calculate button. """
        if (event.key() == QtCore.Qt.Key_Enter 
            or event.key() == QtCore.Qt.Key_Return):
                self.calculate()
        else:
            super().keyPressEvent(event)

    @QtCore.pyqtSlot()
    def calculate(self):
        """ Take input, calculate mass spectrum, display in table. """
        if not (self.check_atoms_input() and
                self.check_charges_input() and
                self.check_mz_input()):
            return

        if self.has_minabundance:
            self.minabundance = self.minabundance_input.value()
        else:
            self.minabundance = None
        self.maxsize = self.maxsize_input.value()
        self.chargesign = self.chargesign_input.currentText()
        self.mzrange = self.mzrange_input.value()

        ms = masstable.mass_interference(self.atoms, self.mz, mzrange=self.mzrange, 
            maxsize=self.maxsize, charge=self.charges, chargesign=self.chargesign,
            minabundance=self.minabundance)
        ms.pop('charge')
        ms.columns = ['molecule', 'mass/charge', 'Δmass/charge', 'mz/Δmz (MRP)', 'probability']
        ms.index = range(1, ms.shape[0] + 1)

        model = TableModel(ms)
        self.table_output.setModel(model)
        self.table_output.horizontalHeader().setSectionResizeMode(widgets.QHeaderView.Stretch)

    @QtCore.pyqtSlot()
    def show_isotopes(self):
        """ Show the standard ratios. """
        if not self.check_atoms_input():
            return

        d = masstable.standard_ratios(self.atoms)
        # recreate DataFrame to get index (isotopes) as a column.
        ratios = DataFrame([d.index, d['isotope_mass'], d['abundance'], d['inverse']]).T
        ratios.columns = ['isotope', 'mass', 'abundance', 'inverse']
        ratios.index = range(1, ratios.shape[0] + 1)

        model = TableModel(ratios)
        self.table_output.setModel(model)
        self.table_output.horizontalHeader().setSectionResizeMode(widgets.QHeaderView.Stretch)


def run():
    """ Run the gui. """
    app = widgets.QApplication([])
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    run()
