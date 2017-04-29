#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" GUI for mass spectrum program. """

from PyQt4 import QtCore, QtGui, QtWebKit
from pandas import set_option
import sys, re
import masstable
import sorttable

# Debug: trace warnings
# import warnings
# warnings.filterwarnings('error')

_float_qrx = QtCore.QRegExp(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$')
_float_validator = QtGui.QRegExpValidator(_float_qrx)


class LogSpinBox(QtGui.QDoubleSpinBox):
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


class MainWindow(QtGui.QMainWindow):
    """ Main window for masstable ui. """
    def __init__(self):
        super().__init__()

        # Setup
        self.setWindowTitle('Mass interference calculator')
        self.resize(600, 720)
        self.setMinimumSize(600, 720)
        self.setCentralWidget(MainWidget(parent=self))

        # Menubar
        # For Mac OSX:
        # - preferences and quit will automatically be moved to Application
        #   menu where they belong.
        # - File menu is then empty and will not be rendered.
        # - On pyqt, name of the application is Python and can only be changed when
        #   installed as an .app with name in Info.plist.
        # - Quit menu always exists in Application menu.
        # - Ctrl is mapped to Cmd on Mac.
        # self.preferences_menuitem = QtGui.QAction('Preferences', self)
        # self.preferences_menuitem.setShortcut('Ctrl+.')
        # self.preferences_menuitem.triggered.connect(self.open_preferences)
        # 
        # self.quit_menuitem = QtGui.QAction('&Quit', self)
        # self.quit_menuitem.setShortcut('Ctrl+Q')
        # self.quit_menuitem.triggered.connect(QtGui.qApp.quit)
        # 
        # self.menubar = self.menuBar()
        # self.file_menu = self.menubar.addMenu('&Program')
        # self.file_menu.addAction(self.preferences_menuitem)
        # self.file_menu.addAction(self.quit_menuitem)


# class PreferencesWindow(QtGui.QWidget):


class MainWidget(QtGui.QWidget):
    """ Central widget class for masstable ui. """
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        # Default values
        self.atoms = []
        self.charges = [1]
        self.mz = 0.0
        self.mzrange = 0.3
        self.maxsize = 5
        self.minabundance = 1.0e-8
        self.has_minabundance = False

        # Inputs
        # self.atoms_label = QtGui.QLabel('sample composition', parent=self)
        self.atoms_input = QtGui.QLineEdit(parent=self)
        self.atoms_input.setPlaceholderText('list of atoms in sample: C H Si Ca')

        self.maxsize_label = QtGui.QLabel('max size', parent=self)
        self.maxsize_input = QtGui.QSpinBox(parent=self)
        self.maxsize_input.setValue(self.maxsize)

        self.charges_label = QtGui.QLabel('target charge(s)', parent=self)
        self.charges_input = QtGui.QLineEdit(' '.join([str(c) for c in self.charges]), parent=self)

        self.chargesign_input = QtGui.QComboBox(parent=self)
        self.chargesign_input.addItem('+', '+')
        self.chargesign_input.addItem('-', '-')

        self.mz_label = QtGui.QLabel('target mass/charge', parent=self)
        self.mz_input = QtGui.QLineEdit(parent=self)
        self.mz_input.setPlaceholderText('number or formula: 26.0 or 12C 14N')

        self.mzrange_label = QtGui.QLabel('±', parent=self)
        self.mzrange_input = QtGui.QDoubleSpinBox(parent=self)
        self.mzrange_input.setValue(self.mzrange)
        self.mzrange_input.setSingleStep(0.1)
        self.mzrange_input.setDecimals(3)

        self.minabundance_label = QtGui.QLabel('missing abundance value', parent=self)
        self.minabundance_input = LogSpinBox(parent=self)
        self.minabundance_input.setValue(self.minabundance)
        self.minabundance_input.setEnabled(self.has_minabundance)
        self.minabundance_select = QtGui.QCheckBox()
        self.minabundance_select.setChecked(self.has_minabundance)

        # Action button
        self.interference_button = QtGui.QPushButton('calculate interference', parent=self)
        self.isotopes_button = QtGui.QPushButton('natural abundance', parent=self)

        # Table output
        self.table_output = QtWebKit.QWebView()

        # Qt4 on Mac has a problem with the focus rectangle around input
        # boxes. Not fixed until Qt5.
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
        pol = QtGui.QSizePolicy()
        pol.setHorizontalPolicy(QtGui.QSizePolicy.Ignored)
        pol.setVerticalPolicy(QtGui.QSizePolicy.Maximum)

        self.atoms_group = QtGui.QGroupBox('sample composition')
        self.atoms_group.setSizePolicy(pol)
        self.atoms_group_layout = QtGui.QHBoxLayout()
        self.atoms_group_layout.addWidget(self.atoms_input)
        self.atoms_group.setLayout(self.atoms_group_layout)

        self.parameters_group = QtGui.QGroupBox('interference parameters')
        self.parameters_group.setSizePolicy(pol)
        self.parameters_group_layout = QtGui.QVBoxLayout()
        self.parameters_group.setLayout(self.parameters_group_layout)

        self.param_top_layout = QtGui.QGridLayout()
        self.param_top_layout.addWidget(self.mz_label, 0, 0)
        self.param_top_layout.addWidget(self.mz_input, 1, 0)
        self.param_top_layout.addWidget(self.mzrange_label, 1, 1)
        self.param_top_layout.addWidget(self.mzrange_input, 1, 2)

        self.param_bottom_layout = QtGui.QGridLayout()
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

        self.button_layout = QtGui.QHBoxLayout()
        self.button_layout.addWidget(self.interference_button)
        self.button_layout.addWidget(self.isotopes_button)

        self.output_layout = QtGui.QHBoxLayout()
        self.output_layout.addWidget(self.table_output)

        self.layout = QtGui.QVBoxLayout()
        self.layout.addWidget(self.atoms_group)
        self.layout.addItem(QtGui.QSpacerItem(10,10))
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
        _isotope_rx = r'(\d*[A-Z][a-z]{0,2})'
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
        _charges_rx = r'(\d+)'
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

    def format_float1f(self, num):
        return '{:.1f}'.format(num)

    def format_float5f(self, num):
        return '{:.5f}'.format(num)

    def format_float4g(self, num):
        return '{:.4g}'.format(num)

    def format_float6g(self, num):
        return '{:.6g}'.format(num)

    @QtCore.pyqtSlot()
    def calculate(self):
        """ Take input, calculate mass spectrum, display in table. """
        # Clear table
        self.table_output.setHtml('')

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

        # Make sure molecules column in DataFrame is not elided.
        # DataFrame.to_html() apparently adheres to max_colwidth.
        set_option('max_colwidth', 1024)

        ms = masstable.mass_interference(self.atoms, self.mz, mzrange=self.mzrange, 
                maxsize=self.maxsize, charge=self.charges, chargesign=self.chargesign,
                minabundance=self.minabundance)

        # Don't need charge column
        ms.pop('charge')

        # Nicer labels
        ms.columns = ['molecule', 'mass/charge', 'Δmass/charge', 'mz/Δmz (MRP)', 'probability']

        table = ms.to_html(escape=False, index=False, classes='sortable',
                           justify='left',
                           formatters={'mass/charge': self.format_float5f,
                                       'Δmass/charge': self.format_float5f,
                                       'mz/Δmz (MRP)': self.format_float1f,
                                       'probability': self.format_float6g})

        # to_html sets border=1 with no option to turn it off
        table = re.sub('border=\"1\" ', '', table)

        self.table_output.setHtml(sorttable.html_head + table + sorttable.html_foot)

    @QtCore.pyqtSlot()
    def show_isotopes(self):
        """ Show the standard ratios. """
        # Clear table
        self.table_output.setHtml('')

        if not self.check_atoms_input():
            return

        # Make sure molecules column in DataFrame is not elided.
        # DataFrame.to_html() apparently adheres to max_colwidth.
        set_option('max_colwidth', 1024)

        ratios = masstable.standard_ratios(self.atoms)

        # Nicer labels
        ratios.columns = ['mass', 'abundance', 'inverse']

        table = ratios.to_html(escape=False, index=True, classes='sortable',
                           justify='left', index_names=True,
                           formatters={'mass': self.format_float5f,
                                       'abundance': self.format_float6g,
                                       'inverse': self.format_float1f})

        # to_html sets border=1 with no option to turn it off
        table = re.sub('border=\"1\" ', '', table)

        self.table_output.setHtml(sorttable.html_head + table + sorttable.html_foot)


def run():
    """ Run the gui. """
    app = QtGui.QApplication([])
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    run()
