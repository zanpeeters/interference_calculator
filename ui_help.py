# -*- coding: utf-8 -*-
""" Help and tooltip text for ui.py """
import sys

atoms_input_tooltip = '''
<html><head/><body>
<p><b>List of atoms to include.</b></p>
<p>Give the composition of the sample. List all atoms separated by space.
All isotopes for each element are included automatically.</p>
</body></html>
'''

charges_input_tooltip = '''
<html><head/><body>
<p><b>List of charges to apply.</b></p>
<p>Give a list of ion charges, separated by space, to be considered for the
target ion. Only use numbers, do not include the charge sign. For example,
\'1 2 3\' will calculate mass-to-charge ratios for +, 2+, and 3+ (or -)
charged ions. Default value is 1.</p>
</body></html>
'''

chargesign_input_tooltip = '''
<html><head/><body>
<p><b>Select sign of ionic charge.</b></p>
<p>The mass of the resulting ion depends on the charge sign and is corrected
for the mass of the extra electron (&ndash;), missing electron (+), or not
corrected (o).</p>
</body></html>
'''

mz_input_tooltip = '''
<html><head/><body>
<p><b>Target mass-to-charge ratio.</b></p>
<p>Give a mass-to-charge (m/z) value to filter the results. Only molecules
with m/z ± range and up to max size atoms will be shown.</p>
<p>If empty, all possible combinations of all isotopes of the selected
atoms will be displayed. <i>This may be a very long list!</i></p>
</body></html>
'''

mzrange_input_tooltip = '''
<html><head/><body>
<p><b>Target mass-to-charge range.</b></p>
<p>Give a range of mass-to-charge ratios to filter the results.
Default value is 0.3</p>
</body></html>
'''

maxsize_input_tooltip = '''
<html><head/><body>
<p><b>Maximum molecule size.</b></p>
<p>Give the maximum number of atoms in a molecule. The number of
possible combinations of <i>n</i> atoms in a molecule grows exponentially
with <i>n</i>, so keep this number low to avoid obscenely large lists.</p>
</body></html>
'''

interference_button_tooltip = '''
<html><head/><body>
<p>(enter)</p>
</body></html>
'''

if sys.platform == 'darwin':
    _modifier = '&#8984;'
else:
    _modifier = 'ctrl'

standard_ratio_button_tooltip = '''
<html><head/><body>
<p>({}-enter)</p>
</body></html>
'''.format(_modifier)

help_button_tooltip = '''
<html><head/><body>
<p>Help ({}-H)</p>
</body></html>
'''.format(_modifier)

spectrum_button_tooltip = '''
<html><head/><body>
<p>Show the mass spectrum ({}-D)</p>
</body></html>
'''.format(_modifier)

mz_warning = '''
<html><head/><body>
<p>If you do not specify a target, <b>ALL</b> combinations up to
<i>max size</i> will be calculated. This can take a long time.</p>
<p>Are you sure?</p>
<p></p>
</body></html>
'''

help_text = '''
<html><head/>
<body style="margin-top: 20px; margin-left: 20px; margin-right: 20px; margin-bottom: 20px">
<table>
    <tr>
    <td><img src="{}" width="128" height="128"></td>
    <td>&nbsp;&nbsp;&nbsp;&nbsp;</td>
    <td align="center">
        <h1>Interference calculator</h1>
        <p><strong>version {}</strong></p>
        <p><strong>&copy; 2017, Zan Peeters</strong></p>
        <p><a href="https://github.com/zanpeeters/interference_calculator">https://github.com/zanpeeters/interference_calculator</a></td>
    </tr>
</table>

<br/>
<p>This program can calculate possible mass interferences (molecules of similar mass)
for any molecule based on the composition of a given sample and a target mass. This
is useful for certain types of mass spectrometry such as SIMS or ICP-MS.</p>

<p>The atomic data used in this program is taken from the International Union for
Physical and Applied Chemistry (IUPAC), Commission on Isotopic Abundances and Atomic
Weights (CIAAW).

<ul>
<li>Atomic masses
    <p><a href=http://ciaaw.org/atomic-masses.htm>http://ciaaw.org/atomic-masses.htm</a></p>
    <p>Wang et al, The Ame2012 atomic mass evaluation, <i>Chinese Physics C</i>, <b>2012</b>, 36, 1603&ndash;2014
     <a href=https://doi.org/10.1088/1674-1137/36/12/003>https://doi.org/10.1088/1674-1137/36/12/003</a></p>

<li>Isotopic abundances
    <p><a href=https://www.degruyter.com/table/j/pac.2016.88.issue-3/pac-2015-0503/pac-2015-0503.xml?id=j_pac-2015-0503_tab_001>https://www.degruyter.com/table/j/pac.2016.88.issue-3/pac-2015-0503/pac-2015-0503.xml?id=j_pac-2015-0503_tab_001</a></p>
     <p>Meija et al, Isotopic compositions of the elements 2013, <i>Pure and Applied Chemistry</i>, <b>2016</b>, 88, 293&ndash;306 (table 1), <a href=https://doi.org/10.1515/pac-2015-0503>https://doi.org/10.1515/pac-2015-0503</a></p>

<li>Electron mass
    <p><a href=http://physics.nist.gov/cgi-bin/cuu/Value?me>http://physics.nist.gov/cgi-bin/cuu/Value?me</a></p>
    <p>Mohr et al, CODATA recommended values of the fundamental physical constants: 2014, <i>arXiv</i>, <b>2015</b>, <a href=https://arxiv.org/pdf/1507.07956.pdf>https://arxiv.org/pdf/1507.07956.pdf</a></p>
</ul>
<br/>
</body></html>
'''
