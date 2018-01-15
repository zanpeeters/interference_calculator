# Interference calculator

![interference calculator logo](icon.svg)

Interference calculator calculates all molecules that can be formed (the interferences) from a combination of a list of atoms (sample composition), given a target mass and range. The calculation considers all isotopes of the sample atoms and build molecules up to a given size. The results are displayed in a table that can be sorted and copied, as well as a mass spectrum (click button ▶︎ to show).

The program can also display the standard ratios of the isotopes for any given element. The results include the natural abundance, standard ratio, inverse ratio, and the standard material in which the isotopic ratios where measured.

## Run the GUI.

To run the GUI, simply call the ui.py script:

```python
>>> python interference_calculator/ui.py
```

## From another script.

To call the calculations from the interactive interpreter, use:

```python
>>> import interference_calculator as ic
>>> ic.interference(['Ca', 'O', 'H', 'Si'], 'Fe')
>>> ic.standard_ratio(['Ca', 'O'])
```

See the docstrings for detailed help and options.
