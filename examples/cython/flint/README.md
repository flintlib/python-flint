# Direct-FLINT Cython example

This extension borrows python-flint's `fmpz` storage and calls the read-only
FLINT function `fmpz_bits`. A compatible system FLINT is an explicit build and
runtime dependency.

```console
python -m pip install .
python test_flint.py
```
