# Capsule-only Cython example

This extension uses `fmpz_add` through python-flint's capsule API. It neither
includes FLINT headers nor links to FLINT.

```console
python -m pip install .
python test_capi.py
```
