# Cython API example

This extension uses python-flint's experimental Cython API to borrow the
underlying value of an `fmpz`, then calls the read-only FLINT function
`fmpz_bits`. The extension may link against a separate FLINT library, provided
its FLINT and GMP data ABIs are compatible with those used by python-flint.

```console
python -m pip install .
python test_cython_api.py
```
