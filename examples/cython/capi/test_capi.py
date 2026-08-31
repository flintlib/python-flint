from flint import fmpz
from _flint_capi_example import add

a = fmpz(2) ** 300
b = fmpz(3) ** 200
assert add(a, b) == a + b
try:
    add(1, b)
except TypeError:
    pass
else:
    raise AssertionError("wrong argument type did not raise TypeError")
