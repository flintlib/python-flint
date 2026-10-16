from flint import fmpz
from _flint_direct_example import bit_length

value = fmpz(2) ** 300 + 1
assert bit_length(value) == 301
try:
    bit_length(1)
except TypeError:
    pass
else:
    raise AssertionError("wrong argument type did not raise TypeError")
