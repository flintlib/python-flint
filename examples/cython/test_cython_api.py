from flint import fmpz
from _flint_cython_example import bit_length, copy


value = fmpz(2) ** 300 + 1
assert bit_length(value) == 301
value_copy = copy(value)
assert value_copy == value
assert value_copy is not value

try:
    bit_length(1)
except TypeError:
    pass
else:
    raise AssertionError("wrong argument type did not raise TypeError")
