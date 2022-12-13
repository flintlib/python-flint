import sys
import doctest
import flint

sys.stdout.write("doctests...");
fail, total = doctest.testmod(flint._flint);
if fail == 0:
    print("OK")
else:
    raise AssertionError("%i of %i doctests failed" % (fail, total))
