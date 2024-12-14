**_gr** -- generic rings (unstable interface)
===============================================================================

.. note::
   This module provides a preliminary experimental interface for FLINT's
   generic rings. This is largely untested and the interface should be
   considered unstable. In future the generic rings code will be integrated
   with the rest of python-flint and this module may be removed. For now it
   provides access to many FLINT types that are not wrapped yet in the rest of
   python-flint.

The generic rings module provides access to the generic rings interface in
FLINT. Usage of this interface consists of creating a context object to
represent a particular domain and then using that context object to create
elements of that domain. For example to create a context for polynomials in two
variables over the Gaussian integers :math:`\mathbb{Z}[i][x,y]` we would do::

    >>> from flint.types._gr import gr_fmpzi_ctx, gr_gr_mpoly_ctx
    >>> ctx = gr_gr_mpoly_ctx.new(gr_fmpzi_ctx, ["x", "y"])
    >>> ctx.gens()
    [x, y]

    # XXX: gens_recursive not available in FLINT < 3.1
    # >>> ctx.gens_recursive()
    # [I, x, y]
    # >>> I, x, y = ctx.gens_recursive()

    >>> x, y = ctx.gens()
    >>> p = (x + y)**2
    >>> p
    x^2 + 2*x*y + y^2

Some domains such as ``gr_fmpzi_ctx`` are global and do not need to be created.
Others such as ``gr_gr_mpoly_ctx`` are created using :meth:`gr_ctx.new`.

.. autoclass :: flint.types._gr.gr_ctx
  :members:
  :undoc-members:

.. autoclass :: flint.types._gr.gr_scalar_ctx
  :members:
  :undoc-members:

.. autoclass :: flint.types._gr.gr_poly_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_mpoly_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr._gr_fmpz_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr._gr_fmpq_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr._gr_fmpzi_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr._gr_fexpr_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_nmod_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_fmpz_mod_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_fq_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_fq_nmod_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_fq_zech_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_nf_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_nf_fmpz_poly_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_real_qqbar_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_complex_qqbar_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_real_ca_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_complex_ca_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_real_algebraic_ca_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_complex_algebraic_ca_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_complex_extended_ca_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_real_float_arf_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_complex_float_acf_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_real_arb_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_complex_acb_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_gr_poly_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_gr_mpoly_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr_series_ctx
   :members:
   :undoc-members:

.. autoclass :: flint.types._gr.gr
   :members:
   :inherited-members:
   :undoc-members:
