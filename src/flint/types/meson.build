thisdir = 'flint/types'

pyfiles = [
  '__init__.py',
  'fq_default.pyi',
  'fmpz_mod.pyi',
  'fmpz_mpoly.pyi',
  'fmpq_mpoly.pyi',
  'fmpq_poly.pyi',
  'fq_default_poly.pyi',
  'nmod_mpoly.pyi',
  'fmpq_series.pyi',
  'fmpz_series.pyi',
  'fmpz.pyi',
  'nmod_poly.pyi',
  'fmpq.pyi',
  'fmpz_mod_mpoly.pyi',
  'nmod.pyi',
  'fmpz_mod_poly.pyi',
  'fmpz_poly.pyi',
]

exts = [
  'fmpz',
  'fmpz_poly',
  'fmpz_mat',
  'fmpz_series',
  'fmpz_vec',

  'fmpq',
  'fmpq_poly',
  'fmpq_mat',
  'fmpq_series',
  'fmpq_vec',

  'nmod',
  'nmod_poly',
  'nmod_mpoly',
  'nmod_mat',
  'nmod_series',

  'fmpz_mod',
  'fmpz_mod_poly',
  'fmpz_mod_mat',

  'fq_default',
  'fq_default_poly',

  'arf',

  'arb',
  'arb_poly',
  'arb_mat',
  'arb_series',

  'acb',
  'acb_poly',
  'acb_mat',
  'acb_series',

  'dirichlet',

  'fmpz_mpoly',
  'fmpz_mod_mpoly',
  'fmpq_mpoly',

  '_gr',
]

if have_acb_theta
  exts += ['acb_theta']
endif

py.install_sources(
  pyfiles,
  pure: false,
  subdir: thisdir,
)

foreach ext : exts
  py.extension_module(
    ext,
    ext + '.pyx',
    dependencies: pyflint_deps,
    install: true,
    subdir: thisdir,
  )
endforeach
