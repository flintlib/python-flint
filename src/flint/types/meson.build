thisdir = 'flint/types'

pyfiles = [
  '__init__.py',

  'fmpz.pyi',
  'fmpz_poly.pyi',
  'fmpz_series.pyi',
  'fmpz_mpoly.pyi',
  # 'fmpz_mat.pyi',
  # 'fmpz_vec.pyi',

  'fmpq.pyi',
  'fmpq_poly.pyi',
  'fmpq_series.pyi',
  'fmpq_mpoly.pyi',
  # 'fmpq_mat.pyi',
  # 'fmpq_vec.pyi',

  'nmod.pyi',
  'nmod_poly.pyi',
  'nmod_mpoly.pyi',
  # 'nmod_mat.pyi',
  # 'nmod_series.pyi',

  'fmpz_mod.pyi',
  'fmpz_mod_poly.pyi',
  'fmpz_mod_mpoly.pyi',
  # 'fmpz_mod_mat.pyi',

  'fq_default.pyi',
  'fq_default_poly.pyi',

  'arf.pyi',
  'arb.pyi',
  'acb.pyi',
  #
  'arb_poly.pyi',
  # 'arb_mat.pyi',
  'arb_series.pyi',
  #
  'acb_poly.pyi',
  # 'acb_mat.pyi',
  'acb_series.pyi',
  #
  # 'dirichlet.pyi',
  #
  # '_gr.pyi',
]

exts = [
  'fmpz',
  'fmpz_poly',
  'fmpz_mpoly',
  'fmpz_mat',
  'fmpz_series',
  'fmpz_vec',

  'fmpq',
  'fmpq_poly',
  'fmpq_mpoly',
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
  'fmpz_mod_mpoly',
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
    limited_api: limited_api_version,
  )
endforeach
