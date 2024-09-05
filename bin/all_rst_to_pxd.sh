#!/usr/bin/env bash

set -e

if [ $# -eq 0 ]
  then
    echo "Usage: bin/all_rst_to_pxd.sh /path/to/flint/doc/source"
    exit 1
fi

FLINT_DOC_DIR=$1

modules=(
    "acb_calc"
    "acb_dft"
    "acb_dirichlet"
    "acb_elliptic"
    "acb_hypgeom"
    "acb_mat"
    "acb_modular"
    "acb_poly"
    "acb"
    "acb_theta"
    "arb_fmpz_poly"
    "arb_hypgeom"
    "arb_mat"
    "arb_poly"
    "arb"
    "arf"
    "arith"
    "bernoulli"
    "dirichlet"
    #"fmpq_mat"
    #"fmpq_mpoly_factor"
    #"fmpq_mpoly"
    #"fmpq_poly"
    #"fmpq"
    #"fmpq_vec"
    "fmpz_factor"
    "fmpz_lll"
    "fmpz_mat"
    #"fmpz_mod_mat"
    #"fmpz_mod_mpoly_factor"
    #"fmpz_mod_mpoly"
    #"fmpz_mod_poly_factor"
    #"fmpz_mod_poly"
    #"fmpz_mod"
    #"fmpz_mod_vec"
    #"fmpz_mpoly_factor"
    #"fmpz_mpoly"
    #"fmpz_mpoly_q"
    "fmpz_poly_factor"
    "fmpz_poly"
    "fmpz"
    "fmpz_vec"
    #"fq_default_mat"
    #"fq_default_poly_factor"
    #"fq_default_poly"
    #"fq_default"
    #"fq_mat"
    #"fq_nmod_mat"
    #"fq_nmod_poly_factor"
    #"fq_nmod_poly"
    #"fq_nmod"
    #"fq_poly_factor"
    #"fq_poly"
    #"fq"
    #"fq_zech_mat"
    #"fq_zech_poly_factor"
    #"fq_zech_poly"
    #"fq_zech"
    "mag"
    #"mpoly"
    #"nmod_mat"
    #"nmod_mpoly_factor"
    #"nmod_mpoly"
    #"nmod_poly_factor"
    #"nmod_poly"
    #"nmod"
    #"nmod_vec"
    "partitions"
    "ulong_extras"
)

for module in ${modules[@]}; do
    echo "Processing $module"
    bin/rst_to_pxd.py flint/$module \
        --flint-doc-dir=$FLINT_DOC_DIR \
        > src/flint/flintlib/functions/$module.pxd
done
