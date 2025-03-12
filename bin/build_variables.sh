#!/usr/bin/env bash
#
# Create a local directory .local to be used as --prefix when building
# local installs of python-flint's dependencies. This also sets the PREFIX
# shell variable and environment variables giving the versions to use for each
# dependency. This script should be sourced rather than executed e.g.:
#
#    $ source bin/build_variables.sh
#
# This is used implicitly by the other build scripts and does not need to be
# executed directly.

PREFIX=$(pwd)/.local
mkdir -p $PREFIX

ARBVER=2.23.0 # Not needed with flint >= 3.0.0 (Arb is included in flint)

YASMVER=1.3.0 # Only needed for MPIR
MPIRVER=3.0.0 # MPIR build no longer works (not clear where to download from)

# These are the actual dependencies used (at least by default):
GMPVER=6.3.0
MPFRVER=4.1.0
FLINTVER='3.2.0'
