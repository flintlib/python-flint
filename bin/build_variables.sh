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

GMPVER=6.2.1
YASMVER=1.3.0
MPIRVER=3.0.0
MPFRVER=4.1.0
FLINTVER=2.9.0
ARBVER=2.23.0
