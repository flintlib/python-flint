#!/usr/bin/env bash
#
# Set the PREFIX shell variable and environment variables giving the versions
# to use for each dependency. This script should be sourced rather than
# executed e.g.:
#
#    $ source bin/build_variables.sh
#
# This is used implicitly by the other build scripts and does not need to be
# executed directly.

PREFIX=$(pwd)/.local

# Versions of the dependencies used by python-flint:
GMPVER=6.3.0
MPFRVER=4.2.2
FLINTVER=3.6.0
