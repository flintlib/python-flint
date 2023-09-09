#!/usr/bin/env python3

from collections import defaultdict
import re
import sys
import os
import fnmatch
import argparse

"""
This is relatively rudimentary, but works for the purpose intended

To use this to parse an arb rst file

python rst_to_pxd basename

or to parse a flint rst file

python rst_to_pxd flint/basename

the pxd file is dumped to stdout

for example
python rst_to_pxd flint/fmpz_poly

will output a skeleton fmpz_poly.pxd to stdout


You will need to configure the location of the flintlib submodule directory and the locations of the flint and
arb doc source directories.

The other useful configuration is the comment_types and set which are two representations of the types that are
not implemented and therefore we want to comment functions that reference them

TODO: DRY with comment types, also should be able to add commented types from the command line.
TODO: don't import self

"""

#  recognize a function definition in rst
is_func = re.compile(r"\.\.( )+(c:)?function( )*::")
# rename types to avoid python -- c name collisions
rename_types = [(re.compile(r"\bfmpz\b"),"fmpz_struct"),(re.compile(r"\bfmpq\b"), "fmpq_struct")]
# comment out functions which use these types
comment_types = re.compile(r"(\bFILE\b)|(\bmpz_t\b)|(\bmpq_t\b)")
comment_set = set(["FILE", "mpz_t", "mpq_t"])
c_types = set(["char", "short", "long", "int", "float", "double"])
type_modifers = re.compile(r"\*|(\bconst\b)|(\bunsigned\b)|(\bsigned\b)")
import_dict = {}

def get_cython_struct_types(file):
    """
    Extract cython types from a pxd file.
    """
    ret = []
    for line in file:
        l = line.strip()
        if l[:8] == "ctypedef":
            if l[-1] == ']':
                l = l[:l.rfind('[')]
            else:
                l = l.strip(':')
            ret.append(l.split()[-1])
    return ret

def fill_import_dict(pyflintlibdir):
    """
    Get a map from cython structs to the pxd that defines them
    """
    with os.scandir(pyflintlibdir) as entry:
        for f in entry:
            if fnmatch.fnmatch(f.name, "*.pxd"):
                with open(f.path) as pxd:
                    for t in get_cython_struct_types(pxd):
                        import_dict[t] = f.name.split('.')[0]

def undecorate(str):
    """
    remove variable name, const, *, etc. to just get types
    """
    ret = str.strip()
    ret = ret[:ret.rfind(' ')]
    ret = re.sub(type_modifers, '', ret)
    return ret.strip()

def get_parameter_types(str):
    params = str[str.find("(") + 1 : str.rfind(")")].split(",")
    return [undecorate(s) for s in params]

def clean_types(function):
    ret = function.strip()
    for old, new in rename_types:
        ret = re.sub(old, new, ret)
    return ret

def get_functions(file):
    """
    Get a list of functions from an rst file
    """
    ret = []
    in_list = False
    for line in file:
        m = is_func.match(line)
        if m:
            ret.append( clean_types(line[m.end():]))
            in_list = True
        else:
            if in_list:
                if line.strip() == '':
                    in_list = False
                else:
                    ret.append(clean_types(line))
    return ret

def get_all_types(function_list):
    ret = set()
    for f in function_list:
        for t in get_parameter_types(f):
            ret.add(t)
    return ret

def gen_imports(function_list):
    """
    Generate import statements for known functions.
    """
    imports = defaultdict(list)
    s = get_all_types(function_list)
    s = s - c_types
    ret = set([])
    for t in s:
        if t in import_dict:
            imports[import_dict[t]].append(t)
        else:
            ret.add(t)
    for k,v in imports.items():
        types = ", ".join(v)
        print("from flint.flintlib." + k + " cimport " + types)
    return ret

def generate_pxd_file(h_name, opts):
    fill_import_dict(opts.flint_lib_dir)
    l=[]
    docdir = opts.arb_doc_dir
    name = h_name 
    if name[:6] == "flint/":
        docdir = opts.flint_doc_dir
        name = name[6:]
    with open(os.path.join(docdir, name + ".rst")) as f:
        l = get_functions(f)
        s = gen_imports(l)
        print()
        print ("\n# unimported types ", s - comment_set)
        print()
        print(r'cdef extern from "' + h_name +r'.h":')
        for f in l:
            if comment_types.search(f):
                print("    # " + f)
            else:
                print("    " + f)


def main(*args):
    usage = """
    $ cd /path/to/python-flint
    $ bin/rst_to_pxd.py flint/fmpz --flint-doc-dir=/path/to/flint/doc/source
    """
    parser = argparse.ArgumentParser(description='Generate a pxd file from an rst file',
                                     usage=usage)
    parser.add_argument('--flint-lib-dir', help='location of the flintlib submodule',
                        default="./src/flint/flintlib")
    parser.add_argument('--arb-doc-dir', help='location of the arb doc source directory',
                        default="/Users/davideinstein/projects/arb/doc/source")
    parser.add_argument('--flint-doc-dir', help='location of the flint doc source directory',
                        default="/Users/davideinstein/projects/flint2/doc/source")
    parser.add_argument('name', help='name of the rst file to parse (e.g. flint/fmpz)')

    args = parser.parse_args()

    generate_pxd_file(args.name, args)


if __name__ == "__main__":
    main(*sys.argv[1:])

