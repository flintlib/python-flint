#!/usr/bin/env python3

from collections import defaultdict
import re
import sys
import os
import argparse
from pathlib import Path

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

For FLINT's Markdown documentation, add --use-flint-parser and point
--flint-doc-dir at the checkout's doc/src directory. The parser is loaded
from the adjacent doc/fdoc package. Without the flag, RST parsing is unchanged.


You will need to configure the location of the flintlib submodule directory and the locations of the flint and
arb doc source directories.

The other useful configuration is the comment_types and set which are two representations of the types that are
not implemented and therefore we want to comment functions that reference them

TODO: DRY with comment types, also should be able to add commented types from the command line.
TODO: don't import self

"""

#  recognize a function definition in rst
is_func = re.compile(r"\.\.( )+(c:)?function( )*::")
is_macro = re.compile(r"\.\.( )+(c:)?macro( )*::")
typed_macro = re.compile(r"^.+[\s*][A-Za-z_]\w*\s*\(")
# rename types to avoid python -- c name collisions
rename_types = [
    (re.compile(r"\bfmpz\b"),"fmpz_struct"),
    (re.compile(r"\bfmpq\b"), "fmpq_struct"),
    (re.compile(r"\bin\b"), "in_"),
    (re.compile(r"\blambda\b"), "lambda_"),
]
# comment out functions which use these types
comment_types = re.compile(r"(\bFILE\b)|(\bmpz_t\b)|(\bmpq_t\b)")
comment_set = set(["FILE", "mpz_t", "mpq_t"])
c_types = set(["void", "char", "short", "long", "int", "float", "double"])
type_modifers = re.compile(r"\*|(\bconst\b)|(\bunsigned\b)|(\bsigned\b)")
import_dict = {}
# gr_domains.rst lists functions that are in gr.h
doc_to_header = {'flint/gr_domains': 'flint/gr'}
void_parameter_list = re.compile(r"\(\s*void\s*\)")
macro_declarations = {
    'flint/mag': {
        'MAG_BITS': 'cdef const int MAG_BITS',
    },
    'flint/arf': {
        'ARF_PREC_EXACT': 'cdef const slong ARF_PREC_EXACT',
    },
    'flint/fmpq': {
        'fmpq_numref': 'fmpz * fmpq_numref(const fmpq_t x)',
        'fmpq_denref': 'fmpz * fmpq_denref(const fmpq_t x)',
    },
    'flint/fmpq_poly': {
        'fmpq_poly_numref': 'fmpz * fmpq_poly_numref(fmpq_poly_t poly)',
        'fmpq_poly_denref': 'fmpz * fmpq_poly_denref(fmpq_poly_t poly)',
    },
}



def get_cython_struct_types(file):
    """
    Extract cython types from a pxd file.
    """
    ret = []
    for line in file:
        l = line.strip()
        if l[:8] == "ctypedef":
            if l[-1] == ']':
                # ctypedef foo foo_t[0]
                l = l[:l.rfind('[')]
            elif '(' in l:
                # ctypedef int (*foo_func)(...)
                l = l[l.index('('):].lstrip('(*')
                l = l[:l.index(')')]
            else:
                # ctypedef foo:
                l = l.strip(':')
            ret.append(l.split()[-1])
    return ret


def fill_import_dict(pyflintlibdir):
    """
    Get a map from cython structs to the pxd that defines them
    """
    import_dict.clear()
    import_dict['fmpq_struct'] = 'types.fmpq'

    type_dir = Path(pyflintlibdir) / 'types'
    for f in sorted(type_dir.glob("*.pxd")):
        with open(f) as pxd:
            for t in get_cython_struct_types(pxd):
                module = 'types.' + f.stem
                if t in import_dict:
                    raise RuntimeError(
                        f"Duplicate Cython type definition for {t}: "
                        f"{import_dict[t]} and {module}"
                    )
                import_dict[t] = module


def undecorate(str):
    """
    remove variable name, const, ``*``, etc. to just get types
    """
    ret = str.strip()
    if ' ' in ret:
        ret = ret[:ret.rfind(' ')]
    ret = re.sub(type_modifers, '', ret).strip()
    return ret


def get_parameter_types(str):
    params = str[str.find("(") + 1 : str.rfind(")")].split(",")
    ret_type = str.split('(')[0].rsplit(' ', 1)[0]
    params.append(ret_type)
    return [undecorate(s) for s in params if s]


def clean_types(function):
    ret = function.strip()
    for old, new in rename_types:
        ret = re.sub(old, new, ret)
    # This also applies to callback signatures nested in parameter lists.
    # Cython expects zero-argument functions to use () rather than (void).
    ret = re.sub(void_parameter_list, "()", ret)
    return ret


def get_functions(file):
    """
    Get a list of functions from an rst file
    """
    ret = []
    macros = []
    in_list = False
    for line in file:
        # Keep track of the macros
        # We want to give them types in cython...
        if is_macro.match(line):
            macros.append(line.strip())
            in_list = False
            continue

        if line.strip().endswith('\\'):
            continue

        m = is_func.match(line)
        if m:
            ret.append(clean_types(line[m.end():]))
            in_list = True
        else:
            if in_list:
                if line.strip() == '':
                    in_list = False
                else:
                    ret.append(clean_types(line))
    return ret, macros


def get_macro_declarations(macros, h_name):
    """Use documented types or explicit overrides; never infer macro types.

    Cython can call a function-like C macro using a function declaration.
    Untyped macros (which may also introduce variables or require lvalues)
    cannot in general be represented this way.
    """
    functions = []
    overrides = macro_declarations.get(h_name, {})
    # Some constants (e.g. MAG_BITS) are not documented in older RST docs.
    constants = [decl for decl in overrides.values() if decl.startswith("cdef ")]
    for macro in macros:
        signature = macro.split('::', 1)[1].strip()
        name = signature.split('(', 1)[0].strip()
        declaration = overrides.get(name)
        if declaration is not None:
            if not declaration.startswith('cdef '):
                functions.append(clean_types(declaration))
        elif typed_macro.match(signature):
            functions.append(clean_types(signature))
    return functions, constants


def get_all_types(function_list):
    ret = set()
    for f in function_list:
        for t in get_parameter_types(f):
            ret.add(t)
    return ret


def get_flint_functions(docdir, name):
    """Read a Markdown module using the parser from its FLINT checkout.

    docdir is FLINT's doc/src directory; fdoc lives alongside it.
    """
    docdir = Path(docdir).resolve()
    sys.path.insert(0, str(docdir.parent))
    try:
        from fdoc.markup import BlockParser, split_front_matter
    finally:
        sys.path.pop(0)

    path = docdir / name / 'index.md'
    metadata, body, line = split_front_matter(path.read_text())
    nodes = BlockParser(path=str(path)).parse(body, line)
    functions, macros = [], []
    for node in nodes:
        if node.kind != 'div':
            continue
        if node['kind'] == 'function':
            functions.extend(clean_types(p) for p in node['protos'])
        elif node['kind'] == 'macro':
            macros.extend('.. macro:: ' + p for p in node['protos'])
    return functions, macros, metadata.get('header')


def has_types(line, types):
    #ts = set(t for t in get_parameter_types(line) if t in types)
    #if ts:
    #    print(ts, file=sys.stderr)
    return any(t in types for t in get_parameter_types(line))


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
    for k, v in sorted(imports.items()):
        types = ", ".join(sorted(v))
        print("from flint.flintlib." + k + " cimport " + types)
    return sorted(ret)




def generate_pxd_file(h_name, opts):
    fill_import_dict(opts.flint_lib_dir)
    l = []
    docdir = opts.arb_doc_dir
    name = h_name
    h_name = doc_to_header.get(h_name, h_name)

    if name[:6] == "flint/":
        docdir = opts.flint_doc_dir
        name = name[6:]
    if opts.use_flint_parser:
        l, macros, header = get_flint_functions(docdir, name)
        if header:
            h_name = 'flint/' + header.removesuffix('.h')
    else:
        with open(os.path.join(docdir, name + ".rst")) as f:
            l, macros = get_functions(f)

    macro_functions, macro_defs = get_macro_declarations(macros, h_name)
    l.extend(macro_functions)
    if opts.use_flint_parser and name == 'gr_domains':
        # Context predicates moved from gr_domains to gr in the Markdown docs.
        # Preserve the import path used by bindings generated from older docs.
        print('from flint.flintlib.functions.gr cimport *')
    unknown_types = gen_imports(l)
    print()
    for t in unknown_types:
        print("# unknown type " + t)
    print()
    for m in macros:
        print("# " + m)
    print()
    print(r'cdef extern from "' + h_name + r'.h":')
    for declaration in macro_defs:
        print("    " + declaration)
    for f in l:
        if has_types(f, unknown_types):
            print("    # " + f)
        else:
            print("    " + f)


def main(*args):
    usage = """
    $ cd /path/to/python-flint
    $ bin/rst_to_pxd.py flint/fmpz --flint-doc-dir=/path/to/flint/doc/source
    $ bin/rst_to_pxd.py flint/fmpz --use-flint-parser --flint-doc-dir=/path/to/flint/doc/src
    """
    parser = argparse.ArgumentParser(description='Generate a pxd file from FLINT documentation',
                                     usage=usage)
    parser.add_argument('--flint-lib-dir', help='location of the flintlib submodule',
                        default="./src/flint/flintlib")
    parser.add_argument('--arb-doc-dir', help='location of the arb doc source directory',
                        default="/Users/davideinstein/projects/arb/doc/source")
    parser.add_argument('--flint-doc-dir', help='location of the flint doc source directory',
                        default="/Users/davideinstein/projects/flint2/doc/source")
    parser.add_argument('--use-flint-parser', action='store_true',
                        help="parse Markdown using FLINT's doc/fdoc parser; "
                             "--flint-doc-dir must point to FLINT's doc/src")
    parser.add_argument('name', help='documentation module to parse (e.g. flint/fmpz)')

    args = parser.parse_args()

    generate_pxd_file(args.name, args)


if __name__ == "__main__":
    main(*sys.argv[1:])
