"""
A Cython plugin for coverage.py suitable for a spin/meson project.

This follows the same general approach as Cython's coverage plugin and uses the
Cython plugin for parsing the C files. The difference here is that files are
laid out very differently in a meson project. Assuming meson makes it a lot
easier to find all the C files because we can just parse the build.ninja file.

https://coverage.readthedocs.io/en/latest/api_plugin.html
https://github.com/cython/cython/blob/master/Cython/Coverage.py
"""
import re
from collections import defaultdict

from coverage.plugin import CoveragePlugin, FileTracer, FileReporter

from functools import cache
from pathlib import Path


# Paths used by spin/meson in a src-layout:
root_dir = Path(__file__).parent
build_dir = root_dir / 'build'
build_install_dir = root_dir / 'build-install'
src_dir = root_dir / 'src'


def get_ninja_build_rules():
    """Read all build rules from build.ninja."""
    rules = []
    with open(build_dir / 'build.ninja') as build_ninja:
        for line in build_ninja:
            line = line.strip()
            if line.startswith('build '):
                line = line[len('build '):]
                target, rule = line.split(': ')
                if target == 'PHONY':
                    continue
                compiler, *srcfiles = rule.split(' ')
                # target is a path relative to the build directory. We will
                # turn that into an absolute path so that all paths in target
                # and srcfiles are absolute.
                target = str(build_dir / target)
                rule = (target, compiler, srcfiles)
                rules.append(rule)
    return rules


def get_cython_build_rules():
    """Get all Cython build rules."""
    cython_rules = []

    for target, compiler, srcfiles in get_ninja_build_rules():
        if compiler == 'cython_COMPILER':
            assert target.endswith('.c')
            assert len(srcfiles) == 1 and srcfiles[0].endswith('.pyx')
            c_file = target
            [cython_file] = srcfiles
            cython_rules.append((c_file, cython_file))

    return cython_rules


@cache
def parse_all_cfile_lines():
    """Parse all generated C files from the build directory."""
    #
    # Each .c file can include code generated from multiple Cython files (e.g.
    # because of .pxd files) being cimported. Each Cython file can contribute
    # to more than one .c file. Here we parse all .c files and then collect
    # together all the executable lines from all of the Cython files into a
    # dict like this:
    #
    #   {filename: {lineno: linestr, ...}, ...}
    #
    # This function is cached because it only needs calling once and is
    # expensive.
    #
    all_code_lines = {}

    for c_file, _ in get_cython_build_rules():

        cfile_lines = parse_cfile_lines(c_file)

        for cython_file, line_map in cfile_lines.items():
            if cython_file == '(tree fragment)':
                continue
            elif cython_file in all_code_lines:
                # Possibly need to merge the lines?
                assert all_code_lines[cython_file] == line_map
            else:
                all_code_lines[cython_file] = line_map

    return all_code_lines


def parse_cfile_lines(c_file):
    """Use Cython's coverage plugin to parse the C code."""
    from Cython.Coverage import Plugin
    return Plugin()._parse_cfile_lines(c_file)


class Plugin(CoveragePlugin):
    """A coverage plugin for a spin/meson project with Cython code."""

    def file_tracer(self, filename):
        """Find a tracer for filename to handle trace events."""
        path = Path(filename)

        if path.suffix in ('.pyx', '.pxd') and root_dir in path.parents:
            # A .pyx file from the src directory. The path has src
            # stripped out and is not a real absolute path but it looks
            # like one. Remove the root prefix and then we have a path
            # relative to src_dir.
            srcpath = path.relative_to(root_dir)
            return CyFileTracer(srcpath)
        else:
            # All sorts of paths come here and we reject them
            return None

    def file_reporter(self, filename):
        """Return a file reporter for filename."""
        srcfile = Path(filename).relative_to(src_dir)
        return CyFileReporter(srcfile)


class CyFileTracer(FileTracer):
    """File tracer for Cython files (.pyx,.pxd)."""

    def __init__(self, srcpath):
        print(srcpath)
        assert (src_dir / srcpath).exists()
        self.srcpath = srcpath

    def source_filename(self):
        return self.srcpath

    def has_dynamic_source_filename(self):
        return True

    def dynamic_source_filename(self, filename, frame):
        """Get filename from frame and return abspath to file."""
        # What is returned here needs to match CyFileReporter.filename
        path = frame.f_code.co_filename
        return self.get_source_filename(path)

    # This is called for every traced line. Cache it:
    @staticmethod
    @cache
    def get_source_filename(filename):
        """Get src-relative path for filename from trace event."""
        path = src_dir / filename
        assert src_dir in path.parents
        assert path.exists()
        return str(path)


class CyFileReporter(FileReporter):
    """File reporter for Cython or Python files (.pyx,.pxd,.py)."""

    def __init__(self, srcpath):
        abspath = (src_dir / srcpath)
        assert abspath.exists()

        # filepath here needs to match dynamic_source_filename
        super().__init__(str(abspath))

        self.srcpath = srcpath

    def relative_filename(self):
        """Path displayed in the coverage reports."""
        return str(self.srcpath)

    def lines(self):
        """Set of line numbers for possibly traceable lines."""
        srcpath = str(self.srcpath)
        all_line_maps = parse_all_cfile_lines()
        line_map = all_line_maps[srcpath]
        return set(line_map)


def coverage_init(reg, options):
    plugin = Plugin()
    reg.add_configurer(plugin)
    reg.add_file_tracer(plugin)
