"""
A Cython plugin for coverage.py suitable for a spin/meson project.

This is derived from Cython's coverage plugin.

https://coverage.readthedocs.io/en/latest/api_plugin.html
https://github.com/cython/cython/blob/master/Cython/Coverage.py
"""
import re
from collections import defaultdict

from coverage.plugin import CoveragePlugin, FileTracer, FileReporter

from functools import cache
from pathlib import Path


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
def parse_all_cfile_lines(exclude_lines):
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

        cfile_lines = parse_cfile_lines(c_file, exclude_lines)

        for cython_file, line_map in cfile_lines.items():
            if cython_file == '(tree fragment)':
                continue
            elif cython_file in all_code_lines:
                # Possibly need to merge the lines?
                assert all_code_lines[cython_file] == line_map
            else:
                all_code_lines[cython_file] = line_map

    return all_code_lines


def parse_cfile_lines(c_file, exclude_lines):
    """Parse a C file and extract all source file lines."""
    #
    # The C code has comments that refer to the Cython source files. We want to
    # parse those comments and match them up with the __Pyx_TraceLine() calls
    # in the C code. The __Pyx_TraceLine calls generate the trace events which
    # coverage feeds through to our plugin. If we can pair them up back to the
    # Cython source files then we measure coverage of the original Cython code.
    #
    match_source_path_line = re.compile(r' */[*] +"(.*)":([0-9]+)$').match
    match_current_code_line = re.compile(r' *[*] (.*) # <<<<<<+$').match
    match_comment_end = re.compile(r' *[*]/$').match
    match_trace_line = re.compile(r' *__Pyx_TraceLine\(([0-9]+),').match
    not_executable = re.compile(
        r'\s*c(?:type)?def\s+'
        r'(?:(?:public|external)\s+)?'
        r'(?:struct|union|enum|class)'
        r'(\s+[^:]+|)\s*:'
    ).match

    # Exclude e.g. # pragma: nocover
    exclude_pats = [f"(?:{regex})" for regex in exclude_lines]
    line_is_excluded = re.compile("|".join(exclude_pats)).search

    code_lines = defaultdict(dict)
    executable_lines = defaultdict(set)
    current_filename = None

    with open(c_file) as lines:
        lines = iter(lines)
        for line in lines:
            match = match_source_path_line(line)
            if not match:
                if '__Pyx_TraceLine(' in line and current_filename is not None:
                    trace_line = match_trace_line(line)
                    if trace_line:
                        executable_lines[current_filename].add(int(trace_line.group(1)))
                continue
            filename, lineno = match.groups()
            current_filename = filename
            lineno = int(lineno)
            for comment_line in lines:
                match = match_current_code_line(comment_line)
                if match:
                    code_line = match.group(1).rstrip()
                    if not_executable(code_line):
                        break
                    if line_is_excluded(code_line):
                        break
                    code_lines[filename][lineno] = code_line
                    break
                elif match_comment_end(comment_line):
                    # unexpected comment format - false positive?
                    break

    exe_code_lines = {}

    for fname in code_lines:
        # Remove lines that generated code but are not traceable.
        exe_lines = set(executable_lines.get(fname, ()))
        line_map = {n: c for n, c in code_lines[fname].items() if n in exe_lines}
        exe_code_lines[fname] = line_map

    return exe_code_lines


class Plugin(CoveragePlugin):
    """
    A Cython coverage plugin for coverage.py suitable for a spin/meson project.
    """
    def configure(self, config):
        """Configure the plugin based on .coveragerc/pyproject.toml."""
        # Read the regular expressions from the coverage config
        self.exclude_lines = tuple(config.get_option("report:exclude_lines"))

    def file_tracer(self, filename):
        """Find a tracer for filename as reported in trace events."""
        # All sorts of paths come here and we discard them if they do not begin
        # with the path to this directory. Otherwise we return a tracer.
        srcfile = self.get_source_file_tracer(filename)

        if srcfile is None:
            return None

        return MyFileTracer(srcfile)

    def file_reporter(self, filename):
        """Return a file reporter for filename."""
        srcfile = self.get_source_file_reporter(filename)

        return MyFileReporter(srcfile, exclude_lines=self.exclude_lines)

    #
    # It is important not to mix up get_source_file_tracer and
    # get_source_file_reporter. On the face of it these two functions do the
    # same thing i.e. you give a path and they return a path relative to src.
    # However the inputs they receive are different. For get_source_file_tracer
    # the inputs are semi-garbage paths from coverage. In particular the Cython
    # trace events use src-relative paths but coverage merges those with CWD to
    # make paths that look absolute but do not really exist. The paths sent to
    # get_source_file_reporter come indirectly from
    # MyFileTracer.dynamic_source_filename which we control and so those paths
    # are real absolute paths to the source files in the src dir.
    #
    # We make sure that get_source_file_tracer is the only place that needs to
    # deal with garbage paths. It also needs to filter out all of the
    # irrelevant paths that coverage sends our way. Once that data cleaning is
    # done we can work with real paths sanely.
    #

    def get_source_file_tracer(self, filename):
        """Map from coverage path to srcpath."""
        path = Path(filename)

        if build_install_dir in path.parents:
            # A .py file in the build-install directory.
            return self.get_source_file_build_install(path)
        elif root_dir in path.parents:
            # A .pyx file from the src directory. The path has src
            # stripped out and is not a real absolute path but it looks
            # like one. Remove the root prefix and then we have a path
            # relative to src_dir.
            return path.relative_to(root_dir)
        else:
            return None

    def get_source_file_reporter(self, filename):
        """Map from coverage path to srcpath."""
        path = Path(filename)

        if build_install_dir in path.parents:
            # A .py file in the build-install directory.
            return self.get_source_file_build_install(path)
        else:
            # An absolute path to a file in src dir.
            return path.relative_to(src_dir)

    def get_source_file_build_install(self, path):
        """Get src-relative path for file in build-install directory."""
        # A .py file in the build-install directory. We want to find its
        # relative path from the src directory. One of path.parents is on
        # sys.path and the relpath from there is also the relpath from src.
        for pkgdir in path.parents:
            init = pkgdir / '__init__.py'
            if not init.exists():
                sys_path_dir = pkgdir
                return path.relative_to(sys_path_dir)
        assert False


class MyFileTracer(FileTracer):
    """File tracer for Cython or Python files (.pyx,.pxd,.py)."""

    def __init__(self, srcpath):
        assert (src_dir / srcpath).exists()
        self.srcpath = srcpath

    def source_filename(self):
        return self.srcpath

    def has_dynamic_source_filename(self):
        return True

    def dynamic_source_filename(self, filename, frame):
        """Get filename from frame and return abspath to file."""
        # What is returned here needs to match MyFileReporter.filename
        srcpath = frame.f_code.co_filename
        return self.srcpath_to_abs(srcpath)

    # This is called for every traced line. Cache it:
    @staticmethod
    @cache
    def srcpath_to_abs(srcpath):
        """Get absolute path string from src-relative path."""
        abspath = (src_dir / srcpath).absolute()
        assert abspath.exists()
        return str(abspath)


class MyFileReporter(FileReporter):
    """File reporter for Cython or Python files (.pyx,.pxd,.py)."""

    def __init__(self, srcpath, *, exclude_lines):
        abspath = (src_dir / srcpath)
        assert abspath.exists()

        # filepath here needs to match dynamic_source_filename
        filepath = str(abspath)
        super().__init__(filepath)

        self.srcpath = srcpath
        self.abspath = abspath
        self.exclude_lines = exclude_lines

    def relative_filename(self):
        """Path displayed in the coverage reports."""
        return str(self.srcpath)

    def lines(self):
        """Set of line numbers for possibly traceable lines."""
        if self.srcpath.suffix == '.py':
            line_map = self.get_pyfile_line_map()
        else:
            line_map = self.get_cyfile_line_map()
        return set(line_map)

    def get_pyfile_line_map(self):
        """Return all lines from .py file."""
        with open(self.abspath) as pyfile:
            line_map = dict(enumerate(pyfile))
        return line_map

    def get_cyfile_line_map(self):
        """Get all traceable code lines for this file."""
        srcpath = str(self.srcpath)
        all_line_maps = parse_all_cfile_lines(self.exclude_lines)
        line_map = all_line_maps[srcpath]
        return line_map


def coverage_init(reg, options):
    plugin = Plugin()
    reg.add_configurer(plugin)
    reg.add_file_tracer(plugin)
