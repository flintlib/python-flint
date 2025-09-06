#!/usr/bin/env python
"""
Update license information in wheels after running auditwheel etc.

Usage:
   cibw_repair_wheel_licenses.py <wheel_file> \
       --license MIT \
       --license-file ./licenses/license_foo_MIT.txt:foo/LICENSE.txt \
       --license-file ./licenses/license_bar_MIT.txt:bar/LICENSE.txt \
       ...

The wheel_file argument should be the path to the wheel file to be repaired.
It will be overwritten in place.

The --license and --license-file arguments are multi-use. The --license
argument should be a SPDX license expression. It will be combined with the
existing License-Expression field in the wheel's METADATA file.

The --license-file argument should be a pair of :-separated paths. The first
path is the path to the license file and the second is the relative path to
put the license file under the wheel's .dist-info/licenses directory. The first
path can be a glob pattern. The second path can also be a glob pattern in which
case the * is replaced in the same way. This makes it possible to do

   --license-file ./src/foo-*/LICENSE:libs/foo-*/LICENSE

which would find the LICENSE file src/foo-1.0/LICENSE and copy it to the
matched .dist-info/licenses/libs/foo-1.0/LICENSE path in the wheel.

PEP 639 says:

Inside the root license directory, packaging tools MUST reproduce the
directory structure under which the source license files are located
relative to the project root.

It is not clear what that means though if the licenses are coming from code
that is not vendored in the repo/sdist. If they are vendored then presumably
the argument here should be like:

    --license-file ./vendored-foo/LICENSE:vendored-foo/LICENSE

"""
import argparse
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from shutil import copyfile
from os import makedirs


def main(*args: str):
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel_file")
    parser.add_argument("--license", action="append", default=[])
    parser.add_argument("--license-file", action="append", default=[])

    parsed = parser.parse_args(args)
    wheel_file = Path(parsed.wheel_file)
    licenses = parsed.license
    license_files = dict(arg_to_paths(f) for f in parsed.license_file)

    update_licenses_wheel(wheel_file, licenses, license_files)


def arg_to_paths(license_file_arg: str) -> tuple[Path, Path]:
    """
    Convert a --license-file argument to a pair of Paths.
    """
    paths = license_file_arg.strip().split(":")
    if len(paths) != 2:
        raise ValueError("license-file argument must be in the form of <path>:<path>"
                         f" but got {license_file_arg}")
    glob1_str, glob2_str = paths
    paths1 = list(Path().glob(glob1_str))
    if len(paths1) != 1:
        raise ValueError(f"Expected one path from glob pattern {glob1_str}"
                         f" but got {paths1}")
    [path1] = paths1

    if '*' not in glob2_str:
        path2_str = glob2_str
    else:
        # Replace * in glob2_str with the part of path1 that matches glob1_str:
        index1 = glob1_str.index('*')
        part1_glob = glob1_str[:index1]
        part2_glob = glob1_str[index1+1:]
        path1_str = str(path1)
        if len(part2_glob) != 0:
            wildcard = path1_str[len(part1_glob):-len(part2_glob)]
        else:
            wildcard = path1_str[len(part1_glob):]
        assert path1_str.startswith(part1_glob)
        assert path1_str.endswith(part2_glob)
        assert path1_str == part1_glob + wildcard + part2_glob
        path2_str = glob2_str.replace('*', wildcard)

    return path1, Path(path2_str)


def update_licenses_wheel(
    wheel_file: Path, licenses: list[str], license_files: dict[Path, Path]
):
    if wheel_file.exists():
        print("Found wheel at", wheel_file)
    else:
        raise ValueError(f"Wheel not found at {wheel_file}")

    for license_file in license_files:
        if license_file.is_absolute():
            raise ValueError("license-file paths must be relative to project root")
        elif not license_file.exists():
            raise ValueError(f"license-file not found: {license_file}")

    # foo/bar-1.0-cp310-cp310-linux_x86_64.whl -> bar-1.0
    name, version = wheel_file.stem.split('-')[:2]
    base = f"{name}-{version}"

    with TemporaryDirectory() as tmpdir:

        print("temp dir:", tmpdir)
        tmpdir = Path(tmpdir)

        run(["wheel", "unpack", "--dest", tmpdir, wheel_file], check=True)

        dist_info = tmpdir / base / f"{base}.dist-info"

        print(f"Adding licenses in {dist_info}")
        update_license_dist_info(dist_info, licenses, license_files)

        run(["wheel", "pack", "--dest-dir", tmpdir, tmpdir / base], check=True)

        # glob for *.whl in tmpdir
        wheels = list(tmpdir.glob(f"{base}-*.whl"))
        if len(wheels) != 1:
            raise ValueError(f"Expected one wheel in {tmpdir}, got {wheels}")
        new_wheel_file = wheels[0]

        print(f"Repaired wheel: {new_wheel_file}")
        print(f"Copying back to: {wheel_file}")
        copyfile(new_wheel_file, wheel_file)


def update_license_dist_info(
    dist_info: Path, licenses: list[str], license_files: dict[Path, Path]
):
    for src, dst in license_files.items():
        wheel_license_path = dist_info / "licenses" / dst
        if wheel_license_path.exists():
            raise ValueError(f"license file already present: {wheel_license_path}")
        #
        # PEP 639 says:
        #
        # Inside the root license directory, packaging tools MUST reproduce the
        # directory structure under which the source license files are located
        # relative to the project root.
        #
        makedirs(wheel_license_path.parent, exist_ok=True)
        copyfile(src, wheel_license_path)
        print(f"Copied license file {src} to {wheel_license_path}")

    metadata_file = dist_info / "METADATA"

    with open(metadata_file, "r") as f:
        lines = f.readlines()

    def brackets(s: str) -> str:
        if ' ' not in s:
            return s
        else:
            return f"({s})"

    for n, line in enumerate(lines):
        if line.startswith("License-Expression: "):
            base_license = line[len("License-Expression: ") :].strip()
            all_licenses = [base_license, *licenses]
            expression = ' AND '.join([brackets(license) for license in all_licenses])
            lines[n] = f"License-Expression: {expression}\n"
            break
    else:
        raise ValueError("Could not find License-Expression in METADATA")

    print("Updated License-Expression from")
    print("    " + base_license)
    print("to")
    print("    " + expression)

    license_files_lines = [line for line in lines if line.startswith("License-File: ")]

    if not license_files_lines:
        raise ValueError("Could not find License-File in METADATA")

    index = lines.index(license_files_lines[-1]) + 1
    new_lines = [f"License-File: {f}\n" for f in license_files]
    lines = lines[:index] + new_lines + lines[index:]

    print("Writing out METADATA with updated License-Expression and License-File fields")
    print("Writing to:", metadata_file)

    with open(metadata_file, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    import sys
    sys.exit(main(*sys.argv[1:]))
