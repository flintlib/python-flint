#!/usr/bin/env python
"""
Update license information in wheels after running auditwheel etc.

Usage:
   cibw_repair_wheel_licenses.py <wheel_file> \
       --license MIT --license-file licenses/license_foo_MIT.txt \
       --license BSD-3-Clause --license-file licenses/license_bar_BSD-3-Clause.txt \
       ...

CWD should be at project root (containing pyproject.toml) and license-file
paths must be relative.

The wheel_file argument should be the path to the wheel file to be repaired.
It will be overwritten in place.
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
    license_files = [Path(f) for f in parsed.license_file]

    update_licenses_wheel(wheel_file, licenses, license_files)


def update_licenses_wheel(
    wheel_file: Path, licenses: list[str], license_files: list[Path]
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
    dist_info: Path, licenses: list[str], license_files: list[Path]
):
    for license_file in license_files:
        wheel_license_path = dist_info / "licenses" / license_file.name
        if wheel_license_path.exists():
            raise ValueError(f"license file already present: {license_file}")
        #
        # PEP 639 says:
        #
        # Inside the root license directory, packaging tools MUST reproduce the
        # directory structure under which the source license files are located
        # relative to the project root.
        #
        makedirs(dist_info / license_file.parent, exist_ok=True)
        copyfile(license_file, wheel_license_path)
        print(f"Added license file {license_file}")

    metadata_file = dist_info / "METADATA"

    with open(metadata_file, "r") as f:
        lines = f.readlines()

    for n, line in enumerate(lines):
        if line.startswith("License-Expression: "):
            base_license = line[len("License-Expression: ") :].strip()
            all_licenses = [base_license, *licenses]
            expression = ' AND '.join([f"({license})" for license in all_licenses])
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
