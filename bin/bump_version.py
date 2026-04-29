#!/usr/bin/env python

filenames_default = [
    "pyproject.toml",
    "src/flint/__init__.py",
    "doc/source/conf.py",
    "src/flint/test/test_all.py",
]


def main(version2=None, *filenames):
    """Bump version number in files.

        $ bin/bump_version.py
        Current version: 0.1.0

        $ bin/bump_version.py 0.1.0 0.1.1
        Set version 0.1.0 to 0.1.1 in:
        pyproject.toml
        src/flint/__init__.py
        doc/source/conf.py
        src/flint/test/test_all.py

    """
    with open("pyproject.toml", "r") as f:
        text = f.read()
    version1 = text.split("version = \"")[1].split("\"")[0]

    if not version2:
        print(f"Current version: {version1}")
        return

    if not filenames:
        filenames = filenames_default

    print(f"Set version {version1} to {version2} in:")
    for filename in filenames:
        print(filename)
        with open(filename, "r") as f:
            text = f.read()
        with open(filename, "w") as f:
            f.write(text.replace(version1, version2))


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
