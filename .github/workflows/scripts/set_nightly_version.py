import argparse
import pathlib
import re


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("build_number")
    args = parser.parse_args()

    path = pathlib.Path("version.py")
    source, replacements = re.subn(
        r'^version = "(.+)"$',
        rf'version = "\g<1>.dev{args.build_number}"',
        path.read_text(encoding="utf-8"),
        flags=re.MULTILINE,
    )
    if replacements != 1:
        raise ValueError("Unable to find the version in version.py")
    path.write_text(source, encoding="utf-8")


if __name__ == "__main__":
    main()
