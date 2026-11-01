import argparse
import pathlib
import subprocess
import sys
import tempfile
import zipfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--dest-dir", required=True)
    args = parser.parse_args()

    # delvewheel has to find the core library, which is inside of the wheel it repairs.
    with tempfile.TemporaryDirectory() as directory:
        zipfile.ZipFile(args.wheel).extract("networkit/networkit_state.dll", directory)
        subprocess.check_call(
            [
                sys.executable,
                "-m",
                "delvewheel",
                "repair",
                "-w",
                args.dest_dir,
                "-v",
                "--add-path",
                str(pathlib.Path(directory) / "networkit"),
                args.wheel,
            ]
        )


if __name__ == "__main__":
    main()
