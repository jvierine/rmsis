"""Compile unmodified, user-provided NRL sources plus our separate C ABI adapter.

Usage: python scripts/build_msis21.py /path/to/original/NRLMSIS2.1
Obtain the original distribution and review its license at
https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.1/
This script neither downloads nor modifies those files.
"""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir",type=Path)
    parser.add_argument("--compiler",default="gfortran")
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    source=args.source_dir.resolve()
    files=[source/f"msis_{name}.F90" for name in ("constants","utils","init","gfn","tfn","dfn","calc")]
    if not all(f.is_file() for f in files): parser.error("Original MSIS 2.1 source files missing")
    compiler=shutil.which(args.compiler)
    if compiler is None: parser.error("Fortran compiler not found")
    suffix=".dll" if sys.platform=="win32" else ".dylib" if sys.platform=="darwin" else ".so"
    output=(args.output or ROOT/"target"/f"librmsis21{suffix}").resolve()
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="rmsis21-build-") as build:
        command=[compiler,"-cpp","-DDBLE","-O2","-fPIC","-shared","-ffree-line-length-none",
                 *map(str,files),str(ROOT/"ffi/msis21_bridge.f90"),"-o",str(output)]
        subprocess.run(command,cwd=build,check=True)
    print(output)

if __name__=="__main__": main()
