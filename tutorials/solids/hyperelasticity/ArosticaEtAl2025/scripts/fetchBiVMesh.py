#!/usr/bin/env python3
"""
Fetch the biventricular mesh and fibres of the Arostica et al. (2025)
benchmark 2 into a case: constant/polyMesh and 0/f0, 0/s0, 0/n0.

They are converted once from the benchmark's own files (see the README.md and
convertBiVMesh.py in the archive) and kept, gzipped, in an archive on the
solids4foam website. The archive is downloaded once into data/, checked
against its SHA-256 sum, and the chosen mesh is copied into the case.
Without internet access, place ArosticaEtAl2025.zip in data/ first, or point
ARO_BIV_ARCHIVE to it.

Usage, from the case directory:
    python3 ../scripts/fetchBiVMesh.py coarse|fine
Needs only the Python 3 standard library.

Created by Claude Opus 5.5, supervised by Chanikya Valeti, UCD.
"""
import hashlib
import os
import shutil
import sys
import urllib.request
import zipfile

URL = "https://www.solids4foam.com/tutorials/archive/ArosticaEtAl2025.zip"
SHA256 = "43bd6f426f5afd0fdd54df347be3101b536ad57fbd01569fdfcca0b7bb07ed1e"
MESHES = ("coarse", "fine")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def archive():
    path = os.environ.get("ARO_BIV_ARCHIVE",
                          os.path.join("data", "ArosticaEtAl2025.zip"))
    if not (os.path.isfile(path) and sha256(path) == SHA256):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        print(f"Downloading {URL}")
        urllib.request.urlretrieve(URL, path)
        if sha256(path) != SHA256:
            os.remove(path)
            sys.exit("ArosticaEtAl2025.zip: SHA-256 mismatch, download rejected")
    return path


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in MESHES:
        sys.exit(__doc__)
    mesh = sys.argv[1]
    prefix = f"ArosticaEtAl2025/benchmark2/{mesh}/"

    with zipfile.ZipFile(archive()) as z:
        names = [n for n in z.namelist()
                 if n.startswith(prefix) and not n.endswith("/")]
        if not names:
            sys.exit(f"No {mesh} mesh in the archive")
        for n in names:
            dest = n[len(prefix):]
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with z.open(n) as src, open(dest, "wb") as dst:
                shutil.copyfileobj(src, dst)
            print(f"    {dest}")


if __name__ == "__main__":
    main()
