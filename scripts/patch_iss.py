#!/usr/bin/env python3
"""Applies the InSilicoSeq patches of ft-db-exp's bin/install_tools.sh to an installed InSilicoSeq.

1. `--mode perfect` crashes in InSilicoSeq 2.0.1 (PerfectErrorModel lacks `store_mutations`), and
   `iss` hides that behind exit status 0 and no output.
2. The `perfect` model's read length (hardcoded 125) becomes ISS_PERFECT_READ_LENGTH, its insert size
   ISS_PERFECT_INSERT_SIZE (default twice the read length, at least 200).
3. The `basic` model's mean phred score (hardcoded 30) becomes ISS_BASIC_PHRED and its read length
   (hardcoded 125) ISS_BASIC_READ_LENGTH. A phred score q gives a per-base substitution rate of
   10^(-q/10).

Without the environment variables the models behave as before. Every patch is idempotent.

Run with the Python of the environment InSilicoSeq is installed in:
    python patch_iss.py
"""
import importlib.util
import os
import sys

# Locate the package without importing its error models, which pull in all of InSilicoSeq's dependencies.
spec = importlib.util.find_spec("iss")
if spec is None:
    sys.exit("InSilicoSeq not found - run this with the Python of the environment it is installed in")
base = os.path.join(os.path.dirname(spec.origin), "error_models")


def patch(name, replacements, marker, import_after):
    path = os.path.join(base, name)
    src = open(path, encoding="utf-8").read()
    if marker in src:
        print(f"{name}: already patched")
        return
    for old, new in replacements:
        if old not in src:
            sys.exit(f"{name}: does not look as expected, not patching (missing {old.strip()!r})")
        src = src.replace(old, new, 1)
    if "import os" not in src:
        src = src.replace(import_after, "import os\n\n" + import_after, 1)
    open(path, "w", encoding="utf-8").write(src)
    print(f"{name}: patched")


perfect = open(os.path.join(base, "perfect.py"), encoding="utf-8").read()
if "store_mutations" not in perfect:
    patch("perfect.py", [
        ("    def __init__(self, fragment_length=None, fragment_sd=None):",
         "    def __init__(self, fragment_length=None, fragment_sd=None, store_mutations=False):"),
        ("        self.fragment_sd = fragment_sd\n        self.quality_forward",
         "        self.fragment_sd = fragment_sd\n        self.store_mutations = store_mutations\n"
         "        self.quality_forward"),
    ], "self.store_mutations = store_mutations", "from iss.error_models import ErrorModel")
else:
    print("perfect.py: store_mutations already present")

patch("perfect.py", [
    ("        self.read_length = 125",
     '        self.read_length = int(os.environ.get("ISS_PERFECT_READ_LENGTH", "125"))'),
    ("        self.insert_size = 200",
     '        self.insert_size = int(os.environ.get("ISS_PERFECT_INSERT_SIZE",\n'
     '                                              str(max(200, 2 * self.read_length))))'),
], "ISS_PERFECT_READ_LENGTH", "from iss.error_models import ErrorModel")

patch("basic.py", [
    ("        self.quality_forward = self.quality_reverse = 30",
     '        self.quality_forward = self.quality_reverse = int(os.environ.get("ISS_BASIC_PHRED", "30"))'),
    ("        self.read_length = 125",
     '        self.read_length = int(os.environ.get("ISS_BASIC_READ_LENGTH", "125"))'),
], "ISS_BASIC_READ_LENGTH", "import numpy as np")
