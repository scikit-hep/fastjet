#!/usr/bin/env python
# Copyright (c) 2021, Aryan Roy
#
# Distributed under the 3-clause BSD license, see accompanying file LICENSE
# or https://github.com/scikit-hep/fastjet for details.

from setuptools import setup  # isort:skip

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension  # isort:skip

import os
import pathlib
import shutil
import subprocess
import sys
import sysconfig
import urllib.request
import zipfile

import setuptools.command.build_ext
import setuptools.command.install

CGAL_ZIP = (
    "https://github.com/CGAL/cgal/releases/download/v5.5.1/CGAL-5.5.1-library.zip"
)

DIR = pathlib.Path(__file__).parent.resolve()
FASTJET = DIR / "fastjet-core"
FASTJET_CONTRIB = DIR / "fastjet-contrib"
PYTHON = DIR / "src/fastjet"
OUTPUT = PYTHON / "_fastjet_core"

LIBS = [
    "fastjet",
    "fastjettools",
    "siscone",
    "siscone_spherical",
    "fastjetplugins",
    "fastjetcontribfragile",
]


def get_version() -> str:
    g = {}
    with open(PYTHON / "version.py") as f:
        exec(f.read(), g)
    return g["__version__"]


class FastJetBuild(setuptools.command.build_ext.build_ext):
    def build_extensions(self):
        if not OUTPUT.exists():
            zip_filename = DIR / pathlib.Path(CGAL_ZIP).parts[-1]

            with urllib.request.urlopen(CGAL_ZIP) as http_obj:
                with open(zip_filename, "wb") as file_obj:
                    shutil.copyfileobj(http_obj, file_obj)

            with zipfile.ZipFile(zip_filename) as zip_obj:
                cgal_dir = DIR / zip_obj.namelist()[0]
                zip_obj.extractall(DIR)

            # Patch for FastJet core version 3.4.0
            # To be removed when https://gitlab.com/fastjet/fastjet/-/merge_requests/1 is merged upstream
            subprocess.run(
                ["patch", "pyinterface/fastjet.i", DIR / "patch_fastjet_i.txt"],
                cwd=FASTJET,
            )

            # Patch for segfault of LimitedWarning
            # For more info see https://github.com/scikit-hep/fastjet/pull/131
            subprocess.run(
                ["patch", "src/ClusterSequence.cc", DIR / "patch_clustersequence.txt"],
                cwd=FASTJET,
            )

            env = os.environ.copy()
            env["PYTHON"] = sys.executable
            env["PYTHON_INCLUDE"] = f'-I{sysconfig.get_path("include")}'
            env["CXXFLAGS"] = "-O3 -Bstatic -lgmp -Bdynamic"
            env["ORIGIN"] = "$ORIGIN"  # if evaluated, it will still be '$ORIGIN'

            args = [
                f"--prefix={OUTPUT}",
                "--enable-allcxxplugins",
                "--enable-cgal-header-only",
                "--enable-cgal",
                f"--with-cgaldir={cgal_dir}",
                "--enable-swig",
                "--enable-pyext",
                "LDFLAGS=-Wl,-rpath=$$ORIGIN/_fastjet_core/lib:$$ORIGIN",
            ]

            try:
                subprocess.run(
                    ["./autogen.sh"] + args, cwd=FASTJET, env=env, check=True
                )
            except Exception:
                subprocess.run(["cat", "config.log"], cwd=FASTJET, check=True)
                raise

            env = os.environ.copy()
            env["ORIGIN"] = "$ORIGIN"  # if evaluated, it will still be '$ORIGIN'
            subprocess.run(["make", "-j"], cwd=FASTJET, env=env, check=True)
            subprocess.run(["make", "install"], cwd=FASTJET, env=env, check=True)

            subprocess.run(
                ["./configure", f"--fastjet-config={FASTJET}/fastjet-config"],
                cwd=FASTJET_CONTRIB,
                env=env,
                check=True,
            )
            subprocess.run(["make", "-j"], cwd=FASTJET_CONTRIB, env=env, check=True)
            subprocess.run(
                ["make", "install"], cwd=FASTJET_CONTRIB, env=env, check=True
            )
            subprocess.run(
                ["make", "fragile-shared"], cwd=FASTJET_CONTRIB, env=env, check=True
            )
            subprocess.run(
                ["make", "fragile-shared-install"],
                cwd=FASTJET_CONTRIB,
                env=env,
                check=True,
            )

        setuptools.command.build_ext.build_ext.build_extensions(self)


class FastJetInstall(setuptools.command.install.install):
    def run(self):
        fastjetdir = pathlib.Path(f"{self.build_lib}/fastjet")

        shutil.copytree(OUTPUT, fastjetdir / "_fastjet_core", symlinks=True)

        pythondir = pathlib.Path(
            subprocess.check_output(
                """make -f pyinterface/Makefile --eval='print-pythondir:
\t@echo $(pythondir)
' print-pythondir""",
                shell=True,
                cwd=FASTJET,
                universal_newlines=True,
            ).strip()
        )

        pyexecdir = pathlib.Path(
            subprocess.check_output(
                """make -f pyinterface/Makefile --eval='print-pyexecdir:
\t@echo $(pyexecdir)
' print-pyexecdir""",
                shell=True,
                cwd=FASTJET,
                universal_newlines=True,
            ).strip()
        )

        shutil.copyfile(pythondir / "fastjet.py", fastjetdir / "_swig.py")
        shutil.copyfile(pyexecdir / "_fastjet.so", fastjetdir / "_fastjet.so")

        setuptools.command.install.install.run(self)


ext_modules = [
    Pybind11Extension(
        "fastjet._ext",
        ["src/_ext.cpp"],
        cxx_std=11,
        include_dirs=[str(OUTPUT / "include")],
        library_dirs=[str(OUTPUT / "lib")],
        runtime_library_dirs=["$ORIGIN/_fastjet_core/lib"],
        libraries=LIBS,
    ),
]


setup(
    version=get_version(),
    ext_modules=ext_modules,
    cmdclass={"build_ext": FastJetBuild, "install": FastJetInstall},
)
