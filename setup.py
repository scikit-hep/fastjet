#!/usr/bin/env python
# Copyright (c) 2021, Aryan Roy
#
# Distributed under the 3-clause BSD license, see accompanying file LICENSE
# or https://github.com/scikit-hep/fastjet for details.

from setuptools import setup  # isort:skip

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension  # isort:skip

# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

import glob
import os.path
import pathlib
import shutil
import subprocess
import sys
import sysconfig

import setuptools.command.build_ext
import setuptools.command.install

# import urllib.request
# import zipfile


CGAL_ZIP = (
    "https://github.com/CGAL/cgal/releases/download/v5.2.1/CGAL-5.2.1-library.zip"
)
cgal_dirname = "CGAL-5.2.1"

DIR = pathlib.Path(__file__).parent.resolve()
FASTJET = DIR / "fastjet-core"
PYTHON = DIR / "src" / "fastjet"


def get_version():
    g = {}
    with open(os.path.join("src", "fastjet", "version.py")) as f:
        exec(f.read(), g)
    return g["__version__"]


class FastJetBuild(setuptools.command.build_ext.build_ext):
    def build_extensions(self):
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        if not (PYTHON / "_fastjet_core").exists():
            zip_filename = DIR / pathlib.Path(CGAL_ZIP).parts[-1]

            # with urllib.request.urlopen(CGAL_ZIP) as http_obj:
            #     print(f"Downloading {CGAL_ZIP} to {zip_filename}")  # noqa T001
            #     with open(zip_filename, "wb") as file_obj:
            #         shutil.copyfileobj(http_obj, file_obj)

            print(f"Downloading {CGAL_ZIP} to {zip_filename}")  # noqa T001
            subprocess.run(["wget", CGAL_ZIP], cwd=DIR, check=True)

            # with zipfile.ZipFile(zip_filename) as zip_obj:
            #     cgal_dirname = zip_obj.namelist()[0]
            #     print(  # noqa T001
            #         f"Unzipping {zip_filename} to {str(DIR / cgal_dirname)}"
            #     )
            #     zip_obj.extractall(DIR)

            print(f"Unzipping {zip_filename} to {str(DIR / cgal_dirname)}")  # noqa T001
            subprocess.run(["unzip", str(zip_filename)], cwd=DIR, check=True)

            subprocess.run(
                [
                    "chmod",
                    "664",
                    str(
                        DIR
                        / cgal_dirname
                        / "include"
                        / "CGAL"
                        / "Exact_predicates_inexact_constructions_kernel.h"
                    ),
                ],
                cwd=DIR,
                check=True,
            )

            subprocess.run(
                [
                    "ls",
                    "-l",
                    str(
                        DIR
                        / cgal_dirname
                        / "include"
                        / "CGAL"
                        / "Exact_predicates_inexact_constructions_kernel.h"
                    ),
                ],
                cwd=DIR,
                check=True,
            )

            env = os.environ.copy()
            env["NOCONFIGURE"] = "1"
            env["PYTHON"] = sys.executable
            env["PYTHON_INCLUDE"] = f'-I{sysconfig.get_path("include")}'
            env["CXXFLAGS"] = "-O3 -Bstatic -lgmp -lgfortran -Bdynamic"
            if sys.platform.startswith("darwin"):
                env["FC"] = "gfortran"

            print("Running autogen.sh")  # noqa T001
            subprocess.run(["./autogen.sh"], cwd=FASTJET, env=env, check=True)

            args = [
                f"--prefix={str(PYTHON / '_fastjet_core')}",
                "--enable-allplugins",
                "--enable-cgal",
                "--enable-cgal-header-only",
                f"--with-cgaldir={str(DIR / cgal_dirname)}",
                "--enable-swig",
                "--enable-pyext",
            ]
            print("Running configure " + " ".join(args))  # noqa T001
            subprocess.run(["./configure"] + args, cwd=FASTJET, check=True, env=env)

            subprocess.run(["make", "-j"], cwd=FASTJET, check=True)
            subprocess.run(["make", "install"], cwd=FASTJET, check=True)

            for pythondir in glob.glob(
                str(PYTHON / "_fastjet_core" / "lib" / "python*")
            ):
                pythondir = pathlib.Path(pythondir)
                shutil.copyfile(
                    pythondir / "site-packages" / "fastjet.py", PYTHON / "_swig.py"
                )
                for sharedobj in glob.glob(str(pythondir / "site-packages" / "*.so*")):
                    sharedobj = pathlib.Path(sharedobj)
                    shutil.copyfile(sharedobj, PYTHON / sharedobj.parts[-1])

        setuptools.command.build_ext.build_ext.build_extensions(self)


class FastJetInstall(setuptools.command.install.install):
    def run(self):
        outerdir = pathlib.Path(
            f"build/lib.{sysconfig.get_platform()}-{'.'.join(map(str, sys.version_info[0:2]))}"
        )

        shutil.copytree(
            PYTHON / "_fastjet_core", outerdir / "fastjet" / "_fastjet_core"
        )

        for pythondir in glob.glob(
            str(outerdir / "fastjet" / "_fastjet_core" / "lib" / "python*")
        ):
            pythondir = pathlib.Path(pythondir)
            (pythondir / "site-packages" / "fastjet.py").rename(
                outerdir / "fastjet" / "_swig.py"
            )
            for sharedobj in glob.glob(str(pythondir / "site-packages" / "*.so*")):
                sharedobj = pathlib.Path(sharedobj)
                sharedobj.rename(outerdir / "fastjet" / sharedobj.parts[-1])
            shutil.rmtree(pythondir)

        setuptools.command.install.install.run(self)


ext_modules = [
    Pybind11Extension(
        "fastjet._ext",
        ["src/_ext.cpp"],
        cxx_std=11,
    ),
]


setup(
    version=get_version(),
    ext_modules=ext_modules,
    cmdclass={"build_ext": FastJetBuild, "install": FastJetInstall},
)
