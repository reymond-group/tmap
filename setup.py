#! /usr/bin/env python3

import multiprocessing
import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from shutil import copyfile, copymode


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
            "-DPYBIND11_CPP_STANDARD=/std:c++17",
        ]

        cfg = "Debug" if self.debug else "Release"
        print(f"Setup.py cfg: {cfg}")

        build_args = ["--config", cfg]
        num_cores = multiprocessing.cpu_count()

        env = os.environ.copy()
        env[
            "CXXFLAGS"
        ] = f'{env.get("CXXFLAGS", "")} -DVERSION_INFO="{self.distribution.get_version()}"'

        # enable post-command args
        build_args += ["--"]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            cmake_args += [
                "-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            cmake_args += ["-G", "Visual Studio 17 2022"]
            cmake_args += ["-A", "x64"]
            cmake_args += ["-T", "ClangCL"]

            # increase job count on windows
            build_args += ["/m"]
        elif platform.system() == "Darwin":
            cmake_args += ["-DOpenMP_C_FLAG=-fopenmp"]
            cmake_args += ["-DOpenMP_CXX_FLAG=-fopenmp"]
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]

            # increase job count on OSX
            build_args += [f"-j{num_cores}"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]

            # increase job count on linux
            build_args += [f"-j{num_cores}"]


        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )
        print()  # Add an empty line for cleaner output


setup(
    name="tmap-viz",
    version="1.0.18",
    author="Daniel Probst",
    author_email="daenuprobst@gmail.com",
    description="A Python package for visualizing large, high-dimensional data sets.",
    long_description="",
    packages=find_packages("src"),
    package_dir={"": "src"},
    ext_modules=[CMakeExtension("_tmap")],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
