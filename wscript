#!/usr/bin/env python
# encoding: utf-8

import os
from traceback import print_tb
from wafbuild.utils import load

VERSION = "1.0.0"
APPNAME = "gp-manifold"

srcdir = "src"
blddir = "build"
libdir = "gp_manifold"
libname = "GpManifold"

compiler = "cxx"
required = ["eigen", "mfem", "ipopt", "nlopt", "kernellib"]
optional = ["utilslib", "controllib", "graphicslib"]


def options(opt):
    # Add build shared library options
    opt.add_option("--shared",
                   action="store_true",
                   help="build shared library")

    # Add build static library options
    opt.add_option("--static",
                   action="store_true",
                   help="build static library")

    # Load library options
    load(opt, compiler, required, optional)

    # Load examples options
    opt.recurse("./src/examples")


def configure(cfg):
    # Load library configurations
    load(cfg, compiler, required, optional)

    # Load examples configurations
    cfg.recurse("./src/examples")

    cfg.env.LIB_MPI += ["mpi_usempif08", "mpi_usempi_ignore_tkr", "mpi_mpifh"]


def build(bld):
    # Library name
    bld.get_env()["libname"] = libname

    # Includes
    includes = []
    for root, _, filenames in os.walk(os.path.join(srcdir, libdir)):
        includes += [os.path.join(root, filename)
                     for filename in filenames if filename.endswith(('.hpp', '.h'))]

    # Sources
    sources = []
    for root, _, filenames in os.walk(os.path.join(srcdir, libdir)):
        sources += [os.path.join(root, filename)
                    for filename in filenames if filename.endswith(('.cpp', '.cc'))]

    # Build library
    bld.shlib(
        features="cxx cxxshlib",
        source=sources,
        target=bld.get_env()["libname"],
        includes=srcdir,
        uselib=bld.get_env()["libs"],
    ) if bld.options.shared else bld.stlib(
        features="cxx cxxstlib",
        source=sources,
        target=bld.get_env()["libname"],
        includes=srcdir,
        uselib=bld.get_env()["libs"],
    )

    # Build executables
    bld.recurse("./src/examples")
    bld.recurse("./src/tests")

    # Install headers
    [bld.install_files("${PREFIX}/include/" + os.path.dirname(f)[4:], f)
     for f in includes]

    # Install libraries
    bld.install_files("${PREFIX}/lib", blddir + "/lib" + bld.get_env()["libname"] + "." + bld.env.SUFFIX) if bld.options.shared else bld.install_files(
        "${PREFIX}/lib", blddir + "/lib" + bld.get_env()["libname"] + ".a")
