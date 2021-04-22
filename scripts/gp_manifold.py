#! /usr/bin/env python
# encoding: utf-8

from waflib.Configure import conf
from utils import check_include, check_lib


def options(opt):
    # Required package options
    opt.load("eigen corrade", tooldir="waf_tools")

    # Options
    opt.add_option(
        "--kernel-path", type="string", help="path to gp-manifold", dest="gpmanifold_path"
    )


@conf
def check_gpmanifold(ctx):
    # Set the search path
    if ctx.options.kernel_path is None:
        path_check = ["/usr/local", "/usr"]
    else:
        path_check = [ctx.options.gpmanifold_path]

    # gp-manifold includes
    check_include(ctx, "GPMANIFOLD", ["gp_manifold"],
                  ["integration.hpp"], path_check)

    # gp-manifold libs
    check_lib(ctx, "GPMANIFOLD", "", ["libGpManifold"], path_check)

    if ctx.env.LIB_GPMANIFOLD or ctx.env.STLIB_GPMANIFOLD:
        # Add dependencies to require libraries
        ctx.get_env()["requires"] = ctx.get_env()[
            "requires"] + ["EIGEN", "CORRADE"]

        # Check for dependencies
        ctx.load("eigen corrade", tooldir="waf_tools")

        # Add library
        ctx.get_env()["libs"] += ["GPMANIFOLD"]


def configure(cfg):
    if not cfg.env.LIB_GPMANIFOLD and not cfg.env.STLIB_GPMANIFOLD:
        cfg.check_gpmanifold()
