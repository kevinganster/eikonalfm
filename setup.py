#!/usr/bin/env python
import sys
from setuptools import setup, Extension
import re


def readme():
    with open("README.md", "r") as f:
        return f.read()

# get metadata from module
with open("eikonalfm/__init__.py") as f:
    metadata = dict(re.findall(r'__(.*)__ ?= ?"(.+)"', f.read()))


pkg_metadata = dict(
    name             = metadata['title'],
    version          = metadata['version'],
    author           = metadata['author'],
    author_email     = metadata['email'],
    description      = "solving the (factored) eikonal equation with the Fast Marching method",
    long_description = readme(),
    long_description_content_type = "text/markdown",
    url              = "https://github.com/kevinganster/eikonalfm",
    license          = "MIT",
    keywords         = "Fast Marching method, factored Fast Marching method, eikonal equation, factored eikonal equation",
    python_requires = ">= 3.0",
    install_requires = [ "setuptools < 60.0",
                         "numpy >= 1.7"],
    # "Development Status :: 4 - Beta",
    classifiers      = [ "License :: OSI Approved :: MIT License",
                         "Operating System :: OS Independent",
                         "Topic :: Scientific/Engineering :: Mathematics",
                         "Programming Language :: C++",
                         "Programming Language :: Python :: 3"],
    packages = ["eikonalfm"]
)


def parse_setuppy_commands():
    """Check the commands and respond appropriately.
    Return a boolean value for whether or not to run the build or not.
    
    Inspired from SciPy's setup.py : https://github.com/scipy/scipy/blob/master/setup.py
    """
    args = sys.argv[1:]

    if not args:
        # User forgot to give an argument probably, let setuptools handle that.
        return True

    info_commands = ['--help-commands', '--name', '--version', '-V',
                     '--fullname', '--author', '--author-email',
                     '--maintainer', '--maintainer-email', '--url',
                     '--license', '--description', '--long-description',
                     '--platforms', '--classifiers', '--keywords',
                     '--provides', '--requires', '--obsoletes',
                     'egg_info', 'install_egg_info', 'rotate']

    for command in info_commands:
        if command in args:
            return False

    good_commands = ('develop', 'sdist', 'build', 'build_ext', 'build_py',
                     'build_clib', 'build_scripts', 'bdist_wheel', 'bdist_rpm',
                     'bdist_wininst', 'bdist_msi', 'bdist_mpkg',
                     'build_sphinx')

    for command in good_commands:
        if command in args:
            return True

    # The following commands are supported, but we need to show more useful messages to the user
    if "install" in args:
        print("Note: if you need reliable uninstall behavior, you should install with pip instead of using `setup.py install`:"
              "  - `pip install .`       (from a git repo or downloaded source release)"
              "  - `pip install eikonalfm`   (last eikonalfm release on PyPI)")
        return True

    return False


if __name__ == "__main__":
    if "--force" in sys.argv:
        run_build = True
        sys.argv.remove("--force")
    else:
        # Raise errors for unsupported commands, improve help output, etc.
        run_build = parse_setuppy_commands()

    # This import is here because it needs to be done before importing setup() from numpy.distutils
    from setuptools import setup

    if run_build:
        # non-build actions should not include the extension module,
        # for example when pip is used to install this when NumPy is not yet present in the system.
        # for wheels numpy has to be already installed
        import numpy as np

        pkg_metadata['ext_modules'] = [
            Extension(
                "eikonalfm.cfm",
                sources=[   "eikonalfm/cfm.cpp",
                            "eikonalfm/marcher.cpp",
                            "eikonalfm/factoredmarcher.cpp",
                            "eikonalfm/heap.cpp"],
                include_dirs=[  "eikonalfm",
                                np.get_include()],
                language="c++",
                extra_compile_args=["-std=c++11"]  # god damn it Mac OS
            )
        ]

    setup(**pkg_metadata)

