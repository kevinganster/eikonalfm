#!/usr/bin/env python
import sys


def readme():
    with open("README.md", "r") as f:
        return f.read()


def configuration(parent_package="", top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True)

    # add the subpackage 'eikonalfm' to the config (to call it's setup.py later)
    config.add_subpackage("eikonalfm")

    return config


#from distutils.command.sdist import sdist

pkg_metadata = dict(
    name             = "eikonalfm",
    version          = "0.9.3",
    description      = "solving the (factored) eikonal equation with the Fast Marching method",
    long_description = readme(),
    long_description_content_type = "text/markdown",
    url              = "https://github.com/kevinganster/eikonalfm",
    author           = "Kevin Ganster",
    author_email     = "kevinganster@gmail.com",
    license          = "MIT",
    keywords         = "Fast Marching method, factored Fast Marching method, eikonal equation, factored eikonal equation",
    configuration    = configuration,
    install_requires = ["numpy >= 1.7"],
    # "Development Status :: 4 - Beta",
    classifiers      = [ "License :: OSI Approved :: MIT License",
                         "Operating System :: OS Independent",
                         "Topic :: Scientific/Engineering :: Mathematics",
                         "Programming Language :: C++",
                         "Programming Language :: Python :: 3"],
#    cmdclass={"sdist": sdist}
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
              "  - `pip install eikonal-fm`   (last eikonal-fm release on PyPI)")
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
        from numpy.distutils.core import setup
    # Don't import numpy in other cases - non-build actions are required to succeed without NumPy,
    # for example when pip is used to install this when NumPy is not yet present in the system.

    setup(**pkg_metadata)

