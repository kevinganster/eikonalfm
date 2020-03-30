# we can safely import numpy here, sice it should be installed by now
# from setuptools import setup
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    config = Configuration("eikonalfm", parent_package, top_path)
    config.add_extension("cfm",
                          sources=["cfm.cpp",
                                   "marcher.cpp",
                                   "factoredmarcher.cpp",
                                   "heap.cpp"],
                          include_dirs=['.'],
                          language="c++",
                          extra_compile_args=["-std=c++11"]) # god damn it Mac OS
    return config


if __name__ == '__main__':
    print("calling setup.py (eikonalfm)")
    setup(**configuration(top_path="").todict(),
          script_args=['build_ext'],
          options={'build_ext': {'inplace': True}}
    )
