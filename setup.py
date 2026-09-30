import sys

import numpy as np
from setuptools import Extension, setup

ext_args = {}
if sys.platform == "darwin":
    ext_args["extra_compile_args"] = ["-std=c++11"]  # god damn it Mac OS

extension = Extension(
    "eikonalfm.cfm",
    sources=[
        "src/eikonalfm/cfm.cpp",
        "src/eikonalfm/marcher.cpp",
        "src/eikonalfm/factoredmarcher.cpp",
        "src/eikonalfm/heap.cpp",
    ],
    include_dirs=["src/eikonalfm", np.get_include()],
    language="c++",
    **ext_args,
)

setup(ext_modules=[extension], license="MIT")
