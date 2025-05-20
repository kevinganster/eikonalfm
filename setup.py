from setuptools import setup, Extension
import numpy as np


pkg_metadata = {}
pkg_metadata['ext_modules'] = [
    Extension(
        "eikonalfm.cfm",
        sources=[   "src/eikonalfm/cfm.cpp",
                    "src/eikonalfm/marcher.cpp",
                    "src/eikonalfm/factoredmarcher.cpp",
                    "src/eikonalfm/heap.cpp"],
        include_dirs=[  "src/eikonalfm",
                        np.get_include()],
        language="c++",
        extra_compile_args=["-std=c++11"]  # god damn it Mac OS
    )
]
setup(**pkg_metadata)