# we can safely import numpy here, sice it should be installed by now
import numpy as np
from setuptools import Extension, setup

ext = Extension(
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


if __name__ == '__main__':
    print("calling setup.py (eikonalfm)")

    setup(
        name = "eikonalfm",
        packages = ['eikonalfm'],
        script_args=['build_ext'],
        options={'build_ext': {'inplace': True}},
        ext_modules=[ext]
    )
