#!/bin/bash
python setup.py build_ext --build-lib src/

# alternatively:
# pip install -e . --force-reinstall --no-deps