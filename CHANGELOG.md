# Changelog

<!--next-version-placeholder-->

## v0.9.8 (22/05/2025)

- Updated legacy setup.py to the pyproject.toml build system using setuptools.
- Adjusted workflows to work with the new build system and also build arm versions of the wheels.
- With this, Python 3.6 is no longer supported because of the changes in setuptools. Instead, we can now build wheels for Python 3.12 and 3.13.
- Removed some compiler warnings.

## v0.9.7 (27/02/2023)

- Multiple Memory Leak fixes

## v0.9.6 (22/07/2022)

- optional output of sequence and orders by use of an optional parameter `output_sensitivities`.
- redone the setup.py script to use pure setuptools (including a pyproject.toml)

## v0.9.5 (21/07/2020)

- Fixed underflow in calculation of the distance field tau0 (factoredmarcher).
- Changed the `eikonal.distance` function to use an index-vector as source x_s instead of a real-valued vector.

## v0.9.4 (10/06/2020)

- Reverted from `size_t` and `ptrdiff_t` back to `unsigned long` and `long`, since it introduced some multi-platform problems.

## v0.9.0 (18/10/2019)

- First release of `eikonalfm`