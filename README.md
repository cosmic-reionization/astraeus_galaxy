# astraeus_galaxy

A lightweight C++ implementation of selected **Astraeus** galaxy evolution routines with merger trees for stellar mass, radiative feedback and metals production exposed through a Python extension module.

This project is intended as a small, self-contained version of the galaxy-evolution part of Astraeus, with the goal of making the relevant routines callable from Python and eventually integrating them with radiative-transfer calculations.

## Features

The current implementation includes:

* Galaxy evolution for a single galaxy over one timestep.
* Continuous halo-mass growth correction.
* Gas accretion from halo growth.
* Metal and dust accretion from the IGM.
* Gas and metal contributions from progenitors.
* Radiative-feedback suppression of the gas reservoir.
* Several radiative-feedback prescriptions:

  * Jeans-mass based feedback
  * Cooling-mass based feedback
  * Sobacchi (2013)
  * Gnedin (2000)
* Python bindings using the Python C API and NumPy C API.

## Project structure

```text
.
├── Makefile
├── python_module.cpp
├── evol_gal.cpp
├── evol_gal.h
├── radiative_feedback.cpp
├── radiative_feedback.h
└── test_astraeus.py
```

### `evol_gal.cpp`

Contains the main galaxy-evolution routines:

```cpp
evolve_gal()
correct_Mvir_to_continously_rise()
get_accreted_and_merged_metals()
get_accreted_and_merged_gas()
apply_radFeedback_and_update_zreion()
```

### `radiative_feedback.cpp`

Contains the radiative-feedback prescriptions and supporting functions, including:

```cpp
calc_radiative_feedback()
get_Mcool()
get_Mjeans()
get_Mjeans_vir()
fg_radfb_Mjeans()
fg_radfb_Mcool()
fg_radfb_Sobacchi()
fg_radfb_Gnedin()
fg_radfb_TempEvol()
```

### `python_module.cpp`

Provides the Python interface.

The resulting module can be imported as:

```python
import astraeus
```

and the main evolution routine is exposed as:

```python
astraeus.evolve_gal(...)
```

## Requirements

The code currently requires:

* C++17 compatible compiler
* `g++`
* Python 3
* NumPy
* Python development headers
* NumPy C headers

For example, on a Linux system:

```bash
python3 --version
g++ --version
python3 -c "import numpy; print(numpy.__version__)"
```

## Compilation

The extension is built using the supplied `Makefile`.

```bash
make
```

This produces a Python extension with a Python-specific suffix, for example:

```text
astraeus.cpython-312-x86_64-linux-gnu.so
```

Clean the build with:

```bash
make clean
```

The extension must be located in a directory visible to Python, such as the current working directory or a directory included in `PYTHONPATH`.

## Python interface

The module is imported with:

```python
import astraeus
```

The main function is:

```python
astraeus.evolve_gal(...)
```

The function accepts scalar galaxy/halo properties together with NumPy arrays containing the IGM metallicity and metal-merger information.

The C++ function modifies several quantities in place.

### Important note about scalar values

Scalar Python values such as:

```python
Mvir = 1.0e10
```

are copied into local C++ variables by `PyArg_ParseTuple()`.

Therefore, changing:

```cpp
float &Mvir
```

inside `evolve_gal()` does **not** modify the original Python variable.

For example:

```python
Mvir = 1.0e10
astraeus.evolve_gal(..., Mvir, ...)
print(Mvir)
```

will still contain the original Python value unless the C++ wrapper explicitly returns the modified value.

NumPy arrays are different: the wrapper obtains their underlying data through `PyArray_DATA()`, so modifications to those arrays can be reflected directly in Python.

## C++ / Python linkage

The Python module is named:

```text
astraeus
```

The corresponding module definition is:

```cpp
"astraeus"
```

and the Python initialization function is:

```cpp
PyInit_astraeus
```

These names must remain consistent with the compiled extension:

```text
astraeus.<Python-extension-suffix>.so
```

For example:

```text
astraeus.cpython-312-x86_64-linux-gnu.so
```

## Development status

This is currently a **development / experimental implementation** rather than a complete replacement for the original Astraeus galaxy-evolution code.

Some parts of the original Astraeus implementation have been removed or simplified, and several routines still contain `TODO` items concerning the exact physical interpretation and integration with the original Astraeus model.

In particular, the current implementation is being adapted for use together with **pyC2Ray**, with radiative feedback evaluated using quantities obtained from the radiative-transfer calculation.

## Physical model

The galaxy evolution currently follows the basic sequence:

```text
Halo properties
      │
      ▼
Continuous Mvir correction
      │
      ▼
Gas accretion
      │
      ├── IGM accretion
      └── progenitor contribution
      │
      ▼
Metal and dust accretion
      │
      ▼
Radiative feedback
      │
      ▼
Updated galaxy gas reservoir
```

The radiative-feedback calculation determines a maximum allowed gas mass through a feedback factor:

```text
Mgas,max = fg × (Ωb / Ωm) × Mvir
```

where `fg` is calculated by the selected radiative-feedback prescription.

## License

Add the appropriate license here according to the licensing terms of the original Astraeus code and the modifications contained in this repository.

## Acknowledgements

This code is based on galaxy-evolution routines originating from **Astraeus** and has been adapted to provide a lightweight C++/Python interface for use in numerical astrophysics and radiative-transfer workflows.

