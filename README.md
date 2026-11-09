# Magnopy

Spin Hamiltonian, ground state, magnons.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/magnopy.svg)](https://badge.fury.io/py/magnopy/)
![Python](https://img.shields.io/pypi/pyversions/magnopy)

[![Documentation Status](https://readthedocs.org/projects/magnopy/badge/?version=latest)](https://magnopy.org/en/latest/?badge=latest)
[![tests (main)](https://img.shields.io/github/actions/workflow/status/magnopy/magnopy/singular-test.yml?branch=main&label=tests%20(main))](https://github.com/magnopy/magnopy/actions/workflows/singular-test.yml?query=branch%3Amain)
[![tests (dev)](https://img.shields.io/github/actions/workflow/status/magnopy/magnopy/singular-test.yml?branch=dev&label=tests%20(dev))](https://github.com/magnopy/magnopy/actions/workflows/singular-test.yml?query=branch%3Adev)


## What is Magnopy?

Magnopy is a Python code that, given a
[spin Hamiltonian](https://docs.magnopy.org/en/latest/user-guide/theory-behind/spin-hamiltonian.html)
in **any**
[convention](https://docs.magnopy.org/en/latest/user-guide/theory-behind/convention.html),
computes bosonic (magnon) Hamiltonian.

* **Any convention** - You declare how your parameters are defined once; Magnopy
  converts when necessary.
* **Read your files** - interfaced with
  [TB2J](https://github.com/mailhexu/TB2J) and
  [GROGU](https://grogupy.readthedocs.io/en/stable/).
* **Interactions up to four spins** - Bilinear exchange, single-ion anisotropy,
  three- and four-spin terms, Zeeman, dipole-dipole. All with full interaction tensors.
* **Ground state and magnons** - Minimization of classical energy and Linear Spin Wave Theory.
* **Two interfaces** - A Python API for scripting. Command-line tools that take you from a file to a plot without writing code.

## Quick example

```python
import magnopy

# Load exchange parameters (TB2J, GROGU, or build your own)
spinham = magnopy.examples.cubic_ferro_nn(J_iso=1.0, S=2.5)

# Find the classical ground state
energy = magnopy.Energy(spinham)
spin_directions = energy.optimize()

# Magnon dispersion from Linear Spin Wave Theory
lswt = magnopy.LSWT(spinham, spin_directions=spin_directions)

for label, k in [("Γ", [0, 0, 0]), ("X", [0.5, 0, 0]), ("R", [0.5, 0.5, 0.5])]:
    print(f"omega({label}) = {lswt.omega(k)[0].real:.4f} meV")
```

```console
omega(Γ) = 0.0000 meV
omega(X) = 0.6121 meV
omega(R) = 1.8363 meV
```

Reading from a file instead:

```python
spinham = magnopy.io.load_tb2j("exchange.out", spin_values=[2.5, 2.5])
```

Or from the command line, no Python scripting required:

```console
magnopy-lswt -ss TB2J -sf exchange.out -sv 2.5 2.5
```

## Documentation

Extensive documentation is available at [magnopy.org](https://magnopy.org).

Magnopy can be used both as a Python library and as a command line tool.

* For details about black box style usage see
  [magnopy-optimize-sd](https://docs.magnopy.org/en/latest/user-guide/cli/magnopy-optimize-sd/index.html)
  or
  [magnopy-lswt](https://docs.magnopy.org/en/latest/user-guide/cli/magnopy-lswt/index.html).
* For code examples see
  [User guide](https://docs.magnopy.org/en/latest/user-guide/index.html).
* For full public API see
  [API](https://docs.magnopy.org/en/latest/api/index.html).
* To get some support and ask questions see
  [User support](https://docs.magnopy.org/en/latest/user-support/index.html).
* For summary of releases see
  [Release notes](https://docs.magnopy.org/en/latest/release-notes/index.html).

## Installation

* Full installation guide is in
  [docs](https://docs.magnopy.org/en/latest/user-guide/installation.html).

To install Magnopy run (you may need to use `pip3`):

```console
pip install magnopy
```

To install with visualization capabilities (recommended) run (you may need to use
`pip3`):

```console
pip install "magnopy[visual]"
```

## User support

If you have a question about Magnopy, have an idea for its improvement, or found a bug,
do not hesitate to contact us via one of the
[support channels](https://docs.magnopy.org/en/latest/user-support/index.html).

## Development

We welcome contributions to Magnopy. Visit the
[development guide](https://docs.magnopy.org/en/latest/development/index.html)
or contact authors of the project to find out more.


## License

The source code of Magnopy is licensed under the  GNU General Public
License (GPL-3.0). See the ["LICENSE" file](https://github.com/magnopy/magnopy/blob/main/LICENSE) in the [Magnopy's repository](https://github.com/magnopy/magnopy).

In addition, if you use magnopy in the scientific publication, cite the preprint
available at `arxiv.org/abs/2611.06833 <https://arxiv.org/abs/2611.06833>`_.

```text
Rybakov, A.  and Marino, M. and Blanter, Y. M. and Coronado, E. and Ferrer, J, 2026.
Quantum Spin Wave Theory with Multi-Spin Interactions for Magnetic Materials Modeling.
Preprint available at https://arxiv.org/abs/2611.06833
```

```latex
@misc{magnopy-preprint,
  doi = {10.48550/ARXIV.2611.06833},
  url = {https://arxiv.org/abs/2611.06833},
  author = {Rybakov,  Andrey and Marino,  Marco and Blanter,  Yaroslav M. and Coronado,  Eugenio and Ferrer,  Jaime},
  keywords = {Other Condensed Matter (cond-mat.other),  Computational Physics (physics.comp-ph),  FOS: Physical sciences},
  title = {Quantum Spin Wave Theory with Multi-Spin Interactions for Magnetic Materials Modeling},
  publisher = {arXiv},
  year = {2026},
}
```

For the detailed guide on how to cite the papers on which Magnopy depends
see [Citation guide](https://docs.magnopy.org/en/latest/cite.html).
