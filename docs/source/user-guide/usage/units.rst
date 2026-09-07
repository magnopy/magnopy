.. _user-guide_usage_units:

*****
Units
*****


Magnopy supports a number of units both for input and output data. On this page we
describe what units are supported. All units keywords are case-insensitive.

.. _user-guide_usage_units_energy:

Energy
======

By default Magnopy outputs energy in the units of "meV", but the user can ask for other
supported units.

======================== =======================
Units                    Keywords
======================== =======================
Electronvolt (eV)        ``"eV"``
Millielectronvolt (meV)  ``"meV"``
Joule                    ``"Joule"``, ``"J"``
Rydberg units of energy  ``"Ry"``, ``"Rydberg"``
Erg                      ``"Erg"``
======================== =======================


.. dropdown:: Conversion factors

    .. doctest::

        >>> import magnopy
        >>> eV = magnopy.si.ELECTRON_VOLT
        >>> meV = magnopy.si.MILLI * magnopy.si.ELECTRON_VOLT
        >>> Joule = magnopy.si.JOULE
        >>> Rydberg = magnopy.si.RYDBERG_ENERGY
        >>> Erg = magnopy.si.ERG

    .. doctest::
        :hide:

        >>> from magnopy._constants._units import _ENERGY_UNITS
        >>> eV == _ENERGY_UNITS["ev"]
        True
        >>> meV == _ENERGY_UNITS["mev"]
        True
        >>> Joule == _ENERGY_UNITS["joule"]
        True
        >>> Rydberg == _ENERGY_UNITS["rydberg"]
        True
        >>> Erg == _ENERGY_UNITS["erg"]
        True

    .. doctest::

        >>> print(f"{eV / meV:.1e} meV in 1 eV")
        1.0e+03 meV in 1 eV
        >>> print(f"{meV / eV:.1e} eV in 1 meV")
        1.0e-03 eV in 1 meV

    .. doctest::

        >>> print(f"{eV / Joule:.8e} Joule in 1 eV")
        1.60217663e-19 Joule in 1 eV
        >>> print(f"{Joule / eV:.8e} eV in 1 Joule")
        6.24150907e+18 eV in 1 Joule

    .. doctest::

        >>> print(f"{eV / Rydberg:.8e} Rydberg in 1 eV")
        7.34986444e-02 Rydberg in 1 eV
        >>> print(f"{Rydberg / eV:.7f} eV in 1 Rydberg")
        13.6056931 eV in 1 Rydberg

    .. doctest::

        >>> print(f"{eV / Erg:.8e} Erg in 1 eV")
        1.60217663e-12 Erg in 1 eV
        >>> print(f"{Erg / eV:.8e} eV in 1 Erg")
        6.24150907e+11 eV in 1 Erg

    .. doctest::

        >>> print(f"{meV / Joule:.8e} Joule in 1 meV")
        1.60217663e-22 Joule in 1 meV
        >>> print(f"{Joule / meV:.8e} meV in 1 Joule")
        6.24150907e+21 meV in 1 Joule

    .. doctest::

        >>> print(f"{meV / Rydberg:.8e} Rydberg in 1 meV")
        7.34986444e-05 Rydberg in 1 meV
        >>> print(f"{Rydberg / meV:.4f} meV in 1 Rydberg")
        13605.6931 meV in 1 Rydberg

    .. doctest::

        >>> print(f"{meV / Erg:.8e} Erg in 1 meV")
        1.60217663e-15 Erg in 1 meV
        >>> print(f"{Erg / meV:.8e} meV in 1 Erg")
        6.24150907e+14 meV in 1 Erg


    .. doctest::

        >>> print(f"{Joule / Rydberg:.8e} Rydberg in 1 Joule")
        4.58742456e+17 Rydberg in 1 Joule
        >>> print(f"{Rydberg / Joule:.8e} Joule in 1 Rydberg")
        2.17987236e-18 Joule in 1 Rydberg


    .. doctest::

        >>> print(f"{Joule / Erg:.1e} Erg in 1 Joule")
        1.0e+07 Erg in 1 Joule
        >>> print(f"{Erg / Joule:.1e} Joule in 1 Erg")
        1.0e-07 Joule in 1 Erg


    .. doctest::

        >>> print(f"{Rydberg / Erg:.8e} Erg in 1 Rydberg")
        2.17987236e-11 Erg in 1 Rydberg
        >>> print(f"{Erg / Rydberg:.8e} Rydberg in 1 Erg")
        4.58742456e+10 Rydberg in 1 Erg


.. _user-guide_usage_units_parameters:

Hamiltonian's parameters
========================

Interaction parameters of the :py:class:`.SpinHamiltonian` are typically stored in the
units of energy (i.e. meV or Joule) or some units that offer direct conversion to an energy
scale (like Kelvin, via Boltzmann constant). Use :py:attr:`.SpinHamiltonian.units` to
check or change the units of the parameters.

For the interaction parameters Magnopy supports all
:ref:`energy units <user-guide_usage_units_energy>` and

====== =====================
Units  Keywords
====== =====================
Kelvin ``"K"``, ``"Kelvin"``
====== =====================

.. note::

    Temperature scale and energy scale are connected via Boltzmann constant
    :math:`E = k_B\cdot T`.


.. dropdown:: Conversion factors

    .. doctest::

        >>> Kelvin = magnopy.si.BOLTZMANN_CONSTANT * magnopy.si.KELVIN

    .. doctest::
        :hide:

        >>> from magnopy._constants._units import _PARAMETER_UNITS
        >>> eV == _PARAMETER_UNITS["ev"]
        True
        >>> meV == _PARAMETER_UNITS["mev"]
        True
        >>> Joule == _PARAMETER_UNITS["joule"]
        True
        >>> Rydberg == _PARAMETER_UNITS["rydberg"]
        True
        >>> Erg == _PARAMETER_UNITS["erg"]
        True
        >>> Kelvin == _PARAMETER_UNITS["kelvin"]
        True

    .. doctest::

        >>> print(f"{Kelvin / eV:.8e} eV for 1 Kelvin")
        8.61733326e-05 eV for 1 Kelvin
        >>> print(f"{eV / Kelvin:.4f} Kelvin for 1 eV")
        11604.5181 Kelvin for 1 eV

    .. doctest::

        >>> print(f"{Kelvin / meV:.10f} meV for 1 Kelvin")
        0.0861733326 meV for 1 Kelvin
        >>> print(f"{meV / Kelvin:.7f} Kelvin for 1 meV")
        11.6045181 Kelvin for 1 meV

    .. doctest::

        >>> print(f"{Kelvin / Joule:.8e} Joule for 1 Kelvin")
        1.38064900e-23 Joule for 1 Kelvin
        >>> print(f"{Joule / Kelvin:.8e} Kelvin for 1 Joule")
        7.24297052e+22 Kelvin for 1 Joule

    .. doctest::

        >>> print(f"{Kelvin / Rydberg:.8e} Rydberg for 1 Kelvin")
        6.33362313e-06 Rydberg for 1 Kelvin
        >>> print(f"{Rydberg / Kelvin:.3f} Kelvin for 1 Rydberg")
        157887.512 Kelvin for 1 Rydberg

    .. doctest::

        >>> print(f"{Kelvin / Erg:.8e} Erg for 1 Kelvin")
        1.38064900e-16 Erg for 1 Kelvin
        >>> print(f"{Erg / Kelvin:.8e} Kelvin for 1 Erg")
        7.24297052e+15 Kelvin for 1 Erg

.. _user-guide_usage_units_magnon-energy:

Magnon energies
===============

Magnon quasiparticles are associated with the oscillatory behavior in the classical
picture. Thus, the list of supported units for magnon energies is extended with the
frequency units.

For the magnon energies Magnopy supports all
:ref:`energy units <user-guide_usage_units_energy>` and

========= ==========================
Units     Keywords
========= ==========================
Hertz     ``"Hertz"``, ``"Hz"``
GigaHertz ``"GigaHertz"``, ``"GHz"``
TeraHertz ``"TeraHertz"``, ``"THz"``
========= ==========================

.. note::

    Frequency scale and energy scale are connected via Planck constant
    :math:`E = h\cdot f`.

.. dropdown:: Conversion factors


    .. doctest::

        >>> Hz = magnopy.si.PLANCK_CONSTANT * magnopy.si.HERTZ
        >>> GHz = magnopy.si.GIGA * Hz
        >>> THz = magnopy.si.TERA * Hz

    .. doctest::
        :hide:

        >>> from magnopy._constants._units import _MAGNON_ENERGY_UNITS
        >>> eV == _MAGNON_ENERGY_UNITS["ev"]
        True
        >>> meV == _MAGNON_ENERGY_UNITS["mev"]
        True
        >>> Joule == _MAGNON_ENERGY_UNITS["joule"]
        True
        >>> Rydberg == _MAGNON_ENERGY_UNITS["rydberg"]
        True
        >>> Erg == _MAGNON_ENERGY_UNITS["erg"]
        True
        >>> Hz == _MAGNON_ENERGY_UNITS["hz"]
        True
        >>> GHz == _MAGNON_ENERGY_UNITS["ghz"]
        True
        >>> THz == _MAGNON_ENERGY_UNITS["thz"]
        True

    .. doctest::

        >>> print(f"{GHz / Hz:.1e} Hz in 1 GHz")
        1.0e+09 Hz in 1 GHz
        >>> print(f"{Hz / GHz:.1e} GHz in 1 Hz")
        1.0e-09 GHz in 1 Hz

    .. doctest::

        >>> print(f"{THz / Hz:.1e} Hz in 1 THz")
        1.0e+12 Hz in 1 THz
        >>> print(f"{Hz / THz:.1e} THz in 1 Hz")
        1.0e-12 THz in 1 Hz

    .. doctest::

        >>> print(f"{THz / GHz:.1e} GHz in 1 THz")
        1.0e+03 GHz in 1 THz
        >>> print(f"{GHz / THz:.1e} THz in 1 GHz")
        1.0e-03 THz in 1 GHz

    .. doctest::

        >>> print(f"{Hz / eV:.8e} eV for 1 Hz")
        4.13566770e-15 eV for 1 Hz
        >>> print(f"{eV / Hz:.8e} Hz for 1 eV")
        2.41798924e+14 Hz for 1 eV

    .. doctest::

        >>> print(f"{Hz / meV:.8e} meV for 1 Hz")
        4.13566770e-12 meV for 1 Hz
        >>> print(f"{meV / Hz:.8e} Hz for 1 meV")
        2.41798924e+11 Hz for 1 meV

    .. doctest::

        >>> print(f"{Hz / Joule:.8e} Joule for 1 Hz")
        6.62607015e-34 Joule for 1 Hz
        >>> print(f"{Joule / Hz:.8e} Hz for 1 Joule")
        1.50919018e+33 Hz for 1 Joule

    .. doctest::

        >>> print(f"{Hz / Rydberg:.8e} Rydberg for 1 Hz")
        3.03965969e-16 Rydberg for 1 Hz
        >>> print(f"{Rydberg / Hz:.8e} Hz for 1 Rydberg")
        3.28984196e+15 Hz for 1 Rydberg

    .. doctest::

        >>> print(f"{Hz / Erg:.8e} Erg for 1 Hz")
        6.62607015e-27 Erg for 1 Hz
        >>> print(f"{Erg / Hz:.8e} Hz for 1 Erg")
        1.50919018e+26 Hz for 1 Erg

    .. doctest::

        >>> print(f"{GHz / eV:.8e} eV for 1 GHz")
        4.13566770e-06 eV for 1 GHz
        >>> print(f"{eV / GHz:.8e} GHz for 1 eV")
        2.41798924e+05 GHz for 1 eV

    .. doctest::

        >>> print(f"{GHz / meV:.8e} meV for 1 GHz")
        4.13566770e-03 meV for 1 GHz
        >>> print(f"{meV / GHz:.6f} GHz for 1 meV")
        241.798924 GHz for 1 meV

    .. doctest::

        >>> print(f"{GHz / Joule:.8e} Joule for 1 GHz")
        6.62607015e-25 Joule for 1 GHz
        >>> print(f"{Joule / GHz:.8e} GHz for 1 Joule")
        1.50919018e+24 GHz for 1 Joule

    .. doctest::

        >>> print(f"{GHz / Rydberg:.8e} Rydberg for 1 GHz")
        3.03965969e-07 Rydberg for 1 GHz
        >>> print(f"{Rydberg / GHz:.8e} GHz for 1 Rydberg")
        3.28984196e+06 GHz for 1 Rydberg

    .. doctest::

        >>> print(f"{GHz / Erg:.8e} Erg for 1 GHz")
        6.62607015e-18 Erg for 1 GHz
        >>> print(f"{Erg / GHz:.8e} GHz for 1 Erg")
        1.50919018e+17 GHz for 1 Erg

    .. doctest::

        >>> print(f"{THz / eV:.8e} eV for 1 THz")
        4.13566770e-03 eV for 1 THz
        >>> print(f"{eV / THz:.6f} THz for 1 eV")
        241.798924 THz for 1 eV

    .. doctest::

        >>> print(f"{THz / meV:.8f} meV for 1 THz")
        4.13566770 meV for 1 THz
        >>> print(f"{meV / THz:.9f} THz for 1 meV")
        0.241798924 THz for 1 meV

    .. doctest::

        >>> print(f"{THz / Joule:.8e} Joule for 1 THz")
        6.62607015e-22 Joule for 1 THz
        >>> print(f"{Joule / THz:.8e} THz for 1 Joule")
        1.50919018e+21 THz for 1 Joule

    .. doctest::

        >>> print(f"{THz / Rydberg:.8e} Rydberg for 1 THz")
        3.03965969e-04 Rydberg for 1 THz
        >>> print(f"{Rydberg / THz:.5f} THz for 1 Rydberg")
        3289.84196 THz for 1 Rydberg

    .. doctest::

        >>> print(f"{THz / Erg:.8e} Erg for 1 THz")
        6.62607015e-15 Erg for 1 THz
        >>> print(f"{Erg / THz:.8e} THz for 1 Erg")
        1.50919018e+14 THz for 1 Erg


.. hint::
    To check what units are hardcoded in Magnopy's source code one can do the following

    .. doctest::

        >>> from magnopy._constants._units import _ENERGY_UNITS
        >>> print(list(_ENERGY_UNITS))
        ['ev', 'mev', 'joule', 'j', 'ry', 'rydberg', 'erg']

    .. doctest::

        >>> from magnopy._constants._units import _PARAMETER_UNITS
        >>> print(list(_PARAMETER_UNITS))
        ['ev', 'mev', 'joule', 'j', 'ry', 'rydberg', 'erg', 'k', 'kelvin']

    .. doctest::

        >>> from magnopy._constants._units import _MAGNON_ENERGY_UNITS
        >>> print(list(_MAGNON_ENERGY_UNITS))
        ['ev', 'mev', 'joule', 'j', 'ry', 'rydberg', 'erg', 'hertz', 'hz', 'gigahertz', 'ghz', 'terahertz', 'thz']
