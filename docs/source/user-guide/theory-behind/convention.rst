.. _user-guide_theory-behind_convention:

*****************************
Convention of the Hamiltonian
*****************************

The convention of the Hamiltonian, i.e. the question of how the Hamiltonian is written,
is very important as the values of the interaction parameters of the Hamiltonian can be
different in different conventions.

We describe the general approach to the convention of the Hamiltonian in the supplementary
note 2 of the |paper-2026|_ and in this page briefly summarize how the convention of the
Hamiltonian is reflected in the code.

The general idea in Magnopy is to support as many conventions as possible. Therefore,
the user is responsible for choosing the convention while creating a spin Hamiltonian.
Magnopy does all the grunt work of converting between different conventions when
necessary.

.. hint::
    If the Hamiltonian is read from a known source (such as |TB2J|_, |GROGU|_), the
    convention is known and set automatically.

Before reading about convention we recommend taking a look at how the Hamiltonian is
written: :ref:`user-guide_theory-behind_spin-hamiltonian`.

Constants before the sum
========================

*   Keywords of the :py:class:`magnopy.Convention`: ``c1``, ``c21``, ``c22``, ``c31``,
    ``c32``, ``c33``, ``c41``, ``c42``, ``c43``, ``c44``, ``c45``.

First property that defines the convention of the Hamiltonian is the set of constants
being written before the summation sign (the popular choices for bilinear term are
:math:`\pm 1` or :math:`\pm \frac{1}{2}`). We expect the user to define either all or some
of the eleven constants placed in front of the summation sign for different types of
Hamiltonian's terms

.. math::
    &C_1
    \\
    &C_{2,1} \qquad
    C_{2,2}
    \\
    &C_{3,1} \qquad
    C_{3,2} \qquad
    C_{3,3}
    \\
    &C_{4,1} \qquad
    C_{4,2} \qquad
    C_{4,3} \qquad
    C_{4,4} \qquad
    C_{4,5}


Normalization of spins
======================

*   Keyword of the :py:class:`magnopy.Convention`: ``spin_normalized``.

In Magnopy we treat spins as operators. Therefore, the question of normalization of spins
is ill defined (one can always write :math:`\hat{S} \rightarrow S \frac{\hat{S}}{S}`,
though). Nevertheless, Magnopy supports the values of the interaction parameters under the
assumption that the spin vectors are normalized to unit length in the classical version of
the quantum spin Hamiltonian.

Let us give an example for the bilinear term of the Hamiltonian:

.. math::
    C_2
    \sum_{\substack{\mu, \nu_2, \\ \alpha_1, \alpha_2, \\ i_1, i_2}}
    J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2}
    S_{\mu, \alpha_1}^{i_1}
    S_{\mu + \nu_2, \alpha_2}^{i_2}

*   If by convention the spins are considered to be "not normalized"
    (``spin_normalized = False``), and user inputs parameters
    :math:`J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2}`, then the Hamiltonian would be
    written as

    .. math::
        C_2
        \sum_{\substack{\mu, \nu_2, \\ \alpha_1, \alpha_2, \\ i_1, i_2}}
        J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2}
        S_{\mu, \alpha_1}^{i_1}
        S_{\mu + \nu_2, \alpha_2}^{i_2}

*   If by convention the spins are considered to be "normalized" (``spin_normalized = True``),
    and user inputs parameters
    :math:`J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2}`, then the Hamiltonian would be
    written as

    .. math::
        C_2
        \sum_{\substack{\mu, \nu_2, \\ \alpha_1, \alpha_2, \\ i_1, i_2}}
        J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2}
        \frac{S_{\mu, \alpha_1}^{i_1}}{S_{\alpha_1}}
        \frac{S_{\mu + \nu_2, \alpha_2}^{i_2}}{S_{\alpha_2}}

where :math:`S_{\alpha_1}` and :math:`S_{\alpha_2}` are the spin quantum numbers of the
spin operators :math:`S_{\mu, \alpha_1}^{i_1}` and :math:`S_{\mu + \nu_2, \alpha_2}^{i_2}`
respectively.

As the Hamiltonian should not change with the change of convention (only the values of the
parameters should), the values of the interaction parameters in these two conventions
are related as

.. math::
    \Biggl(J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2} \Biggr)_{\text{not normalized}}
    =
    \dfrac{1}{S_{\alpha_1} S_{\alpha_2}}
    \Biggl(J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2} \Biggr)_{\text{normalized}}

.. _user-guide_theory-behind_convention_multiple-counting:

Multiple counting
=================

*   Keyword of the :py:class:`magnopy.Convention`: ``multiple_counting``.

This property is known as "double counting" in the case of the bilinear terms
(:ref:`ug_tb_sh_2-2`). It specifies, whether both parameters
:math:`J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2}` and
:math:`J^{i_2, i_1}_{-\nu_2; \alpha_2, \alpha_1}` are included in the summation over the
magnetic sites (``multiple_counting = True``) or only one of them
(``multiple_counting = False``). The two interaction parameters from the example above
form what we call an "equivalent set" of parameters, as only their sum is relevant for the
Hamiltonian and not the individual values of each of them (see also
:ref:`user-guide_theory-behind_equivalent-parameters`).

The generalization of this property for the other terms (with 3 and more components of
spin operators) is not as straightforward as one would expect and is linked with the
:ref:`sets of equivalent parameters <user-guide_theory-behind_equivalent-parameters>`. For
terms with more than two spin operators, the equivalent sets can include more than two
parameters.

*   If by convention the multiple counting is allowed (``multiple_counting = True``)

    The Hamiltonian includes all the parameters from the equivalent set and the user is
    expected to manually input each of them.

*   If by convention the multiple counting is not allowed (``multiple_counting = False``)

    Only one of the parameters from the equivalent set is included in the Hamiltonian
    and user is expected to input that single parameter only.


.. hint::
    When ``multiple_counting = True``, one can pass ``populate_equivalent = True`` to the
    :py:meth:`magnopy.SpinHamiltonian.add` method to automatically populate the equivalent
    parameters and avoid the manual input of each of them.
