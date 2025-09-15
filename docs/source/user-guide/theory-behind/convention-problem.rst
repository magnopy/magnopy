.. _user-guide_theory-behind_convention-problem:

******************
Convention problem
******************


An illustration
===============

To illustrate the problem of the convention we consider an simple example of a
nearest-neighbor model on a cubic lattice with one spin per unit cell

.. math::
    \mathcal{H}
    =
    \sum_{\mu,\nu}
    J_{\nu}
    (\boldsymbol{S}_{\mu}
    \cdot\boldsymbol{S}_{\mu+\nu})

where :math:`J_{\nu} = 1` meV if two spins :math:`\boldsymbol{S}_{\mu}` and
:math:`\boldsymbol{S}_{\mu+\nu}` are nearest-neighbors and zero otherwise.
Moreover, we assume that the value of spin is :math:`|\boldsymbol{S}_{\mu}| = 3/2` and
all spins are oriented along positive :math:`\hat{z}` direction. What is the energy of
the system in this configuration? The answer to this question is, surprisingly,
impossible to give, as we do not know two important facts:

*   Does the sum count both terms :math:`(\mu;\mu+\nu)` and :math:`(\mu+\nu;\mu)`?
    In other words: does the sum include double-counting?

*   Are spin vectors normalized to unity? I.e. does :math:`|\boldsymbol{S}_{\mu}| = 1`
    in the sum and the spin value is accounted for in the exchange parameters or not?

Each of those questions has two answers and together we have four different possible
answers for the energy per unit cell:


(1) :math:`E = 6` meV if bonds are counted twice and spins are normalized;
(2) :math:`E = 13.5` meV if bonds are counted twice and spins are not normalized;
(3) :math:`E = 3` meV if bonds are counted once and spins are normalized;
(4) :math:`E = 6.75` meV if bonds are counted once and spins are not normalized.

If one add the freedom of writing the constant before the sum, then the amount of
different answers will be even larger. Often in the literature the same Hamiltonian is
written with :math:`\pm 1` or :math:`\pm 1/2` in front of the sum.

Therefore, set of parameters, without the *knowledge of the convention* of spin
Hamiltonian is not enough to define it. At the same time is does not matter in what
convention the Hamiltonian is written if the parameters are *correctly converted to it*.
In other words, energy of the system does not depend just on the convention or just
on the parameters. It depend on the both **parameters and convention**.

If the Hamiltonian above is given for the convention (1) with :math:`J_{\nu} = 1` meV. Then
for the remaining three conventions the parameters should be

(2) :math:`J_{\nu} = 6/13.5` meV, then :math:`E = 6` meV;
(3) :math:`J_{\nu} = 2` meV, then :math:`E = 6` meV;
(4) :math:`J_{\nu} = 6/6.75` meV, then :math:`E = 6` meV.

If the pair of the parameters and convention is specified correctly, then the energy of
the system does not depend on it.

Magnopy's approach
==================

Magnopy does not give preference to any convention, supports each one of them and
implements conversion between any pair of conventions. Magnopy was written to
automatically handle everything that concerns convention of the Hamiltonian. User only
have to specify the convention upon the creation of spin Hamiltonian. If the Hamiltonian
is read from the know source (i.e. |TB2J|_, |GROGU|_), then the convention of this source is
known and set automatically.

The Convention in magnopy is implemented in the form of the small :py:class:`.Convention`
class and is fully defined by the following

*   ``c1``, ``c21``, ``c22``, ``c31``, ``c32``, ``c33``, ``c41``, ``c421``, ``c422``,
    ``c43``, ``c44``

    Eleven constants of the expanded form of the spin Hamiltonian

    .. math::
        C_1 \qquad
        C_{2,1} \qquad
        C_{2,2} \qquad
        C_{3,1} \qquad
        C_{3,2} \qquad
        C_{3,3} \\
        C_{4,1} \qquad
        C_{4,2,1} \qquad
        C_{4,2,2} \qquad
        C_{4,3} \qquad
        C_{4,4}

    Each constant is an individual property of the convention and affect the corresponding
    term of the expanded form of the Hamiltonian.

*   ``spin_normalized``

    Whether spins are normalized. This property affects all terms of the expanded form.
*   ``multiple_counting``

    Whether bonds are counted multiple times. This property does not affect first,
    second, fourth and seventh terms of the expanded form.

    This property is responsible for the double counting from the example. However, for
    the terms with more than two spins the equivalent bonds are counter not twice,
    but more times. That is why we do not name this property as "double counting".


Some of the properties above can be undefined if no corresponding parameter are present.
If any of the parameters is affected by some property, then it have to be defined.
