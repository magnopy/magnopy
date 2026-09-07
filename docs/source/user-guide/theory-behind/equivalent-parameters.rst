.. _user-guide_theory-behind_equivalent-parameters:

*********************
Equivalent parameters
*********************

An example
==========

Let one consider the bilinear Hamiltonian that contains (among other) the following term

.. math::

    \mathcal{H}
    =
    ...
    +
    J_{\nu_2; \alpha_1, \alpha_2}^{x,y} S_{\mu, \alpha_1}^x S_{\mu+\nu_2, \alpha_2}^y
    +
    J_{-\nu_2; \alpha_2, \alpha_1}^{y,x} S_{\mu+\nu_2, \alpha_2}^y S_{\mu, \alpha_1}^x
    +
    ...

Now, in the case when
:math:`\boldsymbol{r}_{\mu, \alpha_1} \ne \boldsymbol{r}_{\mu+\nu_2, \alpha_2}`, the
spin operators :math:`S_{\mu, \alpha_1}^x` and :math:`S_{\mu+\nu_2, \alpha_2}^y`
commute, thus the Hamiltonian can be rewritten as

.. math::

    \mathcal{H}
    =
    ...
    +
    \left(
        J_{\nu_2; \alpha_1, \alpha_2}^{x,y}
        +
        J_{-\nu_2; \alpha_2, \alpha_1}^{y,x}
    \right)
    S_{\mu, \alpha_1}^x S_{\mu+\nu_2, \alpha_2}^y
    +
    ...

You can see that the physics of the Hamiltonian is determined only by the sum
:math:`J_{\nu_2; \alpha_1, \alpha_2}^{x,y} + J_{-\nu_2; \alpha_2, \alpha_1}^{y,x}`.
Therefore, that sum can be distributed between the parameters arbitrarily. We call
those two components of the parameters **equivalent**. More generally, a set of
parameters is called a **set of equivalent parameters** (or an **equivalent set**)
if the physics of the Hamiltonian is determined solely by the sum of the (components of
the) parameters from that set.

In the case of the bilinear term, every equivalent set contains two parameters:

.. math::

    J_{\nu_2; \alpha_1, \alpha_2}^{i_1, i_2}
    \quad \text{and} \quad
    J_{-\nu_2; \alpha_2, \alpha_1}^{i_2, i_1}

In the case of the other types of terms, the equivalent sets can contain more than two
parameters. All equivalent sets are listed in the supplementary notes of the
|paper-2026|_ (equations (S.21)—(S.27)).

Symmetrization
==============

As one is free to distribute the sum of the parameters from the equivalent set arbitrarily,
one can choose the parameters from the equivalent set to be equal to each other. For
example, for the bilinear term, one can choose

.. math::

    J_{\nu_2; \alpha_1, \alpha_2}^{i_1, i_2}
    =
    J_{-\nu_2; \alpha_2, \alpha_1}^{i_2, i_1}

We call such choice a case of **symmetrized parameters**.

The user is free to input non-symmetrized parameters to the
:py:class:`magnopy.SpinHamiltonian` class. However, Magnopy will symmetrize them when
it is necessary.

.. note::
    The symmetrization of the parameters **does not** place any restrictions on the
    components of the vector/matrix/tensor of each individual parameter. It only relates
    components of the vector/matrix/tensor of **different** parameters. For instance, it
    **does not** forbid an antisymmetric Dzyaloshinskii-Moriya interaction.


Connection with ``multiple_counting``
=====================================

The choice of the distribution of the sum within the equivalent set is only relevant when
more than one equivalent interaction is included in the summation over the magnetic sites.
Therefore, the choice of the distribution of the sum within the equivalent set is only
relevant when ``multiple_counting = True``.
