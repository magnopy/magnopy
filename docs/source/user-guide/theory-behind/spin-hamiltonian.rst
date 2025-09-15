.. _user-guide_theory-behind_spin-hamiltonian:

****************
Spin Hamiltonian
****************

.. _user-guide_theory-behind_spin-hamiltonian_compact-form:

Compact form
============

Magnopy works with the spin Hamiltonian of the general  **compact** form that includes
at most four linearly coupled spins

.. include:: ../../core-formulas/spin-hamiltonian-compact.inc

.. dropdown:: Meaning of indices

    *   Superscript Latin letters :math:`i, j, u, v` denotes three cartesian components.
        Therefore, the scalar product of two vectors
        :math:`\boldsymbol{S}_{\mu, \alpha}` and :math:`\boldsymbol{J}_{\alpha}`
        is written as

        .. math::
            (\boldsymbol{S}_{\mu, \alpha} \cdot \boldsymbol{J}_{\alpha})
            =
            \sum_k
            S^i_{\mu, \alpha} \cdot J^i_{\alpha}

    *   Subscript Greek letters denote the position of the unit cell in the periodic
        crystal (:math:`\mu, \nu, \lambda, \rho`) or position of the spin within some
        unit cell (:math:`\alpha, \beta, \gamma, \varepsilon`). We understand that the
        indices :math:`\mu` and :math:`\nu` (for an 3-dimensional lattice) can be
        represented by the vector of integer scalars

        .. math::
            \mu &= (\mu_1, \mu_2,, \mu_3)
            \\
            \nu &= (\nu_1, \nu_2,, \nu_3)

        Then the position of the unit cell that contains the spin
        :math:`\boldsymbol{S}_{\mu,\alpha}` is expressed as

        .. math::

            \boldsymbol{r}_{\mu}
            =
            \sum_{n=1}^3 \mu_n\boldsymbol{a}_n

        and the position of the second unit cell that contains the spin
        :math:`\boldsymbol{S}_{\mu+\nu,\beta}` is

        .. math::
            \boldsymbol{r}_{\mu + \nu}
            =
            \sum_{i=n}^3 (\mu_n + \nu_n)\boldsymbol{a}_n
            =
            \boldsymbol{r}_{\mu} + \boldsymbol{r}_{\nu}.

        where :math:`\boldsymbol{a}_n` are the three lattice vectors of the Bravais
        lattice of the underlying crystal. In other word the summation and addition of
        this indices directly translates into the summation and addition of the
        corresponding position vectors. Note, hat the indices :math:`\alpha` and
        :math:`\beta` can be represented by the vectors of three real numbers
        :math:`\in [0, 1]`.

        We hope that three images below will help to understand the Greek indices better
        (click to enlarge)

        .. figure:: ../../../images/positions_cells.png
            :target: ../../_images/positions_cells.png
            :align: center

        .. figure:: ../../../images/positions_sites.png
            :target: ../../_images/positions_sites.png
            :align: center

The parameters :math:`\boldsymbol{J}(\boldsymbol{r}_{\alpha})`,
:math:`\boldsymbol{J}(\boldsymbol{r}_{\nu,\alpha\beta})`,
:math:`\boldsymbol{J}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})`
and
:math:`\boldsymbol{J}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})`
can describe either some external effect or internal interactions between spins. The
Hamiltonian above is discussed and solved in the (TODO) research paper about magnopy.


.. _user-guide_theory-behind_spin-hamiltonian_expanded-form:

Expanded form
=============

However, internally and as a part of public API magnopy stores the Hamiltonian in the
**expanded** form

.. include:: ../../core-formulas/spin-hamiltonian-expanded.inc

Expanded form look monstrous, but it is more convenient as the terms with the same
amount of spins, but different number of unique sites can represent distinct physical
properties or interactions. Moreover, often the numerical constants before the sum for
those terms are different. For example, the on-site quadratic anisotropy is often
written with :math:`C_{2,1} = \pm1` and bilinear exchange interaction with
:math:`C_{2,2} = \pm1/2`. Both of them describe interaction with two spins, however, it
is impossible to write them in the same form with one constant :math:`C_2` without
modification of the parameters themselves.

Magnopy can read this spin Hamiltonian from |TB2J|_ or |GROGU|_ automatically,
taking into account the convention of the Hamiltonian of each source. Real constants
:math:`C_1`, :math:`C_{2,1}`, :math:`C_{2,2}`, :math:`C_{3,1}`, :math:`C_{3,2}`,
:math:`C_{3,3}`, :math:`C_{4,1}`, :math:`C_{4,2,1}`, :math:`C_{4,2,2}`, :math:`C_{4,3}`,
:math:`C_{4,4}` allow magnopy to support **any** convention of the spin Hamiltonian. To
read more about what defines the convention of the spin Hamiltonian go
:ref:`here <user-guide_theory-behind_convention-problem>`.

In the :py:class:`.SpinHamiltonian` class, that is used to store the parameters of the
Hamiltonian, each term is referenced by the numerical indices of the constants. For
example, to have access to the third term (two spins & two sites) one can use
:py:attr:`.SpinHamiltonian.p22` property. To add or remove parameters of this term one
should use :py:meth:`.SpinHamiltonian.add_22` and
:py:meth:`.SpinHamiltonian.remove_22`.

Examples
========

In this part we show how the common terms of the spin Hamiltonian can be written in the
expanded form.

Zeeman interaction
------------------

Linear coupling with the magnetic field is usually written as

.. math::
    \mathcal{H}
    =
    \mu_B\boldsymbol{h}\sum_{\mu,\alpha} g_a \boldsymbol{S}_{\mu, \alpha}

This term can be written as the first term of the expanded form if one defines
:math:`C_1 = 1` and
:math:`\boldsymbol{J}_1(\boldsymbol{r}_{\alpha}) = \mu_B\boldsymbol{h} g_a`.

On-site anisotropy
------------------

We take an example of the triaxial anisotropy, that can be written as

.. math::
    \mathcal{H}
    =
    \sum_{\mu,\alpha}
    \Bigl(
        K^x (S^x_{\mu,\alpha})^2
        +
        K^y (S^y_{\mu,\alpha})^2
        +
        K^z (S^z_{\mu,\alpha})^2
    \Bigr)

This Hamiltonian can be written as the second term of the expanded form if one defines
:math:`C_{2,1} = 1` and

.. math::
    \boldsymbol{J}_{2,1}(\boldsymbol{r}_{\alpha})
    =
    \begin{pmatrix}
        K^x & 0 & 0 \\
        0 & K^y & 0 \\
        0 & 0 & K^z
    \end{pmatrix}

Exchange interaction
--------------------

Bilinear exchange interaction with isotropic and Dzyaloshinsky-Moriya (DM) exchange can
be written as

.. math::
    \mathcal{H}
    =
    \dfrac{1}{2}
    \sum_{\mu,\alpha}
    \sum_{\mu + \nu,\beta}
    \Bigl[
    J^{iso}(\boldsymbol{r}_{\nu,\alpha\beta})
    (\boldsymbol{S}_{\mu,\alpha}
    \cdot
    \boldsymbol{S}_{\mu+\nu,\beta})
    +
    \boldsymbol{D}(\boldsymbol{r}_{\nu,\alpha\beta})
    \cdot
    (\boldsymbol{S}_{\mu,\alpha}
    \times
    \boldsymbol{S}_{\mu+\nu,\beta})
    \Bigr]

where :math:`\boldsymbol{D}(\boldsymbol{r}_{\nu,\alpha\beta})` is a DM vector. This
Hamiltonian can be written as the third term of the expanded form if one defines
:math:`C_{2,2} = 1/2` and

.. math::
    J_{2,2}(\boldsymbol{r}_{\nu,\alpha\beta})
    =
    \begin{pmatrix}
        J^{iso}(\boldsymbol{r}_{\nu,\alpha\beta}) & D^z(\boldsymbol{r}_{\nu,\alpha\beta}) & -D^y(\boldsymbol{r}_{\nu,\alpha\beta}) \\
        -D^z(\boldsymbol{r}_{\nu,\alpha\beta}) & J^{iso}(\boldsymbol{r}_{\nu,\alpha\beta}) & D^x(\boldsymbol{r}_{\nu,\alpha\beta}) \\
        D^y(\boldsymbol{r}_{\nu,\alpha\beta}) & -D^x(\boldsymbol{r}_{\nu,\alpha\beta}) & J^{iso}(\boldsymbol{r}_{\nu,\alpha\beta})
    \end{pmatrix}

Biquadratic exchange
--------------------

Isotropic biquadratic exchange interaction can be written as

.. math::
    \mathcal{H}
    =
    \sum_{\mu,\alpha}
    \sum_{\mu + \nu,\beta}
    J(\boldsymbol{r}_{\nu,\alpha\beta})
    (\boldsymbol{S}_{\mu,\alpha}
    \cdot
    \boldsymbol{S}_{\mu+\nu,\beta})^2

This Hamiltonian can be written as the ninth term of the expanded form is one defines
:math:`C_{4,2,2} = 1` and tensor
:math:`\boldsymbol{J}_{4,2,2}(\boldsymbol{r}_{\nu,\alpha\beta})` such as
:math:`J_{4,2,2}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}) = J(\boldsymbol{r}_{\nu,\alpha\beta})`
if  :math:`(ijuv)` is one of

.. math::
    \begin{matrix}
    (xxxx), & (xyxy), & (xzxz), \\
    (yxyx), & (yyyy), & (yzyz), \\
    (zxzx), & (zyzy), & (yyyy)
    \end{matrix}

and :math:`J_{4,2,2}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}) = 0` otherwise.
