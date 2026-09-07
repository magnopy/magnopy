.. _user-guide_theory-behind_spin-hamiltonian:

****************
Spin Hamiltonian
****************

The basic theory behind Magnopy is described in |paper-2026|_. This page is intended to
serve as a bridge between the paper and the code, by introducing the details of the
formalism and explaining the concepts that are less evident in the paper, but important
for the code.

In |paper-2026|_ we discuss the spin Hamiltonian with arbitrary number :math:`n` of
coupled components of spin angular momentum operator. However, at present time only the
terms with at most four coupled components of spin angular momentum operator (i.e.
:math:`n \le 4`) are implemented in Magnopy. The Hamiltonian then can be written as

.. math::
    \mathcal{H}_{n \le 4}
    =&
    \,C_1
    \sum_{\substack{\mu, \\ \alpha_1, \\ i_1}}
    J^{i_1}_{\alpha_1}
    S_{\mu, \alpha_1}^{i_1}
    +\\&+
    C_2
    \sum_{\substack{\mu, \nu_2, \\ \alpha_1, \alpha_2, \\ i_1, i_2}}
    J^{i_1, i_2}_{\nu_2; \alpha_1, \alpha_2}
    S_{\mu, \alpha_1}^{i_1}
    S_{\mu + \nu_2, \alpha_2}^{i_2}
    +\\&+
    C_3
    \sum_{\substack{\mu, \nu_2, \nu_3, \\ \alpha_1, \alpha_2, \alpha_3, \\ i_1, i_2, i_3}}
    J^{i_1, i_2, i_3}_{\nu_2, \nu_3; \alpha_1, \alpha_2, \alpha_3}
    S_{\mu, \alpha_1}^{i_1}
    S_{\mu + \nu_2, \alpha_1}^{i_2}
    S_{\mu + \nu_3, \alpha_3}^{i_3}
    +\\&+
    C_4
    \sum_{\substack{\mu, \nu_2, \nu_3, \nu_4, \\ \alpha_1, \alpha_2, \alpha_3, \alpha_4, \\ i_1, i_2, i_3, i_4}}
    J^{i_1, i_2, i_3, i_4}_{\nu_2, \nu_3, \nu_4; \alpha_1, \alpha_2, \alpha_3, \alpha_4}
    S_{\mu, \alpha_1}^{i_1}
    S_{\mu + \nu_2, \alpha_2}^{i_2}
    S_{\mu + \nu_3, \alpha_3}^{i_3}
    S_{\mu + \nu_4, \alpha_4}^{i_4}
    :label: eq:user-guide_theory-behind_spin-hamiltonian_the_hamiltonian

where the sum over :math:`n` is expanded.

.. dropdown:: Meaning of indices

    The Hamiltonian is defined on a periodic lattice comprising of
    :math:`N = N_1 \times N_2 \times N_3` unit cells. The unit cell is defined by three
    lattice vectors

    .. math::

        \boldsymbol{a}_1 = (a_1^x, a_1^y, a_1^z)
        \\
        \boldsymbol{a}_2 = (a_2^x, a_2^y, a_2^z)
        \\
        \boldsymbol{a}_3 = (a_3^x, a_3^y, a_3^z)

    that are stored in Magnopy as (:py:attr:`.SpinHamiltonian.cell`)

    .. code-block:: python

        [[a1_x, a1_y, a1_z],
         [a2_x, a2_y, a2_z],
         [a3_x, a3_y, a3_z]]

    *   Superscript indices :math:`i_1, i_2, i_3, i_4` each assume values :math:`x, y, z`
        and refer to the Cartesian components of a vectors, matrix or tensor in the global
        orthonormal reference frame :math:`(\boldsymbol{x}, \boldsymbol{y}, \boldsymbol{z})`.

        For example, the scalar product of two vectors :math:`J^{i_1}_{\alpha_1}` and
        :math:`S_{\mu, \alpha_1}^{i_1}` is written as

        .. math::
            (\boldsymbol{S}_{\mu, \alpha_1} \cdot \boldsymbol{J}_{\alpha_1})
            =
            \sum_{i_1}
            S_{\mu, \alpha_1}^{i_1} \cdot J^{i_1}_{\alpha_1}
            =
            S_{\mu, \alpha_1}^x \cdot J^x_{\alpha_1}
            +
            S_{\mu, \alpha_1}^y \cdot J^y_{\alpha_1}
            +
            S_{\mu, \alpha_1}^z \cdot J^z_{\alpha_1}

    *   Subscript Greek letters :math:`\mu, \nu_m` denote the unit cells in the lattice.
        Each unit cell's index represents a tuple of three integers

        .. math::
            \mu &= (\mu^1, \mu^2, \mu^3)
            \\
            \nu_m &= (\nu_m^1, \nu_m^2, \nu_m^3)

        that define the position of the unit cell in real space relative to the three
        lattice vectors

        .. math::
            \boldsymbol{r}_{\mu+\nu_m}
            =
            (\mu^1 + \nu_m^1) \cdot \boldsymbol{a}_1
            +
            (\mu^2 + \nu_m^2) \cdot \boldsymbol{a}_2
            +
            (\mu^3 + \nu_m^3) \cdot \boldsymbol{a}_3

        The index :math:`\mu` is always taken care of analytically and the code does not
        need to store it (note that interaction parameters do not depend on :math:`\mu`,
        thus highlighting the translational invariance of the Hamiltonian). The indices
        :math:`\nu_m` are stored in Magnopy as tuples of three integers, for example

        .. code-block:: python

            nu_1 = (1, 0, -2)
            nu_2 = (0, 1, 0)

    *   Subscript Greek letters :math:`\alpha_m` denote the magnetic sites within the
        unit cell. Each site's index can be represented by the tuple of three real numbers
        from the half-open interval :math:`[0,1)`

        .. math::
            \alpha_m
            =
            (\alpha_m^1, \alpha_m^2, \alpha_m^3)

        that define the position of the magnetic site within the unit cell

        .. math::
            \boldsymbol{r}_{\alpha_m}
            =
            \alpha_m^1 \cdot \boldsymbol{a}_1
            +
            \alpha_m^2 \cdot \boldsymbol{a}_2
            +
            \alpha_m^3 \cdot \boldsymbol{a}_3

        However, in contrast to the unit cell indices the floating point numbers are ill
        suited for being indices, thus in Magnopy we store indices alpha as integers
        ranging from :math:`0` to :math:`M-1`, where M is the amount of magnetic sites in
        the unit cell, for example

        .. code-block:: python

            alpha_1 = 0
            alpha_2 = 3

    .. hint::
        The position of the magnetic site in real space is defined by the position of the
        unit cell of that site and by the position of the site within unit cell. It is
        denoted by the pair of subscript indices :math:`\mu+\nu_m, \alpha_m` and can be
        expressed as

        .. math::
            \boldsymbol{r}_{\mu+\nu_m, \alpha_m}
            =
            (\mu^1 + \nu_m^1+ \alpha_m^1) \cdot \boldsymbol{a}_1
            +
            (\mu^2 + \nu_m^2+ \alpha_m^2) \cdot \boldsymbol{a}_2
            +
            (\mu^3 + \nu_m^3+ \alpha_m^3) \cdot \boldsymbol{a}_3

For the implementation (and its use) you need to understand one extra concept about this
Hamiltonian that the |paper-2026|_ does not discuss in details: the distinction between
the terms with the same amount of spin operators (same :math:`n`) but different amount
of **unique** magnetic sites.

This distinction is implicit in the code, in particular you will always see an extra
index, that we label as :math:`p_n`, that goes in pair with the index :math:`n` in this
documentation and in the names of some methods implemented in Magnopy.

For the linear, bilinear and trilinear terms of the Hamiltonian the extra index :math:`p_n`
is trivial to construct. For example, there is one case for the linear term,
two cases for the bilinear term and three cases for the trilinear term. Therefore, the
extra index :math:`p_n` for :math:`n = 1, 2, 3` can be defined simply to indicate the
amount of unique magnetic sites in the term

However, for the quadlinear term (:math:`n = 4`) there are five cases that we consider to
be different. In particular the case with two unique magnetic sites splits into two cases:
one where each site has two spin operators associated with it and the other where one site
has three spin operators and the other site has one spin operator. In fact, the amount
of cases for an integer :math:`n` is defined by the amount of its partitions.

There are several strategies for labeling such cases, in Magnopy we do the following:

* Write down each partition of integer :math:`n` as a tuple of integers
* Sort each partition in descending order
* Sort partitions of the same integer lexicographically in descending order.
* Label each partition with an extra index :math:`p_n` that goes from :math:`1` to the
  amount of partitions of :math:`n`.

This strategy has an advantage: for :math:`n < 4` the extra index :math:`p_n` still
indicates the amount of unique magnetic sites. The table below summarizes the labeling
rules that are used in Magnopy. See linked pages for the details of the terms of the
Hamiltonian for each case.

========= =========== ============ =================== =======================================
:math:`n` :math:`p_n` partition    Hamiltonian's terms Convention constant (:math:`C_{n,p_n}`)
========= =========== ============ =================== =======================================
1         1           ( 1 )        :ref:`ug_tb_sh_1-1` :math:`C_{1,1} \equiv C_1`
2         1           (2, 0)       :ref:`ug_tb_sh_2-1` :math:`C_{2,1}`
2         2           (1, 1)       :ref:`ug_tb_sh_2-2` :math:`C_{2,2}`
3         1           (3, 0, 0)    :ref:`ug_tb_sh_3-1` :math:`C_{3,1}`
3         2           (2, 1, 0)    :ref:`ug_tb_sh_3-2` :math:`C_{3,2}`
3         3           (1, 1, 1)    :ref:`ug_tb_sh_3-3` :math:`C_{3,3}`
4         1           (4, 0, 0, 0) :ref:`ug_tb_sh_4-1` :math:`C_{4,1}`
4         2           (3, 1, 0, 0) :ref:`ug_tb_sh_4-2` :math:`C_{4,2}`
4         3           (2, 2, 0, 0) :ref:`ug_tb_sh_4-3` :math:`C_{4,3}`
4         4           (2, 1, 1, 0) :ref:`ug_tb_sh_4-4` :math:`C_{4,4}`
4         5           (1, 1, 1, 1) :ref:`ug_tb_sh_4-5` :math:`C_{4,5}`
========= =========== ============ =================== =======================================

.. toctree::
    :maxdepth: 1
    :hidden:

    spin-hamiltonian-subterms/1-1
    spin-hamiltonian-subterms/2-1
    spin-hamiltonian-subterms/2-2
    spin-hamiltonian-subterms/3-1
    spin-hamiltonian-subterms/3-2
    spin-hamiltonian-subterms/3-3
    spin-hamiltonian-subterms/4-1
    spin-hamiltonian-subterms/4-2
    spin-hamiltonian-subterms/4-3
    spin-hamiltonian-subterms/4-4
    spin-hamiltonian-subterms/4-5

Examples
========

In this part we show how the common terms of the spin Hamiltonian can be written in the
form of equation :eq:`eq:user-guide_theory-behind_spin-hamiltonian_the_hamiltonian`.

Zeeman interaction
------------------

Linear coupling with the magnetic field is usually written as

.. math::
    \mathcal{H}
    =
    \mu_B\boldsymbol{h}\sum_{\mu,\alpha_1} g_{\alpha_1} \boldsymbol{S}_{\mu, \alpha_1}

This term can be written in the form of :ref:`ug_tb_sh_1-1` if one defines
:math:`C_1 = 1` and :math:`\boldsymbol{J}_{\alpha_1} = \mu_B\boldsymbol{h} g_{\alpha_1}`.

On-site anisotropy
------------------

We take an example of the triaxial anisotropy, that can be written as

.. math::
    \mathcal{H}
    =
    \sum_{\mu,\alpha_1}
    \Bigl(
        K^x (S^x_{\mu,\alpha_1})^2
        +
        K^y (S^y_{\mu,\alpha_1})^2
        +
        K^z (S^z_{\mu,\alpha_1})^2
    \Bigr)

This Hamiltonian can be written in the form of :ref:`ug_tb_sh_2-1` if one defines
:math:`C_{2, 1} = 1` and

.. math::
    \boldsymbol{J}_{0;\alpha_1,\alpha_1}
    =
    \begin{pmatrix}
        K^x & 0 & 0 \\
        0 & K^y & 0 \\
        0 & 0 & K^z
    \end{pmatrix}

Exchange interaction
--------------------

Bilinear exchange interaction with isotropic and Dzyaloshinskii-Moriya (DM) exchange can
be written as

.. math::
    \mathcal{H}
    =
    \dfrac{1}{2}
    \sum_{\mu,\nu_2,\alpha_1,\alpha_2}
    \Bigl[
    J^{iso}_{\nu_2;\alpha_1,\alpha_2}
    (\boldsymbol{S}_{\mu,\alpha_1}
    \cdot
    \boldsymbol{S}_{\mu+\nu_2,\alpha_2})
    +
    \boldsymbol{D}_{\nu_2;\alpha_1,\alpha_2}
    \cdot
    (\boldsymbol{S}_{\mu,\alpha_1}
    \times
    \boldsymbol{S}_{\mu+\nu_2,\alpha_2})
    \Bigr]

where :math:`\boldsymbol{D}_{\nu_2;\alpha_1,\alpha_2}` is a DM vector.
This Hamiltonian can be written in the form of :ref:`ug_tb_sh_2-2` if one defines
:math:`C_{2, 2} = 1/2` and

.. math::
    \boldsymbol{J}_{\nu_2;\alpha_1,\alpha_2}
    =
    \begin{pmatrix}
        J^{iso}_{\nu_2;\alpha_1,\alpha_2} & D^z_{\nu_2;\alpha_1,\alpha_2} & -D^y_{\nu_2;\alpha_1,\alpha_2} \\
        -D^z_{\nu_2;\alpha_1,\alpha_2} & J^{iso}_{\nu_2;\alpha_1,\alpha_2} & D^x_{\nu_2;\alpha_1,\alpha_2} \\
        D^y_{\nu_2;\alpha_1,\alpha_2} & -D^x_{\nu_2;\alpha_1,\alpha_2} & J^{iso}_{\nu_2;\alpha_1,\alpha_2}
    \end{pmatrix}

Biquadratic exchange
--------------------

Isotropic biquadratic exchange interaction can be written as

.. math::
    \mathcal{H}
    =
    \sum_{\mu,\nu_2,\alpha_1,\alpha_2}
    J_{\nu_2;\alpha_1,\alpha_2}
    (\boldsymbol{S}_{\mu,\alpha_1}
    \cdot
    \boldsymbol{S}_{\mu+\nu_2,\alpha_2})^2

This Hamiltonian can be written in the form of :ref:`ug_tb_sh_4-3` if one defines
:math:`C_{4,3} = 1` and
:math:`J^{i_1,i_2,i_3,i_4}_{0,\nu_2,\nu_2;\alpha_1, \alpha_1, \alpha_2, \alpha_2} = J_{\nu_2;\alpha_1,\alpha_2}`
when :math:`(i_1,i_2,i_3,i_4) = (xxxx), (xyxy), (xzxz), (yxyx), (yyyy), (yzyz), (zxzx), (zyzy), (zzzz)`
and
:math:`J^{i_1,i_2,i_3,i_4}_{0,\nu_2,\nu_2;\alpha_1, \alpha_1, \alpha_2, \alpha_2} = 0`
otherwise.
