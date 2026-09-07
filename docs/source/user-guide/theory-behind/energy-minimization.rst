.. _user-guide_theory-behind_energy-minimization:

**********************
Minimization of energy
**********************


References
==========

.. [1] Nocedal, J. and Wright, S.J.
       Numerical optimization. New York, NY: Springer New York.
       eds., 1999.
.. [2] Ivanov, A.V., Uzdin, V.M. and Jónsson, H., 2021.
	Fast and robust algorithm for energy minimization of spin systems applied
	in an analysis of high temperature spin configurations in terms of skyrmion
	density.
	Computer Physics Communications, 260, p.107749.


On this page we recall the approach and main formulas from [1]_ and [2]_. The minimization
of the magnetic ground state in Magnopy is an implementation of the method described in
[2]_.

Minimization of the energy function (:math:`E^{(0)}` or :math:`E^{(0)} + E^{corr}`) can be
formulated as a problem of minimizing that function over the :math:`M` vectors of the spin
directions :math:`\boldsymbol{z}_{\alpha}, \alpha = 1, ..., M`

.. math::

	E = F(\boldsymbol{z}_{1}, ..., \boldsymbol{z}_{M})

Directional vectors are unit vectors and vary on the sphere. This fact introduces
complications in the minimization procedure as the optimization space is not a vector
space and the typical (|BFGS|_, for instance) algorithms for linear optimizations cannot
be applied directly. This problem is elegantly solved via parametrization of directional
vectors with the exponents of skew-symmetric matrices [2]_. Given an initial guess
:math:`\boldsymbol{z}_{\alpha}^{(0)}`, any other set of directional vectors
can be obtained by the following formulae:

.. math::

	\boldsymbol{z}_{\alpha}
    =
    e^{\boldsymbol{A}_{\alpha}} \boldsymbol{z}_{\alpha}^{(0)}
    =
    \cos(\theta_{\alpha}) \,\boldsymbol{z}_{\alpha}^{(0)}
    +
    \sin(\theta_{\alpha}) \,(\boldsymbol{r}_{\alpha}\times\boldsymbol{z}_{\alpha}^{(0)})
    +
    (1 - \cos(\theta_{\alpha})) (\boldsymbol{r}_{\alpha} \cdot \boldsymbol{z}_{\alpha}^{(0)})\, \boldsymbol{r}_{\alpha}

where :math:`\boldsymbol{A}_{\alpha}` are skew-symmetric matrices parametrized by three
real numbers as

.. math::

	\boldsymbol{A}_{\alpha}
	=
	\begin{pmatrix}
		0 & -a_{\alpha}^z & a_{\alpha}^y \\
		a_{\alpha}^z & 0 & -a_{\alpha}^x \\
		-a_{\alpha}^y & a_{\alpha}^x & 0
	\end{pmatrix}

and

.. math::

    \theta_{\alpha}
    =
    \sqrt{\left(a_{\alpha}^x\right)^2
    +
    \left(a_{\alpha}^y\right)^2
    +
    \left(a_{\alpha}^z\right)^2}

.. math::

    \boldsymbol{r}_{\alpha}
    =
    \dfrac{(a_{\alpha}^x,
    a_{\alpha}^y,
    a_{\alpha}^z)}{\theta_{\alpha}}


Then energy function can be rewritten as

.. math::

	E
	=
	F(
		e^{\boldsymbol{A}_1} \boldsymbol{z}_{1}^{(0)},
		...,
		e^{\boldsymbol{A}_I} \boldsymbol{z}_{I}^{(0)}
	)

In other words, as a function of vector :math:`\boldsymbol{x}` from the vector space
:math:`\mathbb{R}^{3I}`:

.. math::

	E = F(\boldsymbol{x})
	\qquad
	\boldsymbol{x}
	=(
		a_{1}^x, a_{1}^y, a_{1}^z,
		...,
		a_{I}^x, a_{I}^y, a_{I}^z
	)

Then, energy of the system is minimized with the BFGS algorithm [1]_.


Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm
=================================================

Formula for the inverse Hessian update:

.. math::

	H^{ij}_{k+1}
	=
	\sum_{u,v}(\delta_{i,u} - \rho_ks^i_ky^u_k)
	H^{uv}_k
	(\delta_{v,j} - \rho_ky^v_ks^j_k) + \rho_k s^i_ks^j_k,
	\qquad
	\rho_k = \dfrac{1}{\sum_i y^i_k s^i_k}

Given :ref:`user-guide_theory-behind_energy-minimization_initial-guess`
:math:`\boldsymbol{x}_0` and
:ref:`user-guide_theory-behind_energy-minimization_initial-hessian`
:math:`\boldsymbol{H}_0`,


1.  :math:`k \gets 0`
#.  While convergence is not achieved:

    a)  :ref:`Compute the gradient <user-guide_theory-behind_energy-minimization_gradient>`
        of the function :math:`\boldsymbol{\nabla} F(\boldsymbol{x}_k)`;
    #)  Compute the search direction
        :math:`\boldsymbol{p}_k = -\boldsymbol{H}_k \boldsymbol{\nabla} F(\boldsymbol{x}_k)`;
    #)  Compute length of the step :math:`\alpha_k` via
        :ref:`user-guide_theory-behind_energy-minimization_line-search`;
    #)  Set :math:`\boldsymbol{x}_{k+1} = \boldsymbol{x}_k + \alpha_k \boldsymbol{p}_k`
        and compute gradient :math:`\boldsymbol{\nabla} F(\boldsymbol{x}_{k+1})`;
    #)  Set :math:`\boldsymbol{s}_k = \boldsymbol{x}_{k+1} - \boldsymbol{x}_k` and
        :math:`\boldsymbol{y}_k = \boldsymbol{\nabla} F(\boldsymbol{x}_{k+1}) - \boldsymbol{\nabla} F(\boldsymbol{x}_k)`;
    #)  Update the Hessian matrix :math:`\boldsymbol{H}_{k+1}` by the BFGS formula;
    #)  :math:`k \gets k + 1`.


.. note::
	In our implementation we update the direction vectors at the end of each iteration
	(i.e. at step 2.g). Therefore, the vector :math:`\boldsymbol{x}_k` is always equal to
	:math:`( 0, 0, 0, 0, 0, 0, ..., 0, 0, 0)`.



.. _user-guide_theory-behind_energy-minimization_initial-guess:

Initial guess
=============

Initial guess is provided by the user or randomly generated. The user provides three
components of each directional vector :math:`(z_{\alpha}^x, z_{\alpha}^y, z_{\alpha}^z)`.

.. _user-guide_theory-behind_energy-minimization_initial-hessian:

Initial approximation of the inverse hessian matrix
===================================================

We take an identity matrix as an initial approximation of the hessian matrix and
scale it as

.. math::
    \boldsymbol{H}_0
    =
    \dfrac{\boldsymbol{y}^T_k\boldsymbol{s}_k}{\boldsymbol{y}^T_k\boldsymbol{y}_k}\boldsymbol{I}

before the first update [1]_.


.. _user-guide_theory-behind_energy-minimization_gradient:

Gradient of the function F(x)
=============================

As we choose to update the direction vectors at each step of the BFGS algorithm, then
the gradient with respect to these variables can be computed as [2]_

.. math::
    \dfrac{\partial F}{\partial\boldsymbol{a}_{\alpha}}
    =
	\boldsymbol{t}_{\alpha}
	=
	\boldsymbol{z}_{\alpha} \times \dfrac{\partial E^{(0)}}{\partial\boldsymbol{z}_{\alpha}}

where :math:`\boldsymbol{t}_{\alpha}` is a torque vector and
:math:`\boldsymbol{a}_{\alpha} = (a_{\alpha}^x, a_{\alpha}^y, a_{\alpha}^z)`.

The gradient of the classical energy is computed analytically

.. math::

	\dfrac{\partial E^{(0)}}{\partial z^i_{\alpha}}
	=
    S_{\alpha}\tilde{J}_{\alpha}^i

where :math:`\tilde{J}_{\alpha}^i` is a single-spin renormalized parameter defined by
equation S.68 of |paper-2026-SI|_.

The gradient of the correction energy :math:`E^{corr}` is computed numerically by the two
point formula.

.. _user-guide_theory-behind_energy-minimization_line-search:

Line search
===========

Line search algorithm defines an optimal step length (:math:`\alpha`) for the search
direction :math:`\boldsymbol{p}_k`. It is obtained by minimizing the function

.. math::

	f(\alpha) = F(\boldsymbol{x}_k + \alpha \boldsymbol{p}_k),
	\qquad
	\dfrac{d f(\alpha)}{d \alpha} = \boldsymbol{\nabla} F(\boldsymbol{x}_k + \alpha \boldsymbol{p}_k) \boldsymbol{p}_k

enough to satisfy strong Wolfe conditions:

.. math::

	F(\boldsymbol{x}_k + \alpha\boldsymbol{p}_k)
	&\le
	F(\boldsymbol{x}_k) + c_1 \alpha_k \boldsymbol{\nabla} F(\boldsymbol{x}_k) \boldsymbol{p}_k,
	\\
	\vert\boldsymbol{\nabla} F(\boldsymbol{x}_k + \alpha\boldsymbol{p}_k)\boldsymbol{p}_k\vert
	&\le
	c_2\vert\boldsymbol{\nabla} F(\boldsymbol{x}_k)\boldsymbol{p}_k\vert

Line search algorithm:

Given :math:`\boldsymbol{x}_k` and :math:`\boldsymbol{p}_k`

1.  If :math:`\alpha = 1` satisfies strong Wolfe condition, then return :math:`1`.
#.  Set :math:`\alpha_0 = 0`, :math:`\alpha_{\text{max}} = 2` and choose :math:`\alpha_1`
    via :ref:`user-guide_theory-behind_energy-minimization_cubic-interpolation`;
#.  :math:`i \gets 1`;
#.  While maximum number of iterations is not achieved:

    a)  Compute :math:`f(\alpha_i) = F(\boldsymbol{x}_k + \alpha_i \boldsymbol{p}_k)`;
    #)  If :math:`f(\alpha_i) > f(0) + c_1 \alpha_i f^{\prime}(0)`
        or :math:`f(\alpha_i) \ge f(\alpha_{i-1})`
        and :math:`i > 1`, then return :math:`zoom(\alpha_{i-1}, \alpha_i)`;
    #)  Compute :math:`f^{\prime}(\alpha_i) = \boldsymbol{\nabla} F(\boldsymbol{x}_k + \alpha_i \boldsymbol{p}_k) \boldsymbol{p}_k`;
    #)  If :math:`\vert f^{\prime}(\alpha_i)\vert \le -c_2 f^{\prime}(0)`,
        then return :math:`\alpha_i`;
    #)  If :math:`f^{\prime}(\alpha_i) \ge 0`,
        then return :math:`zoom(\alpha_i, \alpha_{i-1})`;
    #)  Choose :math:`\alpha_{i+1}` via :ref:`user-guide_theory-behind_energy-minimization_cubic-interpolation`;
    #)  :math:`i \gets i + 1`.


:math:`zoom` algorithm:

Given :math:`\alpha_{lo}`, :math:`\alpha_{hi}`

1.  Repeat

    a)  Interpolate :math:`\alpha_j` via :ref:`user-guide_theory-behind_energy-minimization_cubic-interpolation`;
    #)  Compute :math:`f(\alpha_j) = F(\boldsymbol{x}_k + \alpha_j \boldsymbol{p}_k)`;
    #)  Check that value of the function sufficiently decreases.
    #)  If :math:`f(\alpha_j) > f(0) + c_1 \alpha_j f^{\prime}(0)`
        or :math:`f(\alpha_j) \ge f(\alpha_{lo})`,
        then :math:`\alpha_{hi} \gets \alpha_j`
    #)  Else

        i)  If :math:`\vert f^{\prime}(\alpha_j)\vert \le -c_2 f^{\prime}(0)`,
            then return :math:`\alpha_j`;
        #)  If :math:`f^{\prime}(\alpha_j)(\alpha_{hi} - \alpha_{lo}) \ge 0`,
            then :math:`\alpha_{hi} \gets \alpha_{lo}`;
        #) :math:`\alpha_{lo} \gets \alpha_j`.


.. _user-guide_theory-behind_energy-minimization_cubic-interpolation:

Cubic interpolation
-------------------

Given :math:`\alpha_l`, :math:`\alpha_h` and :math:`f(\alpha_l)`, :math:`f(\alpha_h)`
and :math:`f^{\prime}(\alpha_l)`, :math:`f^{\prime}(\alpha_h)` compute new :math:`\alpha_m`
as

.. math::

	\alpha_{min} &= \alpha_h - (\alpha_h - \alpha_l) \dfrac{f^{\prime}(\alpha_h) + d_2 - d_1}{f^{\prime}(\alpha_h) - f^{\prime}(\alpha_l) + 2d_2}
	\\
	d_1 &= f^{\prime}(\alpha_l) + f^{\prime}(\alpha_h) - 3 \dfrac{f(\alpha_l) - f(\alpha_h)}{\alpha_l - \alpha_h}
	\\
	d_2 &= \text{sign}(\alpha_h - \alpha_l) \sqrt{d_1^2 - f^{\prime}(\alpha_l)f^{\prime}(\alpha_h)}

If :math:`d_1^2 - f^{\prime}(\alpha_l)f^{\prime}(\alpha_h) < 0`, then
:math:`\alpha_{min} = \alpha_l` if :math:`f(\alpha_l) \le f(\alpha_h)`, otherwise
:math:`\alpha_{min} = \alpha_h`.
