.. _user-guide_theory-behind_multiple-counting:

*****************
Multiple counting
*****************

Multiple counting occurs in some terms of the Hamiltonian. In those terms the
same set of magnetic centers enters the sum multiple times with different order of
spin operators.

For the sake of the choosing the primary version of the same bond we define the
comparison of two spin operators based on their positions as
(where :math:`\mu = (\mu^1, \mu^2, \mu^3)`)

* :math:`\boldsymbol{S}_{\mu_1,\alpha_1} < \boldsymbol{S}_{\mu_2,\alpha_2}` if one of the
  conditions is met

  * :math:`\mu_2^1 - \mu_1^1 > 0`
  * :math:`\mu_2^1 - \mu_1^1 = 0, \mu_2^2 - \mu_1^2 > 0`
  * :math:`\mu_2^1 - \mu_1^1 = \mu_2^2 - \mu_1^2 = 0, \mu_2^3 - \mu_1^3 > 0`
  * :math:`\mu_2^1 - \mu_1^1 = \mu_2^2 - \mu_1^2 = \mu_2^3 - \mu_1^3 = 0, \alpha_2 - \alpha_1 > 0`

* :math:`\boldsymbol{S}_{\mu_1,\alpha_1} > \boldsymbol{S}_{\mu_2,\alpha_2}` otherwise.

Moreover, we use in this page a short notation for the terms of the Hamiltonian

* With two sites:

  .. math::

    (\alpha, \beta, \nu) &\rightarrow
    J_{2,2}^{ij}(\boldsymbol{r}_{\nu,\alpha\beta})
    S_{\mu,\alpha}^i
    S_{\mu+\nu,\beta}^j
    \\
    (\alpha, \beta, \nu) &\rightarrow
    J_{3, 2}^{iju}(\boldsymbol{r}_{\nu,\alpha\beta})
    S_{\mu,\alpha}^i
    S_{\mu,\alpha}^j
    S_{\mu+\nu,\beta}^u
    \\
    (\alpha, \beta, \nu) &\rightarrow
    J_{4, 2, 1}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})
    S_{\mu,\alpha}^i
    S_{\mu,\alpha}^j
    S_{\mu,\alpha}^u
    S_{\mu+\nu,\beta}^v
    \\
    (\alpha, \beta, \nu) &\rightarrow
    J_{4, 2, 2}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})
    S_{\mu,\alpha}^i
    S_{\mu,\alpha}^j
    S_{\mu+\nu,\beta}^u
    S_{\mu+\nu,\beta}^v

* With three sites

  .. math::

    (\alpha, \beta, \gamma, \nu, \lambda) &\rightarrow
    J_{3, 3}^{iju}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
    S_{\mu,\alpha}^i
    S_{\mu+\nu,\beta}^j
    S_{\mu+\lambda,\gamma}^u
    \\
    (\alpha, \beta, \gamma, \nu, \lambda) &\rightarrow
    J_{4, 3}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
    S_{\mu,\alpha}^i
    S_{\mu,\alpha}^j
    S_{\mu+\nu,\beta}^u
    S_{\mu+\lambda,\gamma}^v

* With four sites

  .. math::

    (\alpha, \beta, \gamma, \varepsilon, \nu, \lambda, \rho) \rightarrow
    J_{4, 3}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})
    S_{\mu,\alpha}^i
    S_{\mu+\nu,\beta}^j
    S_{\mu+\lambda,\gamma}^u
    S_{\mu+\rho,\varepsilon}^v



Two spins & two sites
=====================

For the given term of the Hamiltonian

.. math::

    J_{2,2}^{ij}(\boldsymbol{r}_{\nu,\alpha\beta})
    S_{\mu,\alpha}^i
    S_{\mu+\nu,\beta}^j

two cases are possible

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\alpha, \beta, \nu)` is a primary version with the parameter

    .. math::
      J_{2,2}^{ij}(\boldsymbol{r}_{\nu,\alpha\beta})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} > \boldsymbol{S}_{\mu+\nu,\beta}`, then
    :math:`(\beta, \alpha, -\nu)` is a primary version with the parameter

    .. math::
      J_{2,2}^{ij}(\boldsymbol{r}_{-\nu,\beta\alpha})
      =
      J_{2,2}^{ji}(\boldsymbol{r}_{\nu,\alpha\beta})

Three spins & two sites
=======================

For the given term of the Hamiltonian

.. math::

    J_{3, 2}^{iju}(\boldsymbol{r}_{\nu,\alpha\beta})
    S_{\mu,\alpha}^i
    S_{\mu,\alpha}^j
    S_{\mu+\nu,\beta}^u

two cases are possible

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\alpha, \beta, \nu)` is a primary version with the parameter

    .. math::
      J_{3, 2}^{iju}(\boldsymbol{r}_{\nu,\alpha\beta})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} > \boldsymbol{S}_{\mu+\nu,\beta}`, then
    :math:`(\beta, \alpha, -\nu)` is a primary version with the parameter

    .. math::
      J_{3,2}^{iju}(\boldsymbol{r}_{-\nu,\beta\alpha})
      =
      \dfrac{S_{\alpha}}{S_{\beta}}
      J_{3,2}^{uji}(\boldsymbol{r}_{\nu,\alpha\beta})

Three spins & three sites
=========================

For the given term of the Hamiltonian

.. math::

    J_{3, 3}^{iju}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
    S_{\mu,\alpha}^i
    S_{\mu+\nu,\beta}^j
    S_{\mu+\lambda,\gamma}^u

six cases are possible

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\alpha, \beta, \gamma, \nu, \lambda)` is a primary version with the parameter

    .. math::
      J_{3, 3}^{iju}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\alpha, \gamma, \beta, \lambda, \nu)` is a primary version with the parameter

    .. math::
      J_{3, 3}^{iju}(\boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\nu,\alpha\beta})
      =
      J_{3, 3}^{iuj}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\beta, \alpha, \gamma, -\nu, \lambda - \nu)` is a primary version with the parameter

    .. math::
      J_{3, 3}^{iju}(\boldsymbol{r}_{-\nu,\beta\alpha}, \boldsymbol{r}_{\lambda - \nu,\beta\gamma})
      =
      J_{3, 3}^{jiu}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\beta, \gamma, \alpha, \lambda - \nu, -\nu)` is a primary version with the parameter

    .. math::
      J_{3, 3}^{iju}(\boldsymbol{r}_{\lambda - \nu,\beta\gamma}, \boldsymbol{r}_{-\nu,\beta\alpha})
      =
      J_{3, 3}^{uij}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\gamma, \alpha, \beta, -\lambda, \nu - \lambda)` is a primary version with the parameter

    .. math::
      J_{3, 3}^{iju}(\boldsymbol{r}_{-\lambda, \gamma\alpha}, \boldsymbol{r}_{\nu - \lambda,\gamma\beta})
      =
      J_{3, 3}^{jui}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\gamma, \beta, \alpha, \nu - \lambda, -\lambda)` is a primary version with the parameter

    .. math::
      J_{3, 3}^{iju}(\boldsymbol{r}_{\nu - \lambda,\gamma\beta}, \boldsymbol{r}_{-\lambda, \gamma\alpha})
      =
      J_{3, 3}^{uji}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

Four spins & two sites (1+3)
============================

For the given term of the Hamiltonian

.. math::

    J_{4, 2, 1}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})
    S_{\mu,\alpha}^i
    S_{\mu,\alpha}^j
    S_{\mu,\alpha}^u
    S_{\mu+\nu,\beta}^v

two cases are possible

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\alpha, \beta, \nu)` is a primary version with the parameter

    .. math::
      J_{4, 2, 1}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} > \boldsymbol{S}_{\mu+\nu,\beta}`, then
    :math:`(\beta, \alpha, -\nu)` is a primary version with the parameter

    .. math::
      J_{4, 2, 1}^{ijuv}(\boldsymbol{r}_{-\nu,\beta\alpha})
      =
      \left(\dfrac{S_{\alpha}}{S_{\beta}}\right)^2
      J_{4, 2, 1}^{vjui}(\boldsymbol{r}_{\nu,\alpha\beta})

Four spins & two sites (2+2)
============================

For the given term of the Hamiltonian

.. math::

    J_{4, 2, 2}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})
    S_{\mu,\alpha}^i
    S_{\mu,\alpha}^j
    S_{\mu+\nu,\beta}^u
    S_{\mu+\nu,\beta}^v

two cases are possible

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\alpha, \beta, \nu)` is a primary version with the parameter

    .. math::
      J_{4, 2, 2}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} > \boldsymbol{S}_{\mu+\nu,\beta}`, then
    :math:`(\beta, \alpha, -\nu)` is a primary version with the parameter

    .. math::
      J_{4, 2, 2}^{ijuv}(\boldsymbol{r}_{-\nu,\beta\alpha})
      =
      J_{4, 2, 2}^{uvij}(\boldsymbol{r}_{\nu,\alpha\beta})

Four spins & three sites
========================

For the given term of the Hamiltonian

.. math::

    J_{4, 3}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
    S_{\mu,\alpha}^i
    S_{\mu,\alpha}^j
    S_{\mu+\nu,\beta}^u
    S_{\mu+\lambda,\gamma}^v

six cases are possible

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\alpha, \beta, \gamma, \nu, \lambda)` is a primary version with the parameter

    .. math::
      J_{4, 3}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\alpha, \gamma, \beta, \lambda, \nu)` is a primary version with the parameter

    .. math::
      J_{4, 3}^{ijuv}(\boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\nu,\alpha\beta})
      =
      J_{4, 3}^{ijvu}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\beta, \alpha, \gamma, -\nu, \lambda - \nu)` is a primary version with the parameter

    .. math::
      J_{4, 3}^{ijuv}(\boldsymbol{r}_{-\nu,\beta\alpha}, \boldsymbol{r}_{\lambda - \nu,\beta\gamma})
      =
      \dfrac{S_{\alpha}}{S_{\beta}}
      J_{4, 3}^{ujiv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\beta, \gamma, \alpha, \lambda - \nu, -\nu)` is a primary version with the parameter

    .. math::
      J_{4, 3}^{ijuv}(\boldsymbol{r}_{\lambda - \nu,\beta\gamma}, \boldsymbol{r}_{-\nu,\beta\alpha})
      =
      \dfrac{S_{\alpha}}{S_{\beta}}
      J_{4, 3}^{vjiu}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\gamma, \alpha, \beta, -\lambda, \nu - \lambda)` is a primary version with the parameter

    .. math::
      J_{4, 3}^{ijuv}(\boldsymbol{r}_{-\lambda, \gamma\alpha}, \boldsymbol{r}_{\nu - \lambda,\gamma\beta})
      =
      \dfrac{S_{\alpha}}{S_{\gamma}}
      J_{4, 3}^{ujvi}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\gamma, \beta, \alpha, \nu - \lambda, -\lambda)` is a primary version with the parameter

    .. math::
      J_{4, 3}^{ijuv}(\boldsymbol{r}_{\nu - \lambda,\gamma\beta}, \boldsymbol{r}_{-\lambda, \gamma\alpha})
      =
      \dfrac{S_{\alpha}}{S_{\gamma}}
      J_{4, 3}^{vjui}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})

Four spins & four sites
=======================

For the given term of the Hamiltonian

.. math::

    J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})
    S_{\mu,\alpha}^i
    S_{\mu+\nu,\beta}^j
    S_{\mu+\lambda,\gamma}^u
    S_{\mu+\rho,\varepsilon}^v

twenty four cases are possible

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\rho,\varepsilon}` then
    :math:`(\alpha, \beta, \gamma, \varepsilon, \nu, \lambda, \rho)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\alpha, \beta, \varepsilon, \gamma, \nu, \rho, \lambda)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\rho,\alpha\varepsilon}, \boldsymbol{r}_{\lambda,\alpha\gamma})
      =
      J_{4, 4}^{ijvu}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\rho,\varepsilon}` then
    :math:`(\alpha, \gamma, \beta, \varepsilon, \lambda, \nu, \rho)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\rho,\alpha\varepsilon})
      =
      J_{4, 4}^{iujv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\alpha, \gamma, \varepsilon, \beta, \lambda, \rho, \nu)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon}, \boldsymbol{r}_{\nu,\alpha\beta})
      =
      J_{4, 4}^{ivju}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\alpha, \varepsilon, \beta, \gamma, \rho, \nu, \lambda)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\rho,\alpha\varepsilon}, \boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
      =
      J_{4, 4}^{iuvj}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\alpha, \varepsilon, \gamma, \beta, \rho, \lambda, \nu)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\rho,\alpha\varepsilon}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\nu,\alpha\beta})
      =
      J_{4, 4}^{ivuj}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\rho,\varepsilon}` then
    :math:`(\beta, \alpha, \gamma, \varepsilon, -\nu, \lambda-\nu, \rho-\nu)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{-\nu,\beta\alpha}, \boldsymbol{r}_{\lambda-\nu,\beta\gamma}, \boldsymbol{r}_{\rho-\nu,\beta\varepsilon})
      =
      J_{4, 4}^{jiuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\beta, \alpha, \gamma, \varepsilon, -\nu, \rho-\nu, \lambda-\nu)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{-\nu,\beta\alpha}, \boldsymbol{r}_{\rho-\nu,\beta\varepsilon}, \boldsymbol{r}_{\lambda-\nu,\beta\gamma})
      =
      J_{4, 4}^{jivu}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\rho,\varepsilon}` then
    :math:`(\beta, \gamma, \alpha, \varepsilon, \lambda-\nu, -\nu, \rho-\nu)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\lambda-\nu,\beta\gamma}, \boldsymbol{r}_{-\nu,\beta\alpha}, \boldsymbol{r}_{\rho-\nu,\beta\varepsilon})
      =
      J_{4, 4}^{uijv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\beta, \gamma, \varepsilon, \alpha, \lambda-\nu, \rho-\nu, -\nu)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\lambda-\nu,\beta\gamma}, \boldsymbol{r}_{\rho-\nu,\beta\varepsilon}, \boldsymbol{r}_{-\nu,\beta\alpha})
      =
      J_{4, 4}^{viju}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\beta, \varepsilon, \alpha, \gamma, \rho-\nu, -\nu, \lambda-\nu)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\rho-\nu,\beta\varepsilon}, \boldsymbol{r}_{-\nu,\beta\alpha}, \boldsymbol{r}_{\lambda-\nu,\beta\gamma})
      =
      J_{4, 4}^{uivj}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\beta, \varepsilon, \gamma, \alpha, \rho-\nu, \lambda-\nu, -\nu)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\rho-\nu,\beta\varepsilon}, \boldsymbol{r}_{\lambda-\nu,\beta\gamma}, \boldsymbol{r}_{-\nu,\beta\alpha})
      =
      J_{4, 4}^{viuj}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\rho,\varepsilon}` then
    :math:`(\gamma, \alpha, \beta, \varepsilon, -\lambda, \nu-\lambda, \rho-\lambda)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{-\lambda,\gamma\alpha}, \boldsymbol{r}_{\nu-\lambda,\gamma\beta}, \boldsymbol{r}_{\rho-\lambda,\gamma\varepsilon})
      =
      J_{4, 4}^{juiv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\gamma, \alpha, \varepsilon, \beta, -\lambda, \rho-\lambda, \nu-\lambda)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{-\lambda,\gamma\alpha}, \boldsymbol{r}_{\rho-\lambda,\gamma\varepsilon}, \boldsymbol{r}_{\nu-\lambda,\gamma\beta})
      =
      J_{4, 4}^{jviu}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\rho,\varepsilon}` then
    :math:`(\gamma, \beta, \alpha, \varepsilon, \nu-\lambda, -\lambda, \rho-\lambda)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu-\lambda,\gamma\beta}, \boldsymbol{r}_{-\lambda,\gamma\alpha}, \boldsymbol{r}_{\rho-\lambda,\gamma\varepsilon})
      =
      J_{4, 4}^{ujiv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\gamma, \beta, \varepsilon, \alpha, \nu-\lambda, \rho-\lambda, -\lambda)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu-\lambda,\gamma\beta}, \boldsymbol{r}_{\rho-\lambda,\gamma\varepsilon}, \boldsymbol{r}_{-\lambda,\gamma\alpha})
      =
      J_{4, 4}^{vjiu}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\gamma, \varepsilon, \alpha, \beta, \rho-\lambda, -\lambda, \nu-\lambda)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\rho-\lambda,\gamma\varepsilon}, \boldsymbol{r}_{-\lambda,\gamma\alpha}, \boldsymbol{r}_{\nu-\lambda,\gamma\beta})
      =
      J_{4, 4}^{uvij}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\gamma, \varepsilon, \beta, \alpha, \rho-\lambda, \nu-\lambda, -\lambda)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\rho-\lambda,\gamma\varepsilon}, \boldsymbol{r}_{\nu-\lambda,\gamma\beta}, \boldsymbol{r}_{-\lambda,\gamma\alpha})
      =
      J_{4, 4}^{vuij}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\varepsilon, \alpha, \beta, \gamma, -\rho, \nu-\rho, \lambda-\rho)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{-\rho,\varepsilon,\alpha}, \boldsymbol{r}_{\nu-\rho,\varepsilon,\beta}, \boldsymbol{r}_{\lambda-\rho,\varepsilon,\gamma})
      =
      J_{4, 4}^{juvi}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\varepsilon, \alpha, \gamma, \beta, -\rho, \lambda-\rho, \nu-\rho)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{-\rho,\varepsilon,\alpha}, \boldsymbol{r}_{\lambda-\rho,\varepsilon,\gamma}, \boldsymbol{r}_{\nu-\rho,\varepsilon,\beta})
      =
      J_{4, 4}^{jvui}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\lambda,\gamma}` then
    :math:`(\varepsilon, \beta, \alpha, \gamma, \nu-\rho, -\rho, \lambda-\rho)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu-\rho,\varepsilon,\beta}, \boldsymbol{r}_{-\rho,\varepsilon,\alpha}, \boldsymbol{r}_{\lambda-\rho,\varepsilon,\gamma})
      =
      J_{4, 4}^{ujvi}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\varepsilon, \beta, \gamma, \alpha, \nu-\rho, \lambda-\rho, -\rho)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu-\rho,\varepsilon,\beta}, \boldsymbol{r}_{\lambda-\rho,\varepsilon,\gamma}, \boldsymbol{r}_{-\rho,\varepsilon,\alpha})
      =
      J_{4, 4}^{vjui}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu,\alpha} < \boldsymbol{S}_{\mu+\nu,\beta}` then
    :math:`(\varepsilon, \gamma, \alpha, \beta, \lambda-\rho, -\rho, \nu-\rho)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\lambda-\rho,\varepsilon,\gamma}, \boldsymbol{r}_{-\rho,\varepsilon,\alpha}, \boldsymbol{r}_{\nu-\rho,\varepsilon,\beta})
      =
      J_{4, 4}^{uvji}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})

#.  If :math:`\boldsymbol{S}_{\mu+\rho,\varepsilon} < \boldsymbol{S}_{\mu+\lambda,\gamma} < \boldsymbol{S}_{\mu+\nu,\beta} < \boldsymbol{S}_{\mu,\alpha}` then
    :math:`(\varepsilon, \gamma, \beta, \alpha, \lambda-\rho, \nu-\rho, -\rho)` is a primary version with the parameter

    .. math::
      J_{4, 4}^{ijuv}(\boldsymbol{r}_{\lambda-\rho,\varepsilon,\gamma}, \boldsymbol{r}_{\nu-\rho,\varepsilon,\beta}, \boldsymbol{r}_{-\rho,\varepsilon,\alpha})
      =
      J_{4, 4}^{vuji}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})
