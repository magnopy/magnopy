.. _development_theory-notes:

************
Theory notes
************

In this page we write some formulas that did not make it neither to the publication
about magnopy nor to the user guide, but should be listed somewhere for the convenience
of the developers.


Classical energy
================

The formulas for the classical energy are written for the case of non-normalized spins,
all other convention properties being arbitrary.

.. math::
    E^{(0)}
    =&
    \,C_1
    \sum_{\alpha, i}
    J_1^i(\boldsymbol{r}_{\alpha})
        z^i_{\alpha}
    S_{\alpha}
    +\\&+
    C_{2,1}
    \sum_{\alpha, i,j}
    J_{2,1}^{ij}(\boldsymbol{r}_{\alpha})
        z^i_{\alpha}
        z^j_{\alpha}
        (S_{\alpha})^2
    +\\&+
    C_{2,2}
    \sum_{\substack{\alpha, \beta, \nu, \\ i,j}}
    J_{2,2}^{ij}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^i_{\alpha}
        z^j_{\beta}
        S_{\alpha}
        S_{\beta}
    +\\&+
    C_{3, 1}
    \sum_{\substack{\alpha, \\ i, j, u}}
    J^{iju}_{3, 1}(\boldsymbol{r}_{\alpha})
        z^i_{\alpha}
        z^j_{\alpha}
        z^u_{\alpha}
        (S_{\alpha})^3
    +\\&+
    C_{3, 2}
    \sum_{\substack{\alpha, \beta, \nu, \\ i, j, u}}
    J^{iju}_{3, 2}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^i_{\alpha}
        z^j_{\alpha}
        z^u_{\beta}
        (S_{\alpha})^2
        S_{\beta}
    +\\&+
    C_{3, 3}
    \sum_{\substack{\alpha, \beta, \gamma, \\ \nu, \lambda, \\ i, j, u}}
    J^{iju}_{3, 3}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
        z^i_{\alpha}
        z^j_{\beta}
        z^u_{\gamma}
        S_{\alpha}
        S_{\beta}
        S_{\gamma}
    +\\&+
    C_{4, 1}
    \sum_{\substack{\alpha, \\ i, j, u, v}}
    J_{4, 1}^{ijuv}(\boldsymbol{r}_{\alpha})
        z^i_{\alpha}
        z^j_{\alpha}
        z^u_{\alpha}
        z^v_{\alpha}
        (S_{\alpha})^4
    +\\&+
    C_{4, 2, 1}
    \sum_{\substack{\alpha, \beta, \nu, \\ i, j, u, v}}
    J_{4, 2, 1}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^i_{\alpha}
        z^j_{\alpha}
        z^u_{\alpha}
        z^v_{\beta}
        (S_{\alpha})^3
        S_{\beta}
    +\\&+
    C_{4, 2, 2}
    \sum_{\substack{\alpha, \beta, \nu, \\ i, j, u, v}}
    J_{4, 2, 2}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^i_{\alpha}
        z^j_{\alpha}
        z^u_{\beta}
        z^v_{\beta}
        (S_{\alpha})^2
        (S_{\beta})^2
    +\\&+
    C_{4, 3}
    \sum_{\substack{\alpha, \beta, \gamma, \\ \nu, \lambda, \\ i, j, u, v}}
    J_{4, 3}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
        z^i_{\alpha}
        z^j_{\alpha}
        z^u_{\beta}
        z^v_{\gamma}
        (S_{\alpha})^2
        S_{\beta}
        S_{\gamma}
    +\\&+
    C_{4, 4}
    \sum_{\substack{\alpha, \beta, \gamma, \varepsilon, \nu, \lambda, \rho, \\ \\ i, j, u, v}}
    J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})
        z^i_{\alpha}
        z^j_{\beta}
        z^u_{\gamma}
        z^v_{\varepsilon}
        S_{\alpha}
        S_{\beta}
        S_{\gamma}
        S_{\varepsilon}


Renormalized parameters
=======================

.. math::
    J^{\prime, i}(\boldsymbol{r}_{\alpha})
    =&
    C_1
    J^i_1(\boldsymbol{r}_{\alpha})
    +\\&+
    2C_{2,1}
    \sum_{j}
    J^{ij}_{2,1}(\boldsymbol{r}_{\alpha})
        z^j_{\alpha}S_{\alpha}
    +\\&+
    2C_{2,2}
    \sum_{\beta, \nu, j}
    J^{ij}_{2,2}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^j_{\beta}S_{\beta}
    +\\&+
    3C_{3, 1}
    \sum_{j, u}
    J^{iju}_{3, 1}(\boldsymbol{r}_{\alpha})
        z^j_{\alpha}
        z^u_{\alpha}
        S_{\alpha}
        S_{\alpha}
    +\\&+
    3C_{3, 2}
    \sum_{\substack{\beta, \nu, \\ j, u}}
    J^{iju}_{3, 2}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^j_{\alpha}
        z^u_{\beta}
        S_{\alpha}
        S_{\beta}
    +\\&+
    3C_{3, 3}
    \sum_{\substack{\beta, \gamma, \\ \nu, \lambda, \\ j, u}}
    J^{iju}_{3, 3}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
        z^j_{\beta}
        z^u_{\gamma}
        S_{\beta}
        S_{\gamma}
    +\\&+
    4C_{4, 1}
    \sum_{\substack{j, u, v}}
    J_{4, 1}^{ijuv}(\boldsymbol{r}_{\alpha})
        z^j_{\alpha}
        z^u_{\alpha}
        z^v_{\alpha}
        S_{\alpha}
        S_{\alpha}
        S_{\alpha}
    +\\&+
    4C_{4, 2, 1}
    \sum_{\substack{\beta, \nu, \\ j, u, v}}
    J_{4, 2, 1}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^j_{\alpha}
        z^u_{\alpha}
        z^v_{\beta}
        S_{\alpha}
        S_{\alpha}
        S_{\beta}
    +\\&+
    4C_{4, 2, 2}
    \sum_{\substack{\beta, \nu, \\ j, u, v}}
    J_{4, 2, 2}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^j_{\alpha}
        z^u_{\beta}
        z^v_{\beta}
        S_{\alpha}
        S_{\beta}
        S_{\beta}
    +\\&+
    4C_{4, 3}
    \sum_{\substack{\beta, \gamma \\ \nu, \lambda, \\ j, u, v}}
    J_{4, 3}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
        z^j_{\alpha}
        z^u_{\beta}
        z^v_{\gamma}
        S_{\alpha}
        S_{\beta}
        S_{\gamma}
    +\\&+
    4C_{4, 4}
    \sum_{\substack{\beta, \gamma, \varepsilon, \\ \nu, \lambda, \rho, \\ j, u, v}}
    J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})
        z^j_{\beta}
        z^u_{\gamma}
        z^v_{\varepsilon}
        S_{\beta}
        S_{\gamma}
        S_{\varepsilon}


.. math::
    J^{\prime, ij}(\boldsymbol{r}_{\nu,\alpha\beta})
    =&
    C_{2, 2}
    J^{ij}_{2,2}(\boldsymbol{r}_{\nu,\alpha\beta})+\\&+
    \delta_{\alpha,\beta}
    \Biggl(
        2C_{2,1}
        J^{ij}_{2,1}(\boldsymbol{r}_{\alpha})
        +\\&\phantom{+\delta_{\alpha,\beta}\Biggl(}+
        3C_{3, 1}
        \sum_{u}
        J^{iju}_{3, 1}(\boldsymbol{r}_{\alpha})
            z^u_{\alpha}
            S_{\alpha}
        +\\&\phantom{+\delta_{\alpha,\beta}\Biggl(}+
        6C_{4, 1}
        \sum_{u, v}
        J_{4, 1}^{ijuv}(\boldsymbol{r}_{\alpha})
            z^u_{\alpha}
            z^v_{\alpha}
            S_{\alpha}
            S_{\alpha}
    \Biggr)
    +\\&+
    3C_{3, 2}
    \sum_{\nu, u}
    J^{iuj}_{3, 2}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^u_{\alpha}
        S_{\alpha}
    +\\&+
    3C_{3, 3}
    \sum_{\gamma, \lambda, u}
    J^{iju}_{3, 3}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
        z^u_{\gamma}
        S_{\gamma}
    +\\&+
    6C_{4, 2, 1}
    \sum_{u, v}
    J_{4, 2, 1}^{iuvj}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^u_{\alpha}
        z^v_{\alpha}
        S_{\alpha}
        S_{\alpha}
    +\\&+
    6C_{4, 2, 2}
    \sum_{u, v}
    J_{4, 2, 2}^{iujv}(\boldsymbol{r}_{\nu,\alpha\beta})
        z^u_{\alpha}
        z^v_{\beta}
        S_{\alpha}
        S_{\beta}
    +\\&+
    6C_{4, 3}
    \sum_{\substack{\gamma, \lambda, \\ u, v}}
    J_{4, 3}^{iujv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma})
        z^u_{\alpha}
        z^v_{\gamma}
        S_{\alpha}
        S_{\gamma}
    +\\&+
    6C_{4, 4}
    \sum_{\substack{\gamma, \varepsilon, \\ \lambda, \rho, \\ u, v}}
    J_{4, 4}^{ijuv}(\boldsymbol{r}_{\nu,\alpha\beta}, \boldsymbol{r}_{\lambda,\alpha\gamma}, \boldsymbol{r}_{\rho,\alpha\varepsilon})
        z^u_{\gamma}
        z^v_{\varepsilon}
        S_{\gamma}
        S_{\varepsilon}
