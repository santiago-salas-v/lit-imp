# -*- coding: utf-8 -*-
from _functools import partial
from numerik import nr_ls
import numpy as np


# noinspection PyAugmentAssignment
def calc_xieq(
        n0,
        mm,
        z,
        s_index,
        kc,
        nu_ij,
        neq_0,
        xieq_0,
        method,
        max_it,
        tol,
        method_loops,
        notify_status_func,
        process_func_handle):
    """Newton method for non-linear algebraic system, with line-search
    :return: tuple with neq, xieq, f_0
    :param n0: np.matrix (n x 1) - mol - alimentación
    :param mm: np.matrix (n x 1) - masa molar
    :param z: np.matrix (n x 1) - carga - alimentación
    :param s_index: int - índice de solvente
    :param kc: np.matrix (nr x 1) - "Cte." de equilibrio en reacción j @ T
    :param nu_ij: np.matrix (n x nr) - Coefs. esteq. componente i en reacción j
    :param neq_0: np.matrix (n x 1) - mol - estimado inicial, equilibrio
    :param xieq_0: np.matrix (nr x 1) - avance de reacción j - estimado inicial, equilibrio
    :param method: str - método en uso: 'ideal_solution', 'davies', 'debye-hueckel'
    :param max_it: int - máximo de iteraciones
    :param tol: float - error, tolerancia máxima
    :param method_loops: list (1 x 2): ciclos completados [backtrack, totales]
    :param notify_status_func: función a llamar cada avance de iteración
    :param process_func_handle: función a llamar para cancelación
    """

    n = len(n0)
    nr = nu_ij.shape[1]
    mm_0 = mm[s_index]

    meq_0 = neq_0 / (mm_0 * neq_0[s_index])
    gammaeq_0 = np.matrix(np.ones([n, 1]))
    ionic_str_eq_0 = 1 / 2.0 * np.power(z, 2).T * meq_0
    gammaeq = np.ones([n, 1])
    neq = neq_0
    meq = meq_0
    xieq = xieq_0
    ionic_str_eq = ionic_str_eq_0
    f = partial(
        f_gl_0_ideal,
        n0=n0,
        nu_ij=nu_ij,
        n=n,
        nr=nr,
        kc=kc,
        mm_0=mm_0,
        s_index=s_index)
    j = partial(
        jac_ideal,
        n0=n0,
        nu_ij=nu_ij,
        n=n,
        nr=nr,
        kc=kc,
        mm_0=mm_0,
        s_index=s_index)
    # x is [n_0, n_1, n_2, ..., n_n, xi_1, xi_2, xi_3, ..., xi_{nr}]
    x0 = np.concatenate([neq_0, xieq_0])

    if method == 'ideal_solution':
        # case already contemplated above
        pass
    elif method == 'davies':
        f = partial(
            f_gl_0_ideal,
            n0=n0,
            nu_ij=nu_ij,
            n=n,
            nr=nr,
            kc=kc,
            mm_0=mm_0,
            s_index=s_index)
        j = partial(
            jac_ideal,
            n0=n0,
            nu_ij=nu_ij,
            n=n,
            nr=nr,
            kc=kc,
            mm_0=mm_0,
            s_index=s_index)
        # x is [n_0, m_1, m_2, ..., m_n, xi_1, xi_2, ..., xi_{nr}, \gamma_1,
        #       \gamma_2, ..., \gamma_n, I]
        x0 = np.concatenate(
            [
                np.matrix(neq_0[s_index].item()),
                np.matrix([x.item() for (index, x) in enumerate(meq_0)
                           if index != s_index]).T,
                xieq_0,
                gammaeq_0,
                ionic_str_eq_0
            ])

    progress_k, stop, outer_it_k, outer_it_j, \
        lambda_ls, accum_step, x, \
        diff, f_val, lambda_ls_y, \
        method_loops =    \
        nr_ls(x0=x0,
              f=f,
              j=j,
              tol=tol,
              max_it=max_it,
              inner_loop_condition=lambda x_vec:
              all([item >= 0 for item in x_vec[0:n]]),
              notify_status_func=notify_status_func,
              method_loops=method_loops,
              process_func_handle=process_func_handle)

    if method == 'ideal_solution':
        neq = x[0:n]
        xieq = x[n:n + nr]
        meq = neq / (neq[s_index] * mm_0)
        gammaeq = gammaeq_0
        ionic_str_eq = 1 / 2.0 * np.power(z, 2).T * meq
    elif method == 'davies':
        neq = x[0:n]
        xieq = x[n:n + nr]
    return neq, meq, xieq, gammaeq, ionic_str_eq, method_loops


def f_gl_0_ideal(x, n0, nu_ij, n, nr, kc, mm_0, s_index):
    neq = x[0:n, 0]
    n0_mm0 = neq[s_index] * mm_0
    meq = neq / n0_mm0
    xieq = x[n:n + nr, 0]
    result = np.matrix(np.empty([n + nr, 1], dtype=float))
    result[0:n] = -neq + n0 + nu_ij * xieq
    result[n:n + nr] = -kc + np.prod(np.power(meq, nu_ij), 0).T
    return result


def jac_ideal(x, n0, nu_ij, n, nr, kc, mm_0, s_index):
    neq = x[0:n, 0]
    n0_mm0 = neq[s_index] * mm_0
    meq = neq / n0_mm0
    eins_durch_m = np.diag(np.power(meq, -1).A1, 0)
    quotient = np.diagflat(np.prod(np.power(meq, nu_ij), 0))
    result = np.matrix(np.zeros([n + nr, n + nr], dtype=float))
    result[0:n, 0:n] = -1 * np.eye(n).astype(float)
    result[0:n, n:n + nr] = nu_ij
    # Return Jacobian terms as n calculated from molality (m)
    result[n:n + nr, 0:n] = quotient * nu_ij.T * eins_durch_m * 1 / n0_mm0
    return result
