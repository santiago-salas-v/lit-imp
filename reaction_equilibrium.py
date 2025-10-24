# -*- coding: utf-8 -*-
from _functools import partial
from numerik import nr_ls
from numpy import array,concatenate,ones,zeros,log,power,diagflat,prod,eye,sqrt,multiply,empty,sign,exp

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
        t_abs,
        method,
        max_it,
        tol,
        method_loops,
        notify_status_func,
        process_func_handle):
    """Newton method for non-linear algebraic system, with line-search
    :return: tuple with neq, xieq, f_0
    :param n0: array (n x 1) - mol - alimentación
    :param mm: array (n x 1) - masa molar
    :param z: array (n x 1) - carga - alimentación
    :param s_index: int - índice de solvente
    :param kc: array (nr x 1) - "Cte." de equilibrio en reacción j @ T
    :param nu_ij: array (n x nr) - Coefs. esteq. componente i en reacción j
    :param t_abs: float - temperatura T en Kelvin
    :param neq_0: array (n x 1) - mol - estimado inicial, equilibrio
    :param xieq_0: array (nr x 1) - avance de reacción j - estimado inicial, equilibrio
    :param method: str - método en uso: 'ideal_solution', 'davies', 'debye-hueckel'
    :param max_it: int - máximo de iteraciones
    :param tol: float - error, tolerancia máxima
    :param method_loops: list (1 x 2): ciclos completados [backtrack, totales]
    :param notify_status_func: función a llamar cada avance de iteración
    :param process_func_handle: función a llamar para cancelación
    """

    n = len(n0)
    nr = nu_ij.shape[1]
    mm_0 = mm[s_index].item()
    a_m = a_m_d_h(t_abs)

    meq_0 = neq_0 / (mm_0 * neq_0[s_index])
    ionic_str_eq_0 = 1 / 2.0 * power(z, 2).T.dot(meq_0)
    gammaeq = ones(n)
    gammaeq_0 = gammaeq
    neq = neq_0
    meq = meq_0
    xieq = xieq_0
    ionic_str_eq = ionic_str_eq_0
    m0_ref = 1 / 1000.0  # ref. 1mol/kgsolvent conv. to mol/gsolvent
    ionic_str_eq_adim = ionic_str_eq / m0_ref
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
    x0 = concatenate([neq_0, xieq_0])

    component_order = [s_index]
    component_order.extend(
        [index for (index, x) in enumerate(meq_0) if index != s_index]
    )
    component_order_map = \
        sorted([(origindex, orderedindex) for (orderedindex, origindex)
                in enumerate(component_order)], key=lambda y: y[0])
    return_to_original_indexes = [x[1] for x in component_order_map]

    if method == 'ideal_solution':
        # case already contemplated above
        pass
    elif method == 'davies' or method == 'debye-hueckel':
        # x is [n_0, m_1, m_2, ..., m_n, xi_1, xi_2, ..., xi_{nr}, \gamma_1,
        #       \gamma_2, ..., \gamma_n, I]
        ordered_meq_0 = array(
            [meq_0[index].item() for index in component_order]
        ).T
        ordered_gammaeq_0 = array(
            [gammaeq_0[index].item() for index in component_order]
        ).T
        ordered_nu_ij = concatenate(
            [nu_ij[index] for index in component_order]
        )
        ordered_n0 = array(
            [n0[index].item() for index in component_order]
        ).T
        ordered_z = array(
            [z[index].item() for index in component_order]
        ).T
        ordered_neq_0 = array(
            [neq_0[index].item() for index in component_order]
        ).T
        if method == 'davies':
            ordered_gammaeq_0[0] = gamma_solvent_id(mm_0, ordered_meq_0[1:])
            ordered_gammaeq_0[1:] = multiply(
                gamma_davies(ordered_z[1:], ionic_str_eq_adim, a_m),
                gamma_setchenow(ordered_z[1:], ionic_str_eq_adim, 0.1))
            f = partial(
                f_gl_0_davies,
                n0=ordered_n0,
                nu_ij=ordered_nu_ij,
                n=n,
                nr=nr,
                kc=kc,
                z=ordered_z,
                mm_0=mm_0,
                a_m=a_m)
            j = partial(
                jac_davies,
                n0=ordered_n0,
                nu_ij=ordered_nu_ij,
                n=n,
                nr=nr,
                kc=kc,
                z=ordered_z,
                mm_0=mm_0,
                a_m=a_m)
        elif method == 'debye-hueckel':
            ordered_gammaeq_0[0] = gamma_solvent_id(mm_0, ordered_meq_0[1:])
            ordered_gammaeq_0[1:] = multiply(
                gamma_d_h(ordered_z[1:], ionic_str_eq_adim, a_m),
                gamma_setchenow(ordered_z[1:], ionic_str_eq_adim, 0.1))
            f = partial(
                f_gl_0_d_h,
                n0=ordered_n0,
                nu_ij=ordered_nu_ij,
                n=n,
                nr=nr,
                kc=kc,
                z=ordered_z,
                mm_0=mm_0,
                a_m=a_m)
            j = partial(
                jac_d_h,
                n0=ordered_n0,
                nu_ij=ordered_nu_ij,
                n=n,
                nr=nr,
                kc=kc,
                z=ordered_z,
                mm_0=mm_0,
                a_m=a_m)
        x0 = concatenate(
            [
                ordered_neq_0[0],
                ordered_meq_0[1:],
                xieq_0,
                ordered_gammaeq_0,
                ionic_str_eq_0
            ])
    progress_k, stop, outer_it_k, outer_it_j, \
        lambda_ls, accum_step, x, \
        diff, f_val, lambda_ls_y, \
        method_loops = \
        nr_ls(x0=x0,
              f=f,
              j=j,
              tol=tol,
              max_it=max_it,
              inner_loop_condition=lambda x_vec:
              all([item >= 0 for item in
                   x_vec[0:n]]),
              notify_status_func=notify_status_func,
              method_loops=method_loops,
              process_func_handle=process_func_handle)

    if method == 'ideal_solution':
        neq = x[0:n]
        xieq = x[n:n + nr]
        meq = neq / (neq[s_index] * mm_0)
        gammaeq = gammaeq_0
        ionic_str_eq = 1 / 2.0 * power(z, 2).T.dot(meq)
    elif method == 'davies' or method == 'debye-hueckel':
        neq0 = x[0:n][0]
        meq[1:n] = x[0:n][1:n]
        meq[0] = 1 / mm_0
        neq = meq * neq0 * mm_0
        # Reorder to output original order
        meq = meq[return_to_original_indexes]
        neq = neq[return_to_original_indexes]
        xieq = x[n:n + nr]
        gammaeq = x[n + nr:n + nr + n][return_to_original_indexes]
        ionic_str_eq = x[n + nr + n]
    return neq, meq, xieq, gammaeq, ionic_str_eq, method_loops


def f_gl_0_ideal(x, n0, nu_ij, n, nr, kc, mm_0, s_index):
    neq = x[0:n]
    n0_mm0 = neq[s_index] * mm_0
    m0_ref = 1 / 1000.0  # ref. 1mol/kgsolvent conv. to mol/gsolvent
    meq = neq / n0_mm0  # mol/gsolvent
    xieq = x[n:n + nr]
    result = empty([n + nr], dtype=float)
    result[0:n] = -neq + n0.flatten() + nu_ij.dot(xieq)
    result[n:n + nr] = -kc + prod(array([(meq / m0_ref)**nu_ij[:,j] for j in range(nr)]),axis=1).T
    return result


def jac_ideal(x, n0, nu_ij, n, nr, kc, mm_0, s_index):
    neq = x[0:n]
    n0_mm0 = neq[s_index] * mm_0
    m0_ref = 1 / 1000.0  # ref. 1mol/kgsolvent conv. to mol/gsolvent
    meq = neq / n0_mm0
    diag_1_ov_meq = diagflat(power(meq, -1), 0)
    diag_quotient = diagflat(prod(array([meq**nu_ij[:,j] for j in range(nr)]),axis=1))
    result = array(zeros([n + nr, n + nr], dtype=float))
    result[0:n, 0:n] = -1 * eye(n).astype(float)
    result[0:n, n:n + nr] = nu_ij
    # Return Jacobian terms as n calculated from molality (m)
    result[n:n + nr, 0:n] = \
        diag_quotient.dot(nu_ij.T.dot(diag_1_ov_meq)) * 1 / n0_mm0
    return result


def f_gl_0_davies(x, n0, nu_ij, n, nr, kc, z, mm_0, a_m):
    # f(x) = 0, objective function set for Davies model.
    neq = zeros(n)
    meq = zeros(n)
    # x is [n0 m1 m2 ... m_n xi1 xi2 ... xi_nr gamma1 gamma2 ... gamma_n
    # ionic_str]
    neq[0] = x[0]
    meq[1:n] = x[1:n]
    xieq = x[n:n + nr]
    gammaeq = x[n + nr:n + nr + n]
    ionic_str = x[n + nr + n]

    # calculate neq for all components
    n0_mm0 = neq[0] * mm_0
    meq[0] = 1 / mm_0
    neq = meq * n0_mm0
    m0_ref = 1 / 1000.0  # ref. 1mol/kgsolvent conv. to mol/gsolvent
    ionic_str_adim = ionic_str / m0_ref

    result = array(empty(n + nr + n + 1, dtype=float))
    result[0:n] = -neq + n0 + nu_ij * xieq
    result[n:n + nr] = -kc + multiply(
        prod(power(meq / m0_ref, nu_ij), 0).T,
        prod(power(gammaeq, nu_ij), 0).T
    )
    result[n + nr] = \
        -gammaeq[0] + gamma_solvent_id(mm_0, meq[1:n])
    result[n + nr + 1:n + nr + n] = \
        - gammaeq[1:] + \
        + multiply(
            gamma_davies(z[1:], ionic_str_adim, a_m),
            gamma_setchenow(z[1:], ionic_str_adim, 0.1))
    result[n + nr + n] = \
        -ionic_str + 1 / 2.0 * power(z, 2).T.dot(meq)
    return result


def jac_davies(x, n0, nu_ij, n, nr, kc, z, mm_0, a_m):
    # j(x), Jacobian matrix for Davies model.
    neq = zeros(n)
    meq = zeros(n)
    # x is [n0 m1 m2 ... m_n xi1 xi2 ... xi_nr gamma1 gamma2 ... gamma_n
    # ionic_str]
    neq[0] = x[0].item()
    meq[1:n] = x[1:n]
    xieq = x[n:n + nr]
    gammaeq = x[n + nr:n + nr + n]
    ionic_str = x[n + nr + n].item()

    # calculate neq for all components
    n0_mm0 = (neq[0] * mm_0).item()
    meq[0] = 1 / mm_0
    neq = meq * n0_mm0
    m0_ref = 1 / 1000.0  # ref. 1mol/kgsolvent conv. to mol/gsolvent
    ionic_str_adim = ionic_str / m0_ref
    sqrt_ionic_str_adim = sqrt(ionic_str_adim)

    diag_quotient = diagflat(
        prod(
            power(
                multiply(gammaeq, meq / m0_ref),
                nu_ij),
            0)
    )
    result = array(
        zeros([n + nr + n + 1, n + nr + n + 1], dtype=float)
    )
    result[0:n, 0:n] = \
        diagflat(
            concatenate(
                [-1.0 * array([1]),
                 -n0_mm0 * ones(n - 1)]
            )
    )
    result[1:n] = -meq[1:] * mm_0
    result[0:n, n:n + nr] = nu_ij
    result[n + nr:n + nr + n + 1, n + nr: n + nr + n + 1] = \
        -1.0 * eye(n + 1)
    result[n:n + nr, 0:n] = \
        diag_quotient * nu_ij.T * diagflat(
            concatenate(
                [array(0.0), 1 / meq[1:]]
            )
    )
    result[n:n + nr, n + nr:n + nr + n] = \
        diag_quotient * nu_ij.T * diagflat(
            1 / gammaeq
    )
    gamma0_ov_phi = exp(-1.0 * mm_0 * sum(meq[1:])).item()
    result[n + nr, 1:n] = -1.0 * mm_0 * gamma0_ov_phi
    factor_1 = \
        sqrt_ionic_str_adim / (1 + sqrt_ionic_str_adim) \
        - 0.3 * ionic_str_adim
    dfactor_1_di = \
        (1 / m0_ref) * (-0.3 + 1 / (2 * sqrt_ionic_str_adim *
                                    (1 + sqrt_ionic_str_adim) ** 2))
    factor_2 = power(10,
                        -a_m *
                        power(z[1:], 2) *
                        factor_1 +
                        (1 -
                         power(sign(z[1:]), 2)) *
                        0.1 *
                        ionic_str_adim)
    result[n + nr + 1:n + nr + n, n + nr + n] = \
        multiply(
            log(10.0) * (
                -a_m * power(z[1:], 2) *
                dfactor_1_di +
                (1 - power(sign(z[1:]), 2)) * 0.1 / m0_ref),
            factor_2)
    result[n + nr + n, 1:n] = \
        1 / 2.0 * power(z[1:].T, 2.0)
    return result


def f_gl_0_d_h(x, n0, nu_ij, n, nr, kc, z, mm_0, a_m):
    # f(x) = 0, objective function set for Debye-Hueckel model.
    neq = zeros(n)
    meq = zeros(n)
    # x is [n0 m1 m2 ... m_n xi1 xi2 ... xi_nr gamma1 gamma2 ... gamma_n
    # ionic_str]
    neq[0] = x[0]
    meq[1:n] = x[1:n]
    xieq = x[n:n + nr]
    gammaeq = x[n + nr:n + nr + n]
    ionic_str = x[n + nr + n]

    # calculate neq for all components
    n0_mm0 = neq[0] * mm_0
    meq[0] = 1 / mm_0
    neq = meq * n0_mm0
    m0_ref = 1 / 1000.0  # ref. 1mol/kgsolvent conv. to mol/gsolvent
    ionic_str_adim = ionic_str / m0_ref

    result = array(empty(n + nr + n + 1, dtype=float))
    result[0:n] = -neq + n0 + nu_ij * xieq
    result[n:n + nr] = -kc + multiply(
        prod(power(meq / m0_ref, nu_ij), 0).T,
        prod(power(gammaeq, nu_ij), 0).T
    )
    result[n + nr] = \
        -gammaeq[0] + gamma_solvent_id(mm_0, meq[1:n])
    result[n + nr + 1:n + nr + n] = \
        - gammaeq[1:] + \
        + multiply(
            gamma_d_h(z[1:], ionic_str_adim, a_m),
            gamma_setchenow(z[1:], ionic_str_adim, 0.1))
    result[n + nr + n] = \
        -ionic_str + 1 / 2.0 * power(z, 2).T.dot(meq)
    return result


def jac_d_h(x, n0, nu_ij, n, nr, kc, z, mm_0, a_m):
    # j(x), Jacobian matrix for Debye-Hueckel model.
    neq = zeros(n)
    meq = zeros(n)
    # x is [n0 m1 m2 ... m_n xi1 xi2 ... xi_nr gamma1 gamma2 ... gamma_n
    # ionic_str]
    neq[0] = x[0].item()
    meq[1:n] = x[1:n]
    xieq = x[n:n + nr]
    gammaeq = x[n + nr:n + nr + n]
    ionic_str = x[n + nr + n].item()

    # calculate neq for all components
    n0_mm0 = (neq[0] * mm_0).item()
    meq[0] = 1 / mm_0
    neq = meq * n0_mm0
    m0_ref = 1 / 1000.0  # ref. 1mol/kgsolvent conv. to mol/gsolvent
    ionic_str_adim = ionic_str / m0_ref
    sqrt_ionic_str_adim = sqrt(ionic_str_adim)

    diag_quotient = diagflat(
        prod(
            power(
                multiply(gammaeq, meq / m0_ref),
                nu_ij),
            0)
    )
    result = array(
        zeros([n + nr + n + 1, n + nr + n + 1], dtype=float)
    )
    result[0:n, 0:n] = \
        diagflat(
            concatenate(
                [-1.0 * array([1]),
                 -n0_mm0 * ones(n - 1)]
            )
    )
    result[1:n] = -meq[1:] * mm_0
    result[0:n, n:n + nr] = nu_ij
    result[n + nr:n + nr + n + 1, n + nr: n + nr + n + 1] = \
        -1.0 * eye(n + 1)
    result[n:n + nr, 0:n] = \
        diag_quotient * nu_ij.T * diagflat(
            concatenate(
                [array(0.0), 1 / meq[1:]]
            )
    )
    result[n:n + nr, n + nr:n + nr + n] = \
        diag_quotient * nu_ij.T * diagflat(
            1 / gammaeq
    )
    gamma0_ov_phi = exp(-1.0 * mm_0 * sum(meq[1:])).item()
    result[n + nr, 1:n] = -1.0 * mm_0 * gamma0_ov_phi
    factor_1 = \
        sqrt_ionic_str_adim
    dfactor_1_di = \
        (1 / m0_ref) * (1 / (2 * sqrt_ionic_str_adim))
    factor_2 = power(10,
                        -a_m *
                        power(z[1:], 2) *
                        factor_1 +
                        (1 -
                         power(sign(z[1:]), 2)) *
                        0.1 *
                        ionic_str_adim)
    result[n + nr + 1:n + nr + n, n + nr + n] = \
        multiply(
            log(10.0) * (
                -a_m * power(z[1:], 2) *
                dfactor_1_di +
                (1 - power(sign(z[1:]), 2)) * 0.1 / m0_ref),
            factor_2)
    result[n + nr + n, 1:n] = \
        1 / 2.0 * power(z[1:].T, 2.0)
    return result


def gamma_davies(z, i, a_m):
    # Activity coefficient, Davies model
    sqrt_i = sqrt(i)
    log_gamma = -a_m * power(z, 2) * (sqrt_i / (1 + sqrt_i) - 0.3 * i)
    return power(10, log_gamma)


def gamma_d_h(z, i, a_m):
    # Activity coefficient, Debye-Hueckel model
    sqrt_i = sqrt(i)
    log_gamma = -a_m * power(z, 2) * sqrt_i
    return power(10, log_gamma)


def gamma_solvent_id(mm_0, m):
    phi = 1.0  # ideal
    ln_gamma = - phi * mm_0 * sum(m)
    return exp(ln_gamma)


def gamma_setchenow(z, i, b):
    uncharged_ones = 1 - power(sign(z), 2)
    return power(10, uncharged_ones * b * i)


def a_m_d_h(t_abs=298.15):
    # Parameter a of Debye-Hueckel theory, at T
    epsilon_0 = 8.85418781762e-12  # C^2 * N^-1 * m^-2
    epsilon_r = epsilon_r_water(t_abs)  # C^2 * N^-1 * m^-2
    pi = 3.14159265359
    n_a = 6.02214129e+23  # mol^-1
    e = 1.602176565e-19  # C
    k_b = 1.3806488e-23  # J K^-1
    rho_1 = 0.99714  # kg * L^-1
    a_c = (e ** 2 / (4 * epsilon_r * epsilon_0 * k_b * t_abs)) ** (3 / 2.0) \
        * (2 * n_a / pi ** 2) ** (1 / 2.0) * 1000 ** (1 / 2.0) \
        / log(10.0)  # (mol/L)^(-1/2)
    a_m = a_c * rho_1 ** (1 / 2.0)  # (mol/kg)^(-1/2)
    return a_m


def epsilon_r_water(t_abs=298.15):
    # Static relative permittivity (epsilon_r) of water as func. of T
    # CRC Handbook of Chemistry and Physics, 90th ed., CRC Press, 2009
    a = +0.24921E+03
    b = -0.79069E+00
    c = +0.72997E-03
    d = +0.00000E+00
    return a + b * t_abs + c * t_abs ** 2 + d * t_abs ** 3
