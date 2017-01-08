import numpy as np
from scipy import linalg
import matplotlib
from numerik import nr_ls
import logging
<<<<<<< HEAD
from sympy import solve, nsolve, symbols, pprint
from mpmath import nstr
from _functools import partial
=======
>>>>>>> cba2043dcfd5d7de9f6eccd5102078c2924e7483
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'
# noinspection PyPep8
from matplotlib import pyplot as plt

eps = np.finfo(float).eps
logger = logging.getLogger()
fhandler = logging.FileHandler(filename='./logs/CO2_Equilibria.log')
formatter = logging.Formatter('%(asctime)s;%(message)s')
fhandler.setFormatter(formatter)
logger.addHandler(fhandler)
logger.setLevel(logging.DEBUG)


def eq_set(x, c0, p0co2):
    # x,    [c_h2o, c_h30, c_ho, c_hco3, c_na, c_co3,
    #       c_co2, c_h2co3, pco2,
    #       xi0, xi1, xi2, xi3, xi4]
    c = x[:8]
    xco2 = c[6] / sum(c)
    pco2 = x[8]
    xi = x[9:]
    return np.array([
        1670.0 * 100 * xco2 - pco2,  # Hco2 = pco2/xco2
        10 ** -2.821023053 * c[6] - c[7],  # 10^-pK2 = ch2co3/cco2
        10 ** -3.539912399 * c[7] - c[1] * c[3],  # 10^-pK3 = chco3*ch3o/ch2co3
        10 ** -10.32991986 * c[3] - c[1] * c[5],  # 10^-pK4 = cco3*ch3o/chco3
        10 ** -13.99602524 - c[1] * c[2],  # 10^-pKw = ch3o*cho
        c[0] - c0[0] - (- xi[0] - xi[1] - xi[2] - xi[3] - 2 * xi[4]),
        c[1] - c0[1] - (+ xi[2] + xi[3] + xi[4]),
        c[2] - c0[2] - (+ xi[4]),
        c[3] - c0[3] - (+ xi[2] - xi[3]),
        c[4] - c0[4] - (+ 0),
        c[5] - c0[5] - (+ xi[3]),
        c[6] - c0[6] - (- xi[0] - xi[1]),
        c[7] - c0[7] - (+ xi[1] - xi[2]),
        pco2 - p0co2 - (+ xi[0]) * 8.314 * 298.15
    ]).T


def jac_eq_set(x):
    j = np.zeros([len(x), len(x)])

    c = x[:8]
    # xco2 = c[6] / sum(c)
    # pco2 = x[8]
    # xi = x[9:]

    d_f1_dc = np.zeros([9, ])

    d_f1_dc[0:6] = 1670.0 * 100 * (-1.0 / sum(c) ** 2) * c[6]
    d_f1_dc[6] = 1670.0 * 100 * sum([c_ for i, c_ in enumerate(c) if i != 6]) / sum(c) ** 2
    d_f1_dc[7] = 1670.0 * 100 * (-1.0 / sum(c) ** 2) * c[6]
    d_f1_dc[8] = -1.0

    j[0, :8 + 1] = d_f1_dc

    j[1, 6], j[1, 7] = 10 ** -2.821023053, -1
    j[2, 1], j[2, 3], j[2, 7] = -c[3], -c[1], 10 ** -3.539912399
    j[3, 1], j[3, 3], j[3, 5] = -c[5], 10 ** -10.32991986, -c[1]
    j[4, 1], j[4, 2] = -c[2], -c[1]
    for ind in range(8 + 1):
        j[ind + 5, ind] = +1
    j[5, 9:14] = np.array([+1, +1, +1, +1, +2])
    j[6, 11:14] = np.array([-1, -1, -1])
    j[7, 13] = -1
    j[8, 11:13] = np.array([-1, +1])
    j[10, 12] = -1
    j[11, 9:11] = +1
    j[12, 10:12] = np.array([-1, +1])
    j[13, 9] = -1 * 8.314 * 298.15

    return j


def eq_set_const_p(x, c0, p0co2, p0n2, p0o2):
    # x, [nh2o, nh30, nho, nhco3, nna, nco3, nco2, nh2co3, pco2,
    #     xi1, xi2, xi3, xi4, xi5]
    # c = x[:8]
    # xco2 = c[6] / sum(c)
    pco2 = x[8]
    xi = x[9:]
    f = eq_set(x, c0, p0co2)
    if f.ndim != 1:
        f = f.flatten()
    f[13] = pco2 - p0co2 * (101.325 - pco2) / (p0n2 + p0o2) - \
        (+ xi[0]) * 8.314 * 298.15 * \
        (101.325 - pco2) / (p0n2 + p0o2)
    return f


def jac_eq_set_const_p(x, p0co2, p0n2, p0o2):
    j = jac_eq_set(x)

    # c = x[:8]
    # xco2 = c[6] / sum(c)
    pco2 = x[8]
    xi = x[9:]
    j[13, 8] = +p0co2 / (p0n2 + p0o2) + 1 + \
        xi[0] * 8.314 * 298.15 / (p0n2 + p0o2)
    j[13, 9] = -1 * 8.314 * 298.15 * \
        (101.325 - pco2) / (p0n2 + p0o2)

    return j


def notify_status_func(progress_k, stop_value, k,
                       j_it_backtrack, lambda_ls, accum_step,
                       x, diff, f_val, j_val, lambda_ls_y,
                       method_loops):
    g_min = np.nan
    g1 = np.nan
    y = lambda_ls_y
    pr_str = 'progress=' + str(progress_k) + \
             ';k=' + str(k) + \
             ';backtrack=' + str(j_it_backtrack) + \
             ';lambda_ls=' + str(lambda_ls) + \
             ';method loops=' + str(method_loops) + \
             ';accum_step=' + str(accum_step) + \
             ';stop=' + str(stop_value) + \
             ';X=' + '[' + ','.join(map(str, x.T.A1)) + ']' + \
             ';||X(k)-X(k-1)||=' + str((diff.T * diff).item()) + \
             ';f(X)=' + '[' + ','.join(map(str, f_val.T.A1)) + ']' + \
             ';||f(X)||=' + str(np.sqrt((f_val.T * f_val).item())) + \
             ';j(X)=' + str(j_val.tolist()) + \
             ';Y=' + '[' + ','.join(map(str, y.T.A1)) + ']' + \
             ';||Y||=' + str(np.sqrt((y.T * y).item())) + \
             ';g=' + str(g_min) + \
             ';|g-g1|=' + str(abs(g_min - g1))
    logging.debug(pr_str)


def print_variables_vector(x):
    for index, num in enumerate(x):
        if index < 8:
            print 'C' + str(index) + '=' + '%.20e' % num + ','
        elif index == 8:
            print 'pco2' + '=' + '%.20e' % num + ','
        else:
            print 'x' + str(index - 9) + '=' + '%.20e' % num + ','
    print 'pH: ' + '%g' % -np.log10(x[1])
    print 'pOH: ' + '%g' % -np.log10(x[2])
    print 'pH + pOH: ' + '%g' % sum(-np.log10(x[1:2 + 1]))


def main(xw0nahco3=0.07318, ph0=13.99602524 / 2, x0=None):
    comps = np.array([
        'H2O', 'H3O(+)', 'HO(-)', 'HCO3(-)', 'Na(+)',
        'CO3(2-)', 'CO2', 'H2CO3', 'CO2'
    ])

    mm = np.array([
        18, 19, 17, 61, 23, 60, 44, 62, 44
    ], dtype=float)

    pkw = 13.99602524
    mm0 = mm[0]
    rho0 = 1.0
    # ph0 = pkw / 2
    # xw0nahco3 = 0.07318
    p0n2 = 78.12 / (78.12 + 20.96) * 101.325
    p0o2 = 20.96 / (78.12 + 20.96) * 101.325

    x = np.ones_like(mm) * eps
    x[3] = xw0nahco3 / (mm[3] + mm[4]) / \
        (xw0nahco3 / (mm[3] + mm[4]) + xw0nahco3 /
         (mm[3] + mm[4]) + (1 - xw0nahco3) / mm[0])
    x[4] = xw0nahco3 / (mm[3] + mm[4]) / \
        (xw0nahco3 / (mm[3] + mm[4]) + xw0nahco3 /
         (mm[3] + mm[4]) + (1 - xw0nahco3) / mm[0])

    a = np.array([[1 + mm0 / rho0 * 10 ** (-ph0) / 1000.0,
                   mm0 / rho0 * 10 ** (-ph0) / 1000.0],
                  [mm0 / rho0 * 10 ** (ph0 - pkw) / 1000.0,
                   1 + mm0 / rho0 * 10 ** (ph0 - pkw) / 1000.0]])
    b = np.array([
        [mm0 / rho0 * 10 ** (-ph0) / 1000.0 * (1 - x[3] - x[4])],
        [mm0 / rho0 * 10 ** (ph0 - pkw) / 1000.0 * (1 - x[3] - x[4])]
    ])

    x[1:2 + 1] = linalg.inv(a).dot(b).flatten()
    x[0] = 1 - sum(x[1:])

    print str(comps)
    print 'mole frac.'
    print x
    print 'mass frac.'
    print np.multiply(x, mm) / sum(np.multiply(x, mm))

    c0 = np.array(np.ones([len(mm) - 1, ]) * eps)
    xi = np.array(np.ones(5) * eps)
    c0 = np.array(x[:8] / x[0] * rho0 / mm0 * 1000.0)
    p0co2 = np.array(10 ** -3.5 * 101.325)  # 10^-3.5atm = 10^-3.5*101.325 kPa
    c0[6] = p0co2 / (1670.0 * 100) * (sum(c0) - c0[6]) / \
        (1 - p0co2 / (1670.0 * 100))

    if x0 is None:
        x0 = np.append(np.append(c0, p0co2), xi)
    x = x0

    print 'x0:'
    print_variables_vector(x0)
    print 's(c_i*z_i)=' + '%g' % (x[1] - x[2] - x[3] + x[4] - 2 * x[5])

    print 'Constant V'

    progress_k, stop, outer_it_k, outer_it_j, \
        lambda_ls, accum_step, x, \
        diff, f_val, lambda_ls_y, \
        method_loops = \
        nr_ls(x0=np.matrix(x0).T,
              f=lambda x_v: np.matrix(eq_set(x_v, c0, p0co2)).T,
              j=lambda x_v: np.matrix(jac_eq_set(x_v)),
              tol=1e-14,
              max_it=1000,
              inner_loop_condition=lambda x_vec:
              all([item >= 0 for item in
                   x_vec[0:9]]),
              notify_status_func=notify_status_func,
              method_loops=[0, 0],
              process_func_handle=None)

    x = x.A1
    f_val = f_val.A1

    print '||f||: ' + str(np.sqrt(f_val.T.dot(f_val)))
    print 'x:'
    print_variables_vector(x)
    print 's(c_i*z_i)=' + '%g' % (x[1] - x[2] - x[3] + x[4] - 2 * x[5])

    print '\n'
    print 'Constant P'

    progress_k, stop, outer_it_k, outer_it_j, \
        lambda_ls, accum_step, x, \
        diff, f_val, lambda_ls_y, \
        method_loops = \
        nr_ls(x0=np.matrix(x0).T,
              f=lambda x_v:
              np.matrix(eq_set_const_p(x_v, c0, p0co2, p0n2, p0o2)).T,
              j=lambda x_v:
              np.matrix(jac_eq_set_const_p(x_v, p0co2, p0n2, p0o2)),
              tol=1e-14,
              max_it=1000,
              inner_loop_condition=lambda x_vec:
              all([item >= 0 for item in
                   x_vec[0:9]]),
              notify_status_func=notify_status_func,
              method_loops=[0, 0],
              process_func_handle=None)

    x = x.A1
    f_val = f_val.A1

    print '||f||: ' + str(np.sqrt(f_val.T.dot(f_val)))
    print 'x:'
    print_variables_vector(x)
    print 's(c_i*z_i)=' + '%g' % (x[1] - x[2] - x[3] + x[4] - 2 * x[5])
    print 'Vl/Vg0=' + '%f' % 1.0
    print 'Vl/Vg=' + '%f' % (1.0 * ((101.325 - x[8]) / (p0n2 + p0o2)))

    return x


if __name__ == '__main__':
    x0 = main(xw0nahco3=0.001)
    for xw in [0.07308, 0.07309]:
        print 'xw0nahco3 = ' + '%g' % xw
        x0 = main(xw0nahco3=xw, x0=x0)
    # iterate for xwnahco3 vs. P curve, approaching by
    # valid P < 101.325kPa
    xw = 0.0730905  # beyond sat. (?)
    it = 0
    min_pco2 = 90.0
<<<<<<< HEAD
    p0n2 = 78.12 / (78.12 + 20.96) * 101.325
    p0o2 = 20.96 / (78.12 + 20.96) * 101.325
=======
>>>>>>> cba2043dcfd5d7de9f6eccd5102078c2924e7483

    x_val = np.array(xw)
    y_val = np.array(x0[8])

    c_sym = symbols('c:8')
    p_sym = symbols('pco2')
    xi_sym = symbols('xi:5')

    sym_x = np.append(
        np.append(
            np.array(c_sym), np.array(p_sym)
        ), np.array(xi_sym)
    )

    pprint(eq_set(sym_x, x0[:8], x0[9]))

    f0 = eq_set_const_p(
         sym_x, c0=x0[:8], p0co2=x0[9], p0n2=p0n2, p0o2=p0o2
    )

    # f0 = eq_set(sym_x, c0=x0[:8], p0co2=x0[9])

    solution = nsolve(f0, sym_x, x0.tolist(), verbose=True, tol=1e-14)

    solution = np.fromstring(nstr(solution,20).replace('\n','').replace('[','').replace(']',''), sep='   ')

    print -np.log10(solution[1])
    print -np.log10(solution[2])
    print -np.log10(solution[1]) -np.log10(solution[2])
    print x0
    print jac_eq_set_const_p(x0, x0[9], p0n2, p0o2)


    plt.ion()
    hl, = plt.plot([], [], 'o-')
    ax = plt.gca()
    # plt.show()
    hl.set_xdata(x_val)
    hl.set_ydata(y_val)
    plt.show()
    ax.relim()
    ax.autoscale_view()
    plt.draw()
    plt.pause(0.05)
    while it < 100 and (x0[8] < min_pco2 or x0[8] > 101.325):
        print '\n\nxw0nahco3 = ' + '%1.20g' % xw
        x0_n_m_1 = x0
        x0 = main(xw0nahco3=xw, x0=x0)
        if x0[8] > 101.325:
            xw = (0.07309 + xw) / 2
            x0 = x0_n_m_1
        elif x0[8] < min_pco2:
            x_val = np.append(x_val, [xw])
            y_val = np.append(y_val, [x0[8]])
            hl.set_xdata(x_val)
            hl.set_ydata(y_val)
            ax.relim()
            ax.autoscale_view()
            plt.draw()
            plt.pause(0.05)
<<<<<<< HEAD
            xw *= 1.025
=======
            xw *= 1.05
>>>>>>> cba2043dcfd5d7de9f6eccd5102078c2924e7483
        it += 1
    input('press enter to end')
