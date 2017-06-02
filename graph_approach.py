import numpy as np
from matplotlib import pyplot as plt
dpi_res = 250

npoints = 200
figure = plt.figure(dpi= dpi_res)
ax = plt.axes()
# HA + H2O <<==>> H3O(+) + A(-)
# 2H2O <<==>> H3O(+) + HO(-)
pKa = 6.0
pKw = 14.0
pCT = 3.0
x = np.array([14 * x / float(npoints) for x in range(npoints + 1)])

ax.plot(x, -pCT - x - np.log10(10 ** -pKa + 10 ** -x))
ax.plot(x, -pCT - pKa - np.log10(10 ** -pKa + 10 ** -x))
ax.plot(x, -x, '--')
ax.plot(x, +x - pKw, '--')
ax.set_ylim(top=0)
ax.set_xlabel('pH')
ax.set_ylabel('-log Concentration (molar)')
ax.legend(['HA', 'A-', 'H+', 'OH-'], loc='best')
y_lim = ax.get_ylim()
ax.plot([pKa, pKa], y_lim, '--', color='gray')
plt.text(pKa, y_lim[1], 'pKa = ' + str(pKa))

# HNO3 + H2O <<==>> H3O(+) + NO3(-)
# 2H2O <<==>> H3O(+) + HO(-)
pKa = -1.0
pKw = 14.0
pCT = 3.0
x = np.array([-2 + (14 + 2) * x / float(npoints) for x in range(npoints + 1)])

figure2 = plt.figure(dpi=dpi_res)

ax2 = figure2.add_subplot(121)
ax2.plot(x, -pCT - x - np.log10(10 ** -pKa + 10 ** -x))
ax2.plot(x, -pCT - pKa - np.log10(10 ** -pKa + 10 ** -x))
ax2.plot(x, -x, '--')
ax2.plot(x, +x - pKw, '--')
ax2.set_ylim(top=0)
ax2.set_xlabel('pH')
ax2.set_ylabel('-log Concentration (molar)')
ax2.legend(['HNO3', 'NO3-', 'H3O+', 'HO-'], loc='best')
y_lim = ax2.get_ylim()
ax2.plot([pKa, pKa], y_lim, '--', color='gray')
plt.text(pKa, y_lim[1], 'pKa = ' + str(pKa))

# NH4X + H2O <<==>> NH3 + H3O+ + X-
# 2H2O <<==>> H3O(+) + HO(-)
pKa = +9.3
pKw = 14.0
pCT = 3.0
x = np.array([14 * x / float(npoints) for x in range(npoints + 1)])

ax3 = figure2.add_subplot(122)
ax3.plot(x, -pCT - x - np.log10(10 ** -pKa + 10 ** -x))
ax3.plot(x, -pCT - pKa - np.log10(10 ** -pKa + 10 ** -x))
ax3.plot(x, -x, '--')
ax3.plot(x, +x - pKw, '--')
ax3.set_ylim(top=0)
ax3.set_xlabel('pH')
ax3.set_ylabel('-log Concentration (molar)')
ax3.legend(['NH4+', 'NH3', 'H3O+', 'HO-'], loc='best')
y_lim = ax2.get_ylim()
ax3.plot([pKa, pKa], y_lim, '--', color='gray')
plt.text(pKa, y_lim[1], 'pKa = ' + str(pKa))

# Example 3.10
figure3 = plt.figure(dpi=dpi_res)

pKa = +9.57
pKw = 14.0
pCT = 5.0 - np.log10(3.0)
x = np.array([14 * x / float(npoints) for x in range(npoints + 1)])

ax4 = figure3.add_subplot(111)
ax4.plot(x, -pCT - x - np.log10(10 ** -pKa + 10 ** -x))
ax4.plot(x, -pCT - pKa - np.log10(10 ** -pKa + 10 ** -x))
ax4.plot(x, -x, '--')
ax4.plot(x, +x - pKw, '--')
ax4.set_ylim(top=0)
ax4.set_xlabel('pH')
ax4.set_ylabel('-log Concentration (molar)')
ax4.legend(['NH4+', 'NH3', 'H3O+', 'HO-'], loc='best')
y_lim = ax2.get_ylim()
ax4.plot([pKa, pKa], y_lim, '--', color='gray')
ax4.plot([8.5, 8.5], y_lim, '-', color='black')
plt.text(pKa, y_lim[1], 'pKa = ' + str(pKa))
plt.text(8.5,
         -pCT - pKa - np.log10(10 ** -pKa + 10 ** -8.5),
         '[NH3] = ' + '%1.2g' % \
         (10**(-pCT - pKa - np.log10(10 ** -pKa + 10 ** -8.5)))
         )

plt.show()