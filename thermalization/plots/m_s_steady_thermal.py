#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot the steady and thermal value of m_s to show thermalization in the quantum
critical regime.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('science')

L = 7
filename = f"../data/Z2Thermal/Z2_thermal_obs_L{L}.dat"
datas = np.loadtxt(filename)
xdat, thermal = datas[:, 0], datas[:, 3]
plt.plot(xdat, thermal, 'ro-', label="thermal")

massList = []; steady = []; std_steady = []
for mass in np.linspace(0.0, 1.0, 11):
    filename = f"../data/steady/Z2_steady_L{L}_m={mass:.1f}.dat"
    if os.path.exists(filename):
        datas = np.loadtxt(filename)
        massList.append(mass)
        steady.append(datas[-1, 3])
        std_steady.append(datas[-1, 4])

plt.errorbar(massList, steady, yerr=std_steady, fmt='go-', label="steady")

plt.legend(loc=0)
plt.xlabel(r'$m$')
plt.ylabel(r'$m_s$')
plt.title(r"steady vs thermal $m_s$ plot for L=${" + str(L) + "}$")

# plt.savefig(".pdf")
plt.show()
