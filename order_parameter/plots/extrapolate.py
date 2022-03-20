#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('science')

LList = [11, 13, 15, 17, 19, 21, 23, 25]
xdat = [1.0/L for L in LList[:-1]]

# data from ../extrapolate.jl
ydat_m_s = [0.5841, 0.5983, 0.6084, 0.6159, 0.6217, 0.6264, 0.6301]
ydat_m_s_abs = [0.5932, 0.6040, 0.6123, 0.6187, 0.6238, 0.6280, 0.6314]

#------------------------------------------------------------------------------
# linear fit last 5 data points to get critical point
#------------------------------------------------------------------------------
ydat = ydat_m_s; n = 2
ydat = ydat_m_s_abs; n = 2

plt.plot(xdat, ydat, 'o', color='tab:blue', zorder=10)

p = np.polyfit(xdat[n:], ydat[n:], 1)
func = np.poly1d(p)
xfit = np.linspace(0, 1/LList[0]*1.2)
yfit = func(xfit)
plt.plot(xfit, yfit, 'k-')
# plt.plot(xfit, yfit, 'k-', label=str(round(func(0.0), 4)))

plt.xlabel(r"$1/L$", fontsize=11)
plt.ylabel(r"$g_c$", fontsize=12)
plt.gca().set_xlim(left=0)
plt.legend(loc=0)

# plt.savefig("g_c_m_s_abs.pdf")
plt.show()
