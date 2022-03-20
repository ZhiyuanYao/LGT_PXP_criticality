#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('science')

L = 7; mass=0.6; tmax = 1E4

#------------------------------------------------------------------------------
# plot the long time part from [1, tmax]
#------------------------------------------------------------------------------
filename = f"../data/Z2Powerdynamics/Z2Powerdynamics_L{L}_m={mass:.1f}_tmax=100000.dat"
datas = np.loadtxt(filename)
tList, m_sList = datas[:, 0], datas[:, 2]
idx = (np.abs(tList - tmax)).argmin()
# here 0.5 account for the prefactor difference in the code and paper
plt.plot(tList[:idx], 0.5*m_sList[:idx], '.-', color='tab:blue')

#------------------------------------------------------------------------------
# plot the short time part from [0, 1]
#------------------------------------------------------------------------------
filename = f"../data/Z2dynamics_L7_m={mass:.1f}.dat"
datas = np.loadtxt(filename)
tList, m_sList = datas[:, 0], datas[:, 2]
plt.plot(tList[:idx], 0.5*m_sList[:idx], '.-', color='tab:blue')

# plt.xscale('log')
plt.xscale('symlog', linthreshx=1)
plt.xlabel(r'$t$')
plt.ylabel(r'$\bar{E}(t)$')
plt.legend(loc=1)
plt.gca().set_xlim(right=tmax)
plt.show()
