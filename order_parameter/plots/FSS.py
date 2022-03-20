#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""This program use the theoretical calculations of scaling variables to do
finite-size scaling."""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, mark_inset)
# plt.style.use('science')

LList = [5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]

fig, ax1 = plt.subplots(figsize=(4.8, 3.2))

color_set = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
             '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#972fb3']
ax1.set_prop_cycle(color=color_set)

ylabel = ''
ydat   = []
for L in LList:
    filename = "../data/FSS_L" + str(L) + ".dat"
    datas = np.loadtxt(filename)
    massList, corrList, stag_mList, m_s_abs_List, m_s2List, m_s4List = datas[:,
            0], datas[:, 1], datas[:, 2], datas[:, 3], datas[:, 4], datas[:, 5]
    RList = [1-m_s4List[i]/(3*m_s2List[i]**2) for i in range(len(m_s2List))]
    #--------------------------------------------------------------------------
    ydat = RList;                 ylabel = r'$R = 1 - m_s^4/3m_s^2$'
    ydat = m_s4List*(L)**0.50;    ylabel = r'$m_s^4 \cdot L^{1/2}$'
    ydat = corrList*(L)**0.25;    ylabel = r'$ S_{zz} \cdot L^{1/4}$'
    ydat = m_s2List*(L)**0.25;    ylabel = r'$m_s^2 \cdot L^{1/4}$'
    ydat = m_s_abs_List*(L)**0.125; ylabel = r'$|m_s| \cdot L^{1/8}$'
    #--------------------------------------------------------------------------
    # remember to multiply the correct prefactor since our main.jl program use
    # 1 and -1 for spin basis, while the paper use -1/2 and 1/2.
    #--------------------------------------------------------------------------
    ydat = 0.5*m_s_abs_List*(L)**0.125; ylabel = r'$|m_s| \cdot L^{1/8}$'
    ax1.plot(massList, ydat, 'o-', ms= 3, label=r"$" + str(L) + "$")


#------------------------------------------------------------------------------
# plot zoom data in region [m_min, m_max]
#------------------------------------------------------------------------------
m_min, m_max = 0.56, 0.7
# lower-left corner and its width and heigh, transform defaults to ax.transAxes
# ax2 = ax1.inset_axes([0.55, 0.1, 0.4, 0.4], transform=ax1.transAxes)

ax2 = inset_axes(ax1, width="100%", height="100%", loc='lower left',
        bbox_to_anchor= (0.075, 0.62, 0.35, 0.35), bbox_transform=ax1.transAxes,
        borderpad=0)

color_set = ['#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#972fb3']
ax2.set_prop_cycle(color=color_set)

LList = [11, 13, 15, 17, 19, 21, 23, 25]
for L in LList:
    filename = "../data/FSS_L" + str(L) + ".dat"
    datas = np.loadtxt(filename)
    massList, corrList, stag_mList, m_s_abs_List, m_s2List, m_s4List = datas[:,
            0], datas[:, 1], datas[:, 2], datas[:, 3], datas[:, 4], datas[:, 5]
    RList = [1-m_s4List[i]/(3*m_s2List[i]**2) for i in range(len(m_s2List))]
    #--------------------------------------------------------------------------
    ydat = RList;                 ylabel = r'$R = 1 - m_s^4/3m_s^2$'
    ydat = m_s4List*(L)**0.50;    ylabel = r'$m_s^4 \cdot L^{1/2}$'
    ydat = corrList*(L)**0.25;    ylabel = r'$ S_{zz} \cdot L^{1/4}$'
    ydat = m_s2List*(L)**0.25;    ylabel = r'$m_s^2 \cdot L^{1/4}$'
    ydat = m_s_abs_List*(L)**0.125; ylabel = r'$|m_s| \cdot L^{1/8}$'
    #--------------------------------------------------------------------------
    # remember to multiply the correct prefactor since our main.jl program use
    # 1 and -1 for spin basis, while the paper use -1/2 and 1/2.
    #--------------------------------------------------------------------------
    ydat = 0.5*m_s_abs_List*(L)**0.125; ylabel = r'$|\bar{E}| \cdot L^{1/8}$'
    #--------------------------------------------------------------------------
    index_min, index_max = 0, 0
    for i, mass in enumerate(massList):
        if abs(mass - m_min) < 1E-10:
            index_min = i
        if abs(mass - m_max) < 1E-10:
            index_max = i
    ax2.plot(massList[index_min:index_max+1], ydat[index_min:index_max+1], 'o-', ms=3, zorder=10)


mark_inset(ax1, ax2, loc1=1, loc2=3, fc="none", ec='0.5', zorder=5)

ax1.legend(loc=4, prop={'size': 10}, ncol=3, markerscale=1, frameon=True)
# ax1.legend(loc=4, ncol=2)

ax1.set_xlabel(r'$m/\tilde{t}$', fontsize=14)
ax1.set_ylabel(ylabel, fontsize=13)
ax2.patch.set_facecolor('white')

# Set tick font size
for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(13)

# plt.savefig("FSS_E_abs.pdf")
plt.show()
