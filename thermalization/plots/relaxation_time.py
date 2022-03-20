#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program is to determine the relaxation time for m_s.

The knowledge of relaxation time is necessary for us to do average in a
given time sector to get an unbiased estimate of the steady value.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rcParams["font.family"] = "Times New Roman"

massList = np.linspace(0.0, 1.0, 11)
def check_relaxation_time(L=7, xmax=1E5, tmax=1E5, massList=massList):
    """
    function to determine relaxation time visually: 1E3 safe estimate
    """
    for mass in massList:
        plt.figure(figsize=(9, 4))
        filename = f"../data/Z2Powerdynamics/Z2Powerdynamics_L{L}_m={mass:.1f}_tmax={int(tmax)}.dat"
        datas = np.loadtxt(filename)
        tList, stag_mag = datas[:, 0], datas[:, 2]
        plt.plot(tList, stag_mag, '.-', label=f"m={mass:.1f}")
        plt.xscale('log')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$m_s$')
        plt.legend(loc=1)
        plt.gca().set_xlim(right=xmax)
        plt.show()

if __name__ == '__main__':
    check_relaxation_time(7)
