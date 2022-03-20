#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program plots the calculated steady value in [t0, 2*t0] as a function of
t0 to furthur check the convergence of the steady value.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif')

L = 7
for mass in np.linspace(0.0, 1.0, 11):
    fig, ax = plt.subplots(1, 1)
    filename = f"../data/steady/Z2_steady_L{L}_m={mass:.1f}.dat"
    datas = np.loadtxt(filename)
    t0List, steady, std_steady = datas[:, 0], datas[:, 3], datas[:, 4]
    plt.errorbar(t0List, steady, yerr=std_steady, fmt='go-', label=f"m = {mass:.1f}")
    plt.xlabel(r'$t_0$')
    plt.ylabel(r'$m_s$')
    plt.xscale('log')
    plt.legend(loc=0)
    plt.show()
