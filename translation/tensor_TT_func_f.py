#!/usr/bin/env python
# coding: utf-8

import sys
import matplotlib.pyplot as plt
import numpy as np

import scipy.integrate as integrate
import scipy.special as special

from scipy.interpolate import interp1d
from datetime import datetime
from datetime import timedelta


alpha = 400;
f = 500*10**14.96;
Lambda = 10**-10.36;
m = Lambda*Lambda/f/np.sqrt(2);
GeVtog = 1.783e-24;
GeVtom = 1.973e-16;
G = 6.71e-39;
Mpl = 1/np.sqrt(8*np.pi*G);
H0 = 1.4459e-42;
tau_osc = 5e40;
tau_rei = 7.9941e41;
tau0 = 2.243e42;

kmax = 6e-39;
fk = open("/tomerv3/Sida/b-mode-log-griding/f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat_tau/parameters/klist.dat");
txt = fk.read();
klist = [float(x) for x in txt.split()];
fk.close();

fk = open("/tomerv3/Sida/b-mode-log-griding/f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat_tau/parameters/dk.dat");
txt = fk.read();
dk = [float(x) for x in txt.split()];
fk.close();

fx = open("/tomerv3/Sida/b-mode-log-griding/f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat_tau/parameters/xgrid_HD.dat");
txt = fx.read();
xgrid = [float(x) for x in txt.split()];
fx.close();

fa = open("/tomerv3/Sida/b-mode-log-griding/tab_a_tau.dat");
txt = fa.read();
tab_a = np.reshape([float(x) for x in txt.split()],(4001,2))
fa.close();
func_a = interp1d(tab_a[:,0], tab_a[:,1], kind="cubic")

fa = open("/tomerv3/Sida/b-mode-log-griding/tab_ap_tau.dat");
txt = fa.read();
tab_ap = np.reshape([float(x) for x in txt.split()],(4001,2))
fa.close();
func_ap = interp1d(tab_ap[:,0], tab_ap[:,1], kind="cubic")


fg = open("/tomerv3/Sida/b-mode-log-griding/f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat_tau/GF_new/g1_kk_85.dat");
txt = fg.read();
ndim = int(len(txt.split())/2);
tab_g = np.reshape([float(x) for x in txt.split()],(ndim,2))
fa.close();
func_g1 = interp1d(tab_g[:,0], tab_g[:,1], kind="cubic")

fg = open("/tomerv3/Sida/b-mode-log-griding/f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat_tau/GF_new/g1p_kk_85.dat");
txt = fg.read();
tab_g = np.reshape([float(x) for x in txt.split()],(ndim,2))
fa.close();
func_g1p = interp1d(tab_g[:,0], tab_g[:,1], kind="cubic")

fg = open("/tomerv3/Sida/b-mode-log-griding/f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat_tau/GF_new/g2_kk_85.dat");
txt = fg.read();
tab_g = np.reshape([float(x) for x in txt.split()],(ndim,2))
fa.close();
func_g2 = interp1d(tab_g[:,0], tab_g[:,1], kind="cubic")

fg = open("/tomerv3/Sida/b-mode-log-griding/f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat_tau/GF_new/g2p_kk_85.dat");
txt = fg.read();
tab_g = np.reshape([float(x) for x in txt.split()],(ndim,2))
fa.close();
func_g2p = interp1d(tab_g[:,0], tab_g[:,1], kind="cubic")


fa = open("/tomerv3/Sida/b-mode-log-griding/f_14.96_L_-10.36_tau_1_N_100_kmax_6_nat_tau/test_jl/jl_200_85_d.dat");
txt = fa.read();
ndim = int(len(txt.split())/2);
tab_jl = np.reshape([float(x) for x in txt.split()],(ndim,2))
fa.close();
func_jl = interp1d(tab_jl[:,0], tab_jl[:,1], kind="cubic")


xlist = xgrid[0:5];
ll=200;
int_limit = 200;

def evalfunc(tp):
    ret = integrate.quad(lambda t: 1/func_a(t)*(func_g2(tp)*func_g1p(t)-func_g1(tp)*func_g2p(t)-func_ap(t)/func_a(t)*(func_g2(tp)*func_g1(t)-func_g1(tp)*func_g2(t)))
               *func_jl(t)/((tau0-t)**2 * klist[84]**2), 
               tp, tau0, limit=int_limit);
    return ret[0];

def evalfunc_call(tp):
    ret = integrate.quad(lambda t: 1/func_a(t)*(func_g2(tp)*func_g1p(t)-func_g1(tp)*func_g2p(t)-func_ap(t)/func_a(t)*(func_g2(tp)*func_g1(t)-func_g1(tp)*func_g2(t)))
               *special.spherical_jn(ll, (tau0-t)*klist[84])/((tau0-t)**2 * klist[84]**2), 
               tp, tau0, limit=int_limit);
    return ret[0];

# print(evalfunc(5.*10**40))

#ans = [[tp, evalfunc(tp)] for tp in xlist];
#ans

t_before = datetime.now();
print(evalfunc(5.*10**40))
t_after = datetime.now();
print((t_after - t_before).total_seconds());

t_before = datetime.now();
print(evalfunc_call(5.*10**40))
t_after = datetime.now();
print((t_after - t_before).total_seconds());

