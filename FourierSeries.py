# Fourier Class
# 
# methods: (for regular, sin, and cos)
# plot series with N coefs given A and B (or not) 
# plot series with M-N coefs given A and B (or not)
#    

import numpy as np
from matplotlib import pyplot as plt
import sympy as sp
import seaborn as sns
# import IPython
# sp.init_printing(use_unicode=False, wrap_line=False)

class FourierSeries:
    """" A class to compute and plot Fourier Series of a given length and number of coeficients """""
    global x; x = sp.Symbol('x')
    
    def __init__(self, f, length = sp.pi, series_type = None):
        self.f = f
        self.L = length
        self.series_type = series_type
        self.a = []
        self.b = [0]
        self._series = []
    
    def __str__(self):
        if self.series_type:
            return "Fourier {} series for: ".format(self.series_type) + str(self.f(x))
        else:
            return "Fourier series for: " + str(self.f(x))



    def _calculate_a(self, N):
        n = len(self.a)
        if self.series_type == 'sin':
            self.a = [0 for _ in range(N)]
        elif self.series_type == 'cos':
            if(n == 0):
                self.a.append(sp.nsimplify((1 / (self.L)) * sp.integrate(self.f(x), (x, 0, self.L))))
                self._series.append(self.a[0])
                n = 1
            while (n < N): 
                self.a.append(sp.nsimplify((2 / self.L) * sp.integrate(self.f(x) *
                                                                       sp.cos(n * x * sp.pi / self.L), (x, 0, self.L))))
                n += 1
        else:
            if (n == 0):
                self.a.append(sp.nsimplify((1 / (2*self.L)) *
                                           sp.integrate(self.f(x), (x, -self.L, self.L))))
                self._series.append(self.a[0])
                n = 1
            while (n < N):
                self.a.append(sp.nsimplify((1 / self.L) * sp.integrate(self.f(x) * sp.cos(n * x * sp.pi / self.L), (x, -self.L, self.L))))
                n += 1
            

    def _calculate_b(self, N):
        n = len(self.b)
        if self.series_type == 'cos':
            self.b = [0 for _ in range(N)]
        elif self.series_type == 'sin':
            if(n == 0):
                self.b.append(sp.nsimplify((2 / (self.L)) * sp.integrate(self.f(x), (x, 0, self.L))))
                self._series.append(self.b[0])
                n = 1
            while (n < N): 
                self.b.append(sp.nsimplify((2 / self.L) * sp.integrate(self.f(x) *
                                                          sp.sin(n * x * sp.pi / self.L), (x, 0, self.L))))
                n += 1
        else:
            if (n == 0):
                self.b.append(sp.nsimplify((1 / (self.L)) *
                                           sp.integrate(self.f(x), (x, 0, self.L))))
                self._series.append(self.b[0])
                n = 1
            while (n < N):
                self.b.append(sp.nsimplify((1 / self.L) * sp.integrate(self.f(x) * sp.sin(n * x * sp.pi / self.L), (x, -self.L, self.L))))
                n += 1

    def _get_series(self, N):
        self._calculate_a(N)
        self._calculate_b(N)
        for n in range(len(self._series),N):
            self._series.append(0)
            self._series[n] = self._series[n-1] +\
            self.a[n] * sp.cos(n * x * sp.pi / self.L) + \
            self.b[n] * sp.sin(n * x * sp.pi / self.L)

    def series(self, N):
        self._get_series(N+1)
        return self._series[N]

    def plot_series(self, N, scale=1.5, color = 'r'):
        self._get_series(N)
        return sp.plotting.plot(self._series[N], (x, -self.L * scale, self.L * scale), show=False, title=str(self), line_color=color, steps=self.L / 4, label=str(N) + ' coefficients', legend=True)


    def plot_f_series_range(self, min_coef = 1, max_coef = 5, scale=1.5, cmap = 'autumn_r', f_color = 'b'):
        self._get_series(max_coef+1)
        
        p = sp.plotting.plot(self.f(x), (x, -self.L * scale, self.L * scale), show=False, line_color=f_color,
                            title=str(self), legend=False)
        colors = sns.color_palette(
            cmap, n_colors=max_coef - min_coef + 2)
        for N in range(min_coef, max_coef+1):
            p.extend(sp.plotting.plot(self._series[N], (x, -self.L * scale, self.L * scale), show=False, title=str(
                self), line_color=colors[N], steps=self.L / 4, label=str(N) + ' coefficients', legend=True))
        return p

    def plot_series_range(self, min_coef = 1, max_coef = 5, scale=1.5, cmap = 'autumn_r', f_color = 'b', size = (10,8)):
        self._get_series(max_coef+1)
        
        p = sp.plotting.plot(show=False, line_color=f_color,
                            title=str(self), legend=True)
        colors = sns.color_palette(
            cmap, n_colors=max_coef - min_coef + 2)
        for N in range(min_coef, max_coef+1):
            p.extend(sp.plotting.plot(self._series[N], (x, -self.L * scale, self.L * scale), show=False, title=str(
                self), line_color=colors[N], steps=self.L / 4, label=str(N) + ' coefficients', legend=True, size = size))
            p.rcParams['figure.figsize'] = 10, 8
        return p
