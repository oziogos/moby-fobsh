#!/usr/bin/env python

import json
from scipy.optimize import curve_fit
import numpy as np

with open('/deploy/src/pyAOMlite/pyAOMlite-main/test/dimers.json') as fp:
        reference_data = json.load(fp)

with open('fobsh_Sab.json') as fp:
	fobsh_data = json.load(fp)

x = []
y = []
for mol in reference_data.keys():
	for d in reference_data[mol]['dimers'].keys():
		x.append(float(reference_data[mol]["dimers"][d]["reference_Sab"]))
		y.append(float(fobsh_data[mol]["dimers"][d]["fobsh_Sab"]))

def f(x, a, b):
	return a * x + b

popt, pcov = curve_fit(f, x, y)
residuals = np.array([y_-f(x[c],popt[0],popt[1]) for c,y_ in enumerate(y)])
ss_res = np.sum(residuals**2)
ss_tot = np.sum((np.array(y)-np.mean(y))**2)
R2 = 1 - (ss_res / ss_tot)

print(f'y = a * x + b linear fit:\na = {popt[0]:.6f}; b = {popt[1]:.6f}')
print(f'R2 = {R2:.6f}')
