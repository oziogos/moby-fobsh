#!/usr/bin/env python

import os
import json
import  numpy as np

def make_dir(path_):
	if os.path.exists(path_) == False:
	        os.mkdir(path_)

species_dict = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'S': 16}
inv_species_dict = {i: j for j, i in species_dict.items()}

AOM_dict={
    'C': 1.385600,
    'N': 1.617102,
    'O': 1.505135,
    'F': 1.665190,
    'S': 1.641119,
}

mol_list = [i for i in os.listdir('output') if os.path.isdir(f'output/{i}')]
with open('/deploy/src/pyAOMlite/pyAOMlite-main/test/dimers.json') as fp:
	dimers = json.load(fp)
make_dir('fobsh')
make_dir('fobsh/topologies')

# fetch psf and pot files
for mol in mol_list:
	os.system(f'cp output/{mol}/{mol}.psf fobsh/topologies; cp output/{mol}/{mol}.pot fobsh/topologies')

# read templates
fp = open('/deploy/config/FORCE_EVAL.include.template.neutral')
force_eval_template = fp.readlines()
fp.close()
fp = open('/deploy/config/run.inp.template.neutral')
run_template = fp.readlines()
fp.close()

# AOM, DECOMP, FORCE_EVAL include files + VELOC.init, run.inp
for mol in mol_list:
	make_dir(f'fobsh/{mol}')
	aom = dimers[mol]['AOM_dict']
	dimer_types = list(dimers[mol]['dimers'].keys())
	species = dimers[mol]['dimers'][dimer_types[0]]['species']
	aom_coeffs = list(aom.values())[0]
	with open(f'fobsh/{mol}/AOM_COEFF.include', mode = 'w') as fp:
		for c,i in enumerate(aom_coeffs):
			print(f'{species[c]} {species_dict[species[c]]} 1 0.0 {i}', file = fp)
		for c,i in enumerate(aom_coeffs):
	        	print(f'{species[c]} {species_dict[species[c]]} 1 0.0 {i}', file = fp)
	with open(f'fobsh/{mol}/DECOMP.include', mode = 'w') as fp:
		print('&ENERGY_DECOMP', file = fp)
		print('\tINDEX_MOL_DECOMP\t1 2', file = fp)
		print(f'\tNUM_ACTIVE_ATOMS\t{len(species)}', file = fp)
		print('&END ENERGY_DECOMP', file = fp)
	force_eval = [i for i in force_eval_template]
	for c, i in enumerate(force_eval):
		if i.find('_mol_') != -1:
			force_eval[c] = force_eval[c].replace('_mol_',f'{mol}.pot')
	with open(f'fobsh/{mol}/FORCE_EVAL.include', mode = 'w') as fp:
		for i in force_eval:
			print(i, end = '', file = fp)
	with open(f'fobsh/{mol}/VELOC.init', mode = 'w') as fp:
		for i in range(len(species)):
			print('0.0 0.0 0.0', file = fp)
	run_actual = [i for i in run_template]
	sto_exp = ','.join([ f'{i}:{AOM_dict[inv_species_dict[i]]}' for i in sorted([species_dict[i] for i in set(species) if i != 'H'])])
	for c, i in enumerate(run_actual):
		if i.find('_atoms_per_site_') != -1:
			run_actual[c] = run_actual[c].replace('_atoms_per_site_',f'{int(len(species)/2)}')
		if i.find('_sto_exp_') != -1:
			run_actual[c] = run_actual[c].replace('_sto_exp_',sto_exp)

	with open(f'fobsh/{mol}/run.inp', mode = 'w') as fp:
		for i in run_actual:
			print(i, end = '', file = fp)

# coords
offset = 4.0
for mol in mol_list:
	dimer_types = list(dimers[mol]['dimers'].keys())
	for d in dimer_types:
		species = dimers[mol]['dimers'][d]['species']
		x, y, z = [np.array(dimers[mol]['dimers'][d][i]) for i in ['x','y','z']]
		X, Y, Z = [np.mean(x), np.mean(y), np.mean(z)]	
		xmin , ymin , zmin = [np.min(x), np.min(y), np.min(z)]
		xmax , ymax , zmax = [np.max(x), np.max(y), np.max(z)]
		lx, ly, lz = [xmax - xmin + 2 * offset , ymax - ymin + 2 * offset , zmax - zmin + 2 * offset]
		x = x - X + lx / 2
		y = y - Y + ly / 2
		z = z - Z + lz / 2
		atoms = len(x)
		with open(f'fobsh/{mol}/{mol}_{d}.xyz', mode = 'w') as fp:
			print(atoms, file = fp)
			print(f'{d} {lx} {ly} {lz}', file = fp)
			for i in range(atoms):
				print(f'{species[i]} {x[i]} {y[i]} {z[i]}', file = fp)
