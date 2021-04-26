#!/usr/bin/env python

import os
import json

parent_dir = os.getcwd()

with open('/deploy/src/pyAOMlite/pyAOMlite-main/test/dimers.json') as fp:
        reference_data = json.load(fp)

mol_list = [i for i in os.listdir('fobsh/') if os.path.isdir(f'fobsh/{i}')]

fp = open(f'/deploy/config/TOPOLOGY-NEUTRAL-ONLY.include.template')
topology_template = fp.readlines()
fp.close()

topology_footer = {'97a': ['&KIND nf\n', 'ELEMENT N\n', '&END KIND\n']}

fobsh_results = {i: {'dimers':{}} for i in mol_list}

for mol in mol_list:
	xyz_files = [i for i in os.listdir(f'fobsh/{mol}') if i.endswith('.xyz') and i.find(mol) == 0]
	os.chdir(f'fobsh/{mol}')
	for xyz_file in xyz_files:
		fp = open(f'{xyz_file}')
		xyz = fp.readlines()
		fp.close()
		atoms = xyz[0].strip()
		name = f'{mol}_{xyz[1].split()[0]}.xyz'
		A, B, C = [xyz[1].split()[i] for i in [1,2,3]]
		topology = [i for i in topology_template]
		if mol in topology_footer:
			topology += topology_footer[mol]
		for c, i in enumerate(topology):
			if i.find('_A_') != -1:
				topology[c] = topology[c].replace('_A_',A)
				topology[c] = topology[c].replace('_B_',B)
				topology[c] = topology[c].replace('_C_',C)
			if i.find('_atoms_') != -1:
				topology[c] = topology[c].replace('_atoms_',atoms)
			if i.find('_xyz_file_') != -1:
				topology[c] = topology[c].replace('_xyz_file_',name)
			if i.find('_psf_') != -1:
				topology[c] = topology[c].replace('_psf_',f'{mol}.psf')
		with open(f'TOPOLOGY-NEUTRAL-ONLY.include_{name.split(".xyz")[0]}', mode = 'w') as fp:
			for i in topology:
				print(i, end = '', file = fp)
		os.system(f'cp TOPOLOGY-NEUTRAL-ONLY.include_{name.split(".xyz")[0]} TOPOLOGY-NEUTRAL-ONLY.include')
		os.system('$CP2K_EXEC -i run.inp > run.out')
		fp = open('run-pseudo-hamilt-1.xyz')
		res = fp.readlines()
		fp.close()
		os.system('rm run-pseudo*')
		value = res[3].split()[-1]
		fobsh_results[mol]['dimers'][xyz[1].split()[0]] = {'fobsh_Sab': value}
		print(f"{mol}\t{xyz[1].split()[0]}\t{reference_data[mol]['dimers'][xyz[1].split()[0]]['reference_Sab']}\t{value}")
	os.chdir(parent_dir)

with open('fobsh_Sab.json', mode = 'w') as fp:
	json.dump(fobsh_results,fp)
