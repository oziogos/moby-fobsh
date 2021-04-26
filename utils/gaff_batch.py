#!/usr/bin/env python

import os

parent_dir = os.getcwd()

if os.path.exists('output') == False:
	os.mkdir('output')
for mol2 in [i for i in os.listdir('single_molecules') if i.endswith('.mol2')]:
	os.chdir(parent_dir)
	mol2_actual = mol2.split('.mol2')[0]
	if os.path.exists('output/' + mol2_actual) == False:
		os.mkdir('output/' + mol2_actual)
	os.chdir(parent_dir + '/output/' + mol2_actual)
	os.system(f'cp ../../single_molecules/{mol2} .')
	os.system(f'/deploy/utils/run_ambertools.sh {mol2_actual} > /dev/null')
	os.system('echo writeFrcmod | parmed -p prmtop -c prmcrd > /dev/null')
	os.system(f'/deploy/utils/apply_FF_types.sh {mol2_actual} > /dev/null')
	os.system(f'/deploy/utils/convert_frcmod frcmod > neutral.ff')
	fp = open('neutral.ff')
	neutral_ff = fp.readlines()
	fp.close()
	header = ['begin\n',f'source {mol2_actual}_gafftypes\n',f'name {mol2_actual}_noQ_noImp\n','mol2tolammps molecular\n','supercell ortho -100 100 -100 100 -100 100\n','boundary f f f\n','topo\n']
	with open('instruct.txt', mode='w') as fp:
		for i in header + neutral_ff:
			print(i,file=fp,end='')
		print('end',file=fp)
	os.system('/deploy/utils/worker2 instruct.txt > /dev/null')		
	os.system(f"sed -i -e '2s/.*/{mol2_actual}/' {mol2_actual}_noQ_noImp.mol2")	
	os.system(f'/deploy/utils/lmp2psf2 {mol2_actual}_noQ_noImp.mol2 {mol2_actual}_noQ_noImp.lammps instruct.txt > {mol2_actual}.psf')
	print(mol2_actual)
	print('--------------------------------------------------------------------------------')
	os.system(f'cat {mol2_actual}.pot')
	print('--------------------------------------------------------------------------------')
