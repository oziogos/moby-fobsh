#!/bin/bash
cd $(dirname $0)
gcc -o convert_frcmod convert_frcmod.c -lm -O3
gcc -o lmp2psf2 lmp2psf2.c -lm -O3
gcc -o worker2 worker2.c -lm -O3
gcc -o topo2 topo2.c -lm -O3
