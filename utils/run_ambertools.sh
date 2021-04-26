#!/bin/bash

# Create the prepi file (along with others...):
antechamber -i $1.mol2 -fi mol2 -o $1.prepin -fo prepi -j 4 -at gaff

# Check for missing FF parameters:
parmchk2 -i $1.prepin -o $1.frcmod -f prepi

# Prepare the input file for tleap (leap.in):
echo "source leaprc.gaff" > leap.in
echo "mods = loadAmberParams $1.frcmod" >> leap.in
echo "loadAmberPrep $1.prepin" >> leap.in
echo "saveAmberParm RES prmtop prmcrd" >> leap.in
echo "quit" >> leap.in

# Run tleap:
tleap -s -f leap.in

