#!/bin/bash

# Create mol2 file with gaff atom names
antechamber -i $1.mol2 -fi mol2 -o $1_gafftypes.mol2 -fo mol2
# Use gaff atom types as mol2 atom names
# get number of atoms from the 3rd line of the mol2 file
STR=`sed -n "3,3p" $1_gafftypes.mol2`
atoms=`awk '{print $1}' <<< "$STR"`
# isolate mol2 types from original file
target="@<TRIPOS>ATOM"
STR=`grep -n $target $1.mol2`
IFS=':' read -ra NAMES <<< "$STR"   
line=${NAMES[0]}
sed -n ""$((line+1)),$((line+atoms))""p"" $1.mol2 | awk '{print $6}' > tripos_mol2_types 
# working on the gafftype file; print the header in the temp file
target="@<TRIPOS>ATOM"
STR=`grep -n $target $1_gafftypes.mol2`
IFS=':' read -ra NAMES <<< "$STR"   
line=${NAMES[0]}
head -n $line $1_gafftypes.mol2 > temp
# isolate gaff types
sed -n ""$((line+1)),$((line+atoms))""p"" $1_gafftypes.mol2 | awk '{print $6}' > gaff_types
# join columns
sed -n ""$((line+1)),$((line+atoms))""p"" $1_gafftypes.mol2 | awk '{print $1}' > left
sed -n ""$((line+1)),$((line+atoms))""p"" $1_gafftypes.mol2 | awk '{print $2}' > numeric_types
sed -n ""$((line+1)),$((line+atoms))""p"" $1_gafftypes.mol2 | awk '{print $3"\t"$4"\t"$5}' > middle
sed -n ""$((line+1)),$((line+atoms))""p"" $1_gafftypes.mol2 | awk '{print $7"\t"$8"\t"$9}' > right
paste left gaff_types middle tripos_mol2_types right >> temp
# print the rest...
sed -n "$((line+atoms+1)),\$p" $1_gafftypes.mol2 >> temp
# save mapping
paste numeric_types gaff_types > type_mapping
# clean up
rm left middle right
mv temp $1_gafftypes.mol2
rm gaff_types tripos_mol2_types numeric_types
