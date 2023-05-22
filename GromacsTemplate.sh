#!/bin/bash
# This is a script for the generation of a single replica of protein-excipient interaction

# The I.D. number for the replica 
cnt=BLANKCNT

# The name of the excipient molecule
Molecule=BLANKMOL

# The number of molecules of excipient that are required in the system to make the concentration = 0.1 % w/w
MolNumber=BLANKNUM


# Make directory for each repetition
mkdir rep$cnt

# Copy the template topology into the directory
cp System.top ./rep$cnt/martini.top

# Make truncated octahedral box around the MARTINI protein 
gmx editconf -f Protein.gro -o ./rep$cnt/martininewbox.gro -bt octahedron -box 34.2400979892

# Randomly insert enough molecules of excipient to make the required concentration
gmx insert-molecules -f ./rep$cnt/martininewbox.gro -ci ../$Molecule.gro -o ./rep$cnt/MartiniVac.gro -nmol $MolNumber
echo "$Molecule       $MolNumber" >> ./rep$cnt/martini.top

# Perform preprocessing for minimisation in vacuum
gmx grompp -f ../minimisation.mdp -c ./rep$cnt/MartiniVac.gro -o ./rep$cnt/MinVac  -p ./rep$cnt/martini.top

# run a short minimisation in vacuum 
gmx mdrun -deffnm ./rep$cnt/MinVac -ntmpi 1

# Solvate with MARTINI polarisable water
gmx solvate -cp ./rep$cnt/MinVac.gro -cs polarwater.gro -p ./rep$cnt/martini.top -o ./rep$cnt/MartiniSolv.gro

# Perform preprocessing for the insertion of ions into the box to neutralise the system
gmx grompp -f ../ions.mdp -c ./rep$cnt/MartiniSolv.gro -p ./rep$cnt/martini.top -o ./rep$cnt/ions.tpr 

# Insert Na/Cl ions into system to neutralise
echo "14" | gmx genion -s ./rep$cnt/ions.tpr -o ./rep$cnt/MartiniSolvIons.gro -p ./rep$cnt/martini.top -pname NA -nname CL -neutral

# Preprocess and perform a second minimisation post solvation. This is probably overkill.
gmx grompp -f ../minimisation.mdp -c ./rep$cnt/MartiniSolvIons.gro -p ./rep$cnt/martini.top -o ./rep$cnt/minimisation 
gmx mdrun  -deffnm ./rep$cnt/minimisation -ntmpi 1

# Preprocess and perform relaxation MD, which is in this case in the same ensemble (NPT) and a shorter timestep
gmx grompp -f ../relax.mdp -c ./rep$cnt/minimisation.gro -p ./rep$cnt/martini.top -o ./rep$cnt/relax 
gmx mdrun  -deffnm ./rep$cnt/relax -ntmpi 1

# Preprocess and run production MD
gmx grompp -f ../md.mdp -c ./rep$cnt/relax.gro -p ./rep$cnt/martini.top -o ./rep$cnt/md 
gmx mdrun  -deffnm ./rep$cnt/md -ntmpi 1

