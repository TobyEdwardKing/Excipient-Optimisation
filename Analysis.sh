#!/bin/bash

# These are the APRs of HSA as defined by Aggrescan, used to define APR/non-APR contact
APRs="27 26 25 24 23 22 21 20 19 18 17 16 15 14 347 346 345 344 343 342 341 340 339 338 337 336 329 328 327 326 325 324 323 322 321 320 319 318 157 156 155 154 153 152 151 150 149 148 457 456 455 454 453 452 451 450 449 448 447 446 445 135 134 133 132 131 130 409 408 407 406 405 403 402 401 400 399 398 71 70 69 68 67 66 65 64 551 550 549 548 547 546 545 544 509 508 507 506 505 504 120 119 118 117 116 115 234 233 232 231 230 229 228 227 226 225 373 372 371 370 369 368 367 46 45 44 43 42 41 40 39 529 528 527 526 525 524 523 522 424 423 422 421 420 577 576 575 574 573"

# Make subdirectories
mkdir ./Indices
mkdir ./mindist
mkdir ./Analyses
mkdir ./SASA
# Make an index of the entire system
# This should be the same for every replica if ions are not included
# If ions are included, then an index will have to be made for each replica
# This allows us to remove waters from the trajectory.
	(echo "del 2-12"; echo "del 3-4"; echo "del 0"; echo "0|1"; echo "keep 2"; echo "q") | gmx make_ndx -f ./take1/md.gro -o ./Indices/System.ndx
# Make an index of two groups: all excipient beads (one group) and all HSA beads (another group)
	(echo "del 2-12"; echo "del 3-4"; echo "del 0"; echo "q") | gmx make_ndx -f ./take1/md.gro -o ./Indices/SurfToHSA.ndx

# Make an index of two groups: all excipient beads (one group) and all APR beads (another group)
	(echo "keep 13"; echo "ri $APRs"; echo "q") | gmx make_ndx -f ./take1/md.gro -o ./Indices/SurftoAPR.ndx




 
# cnt and cntmax to run through each replica consecutively
cnt=1
cntmax=5
# while loop to do so
while [ ${cnt} -le ${cntmax} ]; do 	
# make subdirectory for that replica
			mkdir ./mindist/$cnt
			mkdir ./Analyses/$cnt
# use gmx trjconv to remove the waters for both the structure file and the trajector
# removing the PBC at the same time
				gmx trjconv -f take$cnt/md.gro -o ./Analyses/$cnt/waterless.gro -s ./take$cnt/md.tpr -n ./Indices/System.ndx -boxcenter tric -pbc mol
				gmx trjconv -f take$cnt/md.xtc -o ./Analyses/$cnt/waterless.xtc -s ./take$cnt/md.tpr -n ./Indices/System.ndx -boxcenter tric -pbc mol
				# prepare a .tpr for SASA 
				gmx grompp -f minimisation.mdp -c ./Analyses/$cnt/waterless.gro -o ./Analyses/$cnt/analysis.tpr -p ./Analysis.top
# make subdirectory for the specific contact i.e. just APR or all of HSA
	mkdir ./mindist/$cnt/SurftoAPR
	mkdir ./mindist/$cnt/SurfToHSA
	# Find contacts between excipient and APR
	echo "0 1" | gmx mindist -f ./Analyses/$cnt/waterless.xtc -n ./Indices/SurftoAPR.ndx -on ./mindist/$cnt/SurftoAPR/numcount.xvg -d 0.6 -od ./mindist/$cnt/SurftoAPR/mindist.xvg -s ./Analyses/$cnt/analysis.tpr
	# Find contacts between excipient and HSA (all)
	echo "0 1" | gmx mindist -f ./Analyses/$cnt/waterless.xtc -n ./Indices/SurfToHSA.ndx -on ./mindist/$cnt/SurfToHSA/numcount.xvg -d 0.6 -od ./mindist/$cnt/SurfToHSA/mindist.xvg -s ./Analyses/$cnt/analysis.tpr
	# find SASA
	gmx sasa -f Analyses/$cnt/waterless.xtc -s Analyses/$cnt/analysis.tpr -o SASA/$cnt/SASA.xvg -surface 0 -output 1 -probe 0.26 -ndots 4800 

	# Move the counter on 
cnt=$[cnt+1]
done 
