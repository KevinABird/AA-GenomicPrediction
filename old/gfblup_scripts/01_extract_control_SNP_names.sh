#!/bin/bash

for i in {1..1000}
do
    	 plink --noweb --ped ../data/CS_subsets/plinkrandomSubset_${i}.ped --map ../data/CS_subsets/plinkrandomSubset_${i}.map --write-snplist
	 mv plink.snplist ../data/control_markers/plink_control_snplist_${i}.txt
done
	 
