#!/bin/bash 

analisis="all_aligned"
sistemas=("apo" "RNA")

for tipo in {1..4}; do
	for sistema in "${sistemas[@]}"; do
		cd den$tipo/$sistema/$analisis

		rm -rf ./ETA_C1_DATA 
		rm -rf ./ETA_O1_DATA
		rm  eta.out
		rm watsel.pdb
		rm watall.pdb
		pwd

		cp /home/usuario/ns3/inhibidores/md_mixsv/eta_clust.tcl .
		
		`vmd -dispdev text -e eta_clust.tcl > eta.out`
		ls
	#	wait
	
	done

done


