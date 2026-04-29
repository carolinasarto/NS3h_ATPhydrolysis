#!/bin/bash 

#Analisis completo Helicasa
#Esto va en la carpeta analisis/ de cada sistema

if [ -e "traj_walker_otro.in" ]

       then
               rm traj_walker_otro.in
       fi

echo "trajin ../cent.nc" >> traj_walker_otro.in

echo "


distance out gly_glu.dat :235@N :106@CD
distance out arg1_gln.dat :281@CZ :277@OE1
distance out gln_glu.dat :277@NE2 :106@CD
distance out gln_thr.dat :277@NE2 :138@OG1


" >> traj_walker_otro.in
	` /home/usuario/Programas/Amber22/amber22/bin/cpptraj ../x.prmtop traj_walker_otro.in > traj_walker_otro.out`

