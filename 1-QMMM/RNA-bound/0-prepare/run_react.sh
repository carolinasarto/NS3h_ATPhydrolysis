#!/bin/bash

echo 'Minimizacion MM'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c x.rst7 -i min_mm.in -o min_mm.out -r min_mm.rst7  -ref x.rst7

echo 'Minimizacion QM/MM'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c min_mm.rst7 -i min_qmmm.in -o min_qmmm.out -r min_qmmm.rst7 -ref min_mm.rst7

echo 'Termalizacion MM 0 a 100K'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c min_qmmm.rst7 -i term1_mm.in -o term1_mm.out -r term1_mm.rst7 -x term1_mm_traj.nc -ref min_qmmm.rst7

echo 'Termalizacion QM/MM 0 a 100K'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c term1_mm.rst7 -i term1_qmmm.in -o term1_qmmm.out -r term1_qmmm.rst7 -x term1_qmmm_traj.nc -ref term1_mm.rst7

echo 'Termalizacion MM 100 a 200K'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c term1_qmmm.rst7 -i term2_mm.in -o term2_mm.out -r term2_mm.rst7 -x term2_mm_traj.nc -ref term1_qmmm.rst7

echo 'Termalizacion QM/MM 100 a 200K'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c term2_mm.rst7 -i term2_qmmm.in -o term2_qmmm.out -r term2_qmmm.rst7 -x term2_qmmm_traj.nc -ref term2_mm.rst7

echo 'Termalizacion MM 200 a 300K'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c term2_qmmm.rst7 -i term3_mm.in -o term3_mm.out -r term3_mm.rst7 -x term3_mm_traj.nc -ref term2_qmmm.rst7

echo 'Termalizacion QM/MM de 200 a 300K'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c term3_mm.rst7 -i term3_qmmm.in -o term3_qmmm.out -r term3_qmmm.rst7 -x term3_qmmm_traj.nc -ref term3_mm.rst7

echo 'Equilibracion MM a 300K'
$AMBERHOME/bin/pmemd.cuda_SPFP -O -p x.prmtop -c term3_qmmm.rst7 -i eq1_mm.in -o eq1_mm.out -r eq1_mm.rst7 -x eq1_mm_traj.nc -ref term3_qmmm.rst7

echo 'Equilibracion QM/MM a 300K'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c eq1_mm.rst7 -i eq1_qmmm.in -o eq1_qmmm.out -r eq1_qmmm.rst7 -x eq1_qmmm_traj.nc -ref eq1_mm.rst7

echo 'Equilibracion MM a 1 bar'
$AMBERHOME/bin/pmemd.cuda_SPFP -O -p x.prmtop -c eq1_qmmm.rst7 -i eq2_mm.in -o eq2_mm.out -r eq2_mm.rst7 -x eq2_mm_traj.nc -ref eq1_qmmm.rst7

echo 'Equilibracion QM/MM a 1 bar'
/home/usuario/Programas/anaconda3/envs/AmberTools22/bin/sander -O -p x.prmtop -c eq2_mm.rst7 -i eq2_qmmm.in -o eq2_qmmm.out -r eq2_qmmm.rst7 -x eq2_qmmm_traj.nc -ref eq1_mm.rst7

echo 'Dinamica MM NVT'
$AMBERHOME/bin/pmemd.cuda_SPFP -O -p x.prmtop -c eq2_qmmm.rst7 -i md.in -o md1.out -r md1.nc -x md1_traj.nc -ref eq2_qmmm.rst7

echo 'vmd x.prmtop *_traj.nc
same resid as within 5 of resname MG'



