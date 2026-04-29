import os
import subprocess
import shutil
import sys


analisis='/all_aligned'

dir_base='~/NS3h_ATPhydrolysis/3-MixedSV-MD'
for i in range (1,5):
    os.chdir(dir_base+"/den%d" % i)
    for system in ['/apo','/RNA']:

        os.chdir(dir_base+"/den%d" % i+system)
        curr_dir=os.getcwd()

        if not os.path.isdir(curr_dir+analisis):
            os.mkdir(curr_dir+analisis)
        
        os.chdir(curr_dir+analisis)
        print(os.getcwd())

        if os.path.isfile(curr_dir+analisis+"/traj.in"):
                os.remove(curr_dir+analisis+"/traj.in")

        f=open(curr_dir+analisis+"/traj.in","a")

        # La referencia es el primer frame de la replica 1 de apo den1 para todos los sistemas
        f.write("parm "+dir_base+"""/den1/apo/rep1/x.prmtop\nparm ../rep1/x.prmtop\n\n

trajin ../rep1/md1.nc 101 -1 1 parmindex 1\n
trajin ../rep2/md1.nc 101 -1 1 parmindex 1\n
trajin ../rep3/md1.nc 101 -1 1 parmindex 1\n

center :1-439\nimage center familiar\n\n

reference """+dir_base+"""/den1/apo/rep1/md1.nc frame 1 parmindex 0\n\n
rms :1-439@CA,C,O,N reference :1-439@CA,C,O,N\n\n
trajout cent.nc netcdf parmindex 1\n\ntrajout ref.pdb onlyframes 1 1 parmindex 1\n\n""")

        f.close()

        subprocess.check_call(["/home/usuario/Programas/Amber22/amber22/bin/cpptraj" , "-i" , curr_dir+analisis+"/traj.in", "-o", curr_dir+analisis+"/traj.out"])
        print("cent hecho en "+curr_dir)
        
