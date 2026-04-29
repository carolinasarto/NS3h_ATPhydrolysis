import numpy as np
import os
import subprocess
import shutil

ATP=['/home/usuario/ns3/qmmm/paper/md/atp']
RNA7=['/home/usuario/ns3/qmmm/paper/md/rna7A']
RNA16=['/home/usuario/ns3/qmmm/paper/md/rna16A//clustering2/5endDown','/home/usuario/ns3/qmmm/paper/md/rna16A/clustering2/5endUp']

dir_list = ATP+RNA7+RNA16

root_dir=os.getcwd()
for directory in dir_list:
    os.chdir(directory)
    curr_dir=os.getcwd() #curr_dir and directory do not have / at the end
    path=curr_dir+'/analisis'

    if not os.path.isdir(curr_dir+'/analisis'):
        os.mkdir(path)
        print(path+' hecho')

    script = '/runPtraj_allrep.sh'
    if os.path.isfile(path+script):
        os.remove(path+script)
    print(curr_dir)
   
    source='/home/usuario/ns3/qmmm/paper/md'+script
    os.chdir(path)

    dest=path+script
    shutil.copyfile(source,dest)
    subprocess.check_call(['chmod','+x',"."+script])
    print('Corriendo '+script+' en '+directory)
    subprocess.check_call('.'+script)
    

