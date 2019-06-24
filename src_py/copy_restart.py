#Code to copy all restart or any type of files to a common directory
import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

free_chains  = [72] #[16,32,48,64,80,100,150,200]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4]

maindir = os.getcwd()
archivedir = maindir + '/' + 'archive_restart'

if not os.path.isdir(archivedir):
    print("Generating Archive Directory")
    os.mkdir(archivedir)

scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_pework/src_lmp'
   
for ifree in range(len(free_chains)):
    
    print("Free Number of Chains: ", free_chains[ifree])
    workdir1 = scratchdir + 'pework'
    
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdir2 = workdir1 + '/n_' + str(free_chains[ifree])
	
    if not os.path.isdir(workdir2):
        os.mkdir(workdir2)
        
    for iarch in range(len(archarr)):

        free_ntot = free_chains[ifree]*free_mons
        graft_ntot = graft_chains*(graft_mons-tail_mons)

        if archarr[iarch] == 1:
            print("Archval: Block_Block")
            dirstr = 'bl_bl'
            fylstr = 'block_block'
        elif archarr[iarch] == 2:
            print("Archval: Block_Alter")
            dirstr = 'bl_al'
            fylstr = 'block_alter'
        elif archarr[iarch] == 3:
            print("Archval: Alter_Block")
            dirstr = 'al_bl'
            fylstr = 'alter_block'
        elif archarr[iarch] == 4:
            print("Archval: Alter_Alter")
            dirstr = 'al_al'
            fylstr = 'alter_alter'
        else:
            print("Unknown Architecture")
            
        workdir3 = workdir2 + '/' + dirstr
            
        if not os.path.isdir(workdir3):
            os.mkdir(workdir3)
                
        os.chdir(workdir3)
        destdir = os.getcwd()

        print("Copying Files")
            
        # * Heart of the Code

        restart_names = destdir + '/archival*'
        list_of_files = glob.glob(restart_names) 
        
        print("Dataname for LAMMPS: ", restart_names)
        print("Generating Restart Dir")

        restart_dir = destdir + '/restartfiles_n_' + \
            str(free_chains[ifree]) + '_' + dirstr

        if not os.path.exists(restart_dir):
                
            os.mkdir(restart_dir)


        if list_of_files:
	
            for fyl in list_of_files:
                
                subprocess.call(["mv", fyl, restart_dir])
                
        else:
		
            print(restart_names, "not found in", destdir)

        if restart_dir:
                
            print("Zipping Directories ..")
            tar_restart = restart_dir + ".tar.gz"
            print(tar_restart)
            subprocess.call(["tar", "-cvzf", tar_restart, restart_dir])
            desttar = archivedir + '/'
            subprocess.call(["mv", tar_restart, desttar])


