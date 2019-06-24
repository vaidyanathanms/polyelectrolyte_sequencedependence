import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

free_chains  = [16,32,48,64,72,80,100]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4]
py_cutoff    = 1.5

maindir = os.getcwd()
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
        ntot_chains = free_chains[ifree] + graft_chains

        if archarr[iarch] == 1:
            print( "Archval: Block_Block")
            dirstr = 'bl_bl'
            fylstr = 'block_block'
        elif archarr[iarch] == 2:
            print( "Archval: Block_Alter")
            dirstr = 'bl_al'
            fylstr = 'block_alter'
        elif archarr[iarch] == 3:
            print( "Archval: Alter_Block")
            dirstr = 'al_bl'
            fylstr = 'alter_block'
        elif archarr[iarch] == 4:
            print( "Archval: Alter_Alter")
            dirstr = 'al_al'
            fylstr = 'alter_alter'
        else:
            print( "Unknown Architecture")
            
        workdir3 = workdir2 + '/' + dirstr
            
        if not os.path.isdir(workdir3):
            os.mkdir(workdir3)
                
        os.chdir(workdir3)
        destdir = os.getcwd()

        print( "Generating datafiles ...")

        print( "Copying Files")

        # Manipulate Files
            
        dataname = "PEdata_new_"+str(free_chains[ifree])+"_"+fylstr+".dat"
        datapath = destdir + '/' + dataname
            
        restartpath = destdir + '/restart1'
        
        if os.path.exists(restartpath):

            subprocess.call(["mpirun","-np","24", "./lmp_mesabi", "-r",
                             restartpath, dataname])

            if os.path.exists(datapath):
		
                print( "Dataname for LAMMPS: ", dataname)
            
            else:
                print("Something wrong in generating restart datafiles")
                continue
                
        else:
        
            print( dataname, "not found in", datapath)
            continue

        print( "Copy Successful - Manipulating Input Files")
        
