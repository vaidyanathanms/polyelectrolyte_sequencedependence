#-----Automating SMD and US Calculations--------------------------
#-----Version-1: April-23-2018------------------------------------


import numpy
import os
import shutil
import subprocess
import sys
import glob
import os.path
from subprocess import call
import re
import random


#--Chain Details-------------------------------------------------
free_chains  = [32,48,64,72]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4]

#--Main Directories----------------------------------------------
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_pework/src_lmp'

for ifree in range(len(free_chains)):
    
    print("Copying Files: ", free_chains[ifree])
    smdsrcdir = maindir + '/smdinp_all'

    if not os.path.isdir(smdsrcdir):
        os.mkdir(smdsrcdir)
        
    smdsrcdir2 = smdsrcdir + '/n_' + str(free_chains[ifree])
    if not os.path.isdir(smdsrcdir2):
        os.mkdir(smdsrcdir2)

    workdir1 = scratchdir + 'pework'
    
    if not os.path.isdir(workdir1):
        print("Path not found", workdir1)
        continue

    workdir2 = workdir1 + '/n_' + str(free_chains[ifree])
	
    if not os.path.isdir(workdir2):
        print("Path not found", workdir2)
        continue
        
    for iarch in range(len(archarr)):

        ranval = random.random()
        refchain = int(ranval*free_chains[ifree]) + 1
        free_ntot = free_chains[ifree]*free_mons
        graft_ntot = graft_chains*(graft_mons-tail_mons)

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
            
        smdworkdir = workdir2 + '/' + dirstr
            
        if not os.path.isdir(smdworkdir):
            print("Path not found",smdworkdir)
            continue

        os.chdir(smdworkdir)

        smdfyl_work = smdworkdir + '/smd_colfile.inp'

        if not os.path.exists(smdfyl_work):
            print("File not found in", free_chains[ifree],"_",dirstr)
            continue

        smdfyl_backup = smdsrcdir2 + '/smd_colfile_' + \
            str(free_chains[ifree]) + '_' + str(dirstr) + \
            '.inp'
        shutil.copy2(smdfyl_work,smdfyl_backup)
