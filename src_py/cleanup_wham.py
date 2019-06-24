# To clean up after WHAM analysis and save the restart files/colvars------
# file from US sampling---------------------------------------------------
#------------Version: April-28-2018---------------------------------------


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

#--SMD or UmbrellaSampling---------------------------------------

logic_umb = 1 #0 - SMD US; 1 - US

#--Chain Details-------------------------------------------------
free_chains  = [80,100]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4]
rvalarr = []

#--Umbrella Sampling Details-------------------------------------
prefix  = 'US_dir_'
opcolfile = 'out.colvars.traj'
ipcolfile = 'umb_colfile.inp'
outwham   = 'whamout.txt'
force_constant =  0.40
stripterm = 'new' #new or NULL

if stripterm == 'new':
    prefix = prefix + stripterm


#--Main Directories----------------------------------------------
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
            print("Path not found", smdworkdir)
            continue

        os.chdir(smdworkdir)
        destdir = os.getcwd()
        
        # Create backup directory in ${HOME} folder
        backupdir = maindir + '/USbackup_colvar'
        if not os.path.isdir(backupdir):
            os.mkdir(backupdir)

        backdir = backupdir + '/n_' + str(free_chains[ifree])
        if not os.path.isdir(backdir):
            os.mkdir(backdir)        
        
        # Create backup directory in ${SCRATCH} folder
        whamoutdir = workdir1 + '/whamout_all'

        if not os.path.isdir(whamoutdir):
            os.mkdir(whamoutdir)
            print(whamoutdir, "created")

        whaminddir = whamoutdir + '/n_' + str(free_chains[ifree])+\
                     '_'+ dirstr
        if not os.path.isdir(whaminddir):
            os.mkdir(whaminddir)
            print(whaminddir, "created")


        #Copy whamout.txt
        outfyle = whamoutdir + '/whamout_' + dirstr + '_' \
            + str(free_chains[ifree]) + '.txt'
        srcfyle = smdworkdir + '/' + outwham
        
        if not os.path.exists(srcfyle):
            print(srcfyle,"Not found")
        else:
            shutil.copy2(srcfyle,outfyle)

        # Start copying all harmonic centers data

        alldirs = prefix + '*'
        list_of_dirs = glob.glob(alldirs)

        rvalarr = []
        if list_of_dirs:

            for dirval in range(len(list_of_dirs)):

                splice  = list_of_dirs[dirval].split('_')
                rval = splice[len(splice)-1]
                if stripterm:
                    rval = rval.strip(stripterm)

                rvalarr.append(float(rval))

            for dirval in sorted(rvalarr):

                dirname = smdworkdir + '/' + prefix + str(dirval)

                #Copy output colfile
                fylname = dirname + '/' + opcolfile
                if os.path.exists(fylname):
                    desfyl = backdir + '/' + opcolfile + '_' + fylstr \
                             + '_' +str(dirval)
                    shutil.copy2(fylname,desfyl)

                    desfyl2 = whaminddir + '/' + opcolfile + '_'+str(dirval)
                    shutil.copy2(fylname,desfyl2)

                else:
                    print(opcolfile,"not found in", dirname)
                    


                #Copy input colfile
                fylname = dirname + '/' + ipcolfile
                if os.path.exists(fylname):
                    desfyl = backdir + '/' + ipcolfile + '_' + fylstr \
                             + '_' +str(dirval)
                    shutil.copy2(fylname,desfyl)
                    
                    desfyl2 = whaminddir + '/' + ipcolfile + '_'+str(dirval)
                    shutil.copy2(fylname,desfyl2)
                   
                else:
                    print(ipcolfile,"not found in", dirname)
                    



                #Copy Restart files
                restfiles = glob.glob(dirname + '/' + 'smdarchival*')
                if restfiles:

                    latest_fyl = max(restfiles, key=os.path.getctime)
                    backfyl = backdir + '/' + 'USrestart' + '_' + fylstr \
                    + '_' +str(dirval)
                
                    shutil.copy2(latest_fyl, backfyl)
                
                else:

                    print("no restart files found in", dirname)
            


                
