#--For Bias Simulations of Random Chains to Compute Delta U----#

import numpy
import os
import shutil
import subprocess
import sys
import glob
from subprocess import call
import re
import random

restart    = 0 
pybiascons = 5.0
pyregmax   = 40.0
nsamplechains = 5

free_chains  = [80]#[32,48,64,80,100]#,150,200]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4]

maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_pework/src_lmp'
   
for ifree in range(len(free_chains)):
    
    print( "Free Number of Chains: ", free_chains[ifree])
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

        while refchain > free_chains[ifree] - nsamplechains:
            ranval = random.random()
            refchain = int(ranval*free_chains[ifree]) + 1


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
            print("Main path not found", workdir3)
            continue
           
        if restart == 0:
                            
            os.chdir(workdir3)
            workheaddir = os.getcwd()
            biasworkdir = workdir3 + '/biascalc3'

            if not os.path.isdir(biasworkdir):
                os.mkdir(biasworkdir)

            os.chdir(biasworkdir)
            destdir = os.getcwd()

            print( "Starting from beginning")
            
            print( "Copying Files")
            
            # Copy Files
            
            srcfyl = lmpdir + '/in.bias_var'
            desfyl = destdir + '/in.bias_var'
            shutil.copy2(srcfyl, desfyl)
            srcfyl = lmpdir + '/jobbias.sh'
            desfyl = destdir + '/jobbias.sh'
            shutil.copy2(srcfyl, desfyl)
                        
            srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            desfyl = destdir + '/lmp_mesabi'
            shutil.copy2(srcfyl, desfyl)
            
            # Copying Restart files from main work dir

            srcfyl = maindir + '/archive_restart'+'/restartfiles_n_' + \
                str(free_chains[ifree]) + '_' +  dirstr + '.tar.gz'

            archfiles = workheaddir +  '/archival_*.restart'
            
            list_of_files = glob.glob(archfiles)
        
            oldtime = 0

            if list_of_files:

                for fyl in list_of_files:

                    splitdata = fyl.split("/")
                    fylname  = splitdata[len(splitdata)-1]
                    timedata = fylname.split("_")
                    timeval  = float(timedata[1].split(".")[0])

                    if timeval > oldtime:

                        oldtime = timeval

                finarchfyl=workheaddir+'/archival_' + \
                    str(int(oldtime)) + '.restart'

                desfyl = destdir + '/archival_' + \
                    str(int(oldtime)) + '.restart'
                shutil.copy2(finarchfyl,desfyl)

                subprocess.call(["mpirun","-np","24", "./lmp_mesabi", "-r",
                                 desfyl, "databias.data"])


            elif os.path.exists(srcfyl):
		
                print( "Restart file for LAMMPS: ", srcfyl)
                print("Current Dir:", os.getcwd())
                desfyl= destdir+'/restartfiles_n_'+str(free_chains[ifree]) \
                    + '_' +  dirstr + '.tar.gz'
            
                shutil.copy2(srcfyl, desfyl)

# Untar restart files and copy it as restart1
            
                subprocess.call(["tar","-xvzf", desfyl])

                archfiles = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                    str(free_chains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                    str(free_chains[ifree]) + '_' + dirstr +  '/' +  \
                    'archival_*.restart'

                list_of_files = glob.glob(archfiles)
            
                for fyl in list_of_files:

                    splitdata = fyl.split("/")
                    fylname  = splitdata[len(splitdata)-1]
                    timedata = fylname.split("_")
                    timeval  = float(timedata[1].split(".")[0])
                    
                    if timeval > oldtime:

                        oldtime = timeval

                finarchfyl = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                    str(free_chains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                    str(free_chains[ifree]) + '_' + dirstr +  '/' +  \
                    'archival_' + str(int(oldtime)) + '.restart'

                desfyl = destdir + '/archival_' + str(int(oldtime)) + '.restart'
                subprocess.call(["mpirun","-np","24", "./lmp_mesabi", "-r",
                                 desfyl, "databias.data"])

                shutil.copy2(finarchfyl,desfyl)

            else:
		
                print(srcfyl, "not found in", maindir)
                os.chdir(maindir)

            print("Successfully generated input data file")

            # Manipulate Files
            pyinit   = graft_chains*graft_mons+(refchain-1)*free_mons+1
            pyfin    = pyinit + nsamplechains*free_mons - 1
	    
            launch_fyl = 'in.bias'
            fr = open('in.bias_var','r')
            fw = open(launch_fyl,'w')
            fid = fr.read().replace("py_init",str(pyinit)).\
                replace("py_fin",str(pyfin)).\
                replace("py_regmax",str(pyregmax)).\
                replace("py_biascons",str(pybiascons))
            fw.write(fid)
            fw.close()
            fr.close()
            
            # Removing unnecessary files

            print( "Submitting Jobs..")
                
            subprocess.call(["qsub","jobbias.sh"])
	    
            os.chdir(maindir)
