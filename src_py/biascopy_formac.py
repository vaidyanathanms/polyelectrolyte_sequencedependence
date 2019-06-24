#To copy results of bias calculations to Desktop
import numpy
import os
import shutil
import subprocess
import sys
import glob
import os.path
from subprocess import call

#-------------------------Input Data-------------------------------

free_chains  = [32,48,64,72]#,200]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4]
rcutoff      = '1.50'

#----------Flags for copying to avoid redundancy-------------------

archflag = 1
enerout_bias = 0
enerout_nobias = 0
adsorb_chain = 1

#-------------------------Main calculations------------------------

maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_pework/src_lmp'
   
for ifree in range(len(free_chains)):
    
    print("Free Number of Chains: ", free_chains[ifree])
    workdir1 = scratchdir + 'pework'
    
    if not os.path.isdir(workdir1):
        sys.exit(workdir1, "not found")

    workdir2 = workdir1 + '/n_' + str(free_chains[ifree])
	
    if not os.path.isdir(workdir2):
        sys.exit(workdir2, "not found")
        
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
            sys.exit(workdir3,"not found")
                
        os.chdir(workdir3)
        destdir = os.getcwd()

        print("Copying Files")
            
        # * Heart of the Code
        
        resdir = destdir + '/biascalc3'

        if not os.path.exists(resdir):
            print(resdir, "not found")
            continue

        superdir = workdir1 + '/biasall_results'
        superdir_res = maindir + '/bias_restart'

        if not os.path.exists(superdir):
            os.mkdir(superdir)

        if not os.path.exists(superdir_res):
            os.mkdir(superdir_res)

        if archflag:

            f_name = resdir + '/archival_*'
            list_of_files = glob.glob(f_name) 

            if list_of_files:

                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = superdir_res + '/archival_n_' + \
                    str(free_chains[ifree]) + '_' + dirstr + '.lammpstrj' 
                shutil.copy2(latest_fyl, newf_name)

        if enerout_bias:
            
            f_name = resdir + '/jobbias.sh.o*'
            list_of_files = glob.glob(f_name) 

            if list_of_files:
                
                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = superdir + '/ouputBiasEnergy_n_' \
                    + str(free_chains[ifree]) + '_' + dirstr + '.dat' 
                shutil.copy2(latest_fyl, newf_name)
                
        if enerout_nobias:

            f_name = destdir + '/jobmain2.sh.o*'
            list_of_files = glob.glob(f_name) 

            
            if list_of_files:

                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = superdir + '/ouputNoBiasEnergy_n_' \
                    + str(free_chains[ifree]) + '_' + dirstr + '.dat' 
                shutil.copy2(latest_fyl, newf_name)

            else:

                print("did not find biasenergy in ", destdir)

        if adsorb_chain:

            f_name = resdir + '/adsfracchain_rcut_' + rcutoff + '_*'
            list_of_files = glob.glob(f_name) 

            if list_of_files:

                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = superdir + '/bias_adsfracchain_'+ rcutoff \
                    + '_n_'+ str(free_chains[ifree]) + '_' + dirstr + '.dat' 
                shutil.copy2(latest_fyl, newf_name)

            f_name = destdir + '/adsfracchain_rcut_' + rcutoff + '_*'
            list_of_files = glob.glob(f_name) 

            if list_of_files:

                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = superdir + '/Nobias_adsfracchain_'+ rcutoff \
                    + '_n_'+ str(free_chains[ifree]) + '_' + dirstr + '.dat' 
                shutil.copy2(latest_fyl, newf_name)


