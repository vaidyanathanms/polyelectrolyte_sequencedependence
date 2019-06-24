#Select latest files of a given type and remove the coordinate from
#file name. For ease in MATLAB postprocessing
import numpy
import os
import shutil
import subprocess
import sys
import glob
import os.path
from subprocess import call

#-------------Keys for copying------------------------------------

copywham     = 0 #copy whamdirectories
outresults   = 1 #copy dens/grpdens/adsfrac

#-------------Input data------------------------------------------

free_chains  = [16,32,48,64,72,80,100]#,200]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4]
rcutoff      = '2.00'

#------------Main Code begins here-------------------------------

maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_pework/src_lmp'


log_file = scratchdir + '/pework/log_file' #-write log file
flog = open(log_file,'w')
flog.write('%s\t %s\t %s\n' %('n_pa', 'arch', 'filename'))   

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
        
        resdir = destdir + '/results_' + str(free_chains[ifree]) + '_' \
            + dirstr

        if not os.path.exists(resdir):
            print(resdir, "not found")
            continue

        superdir = workdir1 + '/results_' + str(free_chains[ifree]) + '_' \
            + dirstr


        if not os.path.exists(superdir):
            os.mkdir(superdir)



        if outresults:

            f_name = destdir + '/adsfracchain_rcut_' + rcutoff+ '_config*'
            list_of_files = glob.glob(f_name) 

            if list_of_files:

                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = resdir + '/adsfracchain_rcut_' + rcutoff+ '.lammpstrj'
                shutil.copy2(latest_fyl, newf_name)
                finf_name = superdir + '/adsfracchain_rcut_' + rcutoff+ '.lammpstrj'
                shutil.copy2(latest_fyl, finf_name)
                flog.write("%d\t%s\t%s\n" %(free_chains[ifree],dirstr,latest_fyl))
                
            f_name = destdir + '/adsfracmon_rcut_' + rcutoff+ '_config*'
            list_of_files = glob.glob(f_name) 
            

            if list_of_files:

                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = resdir + '/adsfracmon_rcut_' + rcutoff+ '.lammpstrj'
                shutil.copy2(latest_fyl, newf_name)
                finf_name = superdir + '/adsfracmon_rcut_' + rcutoff+ '.lammpstrj'
                shutil.copy2(latest_fyl, finf_name)

            else:

                print("Did not find any adsfracmon_rcut files in \t",destdir)

            f_name = resdir + '/dens_*'
            list_of_files = glob.glob(f_name) 

            if list_of_files:

                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = resdir + '/dens.lammpstrj'
                shutil.copy2(latest_fyl, newf_name)
                finf_name = superdir + '/dens.lammpstrj'
                shutil.copy2(latest_fyl, finf_name)

            else:

                print("Did not find any dens files in \t",resdir)

            f_name = resdir + '/grpdens_*'
            list_of_files = glob.glob(f_name) 

            if list_of_files:

                latest_fyl = max(list_of_files, key=os.path.getctime)
                newf_name = resdir + '/grpdens.lammpstrj'
                shutil.copy2(latest_fyl, newf_name)
                finf_name = superdir + '/grpdens.lammpstrj'
                shutil.copy2(latest_fyl, finf_name)

            else:

                print("Did not find any grpdens files in \t",resdir)


        if copywham:

            superwham = workdir1 + '/wham_' + str(free_chains[ifree]) + '_' \
                        + dirstr

            if not os.path.exists(superwham):
                print(superwham, " created")
                os.mkdir(superwham)


            f_name = destdir + '/whaminp.txt'

            if not os.path.isfile(f_name):
                print("Did not find whaminp.txt in ", destdir)
                continue
            

            desf_name = superwham + '/whaminp.txt'
            shutil.copy2(f_name, desf_name)

            f_name = destdir + '/whamout.txt'
            if not os.path.isfile(f_name):
                print("Did not find whamout.txt in ", destdir)
                continue

            print("Copying whamout.txt from ", destdir) 
            desf_name = superwham + '/whamout.txt'
            shutil.copy2(f_name, desf_name)

            
            whamdir = destdir + '/wham_ana'
            if not os.path.exists(whamdir):
                print(whamdir, "not found in", destdir)
                continue

            f_name = whamdir + '/whaminp_*'
            list_of_files = glob.glob(f_name) 

            for ifile in range(len(list_of_files)):

                splice  = list_of_files[ifile].split('/')
                fylename = splice[len(splice)-1]
                srcf_name = whamdir + '/' + fylename

                if not os.path.isfile(srcf_name):
                    print("Did not find", srcf_name)
                    continue

                desf_name = superwham + '/' + fylename
                shutil.copy2(srcf_name, desf_name)

                


