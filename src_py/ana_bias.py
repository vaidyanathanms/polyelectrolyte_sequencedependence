import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

free_chains  = [32,48,64,72]
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
        ntot_chains = graft_chains + free_chains[ifree]

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
            print("Did not find", workdir3)
            continue

        workdir3 = workdir3 + '/biascalc3'

        if not os.path.isdir(workdir3):
            print("Did not find", workdir3)
            continue
                
        os.chdir(workdir3)
        destdir = os.getcwd()

        print( "Starting from beginning")

        print( "Copying Files")
            
        # Copy Files
        srcfyl = maindir + '/pe_analyze_v2.f90'
        desfyl = destdir + '/pe_analyze.f90'
        shutil.copy2(srcfyl, desfyl)
        
        srcfyl = maindir + '/pe_params.f90'
        desfyl = destdir + '/pe_params.f90'
        shutil.copy2(srcfyl, desfyl)
        
        srcfyl = maindir + '/anainp_var.txt'
        desfyl = destdir + '/anainp_var.txt'
        shutil.copy2(srcfyl, desfyl)

        srcfyl = maindir + '/jobana.sh'
        desfyl = destdir + '/jobana.sh'
        shutil.copy2(srcfyl, desfyl)


        # Manipulate Files
            
        dataname = "databias.data"
        datapath = destdir + '/' + dataname
            
        restartpath = destdir + '/restart1'
            
        if os.path.exists(datapath):
		
            print( "Dataname for LAMMPS: ", dataname)
		
        elif os.path.exists(restartpath):

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
        

        # * means all if need specific format then *.csv
        trajnames = destdir + '/config*.lammpstrj'
        list_of_files = glob.glob(trajnames) 
        traj_str = max(list_of_files, key=os.path.getctime)

        traj_arr = traj_str.split("/")

        print( "trajval", traj_arr[len(traj_arr)-1])

        launch_fyl = 'anainp.txt'
        fr = open('anainp_var.txt','r')
        fw = open(launch_fyl,'w')
        fid = fr.read().replace("py_datafyl",dataname).\
            replace("py_trajfyl",traj_arr[len(traj_arr)-1]).\
            replace("py_nfree",str(free_ntot)).\
            replace("py_ntotchains",str(ntot_chains)).\
            replace("py_ngraft",str(graft_ntot)).\
            replace("py_cutoff",str(py_cutoff)).\
            replace("py_nfrchains",str(free_chains[ifree])).\
            replace("py_ngrchains", str(graft_chains))
        fw.write(fid)
        fw.close()
        fr.close()

        launch_fyl = 'jobana2.sh'
        fr = open('jobana.sh','r')
        fw = open(launch_fyl,'w')
        fid = fr.read().replace("pydir",dirstr).\
            replace("pynfr",str(free_chains[ifree]))
        fw.write(fid)
        fw.close()
        fr.close()

        subprocess.call(["ifort","-r8","-qopenmp","pe_params.f90",
                         "pe_analyze.f90","-o","ana.o"])

        subprocess.call(["qsub","jobana2.sh"])

        os.chdir(maindir)
