import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

#0-initial run 2-restart equilibration 1- production
restart = 1 # For restarting from given configurations

free_chains  = [64]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
f_charge     = 0.5
archarr      = [3]#[1,2,3,4]

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

           
        if restart == 0:
                
            if not os.path.isdir(workdir3):
                os.mkdir(workdir3)
                
            os.chdir(workdir3)
            destdir = os.getcwd()
            
            print( "Starting from beginning")
            
            print( "Copying Files")
            
            # Copy Files
            srcfyl = maindir + '/lmp_params_var.f90'
            desfyl = destdir + '/lmp_params_var.f90'
            shutil.copy2(srcfyl, desfyl)
            
            srcfyl = maindir + '/lammps_inp.f90'
            desfyl = destdir + '/lammps_inp.f90'
            shutil.copy2(srcfyl, desfyl)
            
            srcfyl = maindir + '/ran_numbers.f90'
            desfyl = destdir + '/ran_numbers.f90'
            shutil.copy2(srcfyl, desfyl)
            
            srcfyl = lmpdir  + '/in.init_var'
            desfyl = destdir + '/in.init_var'
            shutil.copy2(srcfyl, desfyl)
            srcfyl = lmpdir + '/in.run1'
            desfyl = destdir + '/in.run1'
            shutil.copy2(srcfyl, desfyl)
            srcfyl = lmpdir + '/in.longrun'
            desfyl = destdir + '/in.longrun'
            shutil.copy2(srcfyl, desfyl)
            srcfyl = lmpdir + '/jobmain.sh'
            desfyl = destdir + '/jobmain.sh'
            shutil.copy2(srcfyl, desfyl)
            
            
            srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            desfyl = destdir + '/lmp_mesabi'
            shutil.copy2(srcfyl, desfyl)
            
            # Manipulate Files
            
            print( "Copy Successful - Manipulating Input Files")
	    
            launch_fyl = 'lmp_params.f90'
            fr = open('lmp_params_var.f90','r')
            fw = open(launch_fyl,'w')
            fid = fr.read().replace("py_free_chains",str(free_chains[ifree])).\
                replace("py_free_mons",str(free_mons)).\
                replace("py_graft_chains",str(graft_chains)).\
                replace("py_graft_mons",str(graft_mons)).\
                replace("py_tail_mons",str(tail_mons)).\
                replace("py_nsalt",str(nsalt)).\
                replace("py_charge_frac",str(f_charge))
            fw.write(fid)
            fw.close()
            fr.close()
            
            # Generate Required Files for LAMMPS
            
            print( "Running FORTRAN script for generating Datafile")
            subprocess.call(["ifort","-r8","-qopenmp","ran_numbers.f90",
                             "lmp_params.f90","lammps_inp.f90","-o","lmpinp.o"])
            
            subprocess.call(["./lmpinp.o",str(archarr[iarch])])
            
            dataname = "PEdata_"+str(free_chains[ifree])+"_"+fylstr+".dat"
            datapath = destdir + '/' + dataname
            
            
            if os.path.exists(datapath):
		
                print( "Dataname for LAMMPS: ", dataname)
		
            else:
		
                print( dataname, "not found in", datapath)
                sys.exit()
		
            launch_fyl = 'in.init'
            fr = open('in.init_var','r')
            fw = open(launch_fyl,'w')
            fid = fr.read().replace("dataname",dataname)
            fw.write(fid)
            fw.close()
            fr.close()
            
            
            # Removing unnecessary files
            
                
            files = glob.glob(destdir +'/*.txt')
            for f in files:
                os.remove(f)
                files = glob.glob(destdir +'/*var*')
            for f in files:
                os.remove(f)
                files = glob.glob(destdir +'/*.mod')
            for f in files:
                os.remove(f)
                files = glob.glob(destdir +'/*.o')
            for f in files:
                os.remove(f)

            print( "Submitting Jobs..")
                
            subprocess.call(["qsub","jobmain.sh"])
	    
            os.chdir(maindir)

        elif restart == 2:

            print ("Running long equilibrium")

            if not os.path.isdir(workdir3):
                print( workdir3, "not found")
                continue

            os.chdir(workdir3)
            destdir = os.getcwd()

            srcfyl = destdir + '/restart1'
            
            if not os.path.exists(srcfyl):

                print( "Restart file for LAMMPS: ", srcfyl,\
                       "not found")
                continue


            srcfyl = lmpdir + '/in.run2'
            desfyl = destdir + '/in.run2'
            shutil.copy2(srcfyl, desfyl)
                
            srcfyl = lmpdir + '/jobmain3.sh'
            desfyl = destdir + '/jobmain3.sh'
            shutil.copy2(srcfyl, desfyl)
	    
	    		
            srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            desfyl = destdir + '/lmp_mesabi'

            if not os.path.exists(desfyl):
                shutil.copy2(srcfyl, desfyl)

                
            print( "Submitting Jobs..")
	 
            subprocess.call(["qsub","jobmain3.sh"])
            os.chdir(maindir)

        else:

            os.chdir(workdir3)
            destdir = os.getcwd()

            if not os.path.isdir(workdir3):
                print( workdir3, "not found")
                continue

            
            archfiles = destdir + '/archival*'
            list_of_files = glob.glob(archfiles)

            if not list_of_files:
                print("No archival files found in ", destdir)
                continue

            srcfyl = lmpdir + '/in.longrun'
            desfyl = destdir + '/in.longrun'
            shutil.copy2(srcfyl, desfyl)
            
            srcfyl = lmpdir + '/jobmain2.sh'
            desfyl = destdir + '/jobmain2.sh'
            shutil.copy2(srcfyl, desfyl)
            
            fylename = destdir + '/lmp_mesabi'

            if not fylename: 
                srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                desfyl = destdir + '/lmp_mesabi'
                shutil.copy2(srcfyl, desfyl)

            print( "Submitting Jobs..")
            
            subprocess.call(["qsub","jobmain2.sh"])
        
            os.chdir(maindir)
	 
