import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

#0-initial run  1- production
restart = 1 # For restarting from given configurations

free_chains  = [32]#,80,32,48]
free_mons    = 30
graft_chains = 64
graft_mons   = 30
tail_mons    = 10
nsalt        = 510
f_charge     = 0.5
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
            srcfyl = maindir + '/switch_charge.f90'
            desfyl = destdir + '/switch_charge.f90'
            shutil.copy2(srcfyl, desfyl)
            
            srcfyl = maindir + '/switch_charge_params.f90'
            desfyl = destdir + '/switch_charge_params.f90'
            shutil.copy2(srcfyl, desfyl)
            
            srcfyl = maindir + '/ran_numbers.f90'
            desfyl = destdir + '/ran_numbers.f90'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = maindir + '/switch_var.txt'
            desfyl = destdir + '/switch_var.txt'
            shutil.copy2(srcfyl, desfyl)
            
            srcfyl = lmpdir  + '/in.init_switch_var'
            desfyl = destdir + '/in.init_switch_var'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = lmpdir + '/in.longrun'
            desfyl = destdir + '/in.longrun'
            shutil.copy2(srcfyl, desfyl)
            
            srcfyl = lmpdir + '/jobswitch.sh'
            desfyl = destdir + '/jobswitch.sh'
            shutil.copy2(srcfyl, desfyl)
            
            
            srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            desfyl = destdir + '/lmp_mesabi'
            shutil.copy2(srcfyl, desfyl)
            
            # Manipulate Files
            
            print( "Copy Successful - Manipulating Input Files")
            print( "Current Dir ", os.getcwd())
            
            tot_chains = free_chains[ifree] + graft_chains
            ntotal = free_chains[ifree]*free_mons + \
                     graft_chains*(graft_mons + tail_mons) + \
                     2.0*nsalt+f_charge*free_chains[ifree]*free_mons \
                     + f_charge*graft_chains*graft_mons
            print( "Total number of monomers ", ntotal)
            
            dataname = "PEdata_new_"+str(free_chains[ifree])+"_"+fylstr+".dat"
            newdataname = "PEdata_recreate_"+str(free_chains[ifree])+"_"+fylstr+".dat"

            
            launch_fyl = 'switch.txt'
            fr = open('switch_var.txt','r')
            fw = open(launch_fyl,'w')
            fid = fr.read().replace("py_nchains",str(tot_chains)).\
                  replace("py_ntotal",str(int(ntotal))).\
                  replace("py_dataname",str(dataname)).\
                  replace("py_outfile",str(newdataname))
            fw.write(fid)
            fw.close()
            fr.close()
            
            # Generate Required Files for LAMMPS
            
            print( "Running FORTRAN script for generating New Datafile")
            subprocess.call(["ifort","-r8","-qopenmp","ran_numbers.f90",
                             "switch_charge_params.f90","switch_charge.f90",
                             "-o","swit.o"])
            
            subprocess.call(["./swit.o",launch_fyl])
            
            
		
            launch_fyl = 'in.switch'
            fr = open('in.init_switch_var','r')
            fw = open(launch_fyl,'w')
            fid = fr.read().replace("py_dataname",newdataname)
            fw.write(fid)
            fw.close()
            fr.close()
            
            
            # Removing unnecessary files
            
                
            #files = glob.glob(destdir +'/*.txt')
            #for f in files:
             #   os.remove(f)
            
            files = glob.glob(destdir +'/*var*')
            for f in files:
                os.remove(f)
            
            files = glob.glob(destdir +'/*.mod')
            for f in files:
                os.remove(f)
            
            files = glob.glob(destdir +'/*.o')
            for f in files:
                os.remove(f)

            if not os.path.isfile(newdataname):
                print(newdataname, "not found")
                continue
                      

            print( "Submitting Jobs..")
                
            subprocess.call(["qsub","jobswitch.sh"])
	    
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
	 
