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

#--SMD or UmbrellaSampling---------------------------------------

logic_umb = 2 #0 - SMD for initializing US; 1 - US; 2-restart US

#--Chain Details-------------------------------------------------
free_chains  = [100]#,150]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4] #1=bl-bl,2=bl-al,3=al-bl,4=al-al
panfs_scratch = 0 #only for restart. See path for destdir

#--Umbrella Sampling Details-------------------------------------
dirval  = 'US_dir_new3'
#umb_centers = [5.0,8.0,42.0,46.0]
#umb_centers = [32.0, 34.0, 38.0, 39.0, 41.0, 42.0, 44.0,46.0,50.0,53.0]
umb_centers = [3.0,5.0,8.0,10.0, 12.0, 14.0, 16.0, 18.0, 22.0, 26.0,28.0,30.0,32.0, 34.0, 38.0, 39.0, 41.0, 42.0, 44.0,46.0,50.0,53.0]
force_constant = 0.4
#--US details for FORTRAN CODE-----------------------------------
tolval  = 1.0
axisval = 3
comval  = 1


#--SMD Details---------------------------------------------------
targinit = 5.0
targfin  = 45.0

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
            os.mkdir(smdworkdir)
                
        os.chdir(smdworkdir)
        destdir = os.getcwd()
          

#--------Prepare for US Simulations--------------------------------------

        if logic_umb == 1:
            
            print("Preparing for US Simulations ...")
            ref_data = smdworkdir + "/datainp.data"
            
            if not os.path.exists(ref_data):
                print(ref_data) 
                sys.exit("Reference Datafile not found")

            for ucen in range(len(umb_centers)):
                
                umbdir = smdworkdir + '/' + dirval + str(umb_centers[ucen])
            
                if not os.path.isdir(umbdir):
                    os.mkdir(umbdir)
                
                os.chdir(umbdir)
                umbworkdir = os.getcwd()
            
                srcfyl = lmpdir  + '/in.umbsample_var'
                desfyl = umbworkdir + '/in.umbsample_var'
                shutil.copy2(srcfyl, desfyl)
                
                srcfyl = lmpdir  + '/umbinp_var'
                desfyl = umbworkdir + '/umbinp_var'
                shutil.copy2(srcfyl, desfyl)
                
                srcfyl = lmpdir  + '/jobumb.sh'
                desfyl = umbworkdir + '/jobumb.sh'
                shutil.copy2(srcfyl, desfyl)
                
                srcfyl = maindir + '/extract_conf.f90'
                desfyl = smdworkdir + '/extract_conf.f90'
                shutil.copy2(srcfyl,desfyl)

                srcfyl = maindir + '/extract_params.f90'
                desfyl = smdworkdir + '/extract_params.f90'
                shutil.copy2(srcfyl,desfyl)

                srcfyl = maindir + '/extract_inp_var.txt'
                desfyl = smdworkdir + '/extract_inp_var.txt'
                shutil.copy2(srcfyl,desfyl)
                              
                srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                desfyl = umbworkdir + '/lmp_mesabi'
                shutil.copy2(srcfyl, desfyl)


#---Now check whether SMD simulations are okay-----------------------------
#---Call FORTRAN CODE - extract_conf.f90recursively with different
#---centers---

#ref_fin and ref_init are same for US simulations
#obtain the reference chain from SMD simulation inputs and feed here

                os.chdir(smdworkdir)
                fyl_smd = smdworkdir + '/smd_colfile.inp'

                with open(fyl_smd) as fp:
                    for lyne in fp:
                        if "main" in lyne:
                            content = re.findall(r'\d+',lyne)
                            
                ref_init = 1
                pyinit = int(content[0])
                pyfin  = int(content[1])

                init_fyl = 'extract_inp.txt'
                fr = open('extract_inp_var.txt','r')
                fw = open(init_fyl,'w')
                fid = fr.read().replace("py_smddata","datainp.data").\
                    replace("py_refinit",str(ref_init)).\
                    replace("py_reffin", str(ref_init)).\
                    replace("py_init",str(pyinit)).\
                    replace("py_fin",str(pyfin)).\
                    replace("py_savedist",str(umb_centers[ucen])).\
                    replace("py_tol", str(tolval)).\
                    replace("py_axis",str(axisval)).\
                    replace("py_comflag", str(comval))

                fw.write(fid)
                fw.close()
                fr.close()


                subprocess.call(["ifort","-mkl","-qopenmp","-r8",\
                                 "extract_params.f90","extract_conf.f90",\
                                 "-o","extract.o"])
                subprocess.call(["./extract.o","extract_inp.txt"])


                checkinitdata = smdworkdir + "/init_datafile"
                desfyl = umbworkdir + '/init_datafile'

                if os.path.exists(checkinitdata):

                    print(checkinitdata, "prepared for US simulations.")
                    shutil.copy2(checkinitdata,desfyl)
                    
                    #Delete & keep a copy under different name for
                    #future copying if necessary

                    prev_US_fyl = smdworkdir + '/previnit_datafile'
                    shutil.copy2(checkinitdata,prev_US_fyl)
                    os.remove(checkinitdata)
                    
                else: #if US could not find a file near that place -
                    #use SMD initialized file or latest US file
                    
                    print(checkinitdata, "not prepared for", \
                              free_chains[ifree],"at",umb_centers[ucen])

                    prev_US_fyl = smdworkdir + '/previnit_datafile'

                    if os.path.exists(prev_US_fyl):

                        print("Copying previous US file for these purposes")
                        shutil.copy2(prev_US_fyl,desfyl)

                    else:

                        print("Copying SMD file for these purposes")
                        shutil.copy2(ref_data,desfyl)
                        
                #---Copy datafiles for backup before analyzing
                
                datadir = maindir + '/USbackupdir'

                if not os.path.isdir(datadir):
                    os.mkdir(datadir)

                datadir2 = datadir + '/n_' + str(free_chains[ifree])

                if not os.path.isdir(datadir2):
                    os.mkdir(datadir2)

                datafyle = datadir2 + '/n_' + str(free_chains[ifree])\
                    + '_' + dirstr + '_' +  str(umb_centers[ucen])

                shutil.copy2(desfyl,datafyle)
                
                os.chdir(umbworkdir)
                umbmain_fyl = 'umb_colfile.inp'
                fr = open('umbinp_var','r')
                fw = open(umbmain_fyl,'w')
                fid = fr.read().replace("py_init",str(pyinit)).\
                    replace("py_fin",str(pyfin)).\
                    replace("py_cval",str(umb_centers[ucen])).\
                    replace("py_fcon", str(force_constant))
                fw.write(fid)
                fw.close()
                fr.close()
                
                umblmp_fyl = 'in.umbsample'
                fr = open('in.umbsample_var','r')
                fw = open(umblmp_fyl,'w')
                fid = fr.read().replace("py_init",str(pyinit)).\
                    replace("py_fin",str(pyfin)).\
                    replace("ref_init",str(ref_init)).\
                    replace("ref_fin", str(ref_init)) 
                fw.write(fid)
                fw.close()
                fr.close()
                
                print("Submitting US jobs for",free_chains[ifree],"at" \
                          ,umb_centers[ucen])
                subprocess.call(["qsub", "jobumb.sh"])
                
        elif logic_umb == 2:

#----------------Restarting US Simulations-----------------------------
            print("Preparing to US Simulations ...")
            ref_data = smdworkdir + "/datainp.data"
            
            for ucen in range(len(umb_centers)):

                umbdir = smdworkdir + '/' + dirval + str(umb_centers[ucen])
                print(smdworkdir)            
                if not os.path.isdir(umbdir):
                    print(umbdir, "not found")
                    continue
                
                os.chdir(umbdir)
                umbworkdir = os.getcwd()
            
                srcfyl = lmpdir  + '/in.umbre'
                desfyl = umbworkdir + '/in.umbre'
                shutil.copy2(srcfyl, desfyl)
                
                srcfyl = lmpdir  + '/jobre_umb.sh'
                desfyl = umbworkdir + '/jobre_umb.sh'
                shutil.copy2(srcfyl, desfyl)
                
                srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                desfyl = umbworkdir + '/lmp_mesabi'
                shutil.copy2(srcfyl, desfyl)

                rfyle1 = 'restart2'
                rfyle2 = 'out.colvars.state'
                          
                if not os.path.exists(rfyle1):
                    print(rfyle1, "not found in", umbdir)
                    continue

                                
                print("Submitting US jobs for",free_chains[ifree],"at" \
                          ,umb_centers[ucen])
                subprocess.call(["qsub", "jobre_umb.sh"])
                
                          

        elif logic_umb == 0:

#----------------Preparing for SMD Simulations-----------------------------

            print("Preparing for SMD simulations ...")
            print("Copying restart files")
            

            srcfyl = maindir + '/archive_restart'+'/restartfiles_n_' + \
                str(free_chains[ifree]) + '_' +  dirstr + '.tar.gz'

            # Check files are already extracted - this takes the most time
            
            if panfs_scratch == 1:

                archfiles = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                            str(free_chains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                            str(free_chains[ifree]) + '_' + dirstr +  '/' +  \
                            'archival_*.restart'

            else:
                
                archfiles = destdir + "/" + 'archival_*.restart'
            
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

                if panfs_scratch == 1:
                    finarchfyl = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                                 str(free_chains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                                 str(free_chains[ifree]) + '_' + dirstr +  '/' +  \
                                 'archival_' + str(int(oldtime)) + '.restart'

                else:
                    finarchfyl = destdir+"/"+'archival_'+ str(int(oldtime)) + '.restart'
                    
                desfyl = destdir + '/archival_' + str(int(oldtime)) + '_restartsmd'
                shutil.copy2(finarchfyl,desfyl)

                subprocess.call(["mpirun","-np","24", "./lmp_mesabi", "-r",
                                 desfyl, "datainp.data"])


            elif os.path.exists(srcfyl):
		
                print( "Restart file for LAMMPS: ", srcfyl)
                print("Current Dir:", os.getcwd())
                desfyl= destdir+'/restartfiles_n_'+str(free_chains[ifree]) \
                    + '_' +  dirstr + '.tar.gz'
            
                shutil.copy2(srcfyl, desfyl)

# Untar restart files and copy it as restart1
            
                subprocess.call(["tar","-xvzf", desfyl])

                if panfs_scratch == 1:
                    archfiles = destdir + str(free_chains[ifree]) \
                                + "/" +dirstr +'/restartfiles_n_'+ \
                                str(free_chains[ifree])+'_'+dirstr \
                                + '/' + 'archival_*.restart'
                else:
                    archfiles = destdir + "/" + '/restartfiles_n_'+\
                                str(free_chains[ifree]) + '_' +dirstr+\
                                '/' + 'archival_*.restart'
                    

                list_of_files = glob.glob(archfiles)
            
                for fyl in list_of_files:

                    splitdata = fyl.split("/")
                    fylname  = splitdata[len(splitdata)-1]
                    timedata = fylname.split("_")
                    timeval  = float(timedata[1].split(".")[0])
                    
                    if timeval > oldtime:

                        oldtime = timeval

                if panfs_scratch == 1:
                    
                    finarchfyl = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                                 str(free_chains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                                 str(free_chains[ifree]) + '_' + dirstr +  '/' +  \
                                 'archival_' + str(int(oldtime)) + '.restart'

                else:
                    
                    finarchfyl = destdir + "/" + '/restartfiles_n_'+ \
                                 str(free_chains[ifree])  + '_' + \
                                 dirstr  +  '/' + 'archival_' \
                                 + str(int(oldtime)) + '.restart'


                desfyl = destdir + '/archival_' + str(int(oldtime)) + '.restart'
                shutil.copy2(finarchfyl,desfyl)

            else:
		
                print(srcfyl, "not found in", maindir)
                os.chdir(maindir)


            pyinit   = graft_chains*graft_mons+(refchain-1)*free_mons+1
            pyfin    = pyinit + free_mons - 1
            ref_init = 1
            ref_fin  = graft_chains*graft_mons

            srcfyl = lmpdir   + '/in.smd_var'
            desfyl = smdworkdir + '/in.smd_var'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = lmpdir   + '/smdcolfile_var.inp'
            desfyl = smdworkdir + '/smdcolfile_var.inp'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = lmpdir   + '/jobsmd.sh'
            desfyl = smdworkdir + '/jobsmd.sh'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            desfyl = smdworkdir + '/lmp_mesabi'
            shutil.copy2(srcfyl, desfyl)

            smdlmp_fyl = 'in.smd'
            fr = open('in.smd_var','r')
            fw = open(smdlmp_fyl,'w')
            fid = fr.read().replace("py_init",str(pyinit)).\
                replace("py_fin",str(pyfin)).\
                replace("ref_init",str(ref_init)).\
                replace("ref_fin", str(ref_fin))
            fw.write(fid)
            fw.close()
            fr.close()
                        
            smdcol_fyl = 'smd_colfile.inp'
            fr = open('smdcolfile_var.inp','r')
            fw = open(smdcol_fyl,'w')
            fid = fr.read().replace("py_init",str(pyinit)).\
                replace("py_fin",str(pyfin)).\
                replace("ref_init",str(ref_init)).\
                replace("ref_fin", str(ref_fin)).\
                replace("py_targinit",str(targinit)).\
                replace("py_targfin",str(targfin)).\
                replace("py_fcon", str(force_constant))
            fw.write(fid)
            fw.close()
            fr.close()

            print("Submitting SMD jobs for", free_chains[ifree])
            subprocess.call(["qsub", "jobsmd.sh"])
