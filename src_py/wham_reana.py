# WHAM Analysis (Recalculate with error bars)
# Use with either COLVARS output or from previously stored output traj 

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

#--------------------------------Inputs-----------------------------

frombackup = 1 # frombackup = 0 (from new simulations), 1 (from backup dir)
arch  = [1,2,3,4] #If frombackup == 0, arch can be any value
nval  = 64
correl_time = 20
skipfr  = 1000
kspr = 0.4

#----------------------------WHAM Analysis--------------------------

#--Some error - not able to figure out with subprocess.call
perform_wham = 1 #1-perform, 0-not perform
minhist = 4.0
maxhist = 53.0
maxbins = 100
maxtol  = 0.001
temperature = 1.0
padding = 0
outfile = "whamout.txt"
maxiter = 1000
stripterm = 'new' #new/none. If new - change before the harmonic cntr
refval_stripterm = 300 #for 'new' keywork ONLY
whamnew_pref = 'whaminp_3' #ONLY required for 'new' keyword

#------------------------------Directories--------------------------

f90src   = '/home/dorfmank/vsethura/allfiles/files_pework/src_f90'
scratchdir = '/scratch.global/vaidya/'
whamexec = '/home/dorfmank/vsethura/wham/wham/wham'

if frombackup == 0:
    maindir = os.getcwd()
else:
    maindir = f90src

#-----------------------------Files-----------------------------------------

whampref = 'whaminp_' #Keep * here
whammain = 'whaminp.txt'
colfile = 'out.colvars.traj'
colfile_name = 'umbsampling'

if frombackup == 0:
    # Keep the prefix without number
    if stripterm:
        prefix  = 'US_dir_' + stripterm 
    else:
        prefix  = 'US_dir_new' 
else:
    prefix = 'USbackup_colvar'

whamcopy = maindir + '/wham'

if not os.path.isfile(whamcopy):
    shutil.copy2(whamexec,whamcopy)

#---------------------------Main Loop Starts--------------------------------

curr_dir = os.getcwd()

print("Starting WHAM Analysis")        
if frombackup == 0:

    rvalarr = []
    wham_dir = curr_dir + '/wham_ana'
    if not os.path.exists(wham_dir):
        os.mkdir(wham_dir)

    dir_names = maindir + '/' + prefix + '*'
    list_of_dirs = glob.glob(dir_names) 

    if list_of_dirs:
        
        whamout = outfile
        fmain = open(whammain,'w')
        for dirval in range(len(list_of_dirs)):
            
            splice  = list_of_dirs[dirval].split('_')
            rval = splice[len(splice)-1]
            if stripterm:
                rval = rval.strip(stripterm)

            rvalarr.append(float(rval))
        
        for dirval in sorted(rvalarr):

            dirname = maindir + '/' + prefix + str(dirval)
            os.chdir(dirname)
            fylname = dirname + '/' + colfile
            print(dirname)

            if os.path.isfile(fylname):
                
                wham_fyl = wham_dir + '/' + 'whaminp_' + str(dirval)
                fr = open(fylname)
                fw = open(wham_fyl,'w')
            
                line = fr.readline()
                header = ' '.join(line.split()).split(' ')

                for nvar in range(len(header)):

                    if header[nvar] == colfile_name:
                        
                        rvar = nvar-1
                    
                    elif header[nvar] == 'step':
                        
                        stepvar = nvar-1 
                    
            
                fr.close()
            
                with open(fylname,'r') as fyl:

                    for sval in range(skipfr):

                        next(fyl)
                        
                    for line in fyl:

                        data = ' '.join(line.split()).split(' ')

                        if not data[rvar].isalpha():

                            fw.write(str(data[stepvar])+'\t'+str(data[rvar])+'\n')
                
                fw.close()
            
                fmain.write(wham_fyl + '\t' +str(dirval) +'\t'+ \
                                str(kspr)+'\t' + str(correl_time)+'\n')

        if perform_wham:
            ranseed = 1223 #int(10000*random.random())   
            whamstring = whamexec+' '+' '+str(minhist)+' '+str(maxhist)\
                +' '+str(maxbins)+' '+str(maxtol)+' '+ str(temperature)+\
                ' '+str(padding)+' '+str(whammain)+' '+str(whamout)+\
                ' '+str(maxiter)+' '+str(ranseed)
            print(whamstring)
            subprocess.call(whamstring)

else:

    for archvals in range(len(arch)):

        wham_dir = scratchdir + 'pework'
        rvalarr = []
        if not os.path.exists(wham_dir):
            print(wham_dir, "not found")
            continue

        if arch[archvals] == 1:
            dirstr = 'bl_bl'
        elif arch[archvals] == 2:
            dirstr = 'bl_al'
        elif arch[archvals] == 3:
            dirstr = 'al_bl'
        elif arch[archvals] == 4:
            dirstr = 'al_al'
        else:
            print("Unknown architecture input", arch[archvals])
            continue

        wham_dir2 = wham_dir + '/wham_' + str(nval) + '_' + str(dirstr) 

        if not os.path.exists(wham_dir2):
            print(wham_dir2, "not found")
            continue
        
        os.chdir(wham_dir2)
        whamcopy = wham_dir2 + '/wham' #WHAM executable
        whamhead = whammain #files containing WHAM paths
        whaminp  = whampref + '*' #individual WHAM files
        whamout  = outfile #output WHAM file


        list_of_files = glob.glob(whaminp) 

        if not os.path.isfile(whamcopy):
            shutil.copy2(whamexec,whamcopy)
        else:
            os.remove(whamcopy)  #remove and copy again
            shutil.copy2(whamexec,whamcopy)

        for fylval in range(len(list_of_files)):
            
            
            splice  = list_of_files[fylval].split('_')
            rval = splice[len(splice)-1]

            if stripterm == 'none':
                rvalarr.append(float(rval))
            elif float(rval) < 100:
                rvalarr.append(float(rval)-refval_stripterm/10)
            else:
                rvalarr.append(float(rval)-refval_stripterm)
            

        fmain = open(whamhead,'w')
        for fylval in sorted(rvalarr):

            if stripterm == 'none':
                fylname = wham_dir2 + '/' + whampref + str(fylval)
            else:
                fylname = wham_dir2 + '/' + whamnew_pref + str(fylval)


            if os.path.isfile(fylname):
                
                fmain.write(fylname + '\t' +str(fylval) +'\t'+ \
                                str(kspr)+'\t' + str(correl_time)+'\n')

        fmain.close()
        print(os.getcwd())
        if perform_wham:

            ranseed = int(10000*random.random())   
            fsub = open('out.txt','w')
            subprocess.call(['./wham',str(minhist),str(maxhist),str(maxbins)\
                             ,str(maxtol),str(temperature),str(padding),\
                             str(whammain),str(whamout),str(maxiter),str(ranseed)]\
                            ,stdout = fsub)

            fsub.close()


