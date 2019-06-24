!---------------To analyze Polyelectrolyte Static Properties---------
!---------------Version 3: Dec-19-2018-------------------------------
!---------------Parameter File: pe_params.f90------------------------
!*******************************************************************

PROGRAM PEMAIN

  USE PARAMETERS_PE
  IMPLICIT NONE

  PRINT *, "Static analysis of polyelectrolyte system .."
  PRINT *, "Starting OMP Threads .."

!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
  PRINT *, "Number of threads are: ", nproc
!$OMP END PARALLEL

  CALL READ_ANA_IP_FILE()
  CALL READ_DATAFILE()
  CALL ANALYZE_TRAJECTORYFILE()
  CALL ALLOUTPUTS()
  CALL DEALLOCATE_ARRAYS()

  PRINT *, "All Calculations Completed Succesfully :)"

END PROGRAM PEMAIN

!--------------------------------------------------------------------


SUBROUTINE READ_ANA_IP_FILE()

  USE PARAMETERS_PE

  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i,j
  CHARACTER(256) :: dumchar

  CALL DEFAULTVALUES()

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"

  logflag = 0

  CALL GETARG(nargs,ana_fname)

  OPEN(unit = anaread,file=trim(ana_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(ana_fname), "not found"
     STOP

  END IF

  DO

     READ(anaread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     IF(dumchar == 'datafile') THEN
        
        READ(anaread,*,iostat=ierr) data_fname

     ELSEIF(dumchar == 'trajectory_file') THEN

        READ(anaread,*,iostat=ierr) traj_fname

     ELSEIF(dumchar == 'nframes') THEN

        READ(anaread,*,iostat=ierr) nframes

     ELSEIF(dumchar == 'skipfr') THEN

        READ(anaread,*,iostat=ierr) skipfr

     ELSEIF(dumchar == 'freqfr') THEN

        READ(anaread,*,iostat=ierr) freqfr
        
     ELSEIF(dumchar == 'nchains') THEN

        READ(anaread,*,iostat=ierr) nchains

     ELSEIF(dumchar == 'nwater') THEN

        READ(anaread,*,iostat=ierr) nwater

     ELSEIF(dumchar == 'atomsperchain') THEN

        READ(anaread,*,iostat=ierr) atperchain

     ELSEIF(dumchar == 'compute_rdf') THEN

        rdfcalc = 1
        READ(anaread,*,iostat=ierr) rdffreq, rmaxbin, rdomcut,npairs
        
        ALLOCATE(pairs_rdf(npairs,3),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate pairs_rdf"
      
        DO i = 1,npairs

           READ(anaread,*,iostat=ierr) pairs_rdf(i,1), pairs_rdf(i,2)

        END DO

     ELSEIF(dumchar == 'compute_ancatrdf') THEN

        acrcalc = 1 !acr=anion_cation_rdf
        READ(anaread,*,iostat=ierr) acrfreq,acrmaxbin,acrdomcut&
             &,ncatgrp,nangrp
        
        ALLOCATE(catgrp(ncatgrp),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate catgrp"
        ALLOCATE(angrp(nangrp),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate angrp"

        READ(anaread,*,iostat=ierr) (catgrp(i),i=1,ncatgrp)
        READ(anaread,*,iostat=ierr) (angrp(i),i=1,nangrp)

     ELSEIF(dumchar == 'compute_rg') THEN

        rgcalc = 1
        READ(anaread,*,iostat=ierr) rgfreq,rgall,rgavg

     ELSEIF(dumchar == 'compute_density') THEN

        denscalc = 1
        READ(anaread,*,iostat=ierr) densfreq, dens_axis, ndentypes&
             &,maxden_bin

        ALLOCATE(dentyp_arr(ndentypes),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate dentyp_arr"

        READ(anaread,*,iostat=ierr) (dentyp_arr(i),i=1,ndentypes)

     ELSEIF(dumchar == 'compute_grpdens') THEN
        
        READ(anaread,*,iostat=ierr) ngroups        

        IF(ngroups /= 0) THEN

           ALLOCATE(grp_data(ngroups),stat = AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate grp_data"
           
           ALLOCATE(dengrp_arr(ngroups,10),stat =&
                & AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate dengrp_arr"
           
           dengrp_arr = -1

           DO j = 1,ngroups

              READ(anaread,*,iostat=ierr) grp_data(j)
              IF(grp_data(j) .GT. 10) STOP "Group exceeded max dim"
              READ(anaread,*,iostat=ierr) (dengrp_arr(j,i),i=1&
                   &,grp_data(j))

           END DO
           
        END IF

     ELSEIF(dumchar == 'compute_fracads') THEN
        
        adscalc = 1
        READ(anaread,*,iostat=ierr) brush_ht,mstart,nfree
        PRINT *, "Number of free molecules: ", nfree
        ALLOCATE(free_mols(nfree),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate free_mols"
        
        DO i = 1,nfree

           free_mols(i) = mstart + i - 1

        END DO

        fracadsavg = 0.0

     ELSEIF(dumchar == 'compute_monfracads') THEN

        IF(chainads) STOP "Chain details before monomer details"
        
        monads = 1
        READ(anaread,*,iostat=ierr) nfreegrp, nfreemons,nadsgrp&
             &,nadsmons,adscut
        
        ALLOCATE(free_ptr(nfreegrp),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate free_ptr"
        ALLOCATE(ads_ptr(nadsgrp),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate ads_ptr"
        ALLOCATE(free_grp(nfreemons),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate free_grp"
        ALLOCATE(ads_grp(nadsmons),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate ads_grp"
        
        READ(anaread,*,iostat=ierr) (free_ptr(i),i=1,nfreegrp)
        READ(anaread,*,iostat=ierr) (ads_ptr(i),i=1,nadsgrp)


     ELSEIF(dumchar == 'compute_chainfracads') THEN

        chainads = 1
        READ(anaread,*,iostat=ierr) nfreechains, ngraftchains,&
             & chadscut

        ALLOCATE(free_chainarr(nfreechains),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate free_chainarr"

        IF(monads == 0) THEN !If it is already there, use those
           ! details

           READ(anaread,*,iostat=ierr) nfreegrp, nfreemons,nadsgrp&
                &,nadsmons,adscut
           
           ALLOCATE(free_ptr(nfreegrp),stat = AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate free_ptr"
           ALLOCATE(ads_ptr(nadsgrp),stat = AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate ads_ptr"
           ALLOCATE(free_grp(nfreemons),stat = AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate free_grp"
           ALLOCATE(ads_grp(nadsmons),stat = AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate ads_grp"
        
           READ(anaread,*,iostat=ierr) (free_ptr(i),i=1,nfreegrp)
           READ(anaread,*,iostat=ierr) (ads_ptr(i),i=1,nadsgrp)

        END IF
  
     ELSEIF(dumchar == 'avg_graft_dist') THEN

        READ(anaread,*,iostat=ierr) graft_type
        avg_rgraft_calc = 1

     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword: ", trim(dumchar)
        STOP

     END IF

  END DO

  IF(logflag == 0) log_fname = "log."//trim(adjustl(traj_fname))
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  PRINT *, "Analysis input file read finished .."

END SUBROUTINE READ_ANA_IP_FILE

!--------------------------------------------------------------------

SUBROUTINE DEFAULTVALUES()

  USE PARAMETERS_PE
  IMPLICIT NONE

  ! Frame, Molecules and Processor Details
  nframes = 0; skipfr = 0; freqfr = 0; nfrcntr = 0
  nwater = 0; nchains = 0
  atperchain = 0

  ! Initialize Flags
  rgall = 0; rgcalc = 0;rdfcalc = 0; denscalc = 0; adscalc = 0;monads=0
  chainads = 0; avg_rgraft_calc = 0

  ! Initialize distributions and frequencies
  rdffreq = 0; rgfreq = 0; densfreq = 0; dens_axis = 0

  ! Initialzie Extra Structural Quantities
  rdomcut = 0;  rmaxbin = 0; ndentypes = 0; ngroups = 0
  ngraftmons = 0

  !Averages
  rvolavg = 0; rgavg = 0; avg_ch_adscnt = 0; acrvolavg = 0

END SUBROUTINE DEFAULTVALUES

!--------------------------------------------------------------------

SUBROUTINE READ_DATAFILE()

  USE PARAMETERS_PE

  IMPLICIT NONE

  INTEGER :: i,j,k,ierr,u,AllocateStatus,imax
  INTEGER :: flag, cntr, nwords
  INTEGER :: aid,molid,atype,ix,iy,iz
  REAL    :: charge,rx,ry,rz
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(256) :: rline,dumchar

  CALL COMPUTE_INIT_NLINES(imax)

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"

  WRITE(logout,*) "Datafile used is :", trim(adjustl(data_fname))

  ntotatoms = 0;ntotbonds=0;ntotangls=0;ntotdihds=0;ntotimprs=0
  atomflag =0;velflag = 0;bondflag=0;anglflag=0;dihdflag=0;imprflag=0

  READ(inpread,*)
  READ(inpread,*)

  DO i = 1,imax-2 !Change here according to convenience
       
     READ(inpread,*) u, dumchar
     
        IF(dumchar == "atoms") THEN
           ntotatoms = u
        ELSEIF(dumchar == "bonds") THEN
           ntotbonds = u
        ELSEIF(dumchar == "angles") THEN
           ntotangls = u
        ELSEIF(dumchar == "dihedrals") THEN
           ntotdihds = u
        ELSEIF(dumchar == "atom" .OR. dumchar == "atomtypes") THEN
           ntotatomtypes = u
        ELSEIF(dumchar == "bond" .OR. dumchar == "bondtypes") THEN
           ntotbondtypes = u
        ELSEIF(dumchar == "angle" .OR. dumchar == "atomtypes") THEN
           ntotangltypes = u
        ELSEIF(dumchar == "dihedral" .OR. dumchar == "dihedraltypes") THEN
           ntotdihdtypes = u
        ELSEIF(dumchar == "improper" .OR. dumchar == "impropertypes") THEN
           ntotimprtypes = u
        ELSEIF(dumchar == "Masses") THEN
           
           ALLOCATE(masses(ntotatomtypes,1),stat = AllocateStatus)
           IF(AllocateStatus/=0) STOP "did not allocate masses"
           
           DO j = 1,ntotatomtypes
              
              READ(inpread,*) u, masses(u,1)
              
           END DO
           
        END IF
        
  END DO

  READ(inpread,*)
  READ(inpread,*) xlo, xhi
  READ(inpread,*) ylo, yhi
  READ(inpread,*) zlo, zhi
  
  box_xl = xhi - xlo
  box_yl = yhi - ylo
  box_zl = zhi - zlo

  PRINT *, "x-box  ", "y-box  ", "z-box  "
  PRINT *, box_xl, box_yl, box_zl

  PRINT *, "STATISTICS"
  PRINT *, "Number of atoms/atomtypes: " , ntotatoms,ntotatomtypes
  PRINT *, "Number of bonds/bondtypes: " , ntotbonds,ntotbondtypes
  PRINT *, "Number of angles/angletypes: " , ntotangls,ntotangltypes
  PRINT *, "Number of diheds/dihedtypes: " , ntotdihds,ntotdihdtypes
  flag = 0; cntr = 0

  CALL ALLOCATE_ARRAYS()

  DO 

     READ(inpread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     !READ DATA HERE FOR CHARGES AND MOLID
     !READ EVERYTHING AND OVERWRITE LATER
     IF(trim(dumchar) == "Atoms") THEN
             
        atomflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms

           READ(inpread,*) aid,molid,atype,charge,rx,ry,rz

           rx = rx - xlo
           ry = ry - ylo
           rz = rz - zlo

           aidvals(aid,1)     = aid
           aidvals(aid,2)     = molid
! To take care of the stupidity in input file
           IF (molid > nchains) THEN
              IF (charge ==  1) aidvals(aid,3) = 7
              IF (charge == -1) aidvals(aid,3) = 8
           ELSE
              aidvals(aid,3)     = atype
           END IF
           charge_lmp(aid,1)  = charge
           rxyz_lmp(aid,1)    = rx
           rxyz_lmp(aid,2)    = ry
           rxyz_lmp(aid,3)    = rz

        END DO

     END IF

     IF(trim(dumchar) == "Masses") THEN

        ALLOCATE(masses(ntotatomtypes,1),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate masses"
         
        DO j = 1,ntotatomtypes

           READ(inpread,*) u, masses(u,1)

        END DO

     END IF

     IF(trim(dumchar) == "Velocities") THEN
             
        velflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms

           READ(inpread,*) vel_xyz(j,1),vel_xyz(j,2),vel_xyz(j,3)&
                &,vel_xyz(j,4)

        END DO


     END IF

     IF(trim(dumchar) == "Bonds") THEN
             
        bondflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotbonds

           READ(inpread,*) bond_lmp(j,1),bond_lmp(j,2),bond_lmp(j,3)&
                &,bond_lmp(j,4)

        END DO

     END IF

     IF(trim(dumchar) == "Angles") THEN
             
        anglflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotangls

           READ(inpread,*) angl_lmp(j,1),angl_lmp(j,2),angl_lmp(j,3)&
                &,angl_lmp(j,4),angl_lmp(j,5)

        END DO

     END IF

     IF(trim(dumchar) == "Dihedrals") THEN
             
        dihdflag = 1
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ntotdihds

           READ(inpread,*) dihd_lmp(j,1),dihd_lmp(j,2),dihd_lmp(j,3)&
                &,dihd_lmp(j,4),dihd_lmp(j,5),dihd_lmp(j,6)

        END DO

     END IF
  
     IF(trim(dumchar) == "Impropers") THEN
             
        imprflag = 1
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ntotimprs

           READ(inpread,*) impr_lmp(j,1),impr_lmp(j,2),impr_lmp(j,3)&
                &,impr_lmp(j,4),impr_lmp(j,5),impr_lmp(j,6)

        END DO

     END IF

  END DO
  
  PRINT *, "Fileread finish .."


END SUBROUTINE READ_DATAFILE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_INIT_NLINES(imax)

  USE PARAMETERS_PE

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: imax
  INTEGER :: init, pos, ipos,u,nwords,lcnt,ierr
  CHARACTER(LEN=120) :: charline

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"
  
  lcnt = 0

  READ(inpread,*)

  DO 

     READ(inpread,'(A)',iostat=ierr) charline     

     lcnt = lcnt + 1
     pos = 1
     nwords = 0

     DO

        ipos = VERIFY(charline(pos:),' ')
        IF(ipos == 0) EXIT
        nwords = nwords + 1
        pos = pos + ipos - 1
        ipos = SCAN(charline(pos:),' ')
        IF(ipos == 0) EXIT
        pos = pos + ipos - 1
        
     END DO

     IF(nwords .GE. 4) THEN

        imax = lcnt - 1
        EXIT
        
     END IF

  END DO

  CLOSE(inpread)

END SUBROUTINE COMPUTE_INIT_NLINES

!--------------------------------------------------------------------

SUBROUTINE ANALYZE_TRAJECTORYFILE()

  USE PARAMETERS_PE

  IMPLICIT NONE

  INTEGER :: i,j,aid,ierr,atchk,atype,jumpfr
  REAL :: xlo,xhi,ylo,yhi,zlo,zhi

  OPEN(unit = 15,file =trim(traj_fname),action="read",status="old"&
       &,iostat=ierr)

  IF(ierr /= 0) STOP "trajectory file not found"

  PRINT *, "trajectory file used is :",trim(adjustl(traj_fname))
  WRITE(logout,*) "trajectory file used is :"&
       &,trim(adjustl(traj_fname))


  CALL STRUCT_INIT()
  CALL OPEN_STRUCT_OUTPUT_FILES()

  DO i = 1,skipfr

     DO j = 1,ntotatoms+9

        READ(15,*) 

     END DO

     IF(mod(i,100) == 0) PRINT *, "Skipped ", i, "frames"

  END DO

  DO i = 1,nframes

     nfrcntr = nfrcntr + 1
     IF(mod(i,100) == 0) PRINT *, "Processing ", i+1,"th frame"

     READ(15,*)
     READ(15,*) timestep

     READ(15,*) 
     READ(15,*) atchk
!!$     IF(atchk /= ntotatoms) STOP "Number of atoms do not match"

     READ(15,*) 
     READ(15,*) xlo, xhi
     READ(15,*) ylo, yhi
     READ(15,*) zlo, zhi

     READ(15,*)

     box_xl = xhi - xlo
     box_yl = yhi - ylo
     box_zl = zhi - zlo
     
     boxx_arr(i)  = box_xl
     boxy_arr(i)  = box_yl
     boxz_arr(i)  = box_zl

     DO j = 1,atchk

        READ(15,*) aid,atype,rxyz_lmp(aid,1),rxyz_lmp(aid,2)&
             &,rxyz_lmp(aid,3)

! To take care of the stupidity in input file
        
        IF (aidvals(aid,2) > nchains) THEN
           IF (charge_lmp(aid,1) ==  1) atype = 7
           IF (charge_lmp(aid,1) == -1) atype = 8
        END IF

        IF(atype .NE. aidvals(aid,3)) THEN

           PRINT *, "Incorrect atom ids"
           PRINT *, i,j,aid,atype,aidvals(aid,3)
           STOP

        END IF

     END DO

     DO j = 1,atchk

        rxyz_lmp(j,1) = rxyz_lmp(j,1) - xlo
        rxyz_lmp(j,2) = rxyz_lmp(j,2) - ylo
        rxyz_lmp(j,3) = rxyz_lmp(j,3) - zlo
        
     END DO

     CALL STRUCT_MAIN(nfrcntr)
     
     DO jumpfr = 1,freqfr

        READ(15,*)
        READ(15,*)        
        READ(15,*)
 
        READ(15,*) atchk

        DO j = 1,atchk+5

           READ(15,*) 

        END DO

     END DO

  END DO

  CLOSE(15)

END SUBROUTINE ANALYZE_TRAJECTORYFILE

!--------------------------------------------------------------------

SUBROUTINE OPEN_STRUCT_OUTPUT_FILES()

  USE PARAMETERS_PE

  IMPLICIT NONE

  CHARACTER(LEN=4) :: rcutchar
  
  WRITE(rcutchar,'(F4.2)') adscut 

  IF(rgcalc) THEN
     
     IF(rgavg) THEN
        dum_fname = "rgavg_"//trim(adjustl(traj_fname))
        OPEN(unit = rgavgwrite,file =trim(dum_fname),action="write"&
             &,status="replace")
     END IF

     IF(rgall) THEN
        dum_fname = "rgall_"//trim(adjustl(traj_fname))
        OPEN(unit = rgwrite,file =trim(dum_fname),action="write"&
             &,status="replace")
     END IF
 
  END IF

  IF(adscalc) THEN
     
     dum_fname = "freeads_"//trim(adjustl(traj_fname))
     OPEN(unit =adswrite,file =trim(dum_fname),action="write"&
          &,status="replace")
  END IF

  IF(monads) THEN
     
     dum_fname = "adsfracmon_rcut_"//rcutchar//"_"&
          &//trim(adjustl(traj_fname))
     OPEN(unit =adsmonwrite,file =trim(dum_fname),action="write"&
          &,status="replace")
  END IF

  IF(chainads) THEN
     dum_fname = "adsfracchain_rcut_"//rcutchar//"_"&
          &//trim(adjustl(traj_fname))
     OPEN(unit =adschwrite,file =trim(dum_fname),action="write"&
          &,status="replace")

     dum_fname = "chainadsval_rcut_"//rcutchar//"_"&
          &//trim(adjustl(traj_fname))
     OPEN(unit =adschwrite2,file =trim(dum_fname),action="write"&
          &,status="replace")

  END IF

  IF(avg_rgraft_calc) THEN

     dum_fname = "tetheravg.dat"
     OPEN(unit =tethwrite,file =trim(dum_fname),action="write"&
          &,status="replace")

  END IF

END SUBROUTINE OPEN_STRUCT_OUTPUT_FILES

!--------------------------------------------------------------------

SUBROUTINE STRUCT_INIT()

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER :: i,j,t1,t2,norm,acnt,fcnt,a1id,molid,flagch,flagpr,jmax
  INTEGER :: AllocateStatus

  IF(rdfcalc) THEN

     rdfarray = 0.0
     rbinval = rdomcut/REAL(rmaxbin)

     DO i = 1, npairs

        t1 = 0; t2 = 0

        DO j = 1,ntotatoms

           IF(aidvals(j,3) == pairs_rdf(i,1)) t1 = t1+1
           IF(aidvals(j,3) == pairs_rdf(i,2)) t2 = t2+1
           
        END DO

        pairs_rdf(i,3) = t1*t2

     END DO

  END IF

  IF(acrcalc) THEN

     acrdfarray = 0.0; ncations=0; nanions = 0
     acrbinval = acrdomcut/REAL(acrmaxbin)
     
     DO i = 1, ntotatoms

        DO j = 1,ncatgrp

           IF(aidvals(i,3) == catgrp(j)) ncations=ncations+1

        END DO

        DO j = 1,nangrp

           IF(aidvals(i,3) == angrp(j)) nanions=nanions+1
           
        END DO

     END DO

     PRINT *, "Total Number of polyanions: ", nanions
     PRINT *, "Total Number of polycations: ", ncations

  END IF



  IF(denscalc) THEN

     avgpolyarr = 0.0

     DO i = 0,maxden_bin-1
        
        DO j = 1,ndentypes

           densarray(i,j) = 0.0
           ch_densarray(i,j) = 0.0

        END DO

     END DO
     
     normdens = 0.0; denbinavg = 0.0
     

     IF(ngroups) THEN

        DO i = 0,maxden_bin-1

           DO j = 1,ngroups

              grparray(i,j) = 0.0
              ch_grparray(i,j) = 0.0
              
           END DO

        END DO
        
     END IF

  END IF

  IF(monads) THEN

     fcnt = 0; acnt = 0; avgadscnt = 0
     DO i = 1,ntotatoms

        DO j = 1,nfreegrp

           IF(aidvals(i,3) == free_ptr(j)) THEN

              fcnt = fcnt + 1
              free_grp(fcnt) = aidvals(i,1)
             
           END IF

        END DO

        DO j = 1,nadsgrp

           IF(aidvals(i,3) == ads_ptr(j)) THEN

              acnt = acnt + 1
              ads_grp(acnt) = aidvals(i,1)

           END IF

        END DO
        
     END DO

     IF(acnt .NE. nadsmons .OR. nfreemons .NE. fcnt) THEN
        
        PRINT *, "Unequal number of free/ads mons"
        PRINT *, acnt, nadsmons, fcnt, nfreemons
        STOP
        
     END IF

  END IF

  IF(chainads) THEN

     IF(monads == 0) THEN

        fcnt = 0; acnt = 0; avgadscnt = 0
        DO i = 1,ntotatoms
           
           DO j = 1,nfreegrp
              
              IF(aidvals(i,3) == free_ptr(j)) THEN
                 
                 fcnt = fcnt + 1
                 free_grp(fcnt) = aidvals(i,1)
                 
              END IF

           END DO

           DO j = 1,nadsgrp
              
              IF(aidvals(i,3) == ads_ptr(j)) THEN
                 
                 acnt = acnt + 1
                 ads_grp(acnt) = aidvals(i,1)
                 
              END IF
              
           END DO
        
        END DO

        IF(acnt .NE. nadsmons .OR. nfreemons .NE. fcnt) THEN
        
           PRINT *, "Unequal number of free/ads mons"
           PRINT *, acnt, nadsmons, fcnt, nfreemons
           STOP
           
        END IF

     ELSE !Now we already know the number and type of adsorbed/free
        !monomers. So just need to cross check whether the number of
        ! free/adsorbed chains are correct

        free_chainarr = -1
        fcnt = 0; 
        DO i = 1,nfreemons

           a1id = free_grp(i)
           molid = aidvals(a1id,2)
           j = 1; flagpr = 0

           ! Find the index to which free_chainarr is filled
           DO WHILE(j .LE. nfreechains)
              
              IF(free_chainarr(j) == -1) THEN
                 
                 jmax = j
                 flagpr = 1
                 j = nfreechains + 1
                 
              ELSE
                 
                 j = j + 1
                 

              END IF
              
           END DO

           flagch = 0
           IF(flagpr == 1) THEN !Array is not fully filled. Now find
              ! if the molid is already present 

              j = 1
              DO WHILE(j .LE. jmax)
                 
                 IF(molid == free_chainarr(j)) THEN
                    
                    flagch = 1
                    j = jmax+1

                 ELSE
                    
                    j = j + 1
                    
                 END IF

              END DO
              
           END IF
 
           IF(flagpr == 1 .AND. flagch == 0) THEN !Molid not found in
              ! the filled array list
              fcnt = fcnt + 1
              free_chainarr(fcnt) = molid

           END IF

        END DO

        WRITE(logout,*) "Free chain details: ", free_chainarr

        IF(fcnt .NE. nfreechains) THEN

           PRINT *, "Unequal number of free chains: "
           PRINT *, fcnt, nfreechains

        END IF

     END IF

  END IF


  IF(avg_rgraft_calc) THEN


     DO i = 1, ntotatoms

        IF(aidvals(i,3) == graft_type) ngraftmons=ngraftmons+1

     END DO

     ALLOCATE(graft_mon_ids(ngraftmons),stat=AllocateStatus)
     IF(AllocateStatus/=0) STOP "graft_mon_ids not allocated"

     j = 0
     DO i = 1, ntotatoms
        
        IF(aidvals(i,3) == graft_type) THEN
   
           j = j + 1          
           graft_mon_ids(j) = aidvals(i,1)
   
        END IF

     END DO

     IF(j .NE. ngraftmons) THEN
        
        PRINT *, "ERROR: Unequal graftmons", j, ngraftmons
        STOP

     END IF

  END IF

END SUBROUTINE STRUCT_INIT

!--------------------------------------------------------------------

SUBROUTINE STRUCT_MAIN(tval)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER, INTENT(IN):: tval

  IF(rgcalc .AND. mod(tval-1,rgfreq)==0) CALL COMPUTE_RADGYR(tval)
  IF(rdfcalc .AND. mod(tval-1,rdffreq)==0) CALL COMPUTE_RDF(tval)
  IF(acrcalc .AND. mod(tval-1,acrfreq)==0) CALL COMPUTE_ANCATRDF(tval)
  IF(denscalc .AND. mod(tval-1,densfreq)==0) CALL COMPUTE_DENS(tval)
  IF(adscalc) CALL COMPUTE_ADSORBEDFREE(tval)
  IF(monads) CALL COMPUTE_FREEPENETRATE_MONS(tval)
  IF(chainads) CALL COMPUTE_FREEPENETRATE_CHAINS(tval)
  IF(avg_rgraft_calc) CALL COMPUTE_INTERGRAFT_DISTANCE(tval)

END SUBROUTINE STRUCT_MAIN

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RDF(iframe)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,paircnt,AllocateStatus
  REAL :: rxval,ryval,rzval,rval
  INTEGER :: a1ref,a2ref
  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: dumrdfarray

  rvolval = box_xl*box_yl*box_zl
  rvolavg = rvolavg + rvolval 

  ALLOCATE(dumrdfarray(0:rmaxbin-1,npairs),stat=AllocateStatus)
  IF(AllocateStatus/=0) STOP "dumrdfarray not allocated"
  dumrdfarray = 0


!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval&
!$OMP& ,ibin,paircnt,a1ref,a2ref) REDUCTION(+:dumrdfarray)
  DO paircnt = 1,npairs

     a1ref = pairs_rdf(paircnt,1); a2ref = pairs_rdf(paircnt,2)
     
     DO i = 1,ntotatoms
        
        a1id   = aidvals(i,1)     
        a1type = aidvals(i,3)
        
        DO j = 1,ntotatoms

           a2id   = aidvals(j,1)        
           a2type = aidvals(j,3)

           IF(a1type == a1ref .AND. a2type == a2ref) THEN        

              rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
              ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
              rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
              
              rxval = rxval - box_xl*ANINT(rxval/box_xl)
              ryval = ryval - box_yl*ANINT(ryval/box_yl)
              rzval = rzval ! No periodicity
              
              rval = sqrt(rxval**2 + ryval**2 + rzval**2)
              ibin = FLOOR(rval/rbinval)
              
              IF(ibin .LT. rmaxbin) THEN
              
                 dumrdfarray(ibin,paircnt) = dumrdfarray(ibin&
                      &,paircnt) + 1
              
              END IF
              
           END IF
           
        END DO

     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(i,j)
  DO j = 1,npairs
     
     DO i = 0,rmaxbin-1

        rdfarray(i,j) = rdfarray(i,j) + REAL(dumrdfarray(i,j))&
             &*rvolval/(REAL(2.0*pairs_rdf(j,3)))
        
     END DO
     
  END DO
!$OMP END DO

!$OMP END PARALLEL

  DEALLOCATE(dumrdfarray)

END SUBROUTINE COMPUTE_RDF

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_ANCATRDF(iframe)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,a1type,a2type,ibin,a1id,a2id,AllocateStatus
  INTEGER :: catflag,catcnt,anflag,ancnt
  REAL :: rxval,ryval,rzval,rval,acrvolval
  INTEGER :: a1ref,a2ref
  INTEGER,ALLOCATABLE, DIMENSION(:) :: dumrdfarray

  acrvolval = box_xl*box_yl*box_zl
  acrvolavg = acrvolavg + acrvolval 

  ALLOCATE(dumrdfarray(0:acrmaxbin-1),stat=AllocateStatus)
  IF(AllocateStatus/=0) STOP "dumrdfarray not allocated"

  dumrdfarray = 0


!$OMP PARALLEL 
!$OMP DO PRIVATE(i,j,a1type,a2type,a1id,a2id,rval,rxval,ryval,rzval&
!$OMP& ,ibin,a1ref,a2ref,catflag,catcnt,anflag,ancnt) REDUCTION(+:dumrdfarray)
  DO i = 1,ntotatoms

     a1id   = aidvals(i,1)     
     a1type = aidvals(i,3)

     DO j = 1,ntotatoms

        a2id   = aidvals(j,1)        
        a2type = aidvals(j,3)

        catcnt = 1; catflag = -1
        ancnt  = 1; anflag  = -1

        DO WHILE(catcnt .LE. ncatgrp)
           
           IF(catgrp(catcnt) == a1type) THEN
              
              catflag = 1
              EXIT
              
           ELSE

              catcnt = catcnt + 1

           END IF
           
        END DO
        
        DO WHILE(ancnt .LE. nangrp)
           
           IF(angrp(ancnt) == a2type) THEN
              
              anflag = 1
              EXIT
              
           ELSE

              ancnt = ancnt + 1

           END IF
           
        END DO
             
        
        IF(anflag == 1 .AND. catflag == 1) THEN        

           rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
           ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
           rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval ! No periodicity
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           ibin = FLOOR(rval/acrbinval)
           
           IF(ibin .LT. acrmaxbin) THEN
              
              dumrdfarray(ibin) = dumrdfarray(ibin) + 1

           END IF
           
        END IF
           
     END DO

  END DO
!$OMP END DO

!$OMP DO PRIVATE(i)
  DO i = 0,acrmaxbin-1

     acrdfarray(i) = acrdfarray(i) + REAL(dumrdfarray(i))&
          &*acrvolval/(REAL(nanions*ncations))
        
  END DO
     
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE COMPUTE_ANCATRDF

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RADGYR(iframe)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER :: i,j,molid,atype
  REAL, DIMENSION(1:nchains) :: rgxx, rgyy, rgzz, rgsq
  REAL, DIMENSION(1:nchains) :: rxcm, rycm, rzcm, totmass
  REAL :: rgsqavg, rgxxavg, rgyyavg, rgzzavg
  INTEGER, INTENT(IN) :: iframe

  rgxx = 0.0; rgyy =0.0; rgzz = 0.0; rgsq = 0.0
  totmass = 0.0
  rgsqavg = 0.0; rgxxavg = 0.0; rgyyavg = 0.0; rgzzavg = 0.0
  
  IF(iframe == 1) PRINT *, "Atoms/molecule: ", atperchain
  IF(iframe == 1) PRINT *, masses

  DO i = 1,nchains*atperchain

     molid = aidvals(i,2)
     atype = aidvals(i,3)
     totmass(molid) = totmass(molid) + masses(atype,1)

     rxcm(molid) = rxcm(molid)+ rxyz_lmp(i,1)*masses(atype,1)
     rycm(molid) = rycm(molid)+ rxyz_lmp(i,2)*masses(atype,1)
     rzcm(molid) = rzcm(molid)+ rxyz_lmp(i,3)*masses(atype,1)

  END DO

  DO i = 1,nchains

     rxcm(i) = rxcm(i)/totmass(i)
     rycm(i) = rycm(i)/totmass(i)
     rzcm(i) = rzcm(i)/totmass(i)

  END DO

  
  IF(iframe == 1) THEN

     OPEN(unit = 98,file ="totmasschk.txt",action="write",status="repl&
          &ace")

     DO i = 1,nchains

        WRITE(98,'(I0,1X,4(F14.9,1X))') i, totmass(i),rxcm(i),&
             & rycm(i), rzcm(i)

     END DO

     CLOSE(98)

     OPEN(unit = 98,file ="molidchk.txt",action="write",status="repl&
          &ace")

     DO i = 1,ntotatoms

        WRITE(98,'(I0,1X,I0)') i, aidvals(i,2)

     END DO

     CLOSE(98)

  END IF
  
  DO i = 1,nchains*atperchain

     molid = aidvals(i,2)
     atype = aidvals(i,3)

     rgxx(molid) = rgxx(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
          &-rxcm(molid))**2)
     rgyy(molid) = rgyy(molid) + masses(atype,1)*((rxyz_lmp(i,2)&
          &-rycm(molid))**2)
     rgzz(molid) = rgzz(molid) + masses(atype,1)*((rxyz_lmp(i,3)&
          &-rzcm(molid))**2)

     rgsq(molid) = rgsq(molid) + masses(atype,1)*((rxyz_lmp(i,1)&
          &-rxcm(molid))**2 + (rxyz_lmp(i,2)-rycm(molid))**2 +&
          & (rxyz_lmp(i,3)-rzcm(molid))**2)

  END DO

  DO i = 1,nchains

     rgsq(i) = rgsq(i)/totmass(i)
     rgxx(i) = rgxx(i)/totmass(i)
     rgyy(i) = rgyy(i)/totmass(i)
     rgzz(i) = rgzz(i)/totmass(i)

  END DO


  DO i = 1,nchains

     rgsqavg = rgsqavg + rgsq(i)
     rgxxavg = rgxxavg + rgxx(i)
     rgyyavg = rgyyavg + rgyy(i)
     rgzzavg = rgzzavg + rgzz(i)
     
  END DO
  
  rgsqavg = rgsqavg/REAL(nchains)
  rgxxavg = rgxxavg/REAL(nchains)
  rgyyavg = rgyyavg/REAL(nchains)
  rgzzavg = rgzzavg/REAL(nchains)
  

  IF(rgavg) THEN

     WRITE(rgavgwrite,'(I0,1X,4(F14.6,1X))') timestep, sqrt(rgsqavg),&
          & sqrt(rgxxavg), sqrt(rgyyavg), sqrt(rgzzavg)

  END IF

  IF(rgall) THEN
     
     WRITE(rgwrite,'(2(I0,1X))') timestep, nchains

     DO i = 1,nchains
     
        WRITE(rgwrite,'(I0,1X,4(F14.6,1X))') i,rgxx(i),rgyy(i),&
             & rgzz(i),sqrt(rgsq(i))

     END DO

  END IF
     
END SUBROUTINE COMPUTE_RADGYR

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_DENS(iframe)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,k, arrindx, dflag, binval,sval
  REAL :: boxdir, rval, rbin
  INTEGER,DIMENSION(0:maxden_bin-1,ndentypes):: inst_array,ch_instarray
  INTEGER,DIMENSION(0:maxden_bin-1,ngroups) :: grp_inst_array,chgrpinst
  INTEGER,DIMENSION(0:maxden_bin-1,2):: polydensarr
  IF(dens_axis == 1) THEN
     boxdir = box_xl
  ELSEIF(dens_axis == 2) THEN
     boxdir = box_yl
  ELSEIF(dens_axis == 3) THEN
     boxdir = box_zl
  ELSE
     STOP "Unknown box direction"
  END IF

  rbin = boxdir/REAL(maxden_bin)
  normdens = normdens + (boxdir/(box_xl*box_yl*box_zl*rbin))
  denbinavg = denbinavg + rbin

  DO i = 0, maxden_bin-1
     
     DO j = 1,ndentypes

        inst_array(i,j) = 0.0
        ch_instarray(i,j) = 0.0

     END DO

     DO j = 1,ngroups

        grp_inst_array(i,j) = 0.0
        chgrpinst(i,j) = 0.0

     END DO

  END DO

  polydensarr = 0.0

  DO i = 1,ntotatoms

     dflag = -1;arrindx = -1
     DO j = 1,ndentypes
        
        dflag = -1
        IF(aidvals(i,3) == dentyp_arr(j)) THEN
           
           arrindx = j
           dflag = 1
           EXIT
           
        END IF

     END DO
     
     IF(dflag == 1) THEN

        IF(dens_axis == 1) THEN
           rval = rxyz_lmp(i,1) - box_xl*FLOOR(rxyz_lmp(i,1)/box_xl)
           binval = FLOOR(rval/rbin)
        ELSEIF(dens_axis == 2) THEN
           rval = rxyz_lmp(i,2) - box_yl*FLOOR(rxyz_lmp(i,2)/box_yl)
           binval = FLOOR(rval/rbin)
        ELSEIF(dens_axis == 3) THEN
           rval = rxyz_lmp(i,3) - box_zl*FLOOR(rxyz_lmp(i,3)/box_zl)
           binval = FLOOR(rval/rbin)
        END IF

        IF(binval .LE. maxden_bin) THEN

           inst_array(binval,arrindx)=inst_array(binval,arrindx)+1
           ch_instarray(binval,arrindx)=ch_instarray(binval,arrindx)&
                &+INT(charge_lmp(dentyp_arr(j),1))

        END IF

     END IF

  END DO

  DO i = 0, maxden_bin - 1
        
     DO j = 1, ndentypes
        
        densarray(i,j) = densarray(i,j) + inst_array(i,j)
        ch_densarray(i,j) = ch_densarray(i,j) + ch_instarray(i,j) 

     END DO
     
  END DO

  DO i = 1,ntotatoms

     dflag = -1; arrindx = 0
     
     IF(aidvals(i,3) == 3 .OR. aidvals(i,3) == 4) THEN
           
        arrindx = 1
        dflag = 1

     ELSEIF(aidvals(i,3) == 5 .OR. aidvals(i,3) == 6) THEN

        arrindx = 2
        dflag = 1
          
     END IF
     
     IF(dflag == 1) THEN

        IF(dens_axis == 1) THEN
           rval = rxyz_lmp(i,1) - box_xl*FLOOR(rxyz_lmp(i,1)/box_xl)
           binval = FLOOR(rval/rbin)
        ELSEIF(dens_axis == 2) THEN
           rval = rxyz_lmp(i,2) - box_yl*FLOOR(rxyz_lmp(i,2)/box_yl)
           binval = FLOOR(rval/rbin)
        ELSEIF(dens_axis == 3) THEN
           rval = rxyz_lmp(i,3) - box_zl*FLOOR(rxyz_lmp(i,3)/box_zl)
           binval = FLOOR(rval/rbin)
        END IF

        IF(binval .LE. maxden_bin) THEN

           polydensarr(binval,arrindx)=polydensarr(binval,arrindx)+1

        END IF

     END IF

  END DO


  DO i = 0, maxden_bin - 1

     DO j = 1,2

        avgpolyarr(i,j) = avgpolyarr(i,j) + polydensarr(i,j)

     END DO
     
  END DO

IF(ngroups /= 0) THEN

   DO i = 1,ntotatoms
      
        DO j = 1,ngroups
         
          DO k = 1,grp_data(j)

            dflag = -1

            IF(aidvals(i,3) == dengrp_arr(j,k)) THEN
            
               arrindx = j
               dflag = 1
               EXIT
            
            END IF

         END DO
     
         IF(dflag == 1) THEN

            IF(dens_axis == 1) THEN
               rval = rxyz_lmp(i,1) - box_xl*FLOOR(rxyz_lmp(i,1)&
                    &/box_xl)
               binval = FLOOR(rval/rbin)
            ELSEIF(dens_axis == 2) THEN
               rval = rxyz_lmp(i,2) - box_yl*FLOOR(rxyz_lmp(i,2)&
                    &/box_yl)
               binval = FLOOR(rval/rbin)
            ELSEIF(dens_axis == 3) THEN
               rval = rxyz_lmp(i,3) - box_zl*FLOOR(rxyz_lmp(i,3)&
                    &/box_zl)
               binval = FLOOR(rval/rbin)
            END IF
            
            IF(binval .LE. maxden_bin) THEN
               grp_inst_array(binval,arrindx)= grp_inst_array(binval,arrindx) + 1
               chgrpinst(binval,arrindx) = chgrpinst(binval,arrindx)&
                    &+charge_lmp(dengrp_arr(j,k),1)

            END IF
            
         END IF
            
      END DO

   END DO

   DO i = 0, maxden_bin - 1

      DO j = 1, ngroups
         
         grparray(i,j) = grparray(i,j) + grp_inst_array(i,j)
         ch_grparray(i,j) = ch_grparray(i,j) + chgrpinst(i,j)

      END DO
      
   END DO

END IF

END SUBROUTINE COMPUTE_DENS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_ADSORBEDFREE(iframe)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iframe
  INTEGER :: i,j,free_flag,molid,molptr,cntr,adsorbflag
  INTEGER,DIMENSION(1:nfree) :: ptr_adsorbed
  REAL :: fracads,rval

  ptr_adsorbed = -1

  DO i = 1,ntotatoms

     molid = aidvals(i,2); adsorbflag = -1 !initialize flag
     
     DO j = 1,nfree

        IF(molid == free_mols(j)) THEN

           molptr = j
           IF(ptr_adsorbed(molptr) == -1) THEN
              adsorbflag = 0 ! Not adsorbed - proceed
           ELSEIF(ptr_adsorbed(molptr) == 1) THEN
              adsorbflag = 1 ! Adsorbed - no need to check
           ELSE
              PRINT *, "Unknown flag for adsorption"
              PRINT *, i, molid, ptr_adsorbed(molptr)
              STOP
           END IF

           EXIT

        END IF

     END DO

     IF(adsorbflag == 0) THEN

        IF(dens_axis == 1) THEN
           rval = rxyz_lmp(i,1) - box_xl*FLOOR(rxyz_lmp(i,1)&
                &/box_xl)
        ELSEIF(dens_axis == 2) THEN
           rval = rxyz_lmp(i,2) - box_yl*FLOOR(rxyz_lmp(i,2)&
                &/box_yl) 
        ELSEIF(dens_axis == 3) THEN
           rval = rxyz_lmp(i,3) ! No periodicity
        END IF

        IF(rval .LE. brush_ht) ptr_adsorbed(molptr) = 1

     END IF

  END DO
     
  cntr = 0

  DO i = 1,nfree

     IF(ptr_adsorbed(i) == 1) cntr = cntr + 1

  END DO

  fracads = REAL(cntr)/REAL(nfree)

  WRITE(adswrite,'(I0,1X,F14.6)') iframe, fracads

  fracadsavg = fracadsavg + fracads

 
END SUBROUTINE COMPUTE_ADSORBEDFREE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_FREEPENETRATE_MONS(iframe)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: iframe
  INTEGER :: i,j,ibin,dumadscnt,a1id,a2id
  REAL :: rxval, ryval, rzval, rval

  dumadscnt = 0

  DO i = 1,nfreemons

     a1id = free_grp(i)

     DO j = 1,nadsmons

        a2id = ads_grp(j)

        rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
        ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
        rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
        
        rxval = rxval - box_xl*ANINT(rxval/box_xl)
        ryval = ryval - box_yl*ANINT(ryval/box_yl)
        rzval = rzval !No periodicity
        
        rval = sqrt(rxval**2 + ryval**2 + rzval**2)
        
        IF(rval .LE. adscut) THEN
           
           dumadscnt = dumadscnt + 1
           EXIT 

        END IF
        
     END DO

  END DO
  
  WRITE(adsmonwrite,*) iframe, dumadscnt, REAL(dumadscnt)&
       &/REAL(nfreemons)
     
  avgadscnt = avgadscnt + REAL(dumadscnt)/REAL(nfreemons)

END SUBROUTINE COMPUTE_FREEPENETRATE_MONS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_INTERGRAFT_DISTANCE(iframe)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: iframe
  INTEGER :: i,j,a1id,a2id
  REAL :: min_rgraftsq, rgraft_glob_min
  REAL :: rxval, ryval, rzval, rvalsq
  REAL, DIMENSION(1:ngraftmons) :: glob_min_graft_distarr

  rgraft_dist_avg = 0.0; glob_min_graft_distarr = 0.0

  DO i = 1,ngraftmons
     
     a1id = graft_mon_ids(i)
     min_rgraftsq = 100000.0 !Arbitrarily large squared value

     DO j = 1,ngraftmons

        IF (i .NE. j) THEN

           a2id = graft_mon_ids(j)
           rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
           ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
           rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval !No periodicity
           
           rvalsq = rxval**2 + ryval**2 + rzval**2
           
           IF(rvalsq .LE. min_rgraftsq) min_rgraftsq = rvalsq
           
        END IF

     END DO
        
     rgraft_dist_avg = rgraft_dist_avg + sqrt(min_rgraftsq)
     glob_min_graft_distarr(i) =  sqrt(min_rgraftsq)

  END DO

  rgraft_dist_avg = REAL(rgraft_dist_avg)/REAL(ngraftmons)
  rgraft_glob_min = MINVAL(glob_min_graft_distarr)
  WRITE(tethwrite,'(I0,1X,2(F16.8,1X))') iframe, rgraft_dist_avg&
       &,rgraft_glob_min


END SUBROUTINE COMPUTE_INTERGRAFT_DISTANCE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_FREEPENETRATE_CHAINS(iframe)

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: iframe
  INTEGER :: i,j,ibin,dumads_ch_cnt,a1id,a2id,findindex,indexval
  INTEGER :: molid
  REAL :: rxval, ryval, rzval, rval
  INTEGER,DIMENSION(1:nfreechains) :: chainptr_adsorbed

  chainptr_adsorbed = -1
  dumads_ch_cnt = 0

  DO i = 1,nfreemons
     
     a1id   = free_grp(i)
     molid  = aidvals(a1id,2)
     j = 1; findindex = -1; indexval = 0

     DO WHILE (findindex == -1 .AND. j .LE. nfreechains)

        IF(molid == free_chainarr(j)) THEN

           indexval  = j
           findindex = 1

        ELSE
           
           j = j + 1

        END IF

     END DO

     IF(findindex == -1) THEN

        PRINT *, "Could not find matching free chain index"
        PRINT *, a1id, molid
        STOP

     END IF

     IF(chainptr_adsorbed(indexval) == -1) THEN

        DO j = 1,nadsmons

           a2id = ads_grp(j)
           
           rxval = rxyz_lmp(a1id,1) - rxyz_lmp(a2id,1) 
           ryval = rxyz_lmp(a1id,2) - rxyz_lmp(a2id,2) 
           rzval = rxyz_lmp(a1id,3) - rxyz_lmp(a2id,3) 
           
           rxval = rxval - box_xl*ANINT(rxval/box_xl)
           ryval = ryval - box_yl*ANINT(ryval/box_yl)
           rzval = rzval  ! No periodicity
           
           rval = sqrt(rxval**2 + ryval**2 + rzval**2)
           
           IF(rval .LE. adscut .AND. chainptr_adsorbed(indexval) == &
                &-1) THEN
              
              dumads_ch_cnt = dumads_ch_cnt + 1
              chainptr_adsorbed(indexval) = 1

              WRITE(adschwrite2,'(4(I0,1X))') iframe,a1id,molid,a2id

              EXIT

           END IF
        
        END DO
        
     END IF

  END DO
  
  WRITE(adschwrite,*) iframe, dumads_ch_cnt, REAL(dumads_ch_cnt)&
       &/REAL(ngraftchains)
     
  avg_ch_adscnt = avg_ch_adscnt + REAL(dumads_ch_cnt)&
       &/REAL(ngraftchains)

END SUBROUTINE COMPUTE_FREEPENETRATE_CHAINS

!--------------------------------------------------------------------

SUBROUTINE ALLOUTPUTS()

  USE PARAMETERS_PE
  IMPLICIT NONE

  PRINT *, "Number of frames from start to end: ", nframes/(freqfr+1)
  PRINT *, "Frequency of Frames: ", freqfr
  PRINT *, "Total number of Frames analyzed: ", nfrcntr

  WRITE(logout,*) "Number of frames from start to end: ", nframes&
       &/(freqfr+1)
  WRITE(logout,*) "Frequency of Frames: ", freqfr+1
  WRITE(logout,*) "Total number of Frames analyzed: ", nfrcntr

  IF(rdfcalc .OR. acrcalc) THEN
     PRINT *, "Writing RDFs .."
     CALL OUTPUT_ALLRDF()
  END IF

  IF(denscalc) THEN
     PRINT *, "Writing Density Data .."
     CALL OUTPUT_DENS()
  END IF

  IF(adscalc) THEN
     PRINT *, "Average fraction adsorbed by chain count and ht criteri&
          &a: ",fracadsavg/REAL(nfrcntr)
  END IF


  IF(monads) THEN
     PRINT *, "Average adsorbed monomer fraction by monomer count: "&
          &,avgadscnt/REAL(nfrcntr)
  END IF

  IF(chainads) THEN
     PRINT *, "Average adsorbed monomer fraction by monomer count and &
          &distance criteria: ",avg_ch_adscnt/REAL(nfrcntr)
  END IF

END SUBROUTINE ALLOUTPUTS

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_ALLRDF()

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER :: i,j,ierr
  REAL, PARAMETER :: vconst = 4.0*pival/3.0
  REAL :: rlower,rupper,nideal,rdffrnorm,acrnorm
  
  IF(rdfcalc) THEN

     rdffrnorm = INT(nfrcntr/rdffreq)
     rvolavg = rvolavg/REAL(rdffrnorm)
     PRINT *, "Average volume of box", rvolavg
     
     dum_fname = "rdf_"//trim(adjustl(traj_fname))
     OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
          &,status="replace",iostat=ierr)
     
     IF(ierr /= 0) THEN
        PRINT *, "Could not open", trim(dum_fname)
     END IF
     
     WRITE(dumwrite,'(A,2X)',advance="no") "r"
     
     DO j = 1,npairs
        
        WRITE(dumwrite,'(2(I0,1X))',advance="no") pairs_rdf(j,1)&
             &,pairs_rdf(j,2)
        
     END DO
     
     WRITE(dumwrite,*)
     
     DO i = 0,rmaxbin-1
        
        rlower = real(i)*rbinval
        rupper = rlower + rbinval
        nideal = vconst*(rupper**3 - rlower**3)
        
        WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*rbinval*(REAL(2*i&
             &+1))
        
        DO j = 1,npairs
           
           WRITE(dumwrite,'(F16.9,1X)',advance="no")rdfarray(i,j)&
                &/(rdffrnorm*nideal)
           
        END DO
        
        WRITE(dumwrite,*)
        
     END DO
     
     CLOSE(dumwrite)

  END IF

  IF(acrcalc) THEN

     acrnorm = INT(nfrcntr/acrfreq)
     acrvolavg = acrvolavg/REAL(acrnorm)
     
     dum_fname = "PErdf_"//trim(adjustl(traj_fname))
     OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
          &,status="replace",iostat=ierr)
     
     IF(ierr /= 0) THEN
        PRINT *, "Could not open", trim(dum_fname)
     END IF
     
     WRITE(dumwrite,'(A,2X)',advance="no") "r     g(r)"
     WRITE(dumwrite,'(A)',advance="no") "PolycationIDs: "


     DO i = 1,ncatgrp

        WRITE(dumwrite,'(I0,1X)',advance="no") catgrp(i)

     END DO

     WRITE(dumwrite,'(A)',advance="no") "PolyAnionIDs: "


     DO i = 1,nangrp

        WRITE(dumwrite,'(I0,1X)',advance="no") angrp(i)

     END DO

     WRITE(dumwrite,*)

     DO i = 0,acrmaxbin-1
        
        rlower = real(i)*acrbinval
        rupper = rlower + acrbinval
        nideal = vconst*(rupper**3 - rlower**3)
        
        WRITE(dumwrite,'(2(F16.9,2X))') 0.5*acrbinval*(REAL(2*i+1))&
             &,acrdfarray(i)/(acrnorm*nideal)
        
     END DO
     
     CLOSE(dumwrite)

  END IF

END SUBROUTINE OUTPUT_ALLRDF

!--------------------------------------------------------------------

SUBROUTINE OUTPUT_DENS()

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER :: i,j,k,ierr,densfrnorm
  REAL :: rlower,rupper,nideal
  REAL, DIMENSION(1:ndentypes) :: frac
  REAL, DIMENSION(ngroups) :: grp_frac
  REAL, DIMENSION(2) :: grp_poly
  INTEGER :: dumchwrite

  dumchwrite = dumwrite + 1
  densfrnorm = INT(nfrcntr/densfreq)
  normdens = normdens/REAL(densfrnorm)
  denbinavg = denbinavg/REAL(densfrnorm)
  frac = 0
  grp_poly = 0

  DO i = 1,ntotatoms
     
     IF(aidvals(i,3) == 3 .OR. aidvals(i,3) == 4) THEN
        
        grp_poly(1) = grp_poly(1) + 1

     ELSEIF(aidvals(i,3) == 5 .OR. aidvals(i,3) ==6) THEN
        
        grp_poly(2) = grp_poly(2) + 1

     END IF

     DO j = 1,ndentypes

        IF(aidvals(i,3) == dentyp_arr(j)) THEN
           
           frac(j) = frac(j) + 1

        END IF

     END DO

  END DO

  DO j = 1,ndentypes
     PRINT *, j, frac(j)
  END DO

  dum_fname = "dens_"//trim(adjustl(traj_fname))
  OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)

  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF

  dum_fname = "chdens_"//trim(adjustl(traj_fname))
  OPEN(unit = dumchwrite,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)

  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF


  WRITE(dumwrite,'(A,2X)',advance="no") "r"
  WRITE(dumchwrite,'(A,2X)',advance="no") "r"

  DO j = 1,ndentypes

     WRITE(dumwrite,'(I0,1X)',advance="no") dentyp_arr(j)
     WRITE(dumchwrite,'(I0,1X)',advance="no") dentyp_arr(j)

  END DO

  WRITE(dumwrite,*)
  WRITE(dumchwrite,*)

  DO i = 0,maxden_bin-1

     WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*denbinavg*(REAL(2&
          &*i+1))
     WRITE(dumchwrite,'(F16.5,2X)',advance="no") 0.5*denbinavg&
          &*(REAL(2*i+1))

     DO j = 1,ndentypes

        WRITE(dumwrite,'(F16.9,1X)',advance="no") densarray(i,j)&
             &*normdens/(REAL(densfrnorm*frac(j)))
        WRITE(dumchwrite,'(F16.9,1X)',advance="no")ch_densarray(i,j)&
             &*normdens/(REAL(densfrnorm*frac(j)))

     END DO

     WRITE(dumwrite,*)
     WRITE(dumchwrite,*)
     
  END DO

  CLOSE(dumwrite)
  CLOSE(dumchwrite)

! If Group is involved

  IF(ngroups /= 0) THEN

     grp_frac = 0

     DO i = 1,ntotatoms
        
        DO j = 1,ngroups
           
           DO k = 1,grp_data(j)
           
              IF(aidvals(i,3) == dengrp_arr(j,k)) THEN
                 
                 grp_frac(j) = grp_frac(j) + 1
                 
              END IF
              
           END DO
           
        END DO
        
     END DO
 
     dum_fname = "grpdens_"//trim(adjustl(traj_fname))
     OPEN(unit = dumwrite,file =trim(dum_fname),action="write"&
          &,status="replace",iostat=ierr)
     
     IF(ierr /= 0) THEN
        PRINT *, "Could not open", trim(dum_fname)
     END IF

     dum_fname = "chgrpdens_"//trim(adjustl(traj_fname))
     OPEN(unit = dumchwrite,file =trim(dum_fname),action="write"&
          &,status="replace",iostat=ierr)
     
     IF(ierr /= 0) THEN
        PRINT *, "Could not open", trim(dum_fname)
     END IF

     
     WRITE(dumwrite,'(A,2X)',advance="no") "r"
     WRITE(dumchwrite,'(A,2X)',advance="no") "r"

     DO j = 1,ngroups

        WRITE(dumwrite,'(A)',advance="no") "GroupID/Ntypes/TypeIDs: "
        WRITE(dumwrite,'(2(I0,1X,A,1X))',advance="no") j,"/"&
             &,grp_data(j),"/"

        WRITE(dumchwrite,'(A)',advance="no") "GroupID/Ntypes/TypeIDs: &
             &"
        WRITE(dumchwrite,'(2(I0,1X,A,1X))',advance="no") j,"/"&
             &,grp_data(j),"/"

        DO k = 1,grp_data(j)
           
           IF(dengrp_arr(j,k) .NE. -1) THEN
              
              WRITE(dumwrite,'(I0,1X)',advance="no") dengrp_arr(j,k)
              WRITE(dumchwrite,'(I0,1X)',advance="no") dengrp_arr(j,k)

           END IF
           
        END DO
        
     END DO
     
     WRITE(dumwrite,*)
     WRITE(dumchwrite,*)
     
     DO i = 0,maxden_bin-1
        
        WRITE(dumwrite,'(F16.5,2X)',advance="no") 0.5*denbinavg&
             &*(REAL(2*i+1))
        WRITE(dumchwrite,'(F16.5,2X)',advance="no") 0.5*denbinavg&
             &*(REAL(2*i+1))
        
        DO j = 1,ngroups
           
           WRITE(dumwrite,'(F16.9,1X)',advance="no") grparray(i,j)&
                &*normdens/REAL(densfrnorm*grp_frac(j))
           WRITE(dumchwrite,'(F16.9,1X)',advance="no") ch_grparray(i&
                &,j)*normdens/REAL(densfrnorm*grp_frac(j))

           
        END DO
        
        WRITE(dumwrite,*)
        WRITE(dumchwrite,*)

     END DO
     
     CLOSE(dumwrite)
     CLOSE(dumchwrite)

  END IF



! As a check for groups

  dum_fname = "polydens_"//trim(adjustl(traj_fname))
  OPEN(unit = 42,file =trim(dum_fname),action="write"&
       &,status="replace",iostat=ierr)

  IF(ierr /= 0) THEN
     PRINT *, "Could not open", trim(dum_fname)
  END IF

  DO i = 0,maxden_bin-1
        
     WRITE(42,'(3(F16.12,2X))') 0.5*denbinavg*(REAL(2*i+1)),&
          & REAL(avgpolyarr(i,1)*normdens)/REAL(densfrnorm&
          &*grp_poly(1)),REAL(avgpolyarr(i,2)*normdens)&
          &/REAL(densfrnorm*grp_poly(2))
     
  END DO
  
  CLOSE(42)

END SUBROUTINE OUTPUT_DENS

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARRAYS()

  USE PARAMETERS_PE
  IMPLICIT NONE

  INTEGER :: AllocateStatus

! Allocate LAMMPS Structure

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate aidvals"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(charge_lmp(ntotatoms,1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate charge_lmp"
  ALLOCATE(vel_xyz(ntotatoms,4),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate vel_xyz"

  IF(ntotbonds /= 0) THEN
     ALLOCATE(bond_lmp(ntotbonds,4),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate bond_lmp"
  ELSE
     PRINT *, "Warning: No bonds - Not correct for bonded systems"
     ALLOCATE(bond_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(bond_lmp)
  END IF
  
  IF(ntotangls /= 0) THEN
     ALLOCATE(angl_lmp(ntotangls,5),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate angl_lmp"
  ELSE
     ALLOCATE(angl_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(angl_lmp)
  END IF
     
  IF(ntotdihds /= 0) THEN
     ALLOCATE(dihd_lmp(ntotdihds,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate dihd_lmp"
  ELSE
     ALLOCATE(dihd_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(dihd_lmp)
  END IF
  
  IF(ntotimprs /= 0) THEN
     ALLOCATE(impr_lmp(ntotimprs,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate zlmp"
  ELSE
     ALLOCATE(impr_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(impr_lmp)
  END IF


! Allocate Box details

  ALLOCATE(boxx_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxx_arr"
  ALLOCATE(boxy_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxy_arr"
  ALLOCATE(boxz_arr(nframes),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate boxz_arr"

  IF(rdfcalc) THEN
     ALLOCATE(rdfarray(0:rmaxbin-1,npairs),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate rdfarray"
  ELSE
     ALLOCATE(rdfarray(1,1),stat = AllocateStatus)
     DEALLOCATE(rdfarray)
  END IF

  IF(acrcalc) THEN
     ALLOCATE(acrdfarray(0:acrmaxbin-1),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate acrdfarray"
  ELSE
     ALLOCATE(acrdfarray(1),stat = AllocateStatus)
     DEALLOCATE(acrdfarray)
  END IF

  
  IF(denscalc) THEN
     ALLOCATE(densarray(0:maxden_bin-1,ndentypes),stat =&
          & AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate densarray"
     ALLOCATE(ch_densarray(0:maxden_bin-1,ndentypes),stat =&
          & AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate ch_densarray"

     ALLOCATE(avgpolyarr(0:maxden_bin-1,2),stat=AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate avgpolyarr"

     IF(ngroups /= 0) THEN
        ALLOCATE(grparray(0:maxden_bin-1,ngroups),stat =&
             & AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate grparray"
     ELSE
        ALLOCATE(grparray(1,1),stat = AllocateStatus)
        DEALLOCATE(grparray)
     END IF

     IF(ngroups /= 0) THEN
        ALLOCATE(ch_grparray(0:maxden_bin-1,ngroups),stat =&
             & AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate ch_grparray"
     ELSE
        ALLOCATE(ch_grparray(1,1),stat = AllocateStatus)
        DEALLOCATE(ch_grparray)
     END IF

  ELSE
     ALLOCATE(densarray(1,1),stat = AllocateStatus)
     DEALLOCATE(densarray)
     ALLOCATE(ch_densarray(1,1),stat = AllocateStatus)
     DEALLOCATE(ch_densarray)
  END IF

  IF(adscalc == 0) THEN
     ALLOCATE(free_mols(1),stat = AllocateStatus)
     DEALLOCATE(free_mols)
  END IF

  IF(monads == 0) THEN
     ALLOCATE(free_ptr(1),stat = AllocateStatus)
     DEALLOCATE(free_ptr)
     ALLOCATE(ads_ptr(1),stat = AllocateStatus)
     DEALLOCATE(ads_ptr)
     ALLOCATE(free_grp(1),stat = AllocateStatus)
     DEALLOCATE(free_grp)
     ALLOCATE(ads_grp(1),stat = AllocateStatus)
     DEALLOCATE(ads_grp)
  END IF
  PRINT *, "Successfully allocated memory"

END SUBROUTINE ALLOCATE_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE PARAMETERS_PE

  IMPLICIT NONE

  DEALLOCATE(aidvals)
  DEALLOCATE(rxyz_lmp)
  DEALLOCATE(vel_xyz)

  IF(ntotbonds /= 0) DEALLOCATE(bond_lmp)
  IF(ntotangls /= 0) DEALLOCATE(angl_lmp)
  IF(ntotdihds /= 0) DEALLOCATE(dihd_lmp)
  IF(ntotimprs /= 0) DEALLOCATE(impr_lmp)

END SUBROUTINE DEALLOCATE_ARRAYS

!--------------------------------------------------------------------


