!---------------To analyze Polyelectrolyte Static Properties---------
!---------------Version 2: May-09-2018-------------------------------
!---------------Main File: pe_analyze.f90----------------------------
!********************************************************************

MODULE PARAMETERS_PE

  USE OMP_LIB
  IMPLICIT NONE

  ! Required Input Variables

  INTEGER :: initdist
  INTEGER :: nframes, skipfr, freqfr, nfrcntr
  INTEGER :: nwater, nchains
  INTEGER :: atperchain
  INTEGER :: nproc

  !Structural analysis input details

  INTEGER :: rdffreq,rmaxbin,npairs,ngroups,mstart,nfree
  REAL    :: rvolavg,rdomcut,rbinval
  INTEGER :: rgfreq, densfreq, dens_axis, ndentypes, maxden_bin
  REAL    :: normdens,denbinavg,brush_ht,fracadsavg
  INTEGER :: nfreegrp,nadsgrp, nfreemons, nadsmons
  REAL    :: adscut,avgadscnt,avg_ch_adscnt
  INTEGER :: ncatgrp, nangrp, nanions, ncations
  REAL    :: acrbinval, acrvolavg, acrdomcut,chadscut
  INTEGER :: acrmaxbin,acrfreq
  INTEGER :: nfreechains, ngraftchains, ngraftmons,graft_type
  REAL    :: rgraft_dist_avg

  ! All flags
  
  INTEGER :: rdfcalc, rgcalc, rgall, rgavg, denscalc
  INTEGER :: adscalc,monads,chainads
  INTEGER :: acrcalc, avg_rgraft_calc

  ! File names and unit Numbers
  
  CHARACTER(LEN = 256) :: ana_fname,data_fname,traj_fname,log_fname
  CHARACTER(LEN = 256) :: rdf_fname, dum_fname
  INTEGER, PARAMETER :: anaread = 2,   logout = 3
  INTEGER, PARAMETER :: inpread = 100, rgwrite = 400,rgavgwrite = 300
  INTEGER, PARAMETER :: dumwrite = 200, rgswrite = 250,adswrite=360
  INTEGER, PARAMETER :: adsmonwrite=370,adschwrite = 380
  INTEGER, PARAMETER :: adschwrite2=390, tethwrite = 410

  !Math Constants

  REAL*8, PARAMETER :: pival  = 3.14159265359
  REAL*8, PARAMETER :: pi2val = 2.0*pival

  !Global analysis variables and arrays

  INTEGER :: atomflag, velflag, bondflag, anglflag, dihdflag,imprflag
  INTEGER :: ntotatoms, ntotbonds, ntotangls,ntotdihds,ntotimprs
  INTEGER :: ntotatomtypes,ntotbondtypes,ntotangltypes,ntotdihdtypes&
       &,ntotimprtypes

  !Lammps trajectory file read details

  REAL :: box_xl,box_yl,box_zl, boxval
  INTEGER*8 :: timestep

  !Structural variables

  INTEGER :: rdfpaircnt
  REAL    :: rvolval
  
  !Structural Average Variables

  REAL :: re2ave, re4ave, rg2ave, rg4ave, b2ave
  
  !Required Arrays - LAMMPS

  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rxyz_lmp, vel_xyz, charge_lmp&
       &,masses
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bond_lmp, angl_lmp,&
       & dihd_lmp, impr_lmp,aidvals
  CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: keywords
  REAL,ALLOCATABLE,DIMENSION(:):: boxx_arr, boxy_arr,boxz_arr

  !Required Arrays - Structural Quantities

  REAL,ALLOCATABLE,DIMENSION(:,:):: rdfarray, densarray, grparray
  REAL,ALLOCATABLE,DIMENSION(:,:):: ch_densarray, ch_grparray
  REAL,ALLOCATABLE,DIMENSION(:,:):: avgpolyarr
  REAL,ALLOCATABLE,DIMENSION(:):: acrdfarray
  INTEGER, ALLOCATABLE, DIMENSION(:) :: graft_mon_ids
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: pairs_rdf
  INTEGER, ALLOCATABLE, DIMENSION(:) :: dentyp_arr,grp_data,free_mols
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: dengrp_arr
  INTEGER, ALLOCATABLE, DIMENSION(:) :: catgrp,angrp
  INTEGER, ALLOCATABLE, DIMENSION(:) :: free_ptr, ads_ptr, free_grp&
       &,ads_grp,free_chainarr
  
END MODULE PARAMETERS_PE

!--------------------------------------------------------------------
