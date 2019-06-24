!---------------To analyze Polyelectrolyte Static Properties---------
!---------------Version 2: May-09-2018-------------------------------
!---------------Main File: pe_analyze.f90----------------------------
!********************************************************************

MODULE SWITCH_PARAMS

  USE OMP_LIB
  IMPLICIT NONE

  ! Required Input Variables

  INTEGER :: nframes, skipfr, freqfr, nfrcntr
  INTEGER :: nwater, nchains
  INTEGER :: atperchain
  INTEGER :: nproc

  ! File names and unit Numbers
  
  CHARACTER(LEN = 512) :: ana_fname,data_fname,out_fname,log_fname
  CHARACTER(LEN = 512) :: dum_fname
  INTEGER, PARAMETER :: anaread = 2,   logout = 3
  INTEGER, PARAMETER :: inpread = 100, outwrite = 400

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

  !Required Arrays - LAMMPS

  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rxyz_lmp, vel_xyz, charge_lmp&
       &,masses
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bond_lmp, angl_lmp,&
       & dihd_lmp, impr_lmp,aidvals
  CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: keywords
  REAL,ALLOCATABLE,DIMENSION(:):: boxx_arr, boxy_arr,boxz_arr

  
END MODULE SWITCH_PARAMS

!--------------------------------------------------------------------
