! Input file to generate LAMMPS file for polyelectrolyte simulations
! Use in conjunction with lammps_inp.f90 and ran_numbers.f90
! Change the required parameters
! arch details: "arch=number:polycation-polyanion 
! arch = 1:block-block;arch=2:block-alt;arch=3:alt-block;arch=4:alt
! -alt.
MODULE PARAMS

  USE RAN_NUMBERS

  IMPLICIT NONE

! Parameter data for creating the data file

  INTEGER, PARAMETER :: N = py_free_chains
  INTEGER, PARAMETER :: M = py_free_mons
  INTEGER, PARAMETER :: N_brush = py_graft_chains
  INTEGER, PARAMETER :: M_brush = py_graft_mons
  INTEGER, PARAMETER :: tail_brush = py_tail_mons
  INTEGER, PARAMETER :: N_salt  = py_nsalt
  INTEGER, PARAMETER :: default_dim = 1 !If 1=>53*53*120.Or else use box.dat
  REAL, PARAMETER    :: charge_frac = py_charge_frac
  INTEGER, PARAMETER :: ncntr_brush  = N_brush*INT((M_brush&
       &-tail_brush)*charge_frac)
  INTEGER, PARAMETER :: ncntr_free   = N*INT(M*charge_frac)
  REAL    :: brush_dist
  INTEGER :: arch

! Box/Particle details

  REAL :: boxl_x, boxl_y, boxl_z
  REAL :: volbox, density
  INTEGER :: totpart
  INTEGER :: n_cntr_ions, nchains, npolyatoms

! Flags for creating the data file

  INTEGER, PARAMETER :: atomic = 0
  INTEGER, PARAMETER :: triblock = 0
  INTEGER, PARAMETER :: stretched = 0
  INTEGER, PARAMETER :: numatomtypes = 8
  INTEGER, PARAMETER :: numbondtypes = 1
  INTEGER, PARAMETER :: numangltypes = 1
  INTEGER, PARAMETER :: numdihdtypes = 0
  INTEGER, PARAMETER :: bondtype = 1
  INTEGER, PARAMETER :: angltype = 1
  INTEGER, PARAMETER :: dihdtype = 0
  INTEGER, PARAMETER :: outfile  = 17

! Global Arrays involved in creating data file
  
  REAL,ALLOCATABLE,DIMENSION(:,:) :: rxyz, uxyz
  REAL,ALLOCATABLE,DIMENSION(:) :: charge
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: aidvals,ixyz

! Character Arrays for creating the data file name

  CHARACTER (LEN = 12) :: f_char
  CHARACTER (LEN = 60 ):: datafile

  TYPE (RAN_SAVE) :: X
  INTEGER(4) :: S
  
END MODULE PARAMS
