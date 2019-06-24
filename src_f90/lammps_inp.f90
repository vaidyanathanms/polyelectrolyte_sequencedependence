! To generate LAMMPS input file for polyelectrolyte simulations
! Use in conjunction with lmp_params.f90 and ran_numbers.f90
PROGRAM LAMMPSINP

  USE PARAMS

  IMPLICIT NONE
  
  LOGICAL :: input_coor = .false.
  REAL :: bondl, bondlsq
  INTEGER :: ierror,narg

  bondlsq = 0

  CALL SYSTEM_CLOCK(S)

  narg = IARGC()
  
  IF(narg == 1) THEN

     OPEN (unit = outfile, file = "lmp_input.dat", status ="replace",&
          & action="write",iostat=ierror)
     
     IF(ierror /= 0) STOP "Cannot open lmp_input.txt"
        
     CALL CREATEFILE(narg)
     CALL COMPUTE_BOX()
     CALL ALLOCATE_ARRAYS()
     CALL INPCOR()
     CALL CREATE_ATTYPE()
     CALL LMP_COORD()
     CALL DEALLOCATE_ARRAYS()
     WRITE(outfile,*) "Successfully Written Datafile .."
     CLOSE (unit = outfile)

  ELSE

     WRITE(outfile,*) "Please specify the architecture"
     WRITE(outfile,*) "1: Block - Block"
     WRITE(outfile,*) "1: Block - Alter"
     WRITE(outfile,*) "1: Alter - Block"
     WRITE(outfile,*) "1: Alter - Alter"

  END IF

END PROGRAM LAMMPSINP

!--------------------------------------------------------------------

SUBROUTINE LMP_COORD()

  USE PARAMS
  
  IMPLICIT NONE
  
  INTEGER :: i,j, k, ierror
  INTEGER ::  bondid, anglid, dihdid
  REAL ::  massval, rx, ry, rz
  
  i = 1
  
  PRINT *, "Writing LAMMPS Datafile .. "

20 FORMAT(5X,I0,2X,A)
22 FORMAT(5X,I0,2X,A)
24 FORMAT(5X,I0,2X,F14.6,2X,A)
  
  OPEN (unit=10, file = datafile, status="replace",action=&
       &"write",iostat = ierror)
  
  IF(ierror /= 0) STOP "Failed to open datafile"
     
  WRITE (10,*) "Data for CG-PE simulations for ",f_char,"architecture"
  WRITE (10,*) 
  WRITE (10,20) totpart, "atoms"

  IF(numbondtypes /= 0) THEN
     WRITE (10,20) N*(M-1)+N_brush*(M_brush-1), "bonds"
  ELSE
     WRITE (10,20) 0, "bonds"
  END IF

  IF(numangltypes /= 0) THEN
     WRITE (10,20) N*(M-2)+N_brush*(M_brush-2), "angles"
  ELSE
     WRITE (10,20) 0, "angles"
  END IF

  IF(numdihdtypes /= 0) THEN
     WRITE (10,20) N*(M-3)+N_brush*(M_brush-3), "dihedrals"
  ELSE
     WRITE (10,20) 0, "dihedrals"
  END IF

  WRITE (10,20) 0, "impropers"
  WRITE (10,20) numatomtypes, "atom types"
  WRITE (10,20) numbondtypes, "bond types"
  WRITE (10,22) numangltypes, "angle types"
  WRITE (10,22) numdihdtypes, "dihedral types"
  WRITE (10,22) 0, "improper types"

  WRITE (10,*)
  WRITE (10,24) 0, boxl_x, "xlo xhi"
  WRITE (10,24) 0, boxl_y, "ylo yhi"
  WRITE (10,24) 0, boxl_z, "zlo zhi"
  WRITE (10,*)
  WRITE (10,*) "Masses"
  WRITE (10,*)
  
!!$Fetching image information - Brush

  ixyz = 0
  
  DO i = 1,N_brush
     
     k = (i-1)*M_brush + 1

     ixyz(k,1) = 0
     ixyz(k,2) = 0
     ixyz(k,3) = 0
     
     DO j = 1,M_brush-1
        
        k = (i-1)*M_brush + j

        rx = rxyz(k,1) - rxyz(k+1,1)
        ry = rxyz(k,2) - rxyz(k+1,2)
        rz = rxyz(k,3) - rxyz(k+1,3)
        
        CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
        CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
        CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))
        
     END DO
     
  END DO

!!$Fetching image information - Free

  DO i = 1,N
     
     k = (i-1)*M + 1 + N_brush*M_brush

     ixyz(k,1) = 0
     ixyz(k,2) = 0
     ixyz(k,3) = 0
     
     DO j = 1,M-1
        
        k = (i-1)*M + j + N_brush*M_brush

        rx = rxyz(k,1) - rxyz(k+1,1)
        ry = rxyz(k,2) - rxyz(k+1,2)
        rz = rxyz(k,3) - rxyz(k+1,3)
        
        CALL IMGFLAGS(rx,ixyz(k,1),boxl_x,ixyz(k+1,1))
        CALL IMGFLAGS(ry,ixyz(k,2),boxl_y,ixyz(k+1,2))
        CALL IMGFLAGS(rz,ixyz(k,3),boxl_z,ixyz(k+1,3))

     END DO
     
  END DO

  ! Writing Masses

  DO i = 1,numatomtypes

     massval = 1.000
     WRITE(10,'(I0,1X,F14.8)') i, massval

  END DO
  
  ! Writing atomic corrdinates
  
  WRITE (10,*) 
  WRITE (10,*) "Atoms"
  WRITE (10,*)

  DO i = 1,totpart
     
     WRITE(10,'(3(I0,1X),4(F14.6,1X),3(I0,1X))') aidvals(i,1),&
          & aidvals(i,2), aidvals(i,3), charge(i), rxyz(i,1), rxyz(i&
          &,2),rxyz(i,3), ixyz(i,1), ixyz(i,2), ixyz(i,3)
        
  END DO

  IF(numbondtypes /= 0) THEN

     ! Writing Bond Details  
     
     bondid = 0
     WRITE (10,*)
     WRITE (10,*) "Bonds"
     WRITE (10,*)
     
     DO i = 1,N_brush
        
        DO j = 1,M_brush-1
           
           k = (i-1)*M_brush + j
           bondid = bondid + 1
           
           WRITE(10,'(4(I0,2X))') bondid, bondtype, aidvals(k,1)&
                &,aidvals(k+1,1)
           
        END DO
        

        
     END DO

     DO i = 1,N
        
        DO j = 1,M-1
           
           bondid = bondid + 1
           k = (i-1)*M + j + N_brush*M_brush
           
           WRITE(10,'(4(I0,2X))') bondid, bondtype, aidvals(k,1)&
                &,aidvals(k+1,1)
           
        END DO

     END DO

  END IF

  IF(numangltypes /= 0) THEN

     ! Writing Angle Details

     anglid = 0
     WRITE (10,*)
     WRITE (10,*) "Angles"
     WRITE (10,*)
     
     DO i = 1,N_brush
        
        DO j = 1,M_brush-2
           
           anglid = anglid + 1
           k = (i-1)*M_brush + j           
           WRITE(10,'(5(I0,2X))') anglid, angltype, aidvals(k,1)&
                &,aidvals(k+1,1),aidvals(k+2,1)
           
        END DO
        
     END DO

     DO i = 1,N
        
        DO j = 1,M-2
           
           anglid = anglid + 1
           k = (i-1)*M + j + N_brush*M_brush
           WRITE(10,'(5(I0,2X))') anglid, angltype, aidvals(k,1)&
                &,aidvals(k+1,1),aidvals(k+2,1)
           
        END DO
        
     END DO

  END IF

  IF(numdihdtypes /= 0) THEN

     ! Writing Dihedral Details
     
     dihdid = 0
     WRITE (10,*)
     WRITE (10,*) "Dihedrals"
     WRITE (10,*)
     
     DO i = 1,N_brush
        
        DO j = 1,M_brush-3
           
           dihdid = dihdid + 1
           k = (i-1)*M_brush + j 

           WRITE(10,'(6(I0,2X))') dihdid, dihdtype, aidvals(k,1)&
                &,aidvals(k+1,1), aidvals(k+2,1), aidvals(k+3,1)
           
        END DO
        
     END DO

     DO i = 1,N
        
        DO j = 1,M-3
           
           dihdid = dihdid + 1
           k = (i-1)*M + j + N_brush*M_brush

           WRITE(10,'(6(I0,2X))') dihdid, dihdtype, aidvals(k,1)&
                &,aidvals(k+1,1), aidvals(k+2,1), aidvals(k+3,1)
           
        END DO
        
     END DO

  END IF

  CLOSE(unit = 10)

END SUBROUTINE LMP_COORD

!--------------------------------------------------------------------  

SUBROUTINE IMGFLAGS(dist, img, boxl, imgout)
    
  USE PARAMS
  IMPLICIT NONE
  
  REAL, INTENT(IN) :: dist,boxl
  INTEGER, INTENT(IN) :: img
  INTEGER, INTENT(OUT) :: imgout
  INTEGER :: nx
  
  IF(dist > boxl/2) THEN
     
     nx = img + 1
     
  ELSEIF(dist < -boxl/2) THEN
     
     nx = img - 1
     
  ELSE
     
     nx = img
     
  END IF
  
  imgout = nx
  
END SUBROUTINE IMGFLAGS

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_BOX()

  USE PARAMS

  IMPLICIT NONE
  INTEGER :: ierr

  IF(default_dim == 1) THEN
     boxl_x = 53.0; boxl_y = 53.0; boxl_z = 120.0
  ELSE
     OPEN(unit = 24,file="box.dat",action="read",status="old"&
          &,iostat=ierr)
     IF(ierr /= 0) STOP "box.dat cannot be found"
     READ(24,*) boxl_x, boxl_y, boxl_z
  END IF

  WRITE(outfile,*) "Simulation Details ...."
  PRINT *, "Simulation Details ...."

  totpart = N*M+N_brush*M_brush+2.0*N_salt+ncntr_brush+ncntr_free
  volbox  = boxl_x*boxl_y*boxl_z
  density = REAL(totpart)/REAL(volbox)
  nchains = N+N_brush
  npolyatoms = N*M + N_brush*M_brush
  brush_dist = SQRT(REAL(boxl_x*boxl_y))/SQRT(REAL(2.0*N_brush))
  ! Extra 2.0 factor above so that diagonal is approximately equal to
  !  sqrt(area/num_chains)

  WRITE(outfile,*) "Total particles: ", totpart
  WRITE(outfile,*) "Number of atomtypes: ", numatomtypes
  WRITE(outfile,*) "# of Polyanions: ", N
  WRITE(outfile,*) "Monomers per polyanion: ", M
  WRITE(outfile,*) "# of Polycations: ", N_brush
  WRITE(outfile,*) "Monomers per polycation: ", M_brush
  WRITE(outfile,*) "# of Salt: ", N_salt
  WRITE(outfile,*) "# of Brush counterions: ", ncntr_brush
  WRITE(outfile,*) "# of Free counterions: ", ncntr_free
  WRITE(outfile,*) "# of Polyelectrolytes: ", npolyatoms
  WRITE(outfile,*) "Total Number of Particles: ", totpart
  WRITE(outfile,*) "LX/LY/LZ: ", boxl_x, boxl_y, boxl_z
  WRITE(outfile,*) "Density: ", density
  WRITE(outfile,*) "Minimum distance between brush base: ", brush_dist

END SUBROUTINE COMPUTE_BOX

!--------------------------------------------------------------------

SUBROUTINE INPCOR()
  
  USE PARAMS

  IMPLICIT NONE
  
  INTEGER :: i,j,k,u,v,ierror,brushcnt, kold
  REAL, PARAMETER :: r0init  = 0.97
  REAL, PARAMETER :: r0sq3   = r0init/sqrt(3.0)
  REAL, PARAMETER :: rmaxsq  = r0init*r0init
  REAL, PARAMETER :: math_pi = 3.14159265359
  REAL :: theta, phi
  REAL :: rx, ry, rz,rval
  LOGICAL :: in_box
  CALL RAN_INIT(S,X)
  
  WRITE(outfile,*) "Random Initial Configuration : NRRW"

  i = 1

! Create Brush

  WRITE(outfile,*) "Generating Brush Configurations .. "

  DO WHILE(i .LE. N_brush)
     
     k = (i-1)*M_brush + 1
        
     rxyz(k,1) = RAN1(X)*boxl_x
     rxyz(k,2) = RAN1(X)*boxl_y
     rxyz(k,3) = 0.50000

     !Non overlapping base

     brushcnt = 1

     DO WHILE(brushcnt .LT. i)
        
        kold = M_brush*(brushcnt-1) + 1

        rx = rxyz(k,1) - rxyz(kold,1)
        ry = rxyz(k,2) - rxyz(kold,2)
        rz = rxyz(k,3) - rxyz(kold,3)

        rx = rx - boxl_x*ANINT(rx/boxl_x)
        ry = ry - boxl_y*ANINT(ry/boxl_y)
        rz = rz - boxl_z*ANINT(rz/boxl_z)

        rval = sqrt(rx**2 + ry**2 + rz**2)

        IF(rval .GE. brush_dist) THEN
           brushcnt = brushcnt + 1
        ELSE
           rxyz(k,1) = RAN1(X)*boxl_x
           rxyz(k,2) = RAN1(X)*boxl_y
           brushcnt = 1
        END IF

     END DO

     j = 2

     DO WHILE (j .LE. M_brush)

        k = (i-1)*M_brush + j
        in_box = .false.

        DO WHILE(in_box == .false.)

           theta     = math_pi*RAN1(X)
           phi       = 2*math_pi*RAN1(X)
           rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
           rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
           rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
           IF(rxyz(k,3) .GE. 0.1 .AND. rxyz(k,3) .LE. boxl_z-0.1)&
                & in_box=.true.

        END DO

        j = j + 1
                
     END DO
     
     i = i + 1
        
  END DO

! Create Free Chains

  WRITE(outfile,*) "Generating Free Chain Configurations .. "

  i = 1

  DO WHILE(i .LE. N)
     
     k = (i-1)*M + 1 + N_brush*M_brush
     
     in_box = .false.
     DO WHILE(in_box == .false.)

        rxyz(k,1) = RAN1(X)*boxl_x
        rxyz(k,2) = RAN1(X)*boxl_y
        rxyz(k,3) = RAN1(X)*boxl_z
        IF(rxyz(k,3) .GT. 0.01 .AND. rxyz(k,3) .LT. boxl_z-0.1)&
             & in_box=.true.

     END DO

     j = 2

     DO WHILE (j .LE. M)
        
        k = (i-1)*M + j + N_brush*M_brush
        in_box = .false.

        DO WHILE(in_box == .false.)

           theta     = math_pi*RAN1(X)
           phi       = 2*math_pi*RAN1(X)
           rxyz(k,1) = rxyz(k-1,1) + r0init*sin(theta)*cos(phi)
           rxyz(k,2) = rxyz(k-1,2) + r0init*sin(theta)*sin(phi)
           rxyz(k,3) = rxyz(k-1,3) + r0init*cos(theta)
           IF(rxyz(k,3) .GE. 0.1 .AND. rxyz(k,3) .LE. boxl_z-0.1)&
                & in_box=.true.

        END DO

        j = j + 1
                
     END DO
     
     i = i + 1
        
  END DO

! Create Salt 

  WRITE(outfile,*) "Generating Salt .. "

  k = 1 + N_brush*M_brush + N*M; i = 1
  DO WHILE(i .LE. 2.0*N_salt)
     
     in_box = .false.
     
     DO WHILE(in_box == .false.)

        rxyz(k,1) = RAN1(X)*boxl_x
        rxyz(k,2) = RAN1(X)*boxl_y
        rxyz(k,3) = RAN1(X)*boxl_z
        IF(rxyz(k,3) .GE. 0.1 .AND. rxyz(k,3) .LE. boxl_z-0.1)&
             & in_box=.true.
        
     END DO
     
     i = i + 1; k = k + 1
        
  END DO


! Create brush counterions

  k = 1 + N_brush*M_brush + N*M +2.0*N_salt; i = 1
  DO WHILE(i .LE. ncntr_brush)
     
     in_box = .false.
     
     DO WHILE(in_box == .false.)

        rxyz(k,1) = RAN1(X)*boxl_x
        rxyz(k,2) = RAN1(X)*boxl_y
        rxyz(k,3) = RAN1(X)*boxl_z
        IF(rxyz(k,3) .GE. 0.1 .AND. rxyz(k,3) .LE. boxl_z-0.1)&
             & in_box=.true.
        
     END DO
     
     i = i + 1; k = k + 1
        
  END DO

! Create free counterions

  k = 1 + N_brush*M_brush + N*M +2.0*N_salt + ncntr_brush
  i = 1
  DO WHILE(i .LE. ncntr_free)
     
     in_box = .false.
     
     DO WHILE(in_box == .false.)

        rxyz(k,1) = RAN1(X)*boxl_x
        rxyz(k,2) = RAN1(X)*boxl_y
        rxyz(k,3) = RAN1(X)*boxl_z
        IF(rxyz(k,3) .GE. 0.1 .AND. rxyz(k,3) .LE. boxl_z-0.1)&
             & in_box=.true.
        
     END DO
     
     i = i + 1; k = k + 1
        
  END DO


  ! PBC
  
  DO i = 1,totpart
     
     rxyz(i,1) = rxyz(i,1) - boxl_x*floor(rxyz(i,1)/boxl_x)
     rxyz(i,2) = rxyz(i,2) - boxl_y*floor(rxyz(i,2)/boxl_y)
     rxyz(i,3) = rxyz(i,3) - boxl_z*floor(rxyz(i,3)/boxl_z)
        
  END DO
  
END SUBROUTINE INPCOR

!--------------------------------------------------------------------

SUBROUTINE CREATE_ATTYPE()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,k
  INTEGER :: fin_neut_brush, fin_neut_free
  REAL :: csum

  csum = 0
  IF(arch .LE. 2) fin_neut_brush = tail_brush + INT(0.5*(M_brush&
       &-tail_brush))

  WRITE(outfile,*) "Generating Charges .. "
! Brush

  DO i = 1, N_brush

     DO j = 1,M_brush

        k = (i-1)*M_brush + j

        aidvals(k,1) = k
        aidvals(k,2) = i

        IF(j == 1) THEN
           aidvals(k,3) = 1
           charge(k)  = 0.0 
        ELSEIF(j .LE. tail_brush) THEN
           aidvals(k,3) = 2
           charge(k) = 0.0
        ELSEIF(arch .LE. 2) THEN
           IF(j .LE. fin_neut_brush) THEN
              aidvals(k,3) = 3
              charge(k) = 0.0
           ELSE
              aidvals(k,3) = 4
              charge(k) = 1.0
           END IF
        ELSE
           IF(mod(j-tail_brush,2) /= 0) THEN
              aidvals(k,3) = 3
              charge(k) = 0.0
           ELSE
              aidvals(k,3) = 4
              charge(k) = 1.0
           END IF
        END IF

     END DO

  END DO

! Free

  IF(arch == 1 .OR. arch == 3) fin_neut_free = INT(0.5*M)

  DO i = 1,N

     DO j = 1,M

        k = (i-1)*M + j + N_brush*M_brush

        aidvals(k,1) = k
        aidvals(k,2) = i + N_brush

        IF(arch == 1 .OR. arch == 3) THEN
           IF(j .LE. fin_neut_free) THEN
              aidvals(k,3) = 5
              charge(k) = 0.0
           ELSE
              aidvals(k,3) = 6
              charge(k) = -1.0
           END IF
        ELSE
           IF(mod(j,2) /= 0) THEN
              aidvals(k,3) = 5
              charge(k) = 0.0
           ELSE
              aidvals(k,3) = 6
              charge(k) = -1.0
           END IF
        END IF

     END DO

  END DO

! Salt

  DO i = 1,2*N_salt

     k = i + N_brush*M_brush + N*M     

     aidvals(k,1) = k
     aidvals(k,2) = 1 + N_brush + N

     IF(mod(i,2) /=  0) THEN
        aidvals(k,3) = 7
        charge(k) = 1.0
     ELSE
        aidvals(k,3) = 8
        charge(k) = -1.0
     END IF
     
  END DO

! Counterions in polymers - will have same type as the salt

  DO i = 1,ncntr_brush

     k = i + N_brush*M_brush + N*M + 2*N_salt
     aidvals(k,1) = k
     aidvals(k,2) = 1 + N_brush + N
     aidvals(k,3) = 7
     charge(k) = -1.0

  END DO

  DO i = 1,ncntr_free

     k = i + N_brush*M_brush + N*M + 2*N_salt + ncntr_brush
     aidvals(k,1) = k
     aidvals(k,2) = 1 + N_brush + N
     aidvals(k,3) = 8
     charge(k) = 1.0

  END DO

! Charge Neutrality

  DO i = 1, totpart

     csum = csum + charge(i)
     
  END DO

  IF(abs(csum) .GT. 0.00000001) THEN

     OPEN(unit = 77,file="charge.txt",action="write",status="replace")

     DO i = 1,totpart
        WRITE(77,*), i,aidvals(i,3), charge(i)
     END DO

     WRITE(outfile,*) "System not charge neutral", csum
     PRINT *, "ERROR: System not charge neutral", csum
     STOP
     
  ELSE

     WRITE(outfile,*) "Good Charge Neutrality ", csum
     PRINT *, "Good Charge Neutrality ", csum

  END IF
     
END SUBROUTINE CREATE_ATTYPE

!--------------------------------------------------------------------

SUBROUTINE CREATEFILE(narg)

  USE PARAMS
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: narg
  CHARACTER (LEN = 6) :: nmon_char
  CHARACTER (LEN = 7) :: prefix = "PEdata_"
  CHARACTER (LEN = 1) :: ch_arch

  CALL GETARG(narg,ch_arch)
  
  WRITE(nmon_char,"(I0)") N
  READ(ch_arch,*) arch

  IF(arch == 1) THEN

     WRITE(outfile,*) "Preparing data file for simulations with block-&
          &block"
     f_char = 'block_block'

  ELSEIF(arch == 2) THEN

     WRITE(outfile,*) "Preparing data file for simulations with block-&
          &alter"
     f_char = 'block_alter'

  ELSEIF(arch == 3) THEN

     WRITE(outfile,*) "Preparing data file for simulations with alter-&
          &block"
     f_char = 'alter_block'

  ELSEIF(arch == 4) THEN

     WRITE(outfile,*) "Preparing data file for simulations with alter-&
          &alter"
     f_char = 'alter_alter'

  ELSE

     WRITE(outfile,*) "Unknown architecture: ", arch
     PRINT *, "Unknown architecture: ", arch
     STOP
     
  END IF
  datafile = prefix//trim(adjustl(nmon_char))//"_"&
       &//trim(adjustl(f_char))//'.dat'
     
  WRITE(outfile,*) "Data file generated for simulation is ",&
       & trim(datafile)

  PRINT *, "Data file generated for simulation is ", trim(datafile)
  
END SUBROUTINE CREATEFILE

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARRAYS()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ALLOCATE(rxyz(1:totpart,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate rxyz"
  ALLOCATE(uxyz(1:totpart,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate uxyz"
  ALLOCATE(charge(1:totpart), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate charge"
  ALLOCATE(aidvals(1:totpart,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate aidvals"
  ALLOCATE(ixyz(1:totpart,3), Stat=AllocateStatus)
  IF(AllocateStatus /= 0) STOP "Could not allocate ixyz"


END SUBROUTINE ALLOCATE_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATE_ARRAYS()

  USE PARAMS
  IMPLICIT NONE

  DEALLOCATE(rxyz)
  DEALLOCATE(uxyz)
  DEALLOCATE(ixyz)
  DEALLOCATE(charge)
  DEALLOCATE(aidvals)

END SUBROUTINE DEALLOCATE_ARRAYS

!--------------------------------------------------------------------
