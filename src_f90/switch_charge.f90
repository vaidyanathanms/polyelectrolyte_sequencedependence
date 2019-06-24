!---------------To switch charge in Polyelectrolytes-----------------
!---------------Version 1: Dec-19-2018-------------------------------
!---------------Parameter File: switch_charge_params.f90-------------
!********************************************************************

PROGRAM SWITCH_CHARGE

  USE SWITCH_PARAMS
  IMPLICIT NONE

  PRINT *, "Rewriting Charges..."
  PRINT *, "Starting OMP Threads .."

!$OMP PARALLEL
  nproc = OMP_GET_NUM_THREADS()
  PRINT *, "Number of threads are: ", nproc
!$OMP END PARALLEL

  CALL READ_ANA_IP_FILE()
  CALL READ_WRITE_DATAFILE()

  PRINT *, "All Calculations Completed Succesfully :)"

END PROGRAM SWITCH_CHARGE

!--------------------------------------------------------------------


SUBROUTINE READ_ANA_IP_FILE()

  USE SWITCH_PARAMS

  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,AllocateStatus,i,j
  CHARACTER(256) :: dumchar

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"
  
  logflag = 0
  nchains = 0
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

     ELSEIF(dumchar == 'out_file') THEN

        READ(anaread,*,iostat=ierr) out_fname

     ELSEIF(dumchar == 'nchains') THEN

        READ(anaread,*,iostat=ierr) nchains

     ELSEIF(dumchar == 'natoms') THEN
        
        READ(anaread,*,iostat=ierr) ntotatoms

     ELSEIF(dumchar == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag  = 1

     ELSE
        
        PRINT *, "unknown keyword: ", trim(dumchar)
        STOP

     END IF

  END DO

  IF(nchains == 0 .OR. ntotatoms == 0) THEN
     PRINT *, "ERROR: Cannot have zero chains"
     STOP
  END IF
  IF(logflag == 0) log_fname = "log.switchcharge"
  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  PRINT *, "Analysis input file read finished .."
  
END SUBROUTINE READ_ANA_IP_FILE

!--------------------------------------------------------------------

SUBROUTINE READ_WRITE_DATAFILE()

  USE SWITCH_PARAMS

  IMPLICIT NONE

  INTEGER :: i,j,k,ierr,u,AllocateStatus,imax
  INTEGER :: flag, cntr, nwords
  INTEGER :: aid,molid,atype,ix,iy,iz
  REAL    :: charge,rx,ry,rz,massval
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(256) :: rline,dumchar

  CALL COMPUTE_INIT_NLINES(imax)

  OPEN(unit=inpread,file = trim(data_fname),action =&
       & 'read', status='old',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "Data file not found"

  OPEN(unit=outwrite,file = trim(out_fname),action =&
       & 'write', status='replace',iostat=ierr) 
  
  IF(ierr .NE. 0) STOP "OutData file not found"


  WRITE(logout,*) "Datafile used is :", trim(adjustl(data_fname))

  WRITE(logout,*) "Outdatafile is ", trim(adjustl(out_fname))

  
  READ(inpread,'(A256)') dumchar
  WRITE(outwrite,*) trim(adjustl(dumchar))
  READ(inpread,'(A256)') dumchar
  WRITE(outwrite,*) trim(adjustl(dumchar))

  DO i = 1,imax-2 !Change here according to convenience
       
     READ(inpread,'(A256)') dumchar
     
     WRITE(outwrite,*) adjustl(trim(adjustl(dumchar)))

  END DO

  READ(inpread,'(A256)') dumchar
  WRITE(outwrite,*) trim(adjustl(dumchar))

  READ(inpread,'(A256)') dumchar
  WRITE(outwrite,*) trim(adjustl(dumchar))
  READ(inpread,'(A256)') dumchar
  WRITE(outwrite,*) trim(adjustl(dumchar))
  READ(inpread,'(A256)') dumchar
  WRITE(outwrite,*) adjustl(trim(adjustl(dumchar)))

  flag = 0; cntr = 0

  DO 

     READ(inpread,'(A)',iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     !READ DATA HERE FOR CHARGES AND MOLID
     !READ EVERYTHING AND OVERWRITE LATER

     IF(trim(dumchar(1:5)) == "Atoms") THEN
        
        WRITE(outwrite,'(A)') trim(dumchar)

        WRITE(outwrite,*)

        atomflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,ntotatoms

           READ(inpread,*) aid,molid,atype,charge,rx,ry,rz,ix,iy,iz

           IF(molid .GT. nchains) THEN
              IF(charge == 1) THEN
                 atype = 7
              ELSEIF(charge == -1) THEN
                 atype = 8
              ELSE
                 PRINT *, "Unknown charge", charge
              END IF
           END IF

           WRITE(outwrite,'(3(I0,1X),4(ES23.16E2,1X),3(I0,1X))') aid&
                &,molid,atype,charge,rx,ry,rz,ix,iy,iz

        END DO

     ELSE

        WRITE(outwrite,*) adjustl(trim(dumchar))

     END IF

  END DO

!!$     IF(trim(dumchar) == "Masses") THEN
!!$
!!$        ALLOCATE(masses(ntotatomtypes,1),stat = AllocateStatus)
!!$        IF(AllocateStatus/=0) STOP "did not allocate masses"
!!$         
!!$        DO j = 1,ntotatomtypes
!!$
!!$           READ(inpread,*) u, masses(u,1)
!!$
!!$        END DO
!!$
!!$     END IF
!!$
!!$     IF(trim(dumchar) == "Velocities") THEN
!!$             
!!$        velflag = 1
!!$        print *, "Reading ", trim(dumchar), " info"
!!$
!!$        DO j = 1,ntotatoms
!!$
!!$           READ(inpread,*) vel_xyz(j,1),vel_xyz(j,2),vel_xyz(j,3)&
!!$                &,vel_xyz(j,4)
!!$
!!$        END DO
!!$
!!$
!!$     END IF
!!$
!!$     IF(trim(dumchar) == "Bonds") THEN
!!$             
!!$        bondflag = 1
!!$        print *, "Reading ", trim(dumchar), " info"
!!$
!!$        DO j = 1,ntotbonds
!!$
!!$           READ(inpread,*) bond_lmp(j,1),bond_lmp(j,2),bond_lmp(j,3)&
!!$                &,bond_lmp(j,4)
!!$
!!$        END DO
!!$
!!$     END IF
!!$
!!$     IF(trim(dumchar) == "Angles") THEN
!!$             
!!$        anglflag = 1
!!$        print *, "Reading ", trim(dumchar), " info"
!!$
!!$        DO j = 1,ntotangls
!!$
!!$           READ(inpread,*) angl_lmp(j,1),angl_lmp(j,2),angl_lmp(j,3)&
!!$                &,angl_lmp(j,4),angl_lmp(j,5)
!!$
!!$        END DO
!!$
!!$     END IF
!!$
!!$     IF(trim(dumchar) == "Dihedrals") THEN
!!$             
!!$        dihdflag = 1
!!$        print *, "Reading", trim(dumchar), "info"
!!$
!!$        DO j = 1,ntotdihds
!!$
!!$           READ(inpread,*) dihd_lmp(j,1),dihd_lmp(j,2),dihd_lmp(j,3)&
!!$                &,dihd_lmp(j,4),dihd_lmp(j,5),dihd_lmp(j,6)
!!$
!!$        END DO
!!$
!!$     END IF
!!$  
!!$     IF(trim(dumchar) == "Impropers") THEN
!!$             
!!$        imprflag = 1
!!$        print *, "Reading", trim(dumchar), "info"
!!$
!!$        DO j = 1,ntotimprs
!!$
!!$           READ(inpread,*) impr_lmp(j,1),impr_lmp(j,2),impr_lmp(j,3)&
!!$                &,impr_lmp(j,4),impr_lmp(j,5),impr_lmp(j,6)
!!$
!!$        END DO
!!$
!!$     END IF
!!$
!!$  END DO
!!$  
!!$  PRINT *, "Fileread finish .."
!!$

END SUBROUTINE READ_WRITE_DATAFILE

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_INIT_NLINES(imax)

  USE SWITCH_PARAMS

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



