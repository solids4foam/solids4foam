      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1    RPL,DDSDDT,DRPLDE,DRPLDT,
     2     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      IMPLICIT NONE

      INCLUDE 'ABA_PARAM.INC'

!     Declare input/output variables
!     Take care: unused variables here are set to doubles
      real*8, INTENT(OUT) :: STRESS(NTENS)
      real*8, INTENT(OUT) :: STATEV(NSTATV)
      real*8, INTENT(OUT) :: DDSDDE(NTENS,NTENS)
      real*8, INTENT(IN) :: SSE,SPD,SCD, RPL,DDSDDT,DRPLDE,DRPLDT
      real*8, INTENT(IN)  ::  STRAN(NTENS), DSTRAN(NTENS)
      real*8, INTENT(IN) :: TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     1    CMNAME
      integer, INTENT(IN) :: NDI,NSHR,NTENS,NSTATV
      real*8, INTENT(IN) :: PROPS(NPROPS)
      integer, INTENT(IN) :: NPROPS
      real*8, INTENT(IN) :: COORDS,DROT,PNEWDT, CELENT,DFGRD0,DFGRD1,
     1    NOEL,NPT,LAYER,KSPT,JSTEP,KINC

!     Declare local variables
      real*8 E, NU, ONE, TWO, ALAMBDA, BLAMBDA, CLAMBDA
      integer I, J
      ONE=1.0D0
      TWO=2.0D0

!     Stiffness matrix
      E=PROPS(1)
      NU=PROPS(2)
      ALAMBDA=E/(ONE+nu)/(ONE-TWO*nu)
      BLAMBDA=(ONE-nu)
      CLAMBDA=(ONE-TWO*nu)
      DO I=1,NTENS
         DO J=1,NTENS
            DDSDDE(I,J)=0.0D0
         ENDDO
      ENDDO
      DDSDDE(1,1)=(ALAMBDA*BLAMBDA)
      DDSDDE(2,2)=(ALAMBDA*BLAMBDA)
      DDSDDE(3,3)=(ALAMBDA*BLAMBDA)
      DDSDDE(4,4)=(ALAMBDA*CLAMBDA)
      DDSDDE(5,5)=(ALAMBDA*CLAMBDA)
      DDSDDE(6,6)=(ALAMBDA*CLAMBDA)
      DDSDDE(1,2)=(ALAMBDA*nu)
      DDSDDE(1,3)=(ALAMBDA*nu)
      DDSDDE(2,3)=(ALAMBDA*nu)
      DDSDDE(2,1)=(ALAMBDA*nu)
      DDSDDE(3,1)=(ALAMBDA*nu)
      DDSDDE(3,2)=(ALAMBDA*nu)

!     Calculate stress
      DO I=1,NTENS
         DO J=1,NTENS
            STRESS(I)=STRESS(I)+DDSDDE(I,J)*STRAN(J)
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE UMAT
