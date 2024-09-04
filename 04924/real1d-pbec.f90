
MODULE COMM_DATA
  INTEGER, PARAMETER :: N = 4000, N2 = N/2, NX = N-1
  INTEGER, PARAMETER :: NSTP = 1000000, NPAS = 1000, NRUN = 2 * 80000
  INTEGER, PARAMETER :: NSTORE = 2 * 160
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
END MODULE COMM_DATA
!
 MODULE GPE_DATA
  USE COMM_DATA, ONLY : N, PI
  INTEGER, PARAMETER :: INPUT = 0  !
  COMPLEX (8), PARAMETER :: CI = (0.0D0,1.0D0) ! Complex i
  REAL (8), PARAMETER :: DX = 0.025D0, DT = DX * DX
  !REAL (8), PARAMETER :: DX = 0.0025D0, DT = 0.00002D0
  REAL (8), PARAMETER :: GAMA = 0.01D0
!
  REAL (8), PARAMETER :: G_1D = -10.0D0 ! Nonlinearity (two-body)
  REAL (8), PARAMETER :: GAEFF = 0.0D0, NLGAM = 0.0D0 ! File 1
  COMPLEX (8), PARAMETER :: CI_GAEFF = CI * GAEFF
  COMPLEX (8), PARAMETER :: CI_NLGAM = CI * NLGAM
  REAL (8), DIMENSION(0:N) :: X, X2, V
  COMPLEX (8), DIMENSION(0:N) :: CP
  REAL (8) :: G, GSTP
  COMPLEX (8) :: G_IGAMA
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : N
  COMPLEX (8), DIMENSION(0:N) :: CAL, CGA
  COMPLEX (8) :: CAIPM
END MODULE CN_DATA 

PROGRAM GROSS_PITAEVSKII_SSCN_1D
   
   USE COMM_DATA
   USE GPE_DATA
   
   IMPLICIT NONE
   
!------------------------ INTERFACE BLOCKS -----------------------
    INTERFACE 
     SUBROUTINE INITIALIZE(CP)
       IMPLICIT NONE
       COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
     END SUBROUTINE INITIALIZE
   END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE COEF()
      IMPLICIT NONE
    END SUBROUTINE COEF
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE CALCNU(CP, DT)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE LU(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
    END SUBROUTINE LU
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE RAD(CP2, RMS)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: CP2
      REAL (8), INTENT(OUT) :: RMS
    END  SUBROUTINE RAD
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END  SUBROUTINE CHEM
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE WRITE_DENSITY(FUNIT, U2)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:), INTENT(IN) :: U2
    END SUBROUTINE WRITE_DENSITY
  END INTERFACE
!------------------------ INTERFACE BLOCKS -----------------------
  INTEGER :: K, I
  REAL (8) :: ZNORM, MU, EN, RMS, T, T1, T2, ITOUT,TAU,F
  REAL (8), DIMENSION(0:N) :: P, P2
  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE
!
  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME(T1)
  CALL INITIALIZE(CP)
  OPEN(7, FILE = 'real-out.txt')
  WRITE(7,900)
  WRITE(7,*)
  WRITE(7,903) G_1D, GAEFF, NLGAM
  WRITE(7,904) GAMA
  WRITE(7,*)
  WRITE(7,905) N
  WRITE(7,906) DX
  WRITE(7,907) NSTP, NPAS, NRUN
  WRITE(7,908) DT
  WRITE(7,*)
  900 FORMAT('  New Real time propagation 1d ')
  903 FORMAT('  Nonlinearity G_1D =',F12.4,', Linear L/G = ',F12.4,', NL L/G = ',F12.4)
  904 FORMAT('  Parameter of trap: GAMMA =',F7.2)
  909 FORMAT('  Transverse trap parameter =',F7.3)
  905 FORMAT('# Space Stp:  N = ', I8)
  906 FORMAT('  Space Step: DX = ', F10.6)
  907 FORMAT('# Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  908 FORMAT('  Time Step:   DT = ', F10.6)
!
  CALL CALCULATE_TRAP()
  CALL COEF()
  CALL NORM(CP, ZNORM)
  CALL CHEM(CP, MU, EN)
  P = ABS(CP); P2 = P * P
  CALL RAD(P2, RMS)
!
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)
  WRITE (7, 1003) ZNORM, MU, EN, RMS, P2(N2)
  1001 FORMAT (19X,'-----------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener/N', 6X, '<x>', 6X, '|Psi(0)|^2')
  1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.6, 2F11.5)
!----------------------------------------------------------------------------------------------------------
  OPEN(1, FILE = 'real-den-ini.txt')
  CALL WRITE_DENSITY(1, P2)
  CLOSE(1)
!-------------------------------------------------------------------------------------------------------------
  TAU=NSTP*DT
  open(101,file="check.txt")
  open(102,file="DYNAMICS1.txt")
  IF (NSTP /= 0) THEN
    !GSTP =  G_1D / DFLOAT(NSTP)
    !G = 0.0D0
    DO K = 0, NSTP
	  	
		T = DFLOAT(K) * DT
		IF (T <= TAU) THEN
			F = T / TAU  ! Linearly scale F from 0 to 1 over time TAU
		ELSE
			F = 1.0D0   ! After TAU, F remains 1
		END IF
		G = F * G_1D
		
		CALL CALCULATE_TRAP()
		DO I = 0, N
			V(I) = (1.0D0 - F) * V(I)  ! V decreases from its initial value to 0
		END DO
      write(101,*)T,"	",G,"	",V(20)
	  G_IGAMA = G - CI_NLGAM
      CALL CALCNU(CP, DT)
      CALL LU(CP)
	  
	  IF (MOD(K,2*NSTORE) == 0) THEN
      P = ABS(CP); P2 = P * P
      DO I = 0, N, 4
        WRITE(102,*) X(I),"								",DFLOAT(K) * DT,"								", P2(I)
      END DO
    END IF
   
	END DO
    
	CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    P = ABS(CP); P2 = P * P
    CALL RAD(P2, RMS)
    WRITE (7, 1005) ZNORM, MU, EN, RMS, P2(N2)
    1005 FORMAT('After NSTP iter.:',F8.4, 2F12.6, 2F11.5)
  ELSE
    G = G_1D
  END IF
close(101)
close(102)
!
   G_IGAMA = G - CI_NLGAM
!
  DO K = 1, NPAS
    CALL CALCNU(CP, DT)
    CALL LU(CP)
  END DO
!
  IF(NPAS /= 0)THEN
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    P = ABS(CP); P2 = P * P
    CALL RAD(P2, RMS)
    WRITE (7, 1006) ZNORM, MU, EN, RMS, P2(N2)
  END IF
  1006 FORMAT('After NPAS iter.:',F8.4, 2F12.6, 2F11.5)
!
  OPEN(4, FILE = 'real-den-int.txt')
  CALL WRITE_DENSITY(4, P2)
  CLOSE(4)
  DO K = 0, N
    WRITE(12, '(F12.6, 2G17.8E3)') X(K), CP(K)
  END DO
!
!	TAU=10
!	T=0.0D0
!
	
OPEN(10, FILE = 'dynamics-output.txt')
	
	DO K = 1, NRUN
    CALL CALCNU(CP, DT)
    CALL LU(CP)
    IF (MOD(K,NSTORE) == 0) THEN
      P = ABS(CP); P2 = P * P
      DO I = 0, N, 4
        !WRITE(11, '(1X, F12.6, $)') P2(I)
		WRITE(10,*) X(I),"								",DFLOAT(K) * DT,"								", P2(I)
      END DO
    END IF
  END DO
!
  IF(NRUN /= 0)THEN
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    P = ABS(CP); P2 = P * P
    CALL RAD(P2, RMS)
    WRITE (7, 1007) ZNORM, MU, EN, RMS, P2(N2)
  END IF
  1007 FORMAT('After NRUN iter.:',F8.4, 2F12.6, 2F11.5)
 
 CLOSE(10)
 !
 
  OPEN(4, FILE = 'real-den-fin.txt')
  CALL WRITE_DENSITY(4, P2)
  CLOSE(4)
  DO K = 0, N
    WRITE(13, '(F12.6, 2G17.8E3)') X(K), CP(K)
  END DO
 !
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7, 1001)
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time:', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time:', INT(T2-T1), ' seconds'
  CLOSE (7)
  
  PRINT*,"COMPLETED"
!

END PROGRAM GROSS_PITAEVSKII_SSCN_1D

SUBROUTINE INITIALIZE(CP)
  USE COMM_DATA, ONLY : N, N2, PI, NSTP
  USE GPE_DATA, ONLY : INPUT, DX, X, X2, GAMA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8) ::  TX, TY, TMP, PI4
  INTEGER :: I
  IF (INPUT /= 0) THEN
    DO I = 0, N
      READ (*,*) X(I), TX, TY
      CP(I) = CMPLX(TX, TY)
    END DO
  ELSE
    PI4 = SQRT(SQRT(PI/GAMA))
    DO I = 0, N
      X(I) = (I-N2) * DX
      X2(I) = X(I) * X(I)
      CP(I) = CMPLX(EXP(-GAMA * X2(I) / 2.0D0) / PI4, 0.0D0)
    END DO
  END IF
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
  USE GPE_DATA, ONLY : X2, V, GAMA
  IMPLICIT NONE
  REAL (8) :: GAM2
  GAM2 = GAMA * GAMA
  V = X2 * GAM2 / 2.0D0
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF()
  USE COMM_DATA, ONLY :  NX
  USE GPE_DATA, ONLY : DT, DX, CI
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: J
  REAL (8) :: DX2
  COMPLEX (8) :: CAI0, CDT
  CDT = CI * DT
  DX2 = DX * DX
  CAIPM = -CDT/(4.0D0 * DX2)
  CAI0 = 1.0D0 + CDT/(2.0D0 * DX2)
  CAL(NX) = 0.0D0
  CGA(NX) = -1.0D0/CAI0
  DO J = NX, 1, -1
     CAL(J-1) = CGA(J)*CAIPM
     CGA(J-1) = -1.0D0/(CAI0 + CAIPM * CAL(J-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU(CP, DT)
  USE COMM_DATA, ONLY : N, NX
  USE GPE_DATA, ONLY : X2, V, G,  CI, CI_GAEFF, G_IGAMA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(0:N) :: P, P2
  COMPLEX (8) :: CTMP
  INTEGER :: I
  DO I = 0, N
    P(I) = ABS(CP(I))
    P2(I) = P(I) * P(I)
    CTMP = DT * (V(I) + G_IGAMA * P2(I) + CI_GAEFF)
    CP(I) = CP(I) * EXP(-CI * CTMP)
  END DO
END SUBROUTINE CALCNU

SUBROUTINE LU(CP)
  USE COMM_DATA, ONLY : N, NX
  USE CN_DATA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  COMPLEX (8), DIMENSION(0:NX) :: CBE
  COMPLEX (8) :: CXX
  INTEGER :: I
!
  CBE(NX) = CP(N)
   DO I = NX, 1, -1
      CXX = CP(I) + (CP(I+1) - 2.0D0 * CP(I) + CP(I-1)) * CAIPM ! correction
      CBE(I-1) = CGA(I) * (CAIPM * CBE(I) - CXX)
   END DO
!-----------------------------
! Boundary condition periodic:
!  DO I = 0, NX
!     CP(I+1) = CAL(I)*CP(I) + CBE(I)
!  END DO
!  CP(0) = CP(N)
!  CP(1) = CP(NX)
!-----------------------------
! Boundary condition reflecting:
   CP(0) = (0.D0, 0.D0)
   DO I = 0, NX-1
      CP(I+1) = CAL(I)*CP(I) + CBE(I)
   END DO
   CP(N) = (0.D0, 0.D0)
!-----------------------------
END SUBROUTINE LU

SUBROUTINE NORM(CP, ZNORM)
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
  !----------------------------------------------------------------------
  INTERFACE
      FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN)  :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
  !----------------------------------------------------------------------
  REAL (8), DIMENSION(0:N) :: P, P2
!
  P = ABS(CP)
  P2 = P * P
  ZNORM = SQRT(SIMP(P2, DX))
END SUBROUTINE NORM
 
SUBROUTINE RAD(CP2, RMS)
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX, X2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: CP2
  REAL (8), INTENT(OUT) :: RMS
!----------------------------------------------------------------------
   INTERFACE
    FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!----------------------------------------------------------------------
 REAL (8), DIMENSION(0:SIZE(CP2)-1) :: P2
!
  P2 = X2 * CP2
  RMS = SQRT(SIMP(P2, DX))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
  USE COMM_DATA, ONLY : N,NX  
  USE GPE_DATA, ONLY : DX, V, G, GAEFF, NLGAM, CI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN
  !----------------------------------------------------------------------
  INTERFACE
     FUNCTION DIFF(P, DX) RESULT (DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF
  END INTERFACE
  !----------------------------------------------------------------------
   INTERFACE
    FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!----------------------------------------------------------------------
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: CPR, CPI, CPR2, CPI2, P2, GP2
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: DCP, DP2, TMP1D
!
  REAL (8) :: ZNORM
!
  CPR = REAL(CP)
  CPI = IMAG(CP)
  CPR2 = CPR * CPR
  CPI2 = CPI * CPI
  P2 = CPR2 + CPI2
  GP2 = G * P2
!  
  DCP = DIFF(CPR, DX)
  DP2 = DCP * DCP
  TMP1D = DP2 + ((V + GP2) * CPR + (NLGAM * P2 - GAEFF) * CPI) * CPR
  ZNORM = SIMP(CPR2, DX)
  MU = SIMP(TMP1D, DX) / ZNORM
  DCP = DIFF(CPI, DX)
  DP2 = DCP * DCP
  TMP1D = DP2 + ((V + GP2) * CPI - (NLGAM * P2 - GAEFF) * CPR) * CPI
  ZNORM = SIMP(CPI2, DX)
  IF (ZNORM /= 0D0) EN = SIMP(TMP1D, DX) / ZNORM
END SUBROUTINE CHEM

 FUNCTION SIMP(F, DX) RESULT (VALUE)
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  REAL (8) :: VALUE
  REAL (8) :: F1, F2
  INTEGER :: I, N
!
  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN 
  F2 = F(2) 
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  VALUE = DX*(F(0) + 4.D0*F1 + 2.D0*F2 + F(N))/3.D0
END FUNCTION SIMP

FUNCTION DIFF(P,DX)  RESULT (DP)
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: P
  REAL (8), INTENT(IN) :: DX
  REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
  INTEGER :: I, N
!
  N = SIZE(P) - 1
  DP(0) = 0.D0
  DP(1) = (P(2) - P(0))/(2.D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.D0*P(I-1)+8.D0*P(I+1)-P(I+2))/(12.D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.D0*DX)
  DP(N) = 0.D0
END FUNCTION DIFF
 
SUBROUTINE WRITE_DENSITY(FUNIT, U2)
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : X
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:), INTENT(IN) :: U2
  INTEGER :: I
!
  DO I = 0, N
    WRITE(FUNIT, 999) X(I), U2(I)
  END DO
  999 FORMAT(F12.4, E17.5E3)
END SUBROUTINE WRITE_DENSITY
!# End File : real1d-pbec.f90
