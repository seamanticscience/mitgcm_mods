#include "DIC_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"

CBOP
C !ROUTINE: DIC_SURFFORCING_INIT

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_SURFFORCING_INIT(
     I          myThid)

C !DESCRIPTION:
C  Calculate first guess of pH

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "DIC_VARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "DIC_LOAD.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER  myThid

#ifdef ALLOW_DIC

C !LOCAL VARIABLES: ====================================================
       INTEGER i,j,k,kLev,pco2_init_iter
       INTEGER intimeP, intime0, intime1
       _RL aWght, bWght
C Number of iterations for pCO2 solvers...
C Solubility relation coefficients
C local variables for carbon chem
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      _RL surfalk(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfphos(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfsi(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surftemp(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfsalt(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL surfdic(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C       INTEGER iprt,jprt
      LOGICAL pH_isLoaded
CEOP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      kLev=1

      CALL DIC_INI_ATMOS( startTime, nIter0, myThid )

ccccccccccccccccccccccccccccccccccccccccc
      IF ( periodicExternalForcing ) THEN
c read in silica field
         CALL LEF_ZERO( silica0,myThid )
         CALL LEF_ZERO( silica1,myThid )

C--   Now calculate whether it is time to update the forcing arrays
       CALL GET_PERIODIC_INTERVAL(
     O                   intimeP, intime0, intime1, bWght, aWght,
     I                   externForcingCycle, externForcingPeriod,
     I                   deltaTClock, startTime, myThid )

       _BARRIER
       _BEGIN_MASTER(myThid)
        WRITE(standardMessageUnit,'(A,I10,A,2(2I5,A))')
     &   ' DIC_SURFFORCING_INIT, it=', nIter0,
     &   ' : Reading new data, i0,i1=', intime0, intime1
       _END_MASTER(myThid)

#ifdef DIC_CARBON_COMPONENTS
       IF ( DIC_silicaFile .NE. ' '  ) THEN
         CALL READ_REC_XYZ_RS( DIC_silicaFile,silica0,intime0,
     &        nIter0,myThid )
         CALL READ_REC_XYZ_RS( DIC_silicaFile,silica1,intime1,
     &        nIter0,myThid )
       ENDIF

       _EXCH_XYZ_RS(silica0, myThid )
       _EXCH_XYZ_RS(silica1, myThid )

       IF ( DIC_silicaFile .NE. ' '  ) THEN
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             SILICA(i,j,k,bi,bj)= bWght*silica0(i,j,k,bi,bj)
     &                          + aWght*silica1(i,j,k,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
#else
       IF ( DIC_silicaFile .NE. ' '  ) THEN
         CALL READ_REC_XY_RS( DIC_silicaFile,silica0,intime0,
     &        nIter0,myThid )
         CALL READ_REC_XY_RS( DIC_silicaFile,silica1,intime1,
     &        nIter0,myThid )
       ENDIF

       _EXCH_XY_RS(silica0, myThid )
       _EXCH_XY_RS(silica1, myThid )

       IF ( DIC_silicaFile .NE. ' '  ) THEN
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
             SILICA(i,j,bi,bj)= bWght*silica0(i,j,bi,bj)
     &                        + aWght*silica1(i,j,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
#endif /* DIC_CARBON_COMPONENTS */
c end periodicExternalForcing
      ENDIF

C =================================================================

      jMin=1
      jMax=sNy
      iMin=1
      iMax=sNx

#ifdef DIC_CARBON_COMPONENTS  
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           DO k=1,Nr
            pH(i,j,k,bi,bj) = 8. _d 0
           ENDDO
          ENDDO
         ENDDO
       ENDDO
      ENDDO
#else
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
            pH(i,j,bi,bj) = 8. _d 0
          ENDDO
         ENDDO
       ENDDO
      ENDDO
#endif

#ifdef DIC_CARBON_COMPONENTS           
C Initiate all coefficients here, as we are called before dic_components_init
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           DO k=1,Nr
            fugf (i,j,k,bi,bj) = 0. _d 0
            ff   (i,j,k,bi,bj) = 0. _d 0
            ak0  (i,j,k,bi,bj) = 0. _d 0
            ak1  (i,j,k,bi,bj) = 0. _d 0
            ak2  (i,j,k,bi,bj) = 0. _d 0
            akb  (i,j,k,bi,bj) = 0. _d 0
            ak1p (i,j,k,bi,bj) = 0. _d 0
            ak2p (i,j,k,bi,bj) = 0. _d 0
            ak3p (i,j,k,bi,bj) = 0. _d 0
            aksi (i,j,k,bi,bj) = 0. _d 0
            akw  (i,j,k,bi,bj) = 0. _d 0
            aks  (i,j,k,bi,bj) = 0. _d 0
            akf  (i,j,k,bi,bj) = 0. _d 0
            bt   (i,j,k,bi,bj) = 0. _d 0
            st   (i,j,k,bi,bj) = 0. _d 0
            ft   (i,j,k,bi,bj) = 0. _d 0
            cat  (i,j,k,bi,bj) = 0. _d 0
            co3  (i,j,k,bi,bj) = 0. _d 0
            hco3 (i,j,k,bi,bj) = 0. _d 0
            co2aq(i,j,k,bi,bj)= 0. _d 0
            Ksp_TP_Calc(i,j,k,bi,bj) = 0. _d 0
            Ksp_TP_Arag(i,j,k,bi,bj) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
#else
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
            ak0  (i,j,k,bi,bj) = 0. _d 0
            ak1  (i,j,k,bi,bj) = 0. _d 0
            ak2  (i,j,k,bi,bj) = 0. _d 0
            akb  (i,j,k,bi,bj) = 0. _d 0
            ak1p (i,j,k,bi,bj) = 0. _d 0
            ak2p (i,j,k,bi,bj) = 0. _d 0
            ak3p (i,j,k,bi,bj) = 0. _d 0
            aksi (i,j,k,bi,bj) = 0. _d 0
            akw  (i,j,k,bi,bj) = 0. _d 0
            aks  (i,j,k,bi,bj) = 0. _d 0
            akf  (i,j,k,bi,bj) = 0. _d 0
            fugf (i,j,k,bi,bj) = 0. _d 0
            ff   (i,j,k,bi,bj) = 0. _d 0
            bt   (i,j,k,bi,bj) = 0. _d 0
            st   (i,j,k,bi,bj) = 0. _d 0
            ft   (i,j,k,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif

      pH_isLoaded = .FALSE.
      IF ( nIter0.GT.PTRACERS_Iter0 .OR.
     &    (nIter0.EQ.PTRACERS_Iter0 .AND. pickupSuff.NE.' ')
     &   ) THEN
C       Read pH from a pickup file if needed
        CALL DIC_READ_PICKUP(
     O                        pH_isLoaded,
     I                        nIter0, myThid )
      ENDIF

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

C determine inorganic carbon chem coefficients
        DO j=jMin,jMax
         DO i=iMin,iMax

#ifdef DIC_BIOTIC
#ifdef DIC_BOUNDS
             surfalk(i,j) = max(0.4 _d 0,
     &                          min(10. _d 0,PTRACER(i,j,kLev,bi,bj,2)))
     &                          * maskC(i,j,kLev,bi,bj)
             surfphos(i,j)  = max(1.0 _d -11,
     &                          min(1. _d -1,PTRACER(i,j,kLev,bi,bj,3)))
     &                          * maskC(i,j,kLev,bi,bj)
#else
             surfalk(i,j) = PTRACER(i,j,kLev,bi,bj,2)
     &                          * maskC(i,j,kLev,bi,bj)
             surfphos(i,j)= PTRACER(i,j,kLev,bi,bj,3)
     &                          * maskC(i,j,kLev,bi,bj)
#endif
#else
             surfalk(i,j) = 2.366595 _d 0 *salt(i,j,kLev,bi,bj)/35. _d 0
     &                          * maskC(i,j,kLev,bi,bj)
             surfphos(i,j)= 5.1225 _d -4 * maskC(i,j,kLev,bi,bj)
#endif
#ifdef DIC_CARBON_COMPONENTS            
C FOR 3d NON-INTERACTIVE Si
             surfsi(i,j)   = SILICA(i,j,kLev,bi,bj)
     &                          * maskC(i,j,kLev,bi,bj)
#else
C FOR 2d NON-INTERACTIVE Si
             surfsi(i,j)   = SILICA(i,j,bi,bj)
     &                          * maskC(i,j,kLev,bi,bj)
#endif

#ifdef DIC_BOUNDS
            surftemp(i,j) = max(-4. _d 0,
     &                          min(50. _d 0, theta(i,j,kLev,bi,bj)))
            surfsalt(i,j) = max(4. _d 0,
     &                          min(50. _d 0, salt(i,j,kLev,bi,bj)))
            surfdic(i,j)  = max(0.4 _d 0,
     &                         min(10. _d 0, PTRACER(i,j,kLev,bi,bj,1)))
#else
            surftemp(i,j) = theta(i,j,kLev,bi,bj)
            surfsalt(i,j) = salt(i,j,kLev,bi,bj)
            surfdic(i,j)  = PTRACER(i,j,kLev,bi,bj,1)
     &                    * maskC(i,j,kLev,bi,bj)
#endif
         ENDDO
        ENDDO

#ifdef DIC_CARBON_COMPONENTS
C Fill 3d values using coeffs from dic_carbon_components_chem
       DO k=1,Nr
        CALL CARBON_COEFFS_SURF(
     I              bi,bj,iMin,iMax,jMin,jMax,
     I              k,myThid)
     
         IF ( k .GT. 1 ) THEN
          CALL CARBON_COEFFS_DEPTH(
     I              bi,bj,iMin,iMax,jMin,jMax,
     I              k,myThid)
         ENDIF
       ENDDO
#else
        CALL CARBON_COEFFS(
     I              surftemp,surfsalt,
     I              bi,bj,iMin,iMax,jMin,jMax,myThid)
#endif /* DIC_CARBON_COMPONENTS */

C====================================================================

        IF ( .NOT.pH_isLoaded ) THEN
C set guess of pH for first step here
          WRITE(standardMessageUnit,*) 'QQ: pCO2 approximation method'
c first approximation
C$TAF LOOP = parallel
          DO j=jMin,jMax
C$TAF LOOP = parallel
           DO i=iMin,iMax
            DO pco2_init_iter=1,10
C$TAF STORE pH(i,j,bi,bj)                              = dic_surf
C$TAF STORE surfalk(i,j), surfphos(i,j), surfsi(i,j)   = dic_surf
#ifdef DIC_CARBON_COMPONENTS
C Use 3d values using coeffs from dic_carbon_components_chem
             DO k=1,Nr
              IF ( maskC(i,j,k,bi,bj) .NE. 0. _d 0) THEN
               CALL CALC_PCO2_APPROX(
     I        theta(i,j,k,bi,bj),salt(i,j,k,bi,bj),
     I        PTRACER(i,j,k,bi,bj,iDIC), PTRACER(i,j,k,bi,bj,iPO4),
     I        SILICA(i,j,k,bi,bj),PTRACER(i,j,k,bi,bj,iALK),
     I        ak1(i,j,k,bi,bj),ak2(i,j,k,bi,bj),
     I        ak1p(i,j,k,bi,bj),ak2p(i,j,k,bi,bj),
     I        ak3p(i,j,k,bi,bj),
     I        aks(i,j,k,bi,bj),akb(i,j,k,bi,bj),
     I        akw(i,j,k,bi,bj),
     I        aksi(i,j,k,bi,bj),akf(i,j,k,bi,bj),
     I        ak0(i,j,k,bi,bj), fugf(i,j,k,bi,bj),
     I        ff(i,j,k,bi,bj),
     I        bt(i,j,k,bi,bj),st(i,j,k,bi,bj),
     I        ft(i,j,k,bi,bj),
     U        pH(i,j,k,bi,bj),pCO2(i,j,k,bi,bj),
     U        co3(i,j,k,bi,bj),hco3(i,j,k,bi,bj),
     U        co2aq(i,j,k,bi,bj),
     I        i,j,k,bi,bj,pco2_init_iter,myThid )
              ENDIF
             ENDDO
#else
             IF ( maskC(i,j,kLev,bi,bj) .NE. 0. _d 0) THEN
C$TAF init dic_surf = static, 10
               CALL CALC_PCO2_APPROX(
     I        surftemp(i,j),surfsalt(i,j),
     I        surfdic(i,j), surfphos(i,j),
     I        surfsi(i,j),surfalk(i,j),
     I        ak1(i,j,bi,bj),ak2(i,j,bi,bj),
     I        ak1p(i,j,bi,bj),ak2p(i,j,bi,bj),ak3p(i,j,bi,bj),
     I        aks(i,j,bi,bj),akb(i,j,bi,bj),akw(i,j,bi,bj),
     I        aksi(i,j,bi,bj),akf(i,j,bi,bj),
     I        ak0(i,j,bi,bj), fugf(i,j,bi,bj),
     I        ff(i,j,bi,bj),
     I        bt(i,j,bi,bj),st(i,j,bi,bj),ft(i,j,bi,bj),
     U        pH(i,j,bi,bj),pCO2(i,j,bi,bj),
     U        co3(i,j,bi,bj),hco3(i,j,bi,bj),co2aq(i,j,bi,bj),
     I        i,j,kLev,bi,bj,pco2_init_iter,myThid )
             ENDIF
#endif /* DIC_CARBON_COMPONENTS */
            ENDDO
           ENDDO
          ENDDO 
C          iprt = MIN(20,sNx)
C          jprt = MIN(20,sNy)
C          WRITE(standardMessageUnit,*) 'QQ first guess pH',
C     &        pH(iprt,jprt,bi,bj),
C     &        theta(iprt,jprt,1,bi,bj), salt(iprt,jprt,1,bi,bj),
C     &        surfdic(iprt,jprt), surfphos(iprt,jprt),
C     &        surfsi(iprt,jprt),surfalk(iprt,jprt)
C          CALL PLOT_FIELD_XYRL(pH, 
C     &         'First Guess pH (dic_surfforcing_init)', 1, myThid)
        ENDIF /* .NOT.pH_isLoaded */
C     end bi,bj loops
       ENDDO
      ENDDO

#endif /* ALLOW_DIC */
      RETURN
      END
