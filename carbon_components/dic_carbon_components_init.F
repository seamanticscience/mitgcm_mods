C $Header: /u/gcmpack/MITgcm/pkg/dic/dic_carbon_components_init.F,v 1.33 2011/10/07 21:36:39 dfer Exp $
C $Name:  $

#include "DIC_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"

CBOP
C !ROUTINE: DIC_SURFFORCING_INIT

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_CARBON_COMPONENTS_INIT(
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
       INTEGER i,j,k,kLev,it
       INTEGER intimeP, intime0, intime1
       _RL aWght, bWght, co3dummy
C Number of iterations for pCO2 solvers...
C Solubility relation coefficients
C local variables for carbon chem
      _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      LOGICAL pH_3disLoaded
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP
#ifdef DIC_CARBON_COMPONENTS

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C PCO2 already initiated in dic_surfforcing_init called just prior
C      CALL DIC_ATMOS(0, startTime, nIter0, myThid )

C Initialize 3d coefficient fields
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          DO k=1,Nr
            fugf3d(i,j,k,bi,bj)= 0. _d 0
            ff3d(i,j,k,bi,bj)  = 0. _d 0
            ak03d(i,j,k,bi,bj) = 0. _d 0
            ak13d(i,j,k,bi,bj) = 0. _d 0
            ak23d(i,j,k,bi,bj) = 0. _d 0
            akb3d(i,j,k,bi,bj) = 0. _d 0
            ak1p3d(i,j,k,bi,bj)= 0. _d 0
            ak2p3d(i,j,k,bi,bj)= 0. _d 0
            ak3p3d(i,j,k,bi,bj)= 0. _d 0
            aksi3d(i,j,k,bi,bj)= 0. _d 0
            akw3d(i,j,k,bi,bj) = 0. _d 0
            aks3d(i,j,k,bi,bj) = 0. _d 0
            akf3d(i,j,k,bi,bj) = 0. _d 0
            bt3d(i,j,k,bi,bj)  = 0. _d 0
            st3d(i,j,k,bi,bj)  = 0. _d 0
            ft3d(i,j,k,bi,bj)  = 0. _d 0
            Ksp_TP_Calc3d(i,j,k,bi,bj) = 0. _d 0
          ENDDO  
         ENDDO
        ENDDO
       ENDDO
      ENDDO

       DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          iMin=1-Olx
          iMax=sNx+Olx
          jMin=1-Oly
          jMax=sNy+Oly

#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL(
     &            'CARBON_COEFFS_SURF',myThid)
#endif
          CALL CARBON_COEFFS_SURF(
     I              bi,bj,iMin,iMax,jMin,jMax,myThid)
       
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL(
     &           'CARBON_COEFFS_DEPTH',myThid)
#endif
          CALL CARBON_COEFFS_DEPTH(
     I                   bi,bj,iMin,iMax,jMin,jMax,
     I                   2,Nr,myThid)
         ENDDO
       ENDDO 

      pH_3disLoaded = .FALSE.
      IF ( nIter0.GT.PTRACERS_Iter0 .OR.
     &    (nIter0.EQ.PTRACERS_Iter0 .AND. pickupSuff.NE.' ')
     &   ) THEN
C       Read pH from a pickup file if needed
        CALL DIC_COMPONENTS_READ_PICKUP(
     O                        pH_3disLoaded,
     I                        nIter0, myThid )
      ENDIF

      IF ( .NOT.pH_3disLoaded ) THEN
C set guess of pH for first step here
         WRITE(standardMessageUnit,*) 
     &     'DIC_CARBON_COMPONENTS_INIT:',
     &     'Setting 3d pH fields to first guess'
C          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
C     &                      SQUEEZE_RIGHT, myThid )
       DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k=1,Nr
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
              IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then
                pH_GLOB(i,j,k,bi,bj) = 8. _d 0
#ifdef DIC_INSITU_COMPONENTS       
                pH_INSIT(i,j,k,bi,bj)= 8. _d 0
#endif
                pH_INIT(i,j,k,bi,bj) = 8. _d 0   
              ELSE
                pH_GLOB(i,j,k,bi,bj) = 0. _d 0
#ifdef DIC_INSITU_COMPONENTS       
                pH_INSIT(i,j,k,bi,bj)= 0. _d 0
#endif
                pH_INIT(i,j,k,bi,bj) = 0. _d 0
              ENDIF     
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF

#endif /* DIC_CARBON_COMPONENTS */
#endif /* ALLOW_DIC */
      RETURN
      END
