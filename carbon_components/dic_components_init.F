#include "DIC_OPTIONS.h"
#include "PTRACERS_OPTIONS.h"

CBOP
C !ROUTINE: DIC_COMPONENTS_INIT

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_COMPONENTS_INIT(
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
       INTEGER i,j,k,bi,bj
      LOGICAL pH_3disLoaded
CEOP
#ifdef DIC_CARBON_COMPONENTS
      INTEGER pco2_init_iter 
      _RL csatdummy
      _RL si_pre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C PCO2 already initiated in dic_surfforcing_init called just prior
C      CALL DIC_ATMOS(0, startTime, nIter0, myThid )
C 3d coefficients also initiated in dic_surfforcing_init

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
                pH_glob(i,j,k,bi,bj) = 8. _d 0
#ifdef DIC_INSITU_COMPONENTS       
                pH_loc(i,j,k,bi,bj)= 8. _d 0
#endif
                pH_ini(i,j,k,bi,bj) = 8. _d 0   
              ELSE
                pH_glob(i,j,k,bi,bj) = 0. _d 0
#ifdef DIC_INSITU_COMPONENTS       
                pH_loc(i,j,k,bi,bj)= 0. _d 0
#endif
                pH_ini(i,j,k,bi,bj) = 0. _d 0
              ENDIF     
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

C Initialize the pH values by running the CSAT solver with the correct atm co2
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k=1,Nr
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
             IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then  
C Silica is not interactive, so assume a SI:P ratio 15 for diatoms.               
              si_pre(i,j,k,bi,bj)=SILICA(i,j,k,bi,bj)
     &                       -(R_SIP*(PTRACER(i,j,k,bi,bj,iPO4)
     &                           -PTRACER(i,j,k,bi,bj,iPPRE)))
              DO pco2_init_iter=1,10
C ph_glob goes with the current atmospheric pCO2 and preformed tracers             
               CALL CALC_CSAT(
     I           theta(i,j,k,bi,bj)  ,salt(i,j,k,bi,bj),
     I           AtmospCO2(i,j,bi,bj),PTRACER(i,j,k,bi,bj,iPPRE),
     I           si_pre(i,j,k,bi,bj) ,PTRACER(i,j,k,bi,bj,iAPRE),
     U           pH_glob(i,j,k,bi,bj),csatdummy,
     I           i,j,k,bi,bj,pco2_init_iter,myThid )

#ifdef DIC_INSITU_COMPONENTS            
C ph_loc goes with the "local" atmospheric pCO2 (tracer) and preformed tracers             
               CALL CALC_CSAT(
     I           theta(i,j,k,bi,bj) ,salt(i,j,k,bi,bj),
     I           PTRACER(i,j,k,bi,bj,iPCO2)*1. _d -6,
     I           PTRACER(i,j,k,bi,bj,iPPRE),
     I           si_pre(i,j,k,bi,bj),PTRACER(i,j,k,bi,bj,iAPRE),
     U           pH_loc(i,j,k,bi,bj),csatdummy,
     I           i,j,k,bi,bj,pco2_init_iter,myThid )
#endif

C ph_ini goes with the initial atmospheric pCO2 (dic_pCO2_ini) and preformed tracers             
               CALL CALC_CSAT(
     I           theta(i,j,k,bi,bj)     ,salt(i,j,k,bi,bj),
     I           dic_pCO2_ini(i,j,bi,bj),PTRACER(i,j,k,bi,bj,iPPRE),
     I           si_pre(i,j,k,bi,bj)    ,PTRACER(i,j,k,bi,bj,iAPRE),
     U           pH_ini(i,j,k,bi,bj)    ,csatdummy,
     I           i,j,k,bi,bj,pco2_init_iter,myThid )
              ENDDO
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
