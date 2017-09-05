C $Header: /u/gcmpack/MITgcm/pkg/dic/dic_carbon_components.F,v 1.11 2011/10/07 21:36:39 dfer Exp $
C $Name:  $

#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: CALC_DIC_COMPONENTS

C !INTERFACE: ==========================================================
      SUBROUTINE CALC_DIC_COMPONENTS(PTR_DIC, PTR_ALK, 
     &           PTR_PO4, PTR_DOP,
     &           PTR_APRE, PTR_PPRE,
#ifdef DIC_INSITU_COMPONENTS       
     &           ATM_CO2_INSIT,
#endif
     &           ATM_CO2_INIT, 
     &           bi,bj,iMin,iMax,jMin,jMax,
     &           myIter,myTime,myThid)
C !DESCRIPTION:
C Calculate Williams and Follows (2011) Carbon component breakdown
C Diagnostic Csat (at intitial pco2, current pco2 and pco2
C    at the last surface exposure), Csoft, Ccarb, Cres, Preg, Areg

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DIC_VARS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
C  bioac                :: biological productivity
       _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL  PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL  PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL  PTR_DOP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL  PTR_APRE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL  PTR_PPRE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef DIC_INSITU_COMPONENTS             
       _RL  ATM_CO2_INSIT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
       _RL  ATM_CO2_INIT(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
       INTEGER iMin, iMax, jMin, jMax, bi, bj
       _RL myTime
       INTEGER myIter, myThid

      LOGICAL     DIAGNOSTICS_IS_ON
      EXTERNAL    DIAGNOSTICS_IS_ON
C !OUTPUT PARAMETERS: ===================================================

#ifdef ALLOW_PTRACERS
#ifdef DIC_BIOTIC
#if ((defined DIC_CARBON_COMPONENTS) && (defined DIC_COMPONENT_DIAGS))
C !LOCAL VARIABLES: ====================================================
C  i,j,k                  :: loop indices
c  ko                     :: loop-within-loop index
c depth_u, depth_l        :: depths of upper and lower interfaces
c flux_u, flux_l          :: flux through upper and lower interfaces
c zbase                   :: depth of bottom of current productive layer
       INTEGER i,j,k
       _RL co3dummy
       _RL pCO2local
       _RL pHlocal
       _RL csatlocal
       _RL CSAT_INIT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL CRES_INIT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef DIC_INSITU_COMPONENTS       
       _RL CSAT_INSIT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL CRES_INSIT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL tmp
#endif
       _RL CSAT (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL CSOFT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL CCARB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL CRES (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL PREG (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL AREG (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
       _RL SIPRE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

C Initialize vars to zero
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
#ifdef DIC_INSITU_COMPONENTS       
          CSAT_INSIT(i,j,k) = 0. _d 0
          CRES_INSIT(i,j,k) = 0. _d 0
#endif
          CSAT_INIT (i,j,k) = 0. _d 0
          CRES_INIT (i,j,k) = 0. _d 0
          CSAT      (i,j,k) = 0. _d 0
          CRES      (i,j,k) = 0. _d 0
          CSOFT     (i,j,k) = 0. _d 0
          CCARB     (i,j,k) = 0. _d 0
          PREG      (i,j,k) = 0. _d 0
          AREG      (i,j,k) = 0. _d 0
          SIPRE     (i,j,k) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO

      DO k=1,Nr
       DO j=jMin,jMax
        DO i=iMin,iMax
         IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then
C Calculate regenerated components
            PREG(i,j,k)=PTR_PO4(i,j,k)-PTR_PPRE(i,j,k)
            AREG(i,j,k)=PTR_ALK(i,j,k)-PTR_APRE(i,j,k)
C Silica is not interactive, so assume a SI:P ratio 15 for diatoms.               
            SIPRE(i,j,k)=Silica(i,j,k,bi,bj)
     &                       -(R_SIP*PREG(i,j,k))

C calculate Csoft and Ccarb - Include DOP in the CSOFT reservoir
            CSOFT(i,j,k)      = R_CP*( PREG(i,j,k) 
     &                     + PTR_DOP(i,j,k) )
            CCARB(i,j,k)      = 0.5 _d 0*( AREG(i,j,k) 
     &                     + R_NP*PREG(i,j,k) )
         ENDIF
        ENDDO
       ENDDO
      ENDDO
       
C Csat referenced to current atmospheric pCO2
      IF ( DIAGNOSTICS_IS_ON( 'CSAT    ' , myThid ) ) THEN
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
          IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then

            CALL CALC_CSAT_APPROX(
     I          theta(i,j,k,bi,bj),salt(i,j,k,bi,bj),
     I          AtmospCO2(i,j,bi,bj),PTR_PPRE(i,j,k),
     I          SIPRE(i,j,k),PTR_APRE(i,j,k),
     I          ak13d(i,j,k,bi,bj),ak23d(i,j,k,bi,bj),
     I          ak1p3d(i,j,k,bi,bj),ak2p3d(i,j,k,bi,bj),
     I          ak3p3d(i,j,k,bi,bj),
     I          aks3d(i,j,k,bi,bj),akb3d(i,j,k,bi,bj),
     I          akw3d(i,j,k,bi,bj),
     I          aksi3d(i,j,k,bi,bj),akf3d(i,j,k,bi,bj),
     I          ak03d(i,j,k,bi,bj), fugf3d(i,j,k,bi,bj),
     I          ff3d(i,j,k,bi,bj),
     I          bt3d(i,j,k,bi,bj),st3d(i,j,k,bi,bj),ft3d(i,j,k,bi,bj),
     U          pH_GLOB(i,j,k,bi,bj),CSAT(i,j,k),co3dummy,
     I          i,j,k,bi,bj,myIter,myThid )

            CRES (i,j,k) = PTR_DIC(i,j,k)-CSAT(i,j,k)
     &                     -CSOFT(i,j,k)-CCARB(i,j,k) 

C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: DIC=',PTR_DIC(i,j,k),
C     &      ' CSAT=',CSAT(i,j,k),' CSOFT=',CSOFT(i,j,k),
C     &      ' CCARB=',CCARB(i,j,k),' CRES=',CRES(i,j,k)  
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDIF

C Csat referenced to starting atmospheric pCO2
      IF ( DIAGNOSTICS_IS_ON( 'CSATINIT' , myThid ) ) THEN
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax       
          IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then      

            CALL CALC_CSAT_APPROX(
     I          theta(i,j,k,bi,bj),salt(i,j,k,bi,bj),
     I          ATM_CO2_INIT(i,j) ,PTR_PPRE(i,j,k),
     I          SIPRE(i,j,k),PTR_APRE(i,j,k),
     I          ak13d(i,j,k,bi,bj),ak23d(i,j,k,bi,bj),
     I          ak1p3d(i,j,k,bi,bj),ak2p3d(i,j,k,bi,bj),
     I          ak3p3d(i,j,k,bi,bj),
     I          aks3d(i,j,k,bi,bj),akb3d(i,j,k,bi,bj),
     I          akw3d(i,j,k,bi,bj),
     I          aksi3d(i,j,k,bi,bj),akf3d(i,j,k,bi,bj),
     I          ak03d(i,j,k,bi,bj), fugf3d(i,j,k,bi,bj),
     I          ff3d(i,j,k,bi,bj),
     I          bt3d(i,j,k,bi,bj),st3d(i,j,k,bi,bj),ft3d(i,j,k,bi,bj),
     U          pH_INIT(i,j,k,bi,bj),CSAT_INIT(i,j,k),co3dummy,
     I          i,j,k,bi,bj,myIter,myThid )

            CRES_INIT (i,j,k) = PTR_DIC(i,j,k)-CSAT_INIT(i,j,k)
     &                          -CSOFT(i,j,k)-CCARB(i,j,k) 

C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: DIC=',PTR_DIC(i,j,k),
C     &      ' CSAT_INIT=',CSAT_INIT(i,j,k),' CSOFT=',CSOFT(i,j,k),
C     &      ' CCARB=',CCARB(i,j,k),' CRES_INIT=',CRES_INIT(i,j,k)  
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDIF

#ifdef DIC_INSITU_COMPONENTS
C Csat referenced to atmospheric pCO2 at last surface exposure (in situ)
      IF ( DIAGNOSTICS_IS_ON( 'CSATSIT ' , myThid ) ) THEN
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax       
          IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
C            IF ( ATM_CO2_INSIT(i,j,k).LT.50. _d 0 ) then
C                CSAT_INSIT(i,j,k)=CSAT_INIT(i,j,k)
C            ELSE
            tmp=ATM_CO2_INSIT(i,j,k)*1. _d -6

            CALL CALC_CSAT_APPROX(
     I          theta(i,j,k,bi,bj),salt(i,j,k,bi,bj),
     I          tmp,PTR_PPRE(i,j,k),
     I          SIPRE(i,j,k),PTR_APRE(i,j,k),
     I          ak13d(i,j,k,bi,bj),ak23d(i,j,k,bi,bj),
     I          ak1p3d(i,j,k,bi,bj),ak2p3d(i,j,k,bi,bj),
     I          ak3p3d(i,j,k,bi,bj),
     I          aks3d(i,j,k,bi,bj),akb3d(i,j,k,bi,bj),
     I          akw3d(i,j,k,bi,bj),
     I          aksi3d(i,j,k,bi,bj),akf3d(i,j,k,bi,bj),
     I          ak03d(i,j,k,bi,bj), fugf3d(i,j,k,bi,bj),
     I          ff3d(i,j,k,bi,bj),
     I          bt3d(i,j,k,bi,bj),st3d(i,j,k,bi,bj),ft3d(i,j,k,bi,bj),
     U          pH_INSIT(i,j,k,bi,bj),CSAT_INSIT(i,j,k),co3dummy,
     I          i,j,k,bi,bj,myIter,myThid )
C            ENDIF

                CRES_INSIT (i,j,k) = PTR_DIC(i,j,k)-CSAT_INSIT(i,j,k)
     &                     -CSOFT(i,j,k)-CCARB(i,j,k) 

C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: DIC=',PTR_DIC(i,j,k),
C     &      ' CSAT_INSIT=',CSAT_INSIT(i,j,k),' CSOFT=',CSOFT(i,j,k),
C     &      ' CCARB=',CCARB(i,j,k),' CRES_INSIT=',CRES_INSIT(i,j,k)  
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDIF
#endif

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       IF ( DIAGNOSTICS_IS_ON( 'CSATINIT' , myThid ) ) THEN
C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: writing CSATINIT diags'
        CALL DIAGNOSTICS_FILL(CSAT_INIT ,'CSATINIT',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(CRES_INIT ,'CRESINIT',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH_INIT   ,'pH_INIT ',0,Nr,2,bi,bj,myThid)
       ENDIF

       IF ( DIAGNOSTICS_IS_ON( 'CSAT    ' , myThid ) ) THEN   
C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: writing CSAT diags'       
        CALL DIAGNOSTICS_FILL(CSAT      ,'CSAT    ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(CRES      ,'CRES    ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH_GLOB   ,'pH_GLOB ',0,Nr,2,bi,bj,myThid)
       ENDIF

#ifdef DIC_INSITU_COMPONENTS
       IF ( DIAGNOSTICS_IS_ON( 'CSATSIT ' , myThid ) ) THEN  
C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: writing CSATSIT diags'         
        CALL DIAGNOSTICS_FILL(CSAT_INSIT,'CSATSIT ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(CRES_INSIT,'CRESSIT ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH_INSIT  ,'pH_INSIT',0,Nr,2,bi,bj,myThid)
       ENDIF
#endif

        CALL DIAGNOSTICS_FILL(CSOFT     ,'CSOFT   ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(CCARB     ,'CCARB   ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(PREG      ,'PREG    ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(AREG      ,'AREG    ',0,Nr,2,bi,bj,myThid)
C Pass out Silica to check everything is OK
        CALL DIAGNOSTICS_FILL(SIPRE     ,'SIPRE   ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(Silica    ,'Silica  ',0,Nr,2,bi,bj,myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
#endif /* DIC_CARBON_COMPONENTS */
#endif /* DIC_BIOTIC */
#endif /* ALLOW_PTRACERS */
       RETURN
       END