#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: CALC_DIC_COMPONENTS

C !INTERFACE: ==========================================================
      SUBROUTINE CALC_DIC_COMPONENTS(ptr_dic, ptr_alk, 
     &           ptr_po4, ptr_dop,
     &           ptr_apre, ptr_ppre,
#ifdef DIC_INSITU_COMPONENTS       
     &           atm_co2_loc,
#endif
     &           atm_co2_ini, 
     &           bi,bj,iMin,iMax,jMin,jMax,
     &           myIter,myTime,myThid)
C !DESCRIPTION:
C Calculate Williams and Follows (2011) Carbon component breakdown
C Diagnostic csat_atm (at intitial pco2, current pco2 and pco2
C    at the last surface exposure), csoft, ccarb, cres_atm, p_reg, a_reg

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
       _RL  ptr_dic(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL  ptr_alk(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL  ptr_po4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL  ptr_dop(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL  ptr_apre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL  ptr_ppre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
#ifdef DIC_INSITU_COMPONENTS             
       _RL  atm_co2_loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
#endif
       _RL  atm_co2_ini(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj)
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
       _RL csat_ini(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL cres_ini(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
#ifdef DIC_INSITU_COMPONENTS       
       _RL csat_loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL cres_loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL tmp
#endif
       _RL csat_atm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL cres_atm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)      
       _RL csoft (1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL ccarb (1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL p_reg (1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL a_reg (1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)
       _RL si_pre(1-OLx:sNx+OLx,1-OLy:sNy+OLy,bi,bj,Nr)

C Initialize vars to zero
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
#ifdef DIC_INSITU_COMPONENTS       
          csat_loc(i,j,k,bi,bj) = 0. _d 0
          cres_loc(i,j,k,bi,bj) = 0. _d 0
#endif
          csat_ini(i,j,k,bi,bj) = 0. _d 0
          cres_ini(i,j,k,bi,bj) = 0. _d 0
          csat_atm(i,j,k,bi,bj) = 0. _d 0
          cres_atm(i,j,k,bi,bj) = 0. _d 0
          csoft   (i,j,k,bi,bj) = 0. _d 0
          ccarb   (i,j,k,bi,bj) = 0. _d 0
          p_reg   (i,j,k,bi,bj) = 0. _d 0
          a_reg   (i,j,k,bi,bj) = 0. _d 0
          si_pre  (i,j,k,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO

      DO k=1,Nr
       DO j=jMin,jMax
        DO i=iMin,iMax
         IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then
C Calculate regenerated components
            p_reg(i,j,k,bi,bj)=ptr_po4(i,j,k,bi,bj)
     &                           -ptr_ppre(i,j,k,bi,bj)
            a_reg(i,j,k,bi,bj)=ptr_alk(i,j,k,bi,bj)
     &                           -ptr_apre(i,j,k,bi,bj)
C Silica is not interactive, so assume a SI:P ratio 15 for diatoms.               
            si_pre(i,j,k,bi,bj)=Silica(i,j,k,bi,bj)
     &                       -(R_SIP*p_reg(i,j,k,bi,bj))

C calculate csoft and ccarb - Include DOP in the csoft reservoir
            csoft(i,j,k,bi,bj) = R_CP*( p_reg(i,j,k,bi,bj) 
     &                     + ptr_dop(i,j,k,bi,bj) )
            ccarb(i,j,k,bi,bj) = 0.5 _d 0*( a_reg(i,j,k,bi,bj) 
     &                     + R_NP*p_reg(i,j,k,bi,bj) )
         ENDIF
        ENDDO
       ENDDO
      ENDDO
       
C csat referenced to current atmospheric pCO2
      IF ( DIAGNOSTICS_IS_ON( 'C_satATM' , myThid ) ) THEN
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
          IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then

            CALL CALC_CSAT(
     I          theta(i,j,k,bi,bj)  ,salt(i,j,k,bi,bj),
     I          AtmospCO2(i,j,bi,bj),ptr_ppre(i,j,k,bi,bj),
     I          si_pre(i,j,k,bi,bj) ,ptr_apre(i,j,k,bi,bj),
     U          pH_glob(i,j,k,bi,bj),csat_atm(i,j,k,bi,bj),
     I          i,j,k,bi,bj,myIter,myThid )

            cres_atm (i,j,k,bi,bj) = ptr_dic(i,j,k,bi,bj)
     &                     -csat_atm(i,j,k,bi,bj)
     &                     -csoft   (i,j,k,bi,bj)
     &                     -ccarb(i,j,k,bi,bj) 

C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: DIC=',ptr_dic(i,j,k,bi,bj),
C     &      ' csat_atm=',csat_atm(i,j,k,bi,bj),' csoft=',csoft(i,j,k,bi,bj),
C     &      ' ccarb=',ccarb(i,j,k,bi,bj),' cres_atm=',cres_atm(i,j,k,bi,bj)  
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDIF

C csat referenced to starting atmospheric pCO2
      IF ( DIAGNOSTICS_IS_ON( 'C_satINI' , myThid ) ) THEN
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax       
          IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then      

            CALL CALC_CSAT(
     I          theta(i,j,k,bi,bj)    ,salt(i,j,k,bi,bj),
     I          atm_co2_ini(i,j,bi,bj),ptr_ppre(i,j,k,bi,bj),
     I          si_pre(i,j,k,bi,bj)   ,ptr_apre(i,j,k,bi,bj),
     U          pH_ini(i,j,k,bi,bj)   ,csat_ini(i,j,k,bi,bj),
     I          i,j,k,bi,bj,myIter,myThid )

            cres_ini (i,j,k,bi,bj) = ptr_dic(i,j,k,bi,bj)
     &                     -csat_ini(i,j,k,bi,bj)
     &                     -csoft(i,j,k,bi,bj)
     &                     -ccarb(i,j,k,bi,bj) 

C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: DIC=',ptr_dic(i,j,k,bi,bj),
C     &      ' csat_ini=',csat_ini(i,j,k,bi,bj),' csoft=',csoft(i,j,k,bi,bj),
C     &      ' ccarb=',ccarb(i,j,k,bi,bj),' cres_ini=',cres_ini(i,j,k,bi,bj)  
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDIF

#ifdef DIC_INSITU_COMPONENTS
C csat referenced to atmospheric pCO2 at last surface exposure (in situ)
      IF ( DIAGNOSTICS_IS_ON( 'C_satLOC' , myThid ) ) THEN
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax       
          IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 

            tmp=atm_co2_loc(i,j,k,bi,bj)*1. _d -6

            CALL CALC_CSAT(
     I          theta(i,j,k,bi,bj) ,salt(i,j,k,bi,bj),
     I          tmp                ,ptr_ppre(i,j,k,bi,bj),
     I          si_pre(i,j,k,bi,bj),ptr_apre(i,j,k,bi,bj),
     U          pH_loc(i,j,k,bi,bj),csat_loc(i,j,k,bi,bj),
     I          i,j,k,bi,bj,myIter,myThid )

            cres_loc (i,j,k,bi,bj) = ptr_dic(i,j,k,bi,bj)
     &                     -csat_loc(i,j,k,bi,bj)
     &                     -csoft(i,j,k,bi,bj)
     &                     -ccarb(i,j,k,bi,bj) 

C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: DIC=',ptr_dic(i,j,k,bi,bj),
C     &      ' csat_loc=',csat_loc(i,j,k,bi,bj),' csoft=',csoft(i,j,k,bi,bj),
C     &      ' ccarb=',ccarb(i,j,k,bi,bj),' cres_loc=',cres_loc(i,j,k,bi,bj)  
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDIF
#endif

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       IF ( DIAGNOSTICS_IS_ON('C_satINI', myThid ) ) THEN
C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: writing C_satINI diags'
        CALL DIAGNOSTICS_FILL(csat_ini,'C_satINI',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(cres_ini,'C_resINI',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH_ini  ,'pH_ini  ',0,Nr,2,bi,bj,myThid)
       ENDIF

       IF ( DIAGNOSTICS_IS_ON('C_satATM', myThid ) ) THEN   
C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: writing csat_atm diags'       
        CALL DIAGNOSTICS_FILL(csat_atm,'C_satATM',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(cres_atm,'C_resATM',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH_GLOB ,'pH_glob ',0,Nr,2,bi,bj,myThid)
       ENDIF

#ifdef DIC_INSITU_COMPONENTS
       IF ( DIAGNOSTICS_IS_ON('C_satLOC', myThid ) ) THEN  
C             WRITE(standardMessageUnit,*) 
C     &      'DIC_CARBON_COMPONENTS: writing CSATSIT diags'         
        CALL DIAGNOSTICS_FILL(csat_loc,'C_satLOC',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(cres_loc,'C_resLOC',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH_loc  ,'pH_loc  ',0,Nr,2,bi,bj,myThid)
       ENDIF
#endif

        CALL DIAGNOSTICS_FILL(csoft  ,'C_soft  ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(ccarb  ,'C_carb  ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(p_reg  ,'P_reg   ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(a_reg  ,'A_reg   ',0,Nr,2,bi,bj,myThid)
C Pass out Silica to check everything is OK
        CALL DIAGNOSTICS_FILL(si_pre ,'Si_pre  ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(Silica ,'Silica  ',0,Nr,2,bi,bj,myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
#endif /* DIC_CARBON_COMPONENTS */
#endif /* DIC_BIOTIC */
#endif /* ALLOW_PTRACERS */
       RETURN
       END