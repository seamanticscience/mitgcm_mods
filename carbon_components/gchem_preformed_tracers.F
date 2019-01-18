# include "GCHEM_OPTIONS.h"
#ifdef ALLOW_DIC
# include "DIC_OPTIONS.h"
#endif
#ifdef ALLOW_BLING
# include "BLING_OPTIONS.h"
#endif
#ifdef ALLOW_DARWIN
# include "DARWIN_OPTIONS.h"
#endif


CBOP
C !ROUTINE: GCHEM_FORCING_SEP
C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_PREFORMED_TRACERS( dic_offset, bi, bj, 
     &          myTime, myIter, myThid )

C !DESCRIPTION:
C     calls subroutine that will update passive preformed tracer values
C     with a separate timestep. 

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#include "PTRACERS_PARAMS.h"
#include "GCHEM.h"
#ifdef ALLOW_DIC
# include "DIC_VARS.h"
#endif /* ALLOW_DIC */
#ifdef ALLOW_BLING
# include "BLING_VARS.h"
#endif /* ALLOW_BLING */
#ifdef ALLOW_DARWIN
# include "DARWIN_FLUX.h"
# include "DARWIN_SIZE.h"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      _RL myTime
      INTEGER myIter, myThid
      INTEGER dic_offset
      INTEGER bi,bj
CEOP

#ifdef ALLOW_GCHEM
#ifdef GCHEM_SEPARATE_FORCING
#if ((defined GCHEM_PREFORMED) || (defined GCHEM_PREFORMED_MLD))
C!LOCAL VARIABLES: ====================================================
C  i,j                  :: loop indices
C  bi,bj                :: tile indices
C  k                    :: vertical level
      INTEGER iMin,iMax,jMin,jMax,i,j,k
      INTEGER nTrac
      INTEGER kmld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

#ifdef DIC_CARBON_COMPONENTS
      _RL csat_atm
#endif

      CHARACTER*(MAX_LEN_FNAM) fldName

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('GCHEM_PREFORMED_TRACERS',myThid)
#endif

      iMin=1-Olx
      iMax=sNx+Olx
      jMin=1-Oly
      jMax=sNy+Oly

C Get arrays of mixed layer depth for preformed tracers
C       DO bj = myByLo(myThid), myByHi(myThid)
C        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=jMin,jMax
          DO i=iMin,iMax
           IF ( usingPCoords ) THEN
            kmld(i,j,bi,bj) = Nr
           ELSE
            kmld(i,j,bi,bj) = 1
           ENDIF
#ifdef GCHEM_PREFORMED_TRACERS_MLD
           DO k=1,Nr
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 
     &           .AND. abs(RF(1+kmld(i,j,bi,bj))) 
     &                .LT.hMixLayer(i,j,bi,bj) 
     &          ) THEN
             kmld(i,j,bi,bj)=kmld(i,j,bi,bj)+1
            ENDIF
           ENDDO
#endif
          ENDDO
         ENDDO
C        ENDDO
C       ENDDO

#ifdef DIC_CARBON_COMPONENTS
C Call once to fill in 3d carbon chemistry values
C    below the surface layer (that's updated in 
C    GCHEM_FORCING_SEP>DIC_BIOTIC_FORCING>DIC_SURFFORCING)
C Treat surface and depth values seperately to be 
C    consistent with dic_surfforcing
C These routines take Theta and Salt from DYNVARS.h and
C    store in the 3d coefficient variables (see DIC_VARS.h)
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_CALL(
     &            'CARBON_COEFFS_SURF',myThid)
#endif
        DO k=2,Nr
          CALL CARBON_COEFFS_SURF(
     I            bi,bj,iMin,iMax,jMin,jMax,
     I            k,myThid)
        ENDDO
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_CALL(
     &            'CARBON_COEFFS_DEPTH',myThid)
#endif
C Now add pressure dependency below the surface
        DO k=2,Nr
          CALL CARBON_COEFFS_DEPTH(
     I            bi,bj,iMin,iMax,jMin,jMax,
     I            k,myThid)
        ENDDO
#endif

#ifdef ALLOW_DIC
C JML Do Preformed Tracers associated with the DIC package 
      IF ( gchem_sepFTr_num .GT. dic_offset ) THEN
       DO nTrac=dic_offset+1,gchem_sepFTr_num
        fldName=PTRACERS_names(nTrac)
                       
        IF ( fldName(1:7) .EQ. 'passive' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &      'EVALUATING PASSIVE TRACER',myThid)
#endif 
C Passive tracer set, or relax to a constant concentration at the surface (template)        
         DO j=jMin,jMax
          DO i=iMin,iMax
           DO k=1,kmld(i,j,bi,bj)
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
             pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &           +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &             *(1. _d 0-pTracer(i,j,k,bi,bj,nTrac))
     &             *_hFacC(i,j,k,bi,bj))
     &             *PTRACERS_dTLev(k)
            ELSE
             pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
            ENDIF
           ENDDO
          ENDDO
         ENDDO
                       
        ELSEIF ( fldName(1:4) .EQ. 'cpre' ) THEN 
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING CPRE TRACER',myThid)
#endif     
C Preformed Carbon, set or relax to surface DIC concentration      
         DO j=jMin,jMax
          DO i=iMin,iMax
           DO k=1,kmld(i,j,bi,bj)
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
             pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &           +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &             *(pTracer(i,j,k,bi,bj,iDIC)
     &                    -pTracer(i,j,k,bi,bj,nTrac))
     &             *_hFacC(i,j,k,bi,bj))
     &             *PTRACERS_dTLev(k)  
            ELSE
		     pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
		    ENDIF               
           ENDDO
          ENDDO
         ENDDO
            
        ELSEIF ( fldName(1:4) .EQ. 'apre' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING APRE TRACER',myThid)
#endif 
C Preformed Alkalinity, set or relax to surface ALK concentration
         DO j=jMin,jMax
          DO i=iMin,iMax
           DO k=1,kmld(i,j,bi,bj)
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
             pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &           +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &             *(pTracer(i,j,k,bi,bj,iALK)
     &                    -pTracer(i,j,k,bi,bj,nTrac))
     &             *_hFacC(i,j,k,bi,bj))
     &             *PTRACERS_dTLev(k)  
            ELSE
		     pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
		    ENDIF   
           ENDDO
          ENDDO
         ENDDO

         iAPRE=nTrac 
            
        ELSEIF ( fldName(1:4) .EQ. 'ppre' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING PPRE TRACER',myThid)
#endif 
C Preformed Phosphate, set or relax to surface PO4 concentration    
         DO j=jMin,jMax
          DO i=iMin,iMax
           DO k=1,kmld(i,j,bi,bj)
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
              pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &            +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &              *(pTracer(i,j,k,bi,bj,iPO4)
     &                     -pTracer(i,j,k,bi,bj,nTrac))
     &              *_hFacC(i,j,k,bi,bj))
     &              *PTRACERS_dTLev(k)   
            ELSE
		     pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
			ENDIF  
           ENDDO
          ENDDO
         ENDDO

         iPPRE=nTrac   
                      
        ELSEIF ( fldName(1:4) .EQ. 'opre' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING OPRE TRACER',myThid)
#endif 
C Preformed Oxygen, set or relax to surface O2 concentration            
         DO j=jMin,jMax
          DO i=iMin,iMax
           DO k=1,kmld(i,j,bi,bj)
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
             pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &           +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &             *(pTracer(i,j,k,bi,bj,iDO2)
     &                      -pTracer(i,j,k,bi,bj,nTrac))
     &             *_hFacC(i,j,k,bi,bj))
     &             *PTRACERS_dTLev(k) 
           ELSE
		    pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
		   ENDIF    
          ENDDO
         ENDDO
        ENDDO

        ELSEIF ( fldName(1:6) .EQ.'wm_age') THEN 
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &          'EVALUATING AGE TRACER',myThid)
#endif   
C Watermass Age Tracer, set or relax to zero seconds at the surface           
         DO j=jMin,jMax
          DO i=iMin,iMax
           DO k=1,kmld(i,j,bi,bj)
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
             pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &           +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &             *(0. _d 0-pTracer(i,j,k,bi,bj,nTrac))
     &             *_hFacC(i,j,k,bi,bj))
     &             *PTRACERS_dTLev(k)  
           ELSE
		    pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
		   ENDIF   
          ENDDO
         ENDDO
        ENDDO

        ELSEIF ( fldName(1:7) .EQ. 'atmpco2' ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &      'EVALUATING ATMPCO2 TRACER',myThid)
#endif 
C Atmospheric CO2 tracer, set surface to current atmpco2
C   Records what the atmosphere was last time the watermass was at the 
C   surface - may be useful for calculating anthropogenic CO2
C Set or relax to surface values
         DO j=jMin,jMax
          DO i=iMin,iMax
           DO k=1,kmld(i,j,bi,bj)
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
             pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &           +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &             *(AtmospCO2(i,j,bi,bj)*1. _d 6
     &                      -pTracer(i,j,k,bi,bj,nTrac))
     &             *_hFacC(i,j,k,bi,bj))
     &             *PTRACERS_dTLev(k) 
            ELSE
			 pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
			ENDIF                 
           ENDDO
          ENDDO
         ENDDO

         iPCO2=nTrac
            
        ELSEIF ( fldName(1:4) .EQ. 'csat' ) THEN
C Csat tracer records saturated dic in the surface ocean 
C Surface/mixed layer nutrients should already be preformed
#ifdef DIC_CARBON_COMPONENTS
#ifdef ALLOW_DEBUG
         IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING CSAT TRACER',myThid)
#endif 
		 DO j=jMin,jMax
		  DO i=iMin,iMax
		   DO k=1,kmld(i,j,bi,bj)
		    IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
             CALL CALC_CSAT(
     I           theta(i,j,k,bi,bj)  ,salt(i,j,k,bi,bj),
     I           AtmospCO2(i,j,bi,bj),pTracer(i,j,k,bi,bj,iPO4),
     I           SILICA(i,j,k,bi,bj) ,pTracer(i,j,k,bi,bj,iALK),
     U           pH_glob(i,j,k,bi,bj),csat_atm,
     I           i,j,k,bi,bj,myIter,myThid )

C Set or relax to surface values 
             pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &           +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &             *(csat_atm-pTracer(i,j,k,bi,bj,nTrac))
     &             *_hFacC(i,j,k,bi,bj))
     &             *PTRACERS_dTLev(k) 
			ELSE
			 pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
			ENDIF
		   ENDDO
		  ENDDO
		 ENDDO  
#endif /* DIC_CARBON_COMPONENTS */

        ELSEIF ( fldName(1:4) .EQ. 'cdis' ) THEN
C Disequilibrium tracer recording difference between saturated dic in the surface ocean 
C   and atmospheric pCO2 last time the watermass was at the surface
#ifdef DIC_CARBON_COMPONENTS
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_MSG(
     &         'EVALUATING CDIS TRACER',myThid)
#endif  
         DO j=jMin,jMax
          DO i=iMin,iMax
           DO k=1,kmld(i,j,bi,bj)
            IF ( maskC(i,j,k,bi,bj).NE.0. _d 0 ) then 
             CALL CALC_CSAT(
     I           theta(i,j,k,bi,bj)  ,salt(i,j,k,bi,bj),
     I           AtmospCO2(i,j,bi,bj),pTracer(i,j,k,bi,bj,iPO4),
     I           SILICA(i,j,k,bi,bj) ,pTracer(i,j,k,bi,bj,iALK),
     U           pH_glob(i,j,k,bi,bj),csat_atm,
     I           i,j,k,bi,bj,myIter,myThid )

C This occurs AFTER dic and atmospheric pco2 is updated with air-sea fluxes, i.e.
C    after anthropogenic CO2 uptake or outgassing has occurred 
C Set or relax to surface values
             pTracer(i,j,k,bi,bj,nTrac)=pTracer(i,j,k,bi,bj,nTrac)
     &           +  (1. _d 0/PTRACERS_prefSurfTauRelax(nTrac)
     &             *((pTracer(i,j,k,bi,bj,iDIC)-csat_atm)
     &                      -pTracer(i,j,k,bi,bj,nTrac))
     &             *_hFacC(i,j,k,bi,bj))
     &             *PTRACERS_dTLev(k) 
            ELSE
             pTracer(i,j,k,bi,bj,nTrac)=0. _d 0
            ENDIF
           ENDDO
          ENDDO
         ENDDO          
#endif /* DIC_CARBON_COMPONENTS */
        ENDIF 
       ENDDO
      ENDIF
      
#if ((defined DIC_CARBON_COMPONENTS) && (defined DIC_COMPONENT_DIAGS))
C          WRITE(standardMessageUnit,*) 
C     &      'GCHEM_PREFORMED_TRACERS: iPPRE,iAPRE,iPCO2=',
C     &      iPPRE,iAPRE,iPCO2
       
C Then do Csat, Csoft, Ccarb and Cdis calculation            
        IF (iAPRE.NE.UNSET_I.AND.iPPRE.NE.UNSET_I
#ifdef DIC_INSITU_COMPONENTS       
     &                      .AND.iPCO2.NE.UNSET_I
#endif /* DIC_INSITU_COMPONENTS */
     &                      ) THEN

C Reference pCO2 can be read in as dic_pCO2_ini set in the first call of DIC_ATMOS
C  OR read in as dic_pCO2_ref in data.dic (e.g. for restarts!)
C Current pCO2 can be directly taken from AtmospCO2
C Using pCO2 tracer allows insitu calculation of Csat (useful for anthropogenic emissions)
C set guess of pH for first step here
C          WRITE(standardMessageUnit,*)  'GCHEM_PREFORMED_TRACERS:',
C     &     'Calling CALC_DIC_COMPONENTS'

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_CALL('CALC_DIC_COMPONENTS',myThid)
#endif          
           CALL CALC_DIC_COMPONENTS(
     &         pTracer(1-OLx,1-OLy,Nr,bi,bj,iDIC), 
     &         pTracer(1-OLx,1-OLy,Nr,bi,bj,iALK), 
     &         pTracer(1-OLx,1-OLy,Nr,bi,bj,iPO4),
     &         pTracer(1-OLx,1-OLy,Nr,bi,bj,iDOP),
     &         pTracer(1-OLx,1-OLy,Nr,bi,bj,iAPRE),
     &         pTracer(1-OLx,1-OLy,Nr,bi,bj,iPPRE),
#ifdef DIC_INSITU_COMPONENTS       
     &         pTracer(1-OLx,1-OLy,Nr,bi,bj,iPCO2),
#endif /* DIC_INSITU_COMPONENTS */
     &         dic_pCO2_ini(1-OLx,1-OLy,bi,bj),
     &         bi,bj,iMin,iMax,jMin,jMax,
     &         myIter,myTime,myThid)
        ENDIF
#endif /* DIC_CARBON_COMPONENTS */
#endif /* ALLOW_DIC */        
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('GCHEM_PREFORMED_TRACERS',myThid)
#endif
#endif /* GCHEM_PREFORMED */
#endif /* GCHEM_SEPARATE_FORCING */
#endif /* ALLOW_GCHEM */

      RETURN
      END
