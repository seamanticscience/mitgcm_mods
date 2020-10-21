C $Header: /u/gcmpack/MITgcm/pkg/dic/dic_biotic_forcing.F,v 1.32 2016/01/16 21:57:12 jmc Exp $
C $Name:  $

#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: DIC_BIOTIC_FORCING

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_BIOTIC_FORCING(
     U                      PTR_DIC, PTR_ALK, PTR_PO4, PTR_DOP,
#ifdef ALLOW_O2
     U                      PTR_O2,
#endif
#ifdef ALLOW_FE
     &                      PTR_FE,
#ifdef ALLOW_VARIABLE_LIGANDS
     &                      PTR_LIG,
#endif
#endif
     I                      bi, bj, iMin, iMax, jMin, jMax,
     I                      myIter, myTime, myThid )

C !DESCRIPTION:
C updates all the tracers for the effects of air-sea exchange, biological
c activity and remineralization

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DIC_VARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"

C !INPUT/OUTPUT PARAMETERS: ===================================================
C  PTR_DIC              :: dissolced inorganic carbon
C  PTR_ALK              :: alkalinity
C  PTR_PO4              :: phosphate
c  PTR_DOP              :: dissolve organic phosphurous
c  PTR_O2               :: oxygen
C  PTR_FE               :: iron
c  bi, bj               :: current tile indices
C  myIter               :: current timestep
C  myTime               :: current time
C  myThid               :: thread number
      _RL  PTR_DIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_ALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_PO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  PTR_DOP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef ALLOW_O2
      _RL  PTR_O2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
#ifdef ALLOW_FE
      _RL  PTR_FE (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef ALLOW_VARIABLE_LIGANDS
      _RL  PTR_LIG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
#endif
      INTEGER bi, bj, iMin, iMax, jMin, jMax
      INTEGER myIter
      _RL myTime
      INTEGER myThid

#ifdef ALLOW_PTRACERS
#ifdef DIC_BIOTIC

C !LOCAL VARIABLES: ====================================================
C  i,j,k                  :: loop indices
C  G*                     :: tendency term for the tracers
C  SURA                   :: tendency of alkalinity due to freshwater
C  SURC                   :: tendency of DIC due to air-sea exchange
C                            and virtual flux
C  SURO                   :: tendency of O2 due to air-sea exchange
C  GPO4                   :: tendency of PO4 due to biological productivity,
C                            exchange with DOP pool and reminerization
C  CAR                    :: carbonate changes due to biological
C                             productivity and remineralization
C  BIOac                  :: biological productivity
C  RDOP                   :: DOP sink due to remineralization
C  pflux                  :: changes to PO4 due to flux and remineralization
C  CAR_S                  :: carbonate sink
C  cflux                  :: carbonate changes due to flux and remineralization
C  freefe                 :: iron not bound to ligand
      _RL  GDIC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GALK(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GPO4(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  GDOP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  SURA(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURC(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SURO(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  CAR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  BIOac(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  RDOP(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  pflux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  exportflux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  CAR_S(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  cflux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#ifdef ALLOW_O2
      _RL  GO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
#ifdef ALLOW_FE
      _RL  GFE(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  SURFE(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  SEDFE(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  HYDROFE(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  freefe(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  scav_poc

#ifdef SUBSURFACE_FE_SOURCES
      _RL sed_poc
#endif

#ifdef ALLOW_VARIABLE_LIGANDS
      _RL  GLIG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  lambda_lig(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  lig_incr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL  lig_decr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
#endif
#endif
       INTEGER i,j,k,kBottom,kSurface
#ifdef CAR_DISS
       INTEGER nCALCITEstep
#endif
CEOP

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('DIC_BIOTIC_FORCING',myThid)
#endif

      IF ( useThSIce .OR. useSEAICE .OR. useCoupler ) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('DIC_FIELDS_UPDATE',myThid)
#endif
        CALL DIC_FIELDS_UPDATE(
     I                    bi, bj, myTime, myIter, myThid )
      ENDIF

C kBottom is updated later on
      IF ( usingPCoords ) THEN
       kSurface = Nr
       kBottom = 1
      ELSE
       kSurface = 1
       kBottom = Nr
      ENDIF

       DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           RDOP(i,j,k)     =0. _d 0
           GDIC(i,j,k)     =0. _d 0
           GALK(i,j,k)     =0. _d 0
           GPO4(i,j,k)     =0. _d 0
           GDOP(i,j,k)     =0. _d 0
           CAR(i,j,k)      =0. _d 0
           BIOac(i,j,k)    =0. _d 0
           pflux(i,j,k)    =0. _d 0
           exportflux(i,j,k)=0. _d 0
           cflux(i,j,k)    =0. _d 0
           CAR_S(i,j,k)    =0. _d 0
#ifdef ALLOW_O2
           GO2(i,j,k)      =0. _d 0
#endif
#ifdef ALLOW_FE
           GFE (i,j,k)     =0. _d 0
           freefe(i,j,k)   =0. _d 0
#ifdef ALLOW_VARIABLE_LIGANDS
           GLIG(i,j,k)     =0. _d 0
C  no longer needed after adding full initialisation of freefe in S/R FE_CHEM
c          freefe(i,j,k)  =0. _d 0
           lig_incr(i,j,k) =0. _d 0
           lig_decr(i,j,k) =0. _d 0
#endif
#endif
          ENDDO
         ENDDO
       ENDDO
       DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           SURA(i,j)   =0. _d 0
           SURC(i,j)   =0. _d 0
           SURO(i,j)   =0. _d 0
#ifdef ALLOW_FE
           SURFE(i,j)  =0. _d 0
           SEDFE(i,j)  =0. _d 0
           HYDROFE(i,j)=0. _d 0
#endif           
          ENDDO
         ENDDO

C carbon air-sea interaction
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('DIC_SURFFORCING',myThid)
#endif
       CALL DIC_SURFFORCING(
     I                    PTR_DIC, PTR_ALK, PTR_PO4,
     O                    SURC,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )

C alkalinity air-sea interaction
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('ALK_SURFFORCING',myThid)
#endif
       CALL ALK_SURFFORCING(
     I                    PTR_ALK,
     O                    SURA,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )

#ifdef ALLOW_O2
C oxygen air-sea interaction
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('O2_SURFFORCING',myThid)
#endif
       CALL O2_SURFFORCING(
     I                    PTR_O2,
     O                    SURO,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
#endif

#ifdef ALLOW_FE
C find free iron
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('FE_CHEM',myThid)
#endif
       CALL FE_CHEM( bi, bj, iMin, iMax, jMin, jMax,
     U               PTR_FE,
#ifdef ALLOW_VARIABLE_LIGANDS
     I               PTR_LIG, 
#endif
     O               freefe,
     I               myIter, myThid )
#endif

C biological activity
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('BIO_EXPORT',myThid)
#endif
       CALL BIO_EXPORT(
     I               PTR_PO4,
#ifdef ALLOW_FE
     I               PTR_FE,
#endif
     O               BIOac,
     I               bi, bj, iMin, iMax, jMin, jMax,
     I               myIter, myTime, myThid )

C flux of po4 from layers with biological activity
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('PHOS_FLUX',myThid)
#endif
       CALL PHOS_FLUX(
     I               BIOac,
     U               pflux, exportflux,
     I               bi, bj, iMin, iMax, jMin, jMax,
     I               myIter, myTime, myThid )

C- Carbonate sink
       DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
             CAR_S(i,j,k)=BIOac(i,j,k)*R_CP*rain_ratio(i,j,bi,bj)*
     &                    (1. _d 0-DOPfraction)
          ENDDO
         ENDDO
       ENDDO

C carbonate
#ifdef CAR_DISS
C dissolution only below saturation horizon
C code following method by Karsten Friis
        nCALCITEstep = 3600
        IF(myIter .lt. (nIter0+5) .or.
     &               mod(myIter,nCALCITEstep) .eq. 0)THEN
#ifdef ALLOW_DEBUG
          IF (debugMode) CALL DEBUG_CALL('CALCITE_SATURATION',myThid)
#endif
          CALL CALCITE_SATURATION(
     I                    PTR_DIC, PTR_ALK, PTR_PO4,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
        ENDIF

#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('CAR_FLUX_OMEGA_TOP',myThid)
#endif
        CALL CAR_FLUX_OMEGA_TOP(
     I                    BIOac,
     O                    cflux,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
#else
C old OCMIP way
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('CAR_FLUX',myThid)
#endif
        CALL CAR_FLUX(
     I                    CAR_S,
     U                    cflux,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myIter, myTime, myThid )
#endif

C add all tendencies for PO4, DOP, ALK, DIC
       DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
#ifdef DIC_NO_NEG
           RDOP(i,j,k)= MAX(maskC(i,j,k,bi,bj)*KDOPRemin*PTR_DOP(i,j,k)
     &                     ,0. _d 0)
#else
           RDOP(i,j,k)= maskC(i,j,k,bi,bj)*KDOPRemin*PTR_DOP(i,j,k)
#endif
           GPO4(i,j,k)=-BIOac(i,j,k)+pflux(i,j,k) + RDOP(i,j,k)

           car(i,j,k) = cflux(i,j,k) - CAR_S(i,j,k)

           GDOP(i,j,k)=+BIOac(i,j,k)*DOPfraction - RDOP(i,j,k)

           GALK(i,j,k)=+2. _d 0 *car(i,j,k)-R_NP*GPO4(i,j,k)

           GDIC(i,j,k)=car(i,j,k)+R_CP*GPO4(i,j,k)

#ifdef ALLOW_O2
           if (PTR_O2(i,j,k).GT.O2crit) then
             GO2(i,j,k)= R_OP*GPO4(i,j,k)
           else
             GO2(i,j,k)= 0. _d 0
           endif
#endif
#ifdef ALLOW_FE
#ifdef ALLOW_VARIABLE_LIGANDS
C Prognostic Ligand code...
C Here Ligand is produced by excretion during production (proportional to bioac)
C     and release of iron binding sites during remineralization.           
           lig_incr(i,j,k)  = (pflux(i,j,k)+BIOac(i,j,k)*DOPfraction)
     &                             *R_CP*gamma_lig

C Now calculate timescale for degredation 
C  Initially, the first order way is kind of cheating as it fixes lambda to go exactly
C towards the ligand concentration that we want 
C           lambda_lig(i,j,k)= lig_incr(i,j,k)/ligand_tot
C 
C Calculate depth-dependent degradation timescale
C     constrained by gamma_over_lambda, global productivity/surface area  
C     and nominal ligand concentration
C           lambda_lig(i,j,k)= gamma_lig/(-rC(k)*goverl_lig)
C
C So try to create a more dynamic consumption timescales using Q10
C  lambda_ref calculated in dic_init_fixed from lambda_over_gamma
C  and the value of gamma_lig
           lambda_lig(i,j,k) = lambda_ref*
     &        (lig_q10**((theta(i,j,k,bi,bj)-lig_tref)/10. _d 0))

C Ligand concentration decreases according to the Q10 rate
C Make sure this is always negative        
           lig_decr(i,j,k) = MIN(-lambda_lig(i,j,k)
     &                        *PTR_LIG(i,j,k)*maskC(i,j,k,bi,bj)
C                             may also decrease as DOP is consumed  
C     &                        -gamma_lig*R_CP*RDOP(i,j,k)
     &                        ,0. _d 0)
                
C Now do ligand tendency
           GLIG(i,j,k)=lig_incr(i,j,k)+lig_decr(i,j,k)
#endif /* ALLOW_VARIABLE_LIGANDS */

           IF ( KScav_poc.GT.0. _d 0 ) THEN
CJML       Use particle depended scavenging rate, exploiting exportflux
C          and constant sinking rate (See Parekh etal 2005, GBC equation 11)
               IF ( k.EQ.kSurface ) THEN
                   KScav=KScav_surf
               ELSE
C          Calculate POC concentration from sinking particles
                   IF ( poc_wsp_factor.GT.0. _d 0 ) THEN
C Constant sinking speed if already set, else use a linear 
C   increase may be consistent with the Martin Curve
                       poc_wsp = -rC(k) * poc_wsp_factor
     &                                  + poc_wsp_offset
                   ENDIF
C          Export flux is in mol P.m-3.s-1 and we need mg C.l-1, R_pop2poc does this               
                   scav_poc=MAX(0. _d 0,exportflux(i,j,k)/poc_wsp)
     &                          *R_pop2poc
C          Scavenging rate per second     
                   KScav=KScav_poc*scav_ratio*(scav_poc**scav_exp)            
               ENDIF
           ENDIF
           
C          Modify the scavenging rate with other mechanisms (bkgrnd, dust) 
           KScav=KScav
     &             +KScav_background
     &             +KScav_dust*(
     &             (InputFe(i,j,bi,bj)*recip_hFacC(i,j,kSurface,bi,bj)
     &             /(R_dust2fe*dust_wsp)))
           
           GFE(i,j,k) = R_FeP*GPO4(i,j,k)
     &                        -KScav*freefe(i,j,k)
#endif /* ALLOW_FE */
          ENDDO
         ENDDO
       ENDDO

C Surface Sources and sinks
       DO j=jMin,jMax
         DO i=iMin,iMax         
               GALK(i,j,kSurface)=GALK(i,j,kSurface)+SURA(i,j)
               GDIC(i,j,kSurface)=GDIC(i,j,kSurface)+SURC(i,j)
#ifdef ALLOW_O2
               GO2(i,j,kSurface) =GO2(i,j,kSurface)+SURO(i,j)
#endif
#ifdef ALLOW_FE
       IF ( alpfe.NE.0. _d 0 ) THEN
               SURFE(i,j)=alpfe*InputFe(i,j,bi,bj)
     &                   *recip_drF(kSurface)
     &                   *recip_hFacC(i,j,kSurface,bi,bj)
       ELSE
               SURFE(i,j)=0. _d 0
       ENDIF
       
               GFE(i,j,kSurface)=GFE(i,j,kSurface)+SURFE(i,j)

#ifdef SUBSURFACE_FE_SOURCES
       kBottom   = MAX(kLowC(i,j,bi,bj),1)
       
       IF ( -rC(kBottom).LE.fe_sed_depth_max ) THEN         
          IF ( selectFeSedSource.EQ.1 ) THEN
C Elrod et al (2004) formulation 
C include sediment iron source using the flux of po4 into bottom layer            
C JML pflux*dz*RCP is mol C.m-2.s-1, sedfe_factor converts to mol Fe.m-2.s-1   
               sed_poc=MAX(0. _d 0,pflux(i,j,kBottom))
     &                          *R_CP*drF(kBottom)

               SEDFE(i,j)  = MAX( 
     &              0. _d 0, sedfe_factor * sed_poc + sedfe_offset
     &            )*recip_drF(kBottom)*recip_hFacC(i,j,kBottom,bi,bj)
          ELSEIF ( selectFeSedSource.EQ.2 ) THEN
C Dale et al (2015) formulation. 
C In this paper, gamma = 170 umol Fe m-2 d-1 (170/(86400*1e6) mol Fe m-2 s-1)
C But that leads to HUGE sediment fluxes (>100 Gmol/yr) so here we reduce the
C gamma to ~20 umol Fe m-2 d-1 (20/(86400*1e6) mol Fe m-2 s-1)
C NOTE UNITS: sed_poc is mmol m-2 d-1 (from mol C.m-2.s-1) and O2 is umol l-1
               sed_poc=MAX(0. _d 0,pflux(i,j,kBottom))
     &                          *R_CP*drF(kBottom)
     &                          *1000*86400 

               SEDFE(i,j)=sedfe_factor
     &            * tanh( sed_poc / PTR_O2(i,j,kBottom)*1000 )
     &            *recip_drF(kBottom)*recip_hFacC(i,j,kBottom,bi,bj)
          ENDIF
       ENDIF    

       IF ( selectFeHTVSource.EQ.1 .AND.       
     &           -rC(kBottom).GE.fe_vent_depth_min ) THEN      
C JML include iron source due to hydrothermal input into bottom layer   
C                                         He3 flux is mol.m-2.s-1
               HYDROFE(i,j)=R_FeHe3*solfe*HydroInputHe3(i,j,bi,bj)
     &             *recip_drF(kBottom)
     &             *recip_hFacC(i,j,kBottom,bi,bj)
       ENDIF

               GFE(i,j,kBottom)=GFE(i,j,kBottom)
     &            +SEDFE(i,j)+HYDROFE(i,j)
#endif /* SUBSURFACE_FE_SOURCES */     
#endif /* ALLOW_FE */
          ENDDO
         ENDDO

       IF ( useOBCS ) THEN
        DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
            GDIC(i,j,k) = GDIC(i,j,k)*maskInC(i,j,bi,bj)
            GALK(i,j,k) = GALK(i,j,k)*maskInC(i,j,bi,bj)
            GPO4(i,j,k) = GPO4(i,j,k)*maskInC(i,j,bi,bj)
            GDOP(i,j,k) = GDOP(i,j,k)*maskInC(i,j,bi,bj)
#ifdef ALLOW_O2
            GO2(i,j,k)  = GO2(i,j,k)*maskInC(i,j,bi,bj)
#endif
#ifdef ALLOW_FE
            GFE(i,j,k)  = GFE(i,j,k)*maskInC(i,j,bi,bj)
#ifdef ALLOW_VARIABLE_LIGANDS
            GLIG(i,j,k) = GLIG(i,j,k)*maskInC(i,j,bi,bj)
#endif
#endif
          ENDDO
         ENDDO
        ENDDO
       ENDIF

C update
       DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
           PTR_DIC(i,j,k)=
     &      PTR_DIC(i,j,k)+GDIC(i,j,k)*PTRACERS_dTLev(k)
           PTR_ALK(i,j,k)=
     &      PTR_ALK(i,j,k)+GALK(i,j,k)*PTRACERS_dTLev(k)
           PTR_PO4(i,j,k)=
     &      PTR_PO4(i,j,k)+GPO4(i,j,k)*PTRACERS_dTLev(k)
           PTR_DOP(i,j,k)=
     &      PTR_DOP(i,j,k)+GDOP(i,j,k)*PTRACERS_dTLev(k)
#ifdef ALLOW_O2
           PTR_O2(i,j,k)=
     &      PTR_O2(i,j,k)+GO2(i,j,k)*PTRACERS_dTLev(k)
#endif
#ifdef ALLOW_FE
           PTR_FE(i,j,k)=
     &      PTR_FE(i,j,k)+GFE(i,j,k)*PTRACERS_dTLev(k)
     
#ifdef ALLOW_VARIABLE_LIGANDS
C JML make sure ligand remains positive
           PTR_LIG(i,j,k)=MAX(0. _d 0,
     &      PTR_LIG(i,j,k)+GLIG(i,j,k)*PTRACERS_dTLev(k)
     &      )
#endif
#endif
          ENDDO
         ENDDO
       ENDDO

#ifdef ALLOW_FE
#ifdef MINFE
c find free iron and get rid of insoluble part
#ifdef ALLOW_DEBUG
       IF (debugMode) CALL DEBUG_CALL('FE_CHEM',myThid)
#endif
       CALL FE_CHEM( bi, bj, iMin, iMax, jMin, jMax,
     U               PTR_FE,
#ifdef ALLOW_VARIABLE_LIGANDS
     I               PTR_LIG,
#endif
     O               freefe,
     I               myIter, myThid )
#endif
#endif

#ifdef ALLOW_TIMEAVE
C save averages
       IF ( PTRACERS_taveFreq.GT.0. ) THEN
        DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
            BIOave(i,j,k,bi,bj)   =BIOave(i,j,k,bi,bj)+
     &                             BIOac(i,j,k)*deltaTClock
            CARave(i,j,k,bi,bj)   =CARave(i,j,k,bi,bj)+
     &                             CAR(i,j,k)*deltaTClock
            OmegaCave(i,j,k,bi,bj)=OmegaCave(i,j,k,bi,bj)+
     &                             OmegaC(i,j,k,bi,bj)*deltaTClock
            pfluxave(i,j,k,bi,bj) =pfluxave(i,j,k,bi,bj) +
     &                             pflux(i,j,k)*deltaTClock
            epfluxave(i,j,k,bi,bj)=epfluxave(i,j,k,bi,bj) +
     &                             exportflux(i,j,k)*deltaTClock
            cfluxave(i,j,k,bi,bj) =cfluxave(i,j,k,bi,bj) +
     &                             cflux(i,j,k)*deltaTClock
          ENDDO
         ENDDO
        ENDDO
         DO j=jMin,jMax
          DO i=iMin,iMax
              SURave(i,j,bi,bj)    =SURave(i,j,bi,bj)+
     &                              SURC(i,j)*deltaTClock
#ifdef ALLOW_O2
              SUROave(i,j,bi,bj)   =SUROave(i,j,bi,bj)+
     &                              SURO(i,j)*deltaTClock
#endif
              pCO2ave(i,j,bi,bj)   =pCO2ave(i,j,bi,bj)+
     &                              pCO2(i,j,bi,bj)*deltaTClock
              pHave(i,j,bi,bj)     =pHave(i,j,bi,bj)+
     &                              pH(i,j,bi,bj)*deltaTClock
              fluxCO2ave(i,j,bi,bj)=fluxCO2ave(i,j,bi,bj)+
     &                           fluxCO2(i,j,bi,bj)*deltaTClock
          ENDDO
         ENDDO
         DIC_timeAve(bi,bj) = DIC_timeAve(bi,bj)+deltaTClock
       ENDIF
#endif /* ALLOW_TIMEAVE*/

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN

        CALL DIAGNOSTICS_FILL(BIOac     ,'DICBIOA ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(CAR       ,'DICCARB ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(OmegaC    ,'OMEGAC  ',0,Nr,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pflux     ,'DICPFLUX',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(exportflux,'CEXPFLX ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(cflux     ,'CARBFLUX',0,Nr,2,bi,bj,myThid)

        CALL DIAGNOSTICS_FILL(SURC      ,'DICTFLX ',0,1 ,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pCO2      ,'DICPCO2 ',0,1 ,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(pH        ,'DICPHAV ',0,1 ,1,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(fluxCO2   ,'DICCFLX ',0,1 ,1,bi,bj,myThid)

        CALL DIAGNOSTICS_FILL(SURA      ,'DICAFLX ',0,1 ,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(RDOP      ,'DICRDOP ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(GDIC      ,'DICGDIC ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(GDOP      ,'DICGDOP ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(GPO4      ,'DICGPO4 ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(GALK      ,'DICGALK ',0,Nr,2,bi,bj,myThid)
       
#ifdef ALLOW_O2
        CALL DIAGNOSTICS_FILL(SURO      ,'DICOFLX ',0,1 ,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(GO2       ,'DICGO2  ',0,Nr,2,bi,bj,myThid)
#endif
#ifdef ALLOW_FE
        CALL DIAGNOSTICS_FILL(SURFE     ,'DICFEFLX',0,1 ,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(GFE       ,'DICGFE  ',0,Nr,2,bi,bj,myThid)
#if SUBSURFACE_FE_SOURCES  
        CALL DIAGNOSTICS_FILL(SEDFE     ,'DICSEDFE',0,1 ,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(HYDROFE   ,'DICHTVFE',0,1 ,2,bi,bj,myThid)
#endif /* SUBSURFACE_FE_SOURCES */
#ifdef ALLOW_VARIABLE_LIGANDS
        CALL DIAGNOSTICS_FILL(GLIG      ,'DICGLIG ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(lambda_lig,'LIGLAMDA',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(lig_incr  ,'LIGINCR ',0,Nr,2,bi,bj,myThid)
        CALL DIAGNOSTICS_FILL(lig_decr  ,'LIGDECR ',0,Nr,2,bi,bj,myThid)
#endif /* ALLOW_VARIABLE_LIGANDS */
#endif /* ALLOW_FE */
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('DIC_BIOTIC_FORCING',myThid)
#endif

#endif /* DIC_BIOTIC */
#endif /* ALLOW_PTRACERS */

       RETURN
       END
