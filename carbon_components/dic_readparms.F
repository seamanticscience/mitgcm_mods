C $Header: /u/gcmpack/MITgcm/pkg/dic/dic_readparms.F,v 1.10 2010/04/11 20:54:49 jmc Exp $
C $Name:  $

#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: DIC_READPARMS
C !INTERFACE: ==========================================================
      SUBROUTINE DIC_READPARMS( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R DIC_READPARMS
C     | o Initialise and read dic package parameters
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIC_VARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid    :: My Thread Id. number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIC

C     === Local variables ===
C     msgBuf    :: Informational/error message buffer
C     iUnit     :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C-- Abiotic dic parameters:
C permil   :: set carbon mol/m3 <---> mol/kg conversion factor
C             default permil = 1024.5 kg/m3
C Pa2Atm   :: Conversion factor for atmospheric pressure pLoad (when coupled
C             to atmospheric model) into Atm. Default assumes pLoad in Pascal
C             1 Atm = 1.01325e5 Pa = 1013.25 mb

      NAMELIST /ABIOTIC_PARMS/ permil, Pa2Atm

#ifdef DIC_BIOTIC

C-- Biotic dic parameters:
C   DOPfraction :: fraction of new production going to DOP
C   KDOPRemin   :: DOP remineralization rate (1/s) = 1/(6 month)
C   KRemin      :: remin power law coeff
C   zcrit       :: Minimum Depth (m) over which biological activity
C                  is computed --> determines nlev as the indice of the
C                  first layer deeper than -zcrit
C   O2crit      :: critical oxygen level (mol/m3)
C   R_OP, R_CP  :: stochiometric ratios
C   R_NP, R_FeP, R_pop2poc
C   zca         :: scale depth for CaCO3 remineralization (m)
CC Parameters for light/nutrient limited bioac
C   parfrac     :: fraction of Qsw that is PAR
C   k0          :: light attentuation coefficient (1/m)
C   lit0        :: half saturation light constant (W/m2)
C   KPO4        :: half saturation phosphate constant (mol/m3)
C   KFE         :: half saturation fe constant (mol/m3)
CC Iron chemisty values
C   alpfe       :: solubility of aeolian fe
C   fesedflux_pcm :: ratio of sediment iron to sinking organic matter
C   FeIntSec    :: y-axis crossing for Fe_flux = fesedflux_pcm*pflux + FeIntSec
C   freefemax   :: max solubility of free iron (mol/m3)
CC Control variables
C   KScav       :: iron scavenging rate (s-1)
C   ligand_stab :: ligand-free iron stability constant (m3/mol)
C   ligand_tot  :: total free ligand  (mol/m3)
CC JML Added some extras for particle dependent Fe scavenging rates
CC See Parekh etal 2005, GBC equation 11
C   scav_surf_min :: surface scavenging rate (s-1)
C   scav_wsp    :: Constant particulate sinking rate (m/s)
C   scav_inter  :: ????
C   scav_exp    :: Honeyman exponent
C   alpha       :: timescale for biological activity
C                  read in alphaUniform and filled in 2d array alpha
C   rain_ratio  :: inorganic/organic carbon rain ratio
C                  read in rainRatioUniform and filled in 2d array rain_ratio

      NAMELIST /BIOTIC_PARMS/
     & DOPfraction, KDOPRemin, KRemin, zcrit, zca,
     & O2crit, R_OP, R_CP, R_NP, R_FeP, R_pop2poc, R_SIP,
     & parfrac, k0, lit0, KPO4, KFE, scav_exp,
     & alpfe, freefemax, scav_wsp, scav_inter, scav_ratio,
     & KScav, scav_surf_min, ligand_stab, ligand_tot, 
     & fesedflux_pcm, FeIntSec, solfe, R_FeHe3,
     & fe_vent_depth_min, fe_sed_depth_max,
#ifdef ALLOW_VARIABLE_LIGANDS     
     & gamma_lig, bio_tot_lig,
#endif
     & alphaUniform, rainRatioUniform 
#endif

      NAMELIST /DIC_FORCING/
     &          DIC_windFile, DIC_atmospFile, DIC_iceFile,
     &          DIC_ironFile, DIC_silicaFile, DIC_hydroventFile,
     &          DIC_forcingPeriod, DIC_forcingCycle,
     &          dic_int1, dic_int2, dic_int3, dic_int4, dic_pCO2
#ifdef ALLOW_C13
     &          , dic_p13co2, c13airseafac, c13biofac
#endif
#ifdef DIC_CARBON_COMPONENTS
     &          , dic_pCO2_ref
#endif
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _BEGIN_MASTER(myThid)

       permil      = 1. _d 0 / 1024.5 _d 0
       Pa2Atm      = 1.01325 _d 5

#ifdef DIC_BIOTIC
       DOPfraction = 0.67 _d 0
       KDOPRemin   = 1. _d 0/(6. _d 0*30. _d 0*86400. _d 0)
       KRemin      = 0.9 _d 0
       zcrit       = 500. _d 0
       O2crit      = 4. _d -3
       R_OP        =-170. _d 0
       R_CP        = 117. _d 0
       R_NP        = 16. _d 0
       R_FeP       = 0.000468 _d 0
       R_SIP       = 15. _d 0
       zca         = 3500. _d 0
       parfrac     = 0.4 _d 0
       k0          = 0.02 _d 0
       lit0        = 30. _d 0
       KPO4        = 5. _d -4
       KFE         = 1.2 _d -7
       alpfe       = 0.01 _d 0
       freefemax   = 3. _d -7
       alphaUniform     = 2. _d -3/(360. _d 0 * 86400. _d 0)
       rainRatioUniform = 7. _d -2
#ifdef ALLOW_VARIABLE_LIGANDS 
C JML Added vars for prognostic Ligands       
       gamma_lig   = 0. _d 0
       bio_tot_lig = 0. _d 0
#endif
       KScav       = 0.19 _d 0/(360. _d 0*86400. _d 0)
C JML Added var for particle dependent Fe scavenging
C JML R_pop2poc converts POP in mol P.m-3.s-1 to POC in mg C.l-1
       R_pop2poc   = R_CP * 12
       scav_surf_min = KScav
       scav_exp    = 0. _d 0
       scav_wsp    = 2900. _d 0/(360. _d 0*86400. _d 0)
C scav_inter is in l^scav_exp.mg^-scav_exp.d-1       
C (See Parekh etal 2005, GBC equation 11) covert to s-1 
       scav_ratio  = 1. _d 0   
       scav_inter  = 0.079 _d 0/86400. _d 0
       ligand_stab = 1. _d 8
       ligand_tot  = 1. _d -6
C These values for hydrothermal input of iron from Tagliabue et al (2014)       
       R_FeHe3     = 4.5 _d 8
       solfe       = 2. _d -3
       fe_vent_depth_min = 750. _d 0
C Hajoon's sedimentary iron source (JML different coeffs tho)
       fe_sed_depth_max = 6000. _d 0
C    For Elrod et al (2004) use:
       fesedflux_pcm = 6.8 _d -4       
       FeIntSec      = -0.5 _d -6 / 86400. _d 0
C    For Dale et al (2015) use:
C       fesedflux_pcm = 170. _d 0 / (86400. _d 0 * 1. _d 6)      
C       FeIntSec      = 0
#endif

#ifdef ALLOW_C13
       dic_p13co2 = 0. _d 0
       c13airseafac = 0. _d 0
       c13biofac = 0. _d 0
#endif

       DIC_windFile  = ' '
       DIC_atmospFile= ' '
       DIC_iceFile   = ' '
       DIC_ironFile  = ' '
       DIC_silicaFile= ' '
C JML include iron source due to hydrothermal input into bottom layer            
       DIC_hydroventFile=' '       
       dic_int1    = 0
       dic_int2    = 0
       dic_int3    = 0
       dic_int4    = 0
       dic_pCO2    = 278. _d -6
C Reference for saturated carbon diagnostic CSATINIT
#ifdef DIC_CARBON_COMPONENTS
       dic_pCO2_ref = UNSET_RL
#endif

c default periodic forcing to same as for physics
       DIC_forcingPeriod = externForcingPeriod
       DIC_forcingCycle  = externForcingCycle

      WRITE(msgBuf,'(A)') ' DIC_READPARMS: opening data.dic'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

      CALL OPEN_COPY_DATA_FILE( 'data.dic', 'DIC_READPARMS',
     O                          iUnit, myThid )

C--   Read parameters from open data file:

C-    Abiotic parameters
      READ(UNIT=iUnit,NML=ABIOTIC_PARMS)

#ifdef DIC_BIOTIC
C-    Biotic parameters
      READ(UNIT=iUnit,NML=BIOTIC_PARMS)
#endif

C-    forcing filenames and parameters
      READ(UNIT=iUnit,NML=DIC_FORCING)

      WRITE(msgBuf,'(A)')
     &   ' DIC_READPARMS: finished reading data.dic'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     I                    SQUEEZE_RIGHT, myThid )

C--   Close the open data file
      CLOSE(iUnit)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C-    derive other parameters:

#ifdef DIC_BIOTIC
       QSW_underice = .FALSE.
#ifdef USE_QSW_UNDERICE
       QSW_underice = .TRUE.
#elif (defined (USE_QSW))
C if using Qsw and seaice, then ice fraction is already
C taken into account
       IF ( useSEAICE ) QSW_underice = .TRUE.
       IF ( useThSIce ) QSW_underice = .TRUE.
#endif
#endif /* DIC_BIOTIC */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Print out parameter values :

      iUnit = standardMessageUnit
      WRITE(msgBuf,'(A)') ' '
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)') '// ==================================='
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)') '// DIC package parameters :'
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)
      WRITE(msgBuf,'(A)') '// ==================================='
      CALL PRINT_MESSAGE(msgBuf,iUnit,SQUEEZE_RIGHT,myThid)

C- namelist ABIOTIC_PARMS
       CALL WRITE_0D_RL( permil, INDEX_NONE,'permil =',
     &  ' /* Ref. density to convert mol/m3 to mol/kg */')
       CALL WRITE_0D_RL( Pa2Atm, INDEX_NONE,'Pa2Atm =',
     &  ' /* Atmosph. pressure conversion coeff (to Atm) */')

#ifdef DIC_BIOTIC
C- namelist BIOTIC_PARMS
       CALL WRITE_0D_RL( DOPfraction, INDEX_NONE,'DOPfraction =',
     &  ' /* Fraction of new production going to DOP */')
       CALL WRITE_0D_RL( KDOPRemin, INDEX_NONE,'KDOPRemin =',
     &  ' /* DOP remineralization rate (1/s) */')
       CALL WRITE_0D_RL( KRemin, INDEX_NONE,'KRemin =',
     &  ' /* Remin power law coeff. */')
       CALL WRITE_0D_RL( zcrit, INDEX_NONE,'zcrit =',
     &  ' /* Minimum depth for biological activity (m) */')
       CALL WRITE_0D_RL( O2crit, INDEX_NONE,'O2crit =',
     &  ' /* Critical oxygen level (mol/m3) */')
       CALL WRITE_0D_RL( R_OP, INDEX_NONE,'R_OP =',
     &  ' /* Stochiometric ratio R_OP */')
       CALL WRITE_0D_RL( R_CP, INDEX_NONE,'R_CP =',
     &  ' /* Stochiometric ratio R_CP */')
       CALL WRITE_0D_RL( R_NP, INDEX_NONE,'R_NP =',
     &  ' /* Stochiometric ratio R_NP */')
       CALL WRITE_0D_RL( R_FeP, INDEX_NONE,'R_FeP =',
     &  ' /* Stochiometric ratio R_FeP */')
       CALL WRITE_0D_RL( zca, INDEX_NONE,'zca =',
     &  ' /* Scale depth for CaCO3 remineralization (m) */')
       CALL WRITE_0D_RL( parfrac, INDEX_NONE,'parfrac =',
     &  ' /* Fraction of Qsw that is PAR */')
       CALL WRITE_0D_RL( k0, INDEX_NONE,'k0 =',
     &  ' /* Light attentuation coefficient (1/m) */')
       CALL WRITE_0D_RL( lit0, INDEX_NONE,'lit0 =',
     &  ' /* Half saturation light constant (W/m2) */')
       CALL WRITE_0D_RL( KPO4, INDEX_NONE,'KPO4 =',
     &  ' /* Half saturation phosphate constant (mol/m3) */')
       CALL WRITE_0D_RL( KFE, INDEX_NONE,'KFE =',
     &  ' /* Half saturation fe constant (mol/m3) */')
       CALL WRITE_0D_RL( alpfe, INDEX_NONE,'alpfe =',
     &  ' /* Solubility of aeolian fe */')
       CALL WRITE_0D_RL( freefemax, INDEX_NONE,'freefemax =',
     &  ' /* Max solubility of free iron (mol/m3) */')
       CALL WRITE_0D_RL( KScav, INDEX_NONE,'KScav =',
     &  ' /* Iron scavenging rate */')
#ifdef ALLOW_VARIABLE_LIGANDS
C JML Prognostic Ligand parameters
      IF ( gamma_lig.GT.0. _d 0 ) THEN
            CALL WRITE_0D_RL(gamma_lig, INDEX_NONE,'gamma_lig =',
     &  ' /* Fraction of biomass that produces ligands */')   
C            CALL WRITE_0D_RL(goverl_lig, INDEX_NONE,'goverl_lig =',
C     &  ' /* Gamma over Lambda, ratio of ligand coefficients */')      
      ENDIF
#endif
C JML write out particle depended scavenging parameters
      IF ( scav_exp.GT.0. _d 0 ) THEN
            CALL WRITE_0D_RL(scav_surf_min, INDEX_NONE,'scav_exp =',
     &  ' /* Honeyman exponent */')   
            CALL WRITE_0D_RL( scav_exp, INDEX_NONE,'scav_exp =',
     &  ' /* Honeyman exponent */')
            CALL WRITE_0D_RL( scav_inter, INDEX_NONE,'scav_inter =',
     &  ' /* Iron scavenging rate constant */')
            CALL WRITE_0D_RL( scav_wsp, INDEX_NONE,'scav_wsp',
     &  ' /* Inverse particle sinking speed (s/m) */')
       CALL WRITE_0D_RL( R_pop2poc, INDEX_NONE,'R_pop2poc =',
     &  ' /* Stochiometric ratio R_pop2poc */')
      ENDIF
       CALL WRITE_0D_RL( ligand_stab, INDEX_NONE,'ligand_stab =',
     &  ' /* Ligand-free iron stability constant (m3/mol) */')
       CALL WRITE_0D_RL( ligand_tot, INDEX_NONE,'ligand_tot =',
     &  ' /* Total free ligand  (mol/m3) */')
       CALL WRITE_0D_RL( alphaUniform, INDEX_NONE,'alphaUniform =',
     &  ' /* Timescale for biological activity */')
       CALL WRITE_0D_RL(rainRatioUniform,INDEX_NONE,'rainRatioUniform=',
     &  ' /* Inorganic/organic carbon rain ratio */')

       CALL WRITE_0D_L( QSW_underice, INDEX_NONE, 'QSW_underice  =',
     &  '  /* Flag for Qsw under Sea-Ice (i.e. SI fract included) */')
#endif

C- namelist DIC_FORCING
       CALL WRITE_0D_C( DIC_windFile, -1, INDEX_NONE, 'DIC_windFile =',
     & '  /* File name of wind speeds */')
       CALL WRITE_0D_C( DIC_atmospFile, -1,INDEX_NONE,'DIC_atmospFile=',
     & '  /* File name of atmospheric pressure*/')
       CALL WRITE_0D_C( DIC_iceFile, -1, INDEX_NONE, 'DIC_iceFile =',
     & '  /* File name of seaice fraction */')
       CALL WRITE_0D_C( DIC_ironFile, -1, INDEX_NONE, 'DIC_ironFile =',
     & '  /* File name of aeolian iron flux */')
       CALL WRITE_0D_C( DIC_silicaFile, -1,INDEX_NONE,'DIC_silicaFile=',
     & '  /* File name of surface silica */')
       CALL WRITE_0D_RL( DIC_forcingPeriod,
     &   INDEX_NONE,'DIC_forcingPeriod =',
     &  ' /* Periodic forcing parameter specific for DIC (s) */')
       CALL WRITE_0D_RL( DIC_forcingCycle,
     &   INDEX_NONE,'DIC_forcingCycle =',
     &  ' /* Periodic forcing parameter specific for DIC (s) */')
       CALL WRITE_0D_I( dic_int1, INDEX_NONE, 'dic_int1 =',
     &  '  /*  */')
       CALL WRITE_0D_I( dic_int2, INDEX_NONE, 'dic_int2 =',
     &  '  /*  */')
       CALL WRITE_0D_I( dic_int3, INDEX_NONE, 'dic_int3 =',
     &  '  /*  */')
       CALL WRITE_0D_I( dic_int4, INDEX_NONE, 'dic_int4 =',
     &  '  /*  */')
       CALL WRITE_0D_RL( dic_pCO2, INDEX_NONE,'dic_pCO2 =',
     &  ' /* Atmospheric pCO2 to be read in data.dic */')
#ifdef ALLOW_C13
       CALL WRITE_0D_RL( dic_p13co2, INDEX_NONE,'dic_p13co2 =',
     &  ' /* Atmospheric p13CO2 to be read in data.dic */')
#endif
#ifdef DIC_CARBON_COMPONENTS
       CALL WRITE_0D_RL( dic_pCO2_ref, INDEX_NONE,'dic_pCO2_ref =',
     &  ' /* Atmospheric pCO2 reference for Csat_init component */')
#endif
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( dic_int1.EQ.0 .AND. dic_pCO2.NE.278. _d -6 ) THEN
        WRITE(msgBuf,'(A)')
     &    'DIC_READPARMS: cannot change default dic_pCO2 if dic_int1=0'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R DIC_READPARMS: dic_pCO2 error'
      ENDIF

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_DIC */

      RETURN
      END