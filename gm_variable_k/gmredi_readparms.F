C $Header: /u/gcmpack/MITgcm/pkg/gmredi/gmredi_readparms.F,v 1.22 2011/07/13 22:59:53 jmc Exp $
C $Name:  $

#include "GMREDI_OPTIONS.h"

CBOP
C     !ROUTINE: GMREDI_READPARMS
C     !INTERFACE:
      SUBROUTINE GMREDI_READPARMS( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE GMREDI_READPARMS
C     | o Routine to initialize GM/Redi variables and constants.
C     *==========================================================*
C     | Initialize GM/Redi parameters, read in data.gmredi
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GMREDI.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
      INTEGER myThid

#ifdef ALLOW_GMREDI
C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit
CEOP

C--   GM/Redi parameter
C     GM_Small_Number  :: epsilon used in computing the slope
C     GM_slopeSqCutoff :: slope^2 cut-off value
      NAMELIST /GM_PARM01/
     &          GM_AdvForm, GM_AdvSeparate,
     &          GM_InMomAsStress,
     &          GM_isopycK,
     &          GM_background_K,
     &          GM_iso2dFile, GM_iso1dFile,
     &          GM_bol2dFile, GM_bol1dFile,
     &          GM_taper_scheme,
     &          GM_maxSlope,
     &          GM_Kmin_horiz,
     &          GM_Small_Number, GM_slopeSqCutoff,
     &          GM_Visbeck_alpha, GM_Visbeck_length,
     &          GM_Visbeck_depth,
     &          GM_useHMM2011, GM_useFM2006N2,
     &          GM_Visbeck_minDepth, GM_Visbeck_maxSlope,
     &          GM_Visbeck_minVal_K, GM_Visbeck_maxVal_K,
     &          GM_facTrL2dz, GM_facTrL2ML, GM_maxTransLay,
     &          GM_Scrit, GM_Sd,
     &          GM_MNC,
     &          GM_UseBVP, GM_BVP_cMin, GM_BVP_ModeNumber,
     &          GM_useSubMeso, subMeso_Ceff, subMeso_invTau,
     &          subMeso_LfMin, subMeso_Lmax

      _BEGIN_MASTER(myThid)

C--   Default values GM/Redi
      GM_AdvForm          = .FALSE.
      GM_AdvSeparate      = .FALSE.
      GM_InMomAsStress    = .FALSE.
      GM_isopycK          = -999.
      GM_background_K     = 0. _d 0
      GM_maxSlope         = 1. _d -2
      GM_Kmin_horiz       = 0. _d 0
      GM_Small_Number     = 1. _d -20
      GM_slopeSqCutoff    = 1. _d +48
      GM_taper_scheme     = ' '
      GM_facTrL2dz        = 1.
      GM_facTrL2ML        = 5.
      GM_maxTransLay      = 500.
      GM_Scrit            = 0.004 _d 0
      GM_Sd               = 0.001 _d 0
      GM_MNC              = useMNC
      GM_iso2dFile        = ' '
      GM_iso1dFile        = ' '
      GM_bol2dFile        = ' '
      GM_bol1dFile        = ' '

C--   Default values GM/Redi I/O control
c     GM_dumpFreq         = -1.
c     GM_taveFreq         = -1.

C--   Default values Visbeck
      GM_useHMM2011       =  .FALSE.
      GM_useFM2006N2      =  .FALSE.
      GM_Visbeck_alpha    =    0. _d 0
      GM_Visbeck_length   =  200. _d 3
      GM_Visbeck_depth    = 1000. _d 0
      GM_Visbeck_minDepth =    0. _d 0
      GM_Visbeck_maxSlope = UNSET_RL
      GM_Visbeck_minVal_K =    0. _d 0
      GM_Visbeck_maxVal_K = 2500. _d 0

C--   Default values BVP
      GM_UseBVP           = .FALSE.
      GM_BVP_ModeNumber   = 1
      GM_BVP_cMin         = 1. _d -1

C--   Default values Sub-Meso (corresponding tau = 5.8 days):
      GM_useSubMeso = .FALSE.
      subMeso_invTau =  2.0 _d -6
      subMeso_LfMin  =  1.0 _d +3
      subMeso_Ceff   =  7.0 _d -2
      subMeso_Lmax   = 110. _d +3

      WRITE(msgBuf,'(A)') ' GM_READPARMS: opening data.gmredi'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
      CALL OPEN_COPY_DATA_FILE(
     I                          'data.gmredi', 'GM_READPARMS',
     O                          iUnit,
     I                          myThid )

C     Read parameters from open data file
      READ(UNIT=iUnit,NML=GM_PARM01)
      WRITE(msgBuf,'(A)') ' GM_READPARMS: finished reading data.gmredi'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)
C     Close the open data file
      CLOSE(iUnit)

C     Default value for GM_isopycK is equal to GM_background_K :
      IF (GM_isopycK.EQ.-999.) GM_isopycK = GM_background_K

C     Default value for GM_Visbeck_maxSlope is equal to GM_maxSlope :
      IF ( GM_Visbeck_maxSlope .EQ. UNSET_RL )
     &     GM_Visbeck_maxSlope = GM_maxSlope

C     Some constants
      GM_rMaxSlope = 0.
      if (GM_maxSlope.NE.0.) GM_rMaxSlope = 1. _d 0 / GM_maxSlope

      IF (GM_AdvForm) THEN
        GM_skewflx = 0.
        GM_advect  = 1.
        GM_ExtraDiag = GM_Visbeck_alpha.NE.0. .OR. GM_isopycK.NE.0.
      ELSE
        GM_skewflx = 1.
        GM_advect  = 0.
        GM_ExtraDiag = GM_isopycK.NE.GM_background_K
      ENDIF
      IF ( GM_iso2dFile .NE. GM_bol2dFile .OR.
     &     GM_iso1dFile .NE. GM_bol1dFile ) THEN
        GM_ExtraDiag = .TRUE.
      ENDIF

      IF (GM_UseBVP) THEN
        GM_BVP_rModeNumber = 1. _d 0/DBLE(GM_BVP_ModeNumber)
        GM_BVP_cHat2Min = GM_BVP_cMin**2*rhoConst/gravity
      ENDIF

C     Make sure that we locally honor the global MNC on/off flag
      GM_MNC = GM_MNC .AND. useMNC
#ifndef ALLOW_MNC
C     Fix to avoid running without getting any output:
      GM_MNC = .FALSE.
#endif
      GM_MDSIO = (.NOT. GM_MNC) .OR. outputTypesInclusive

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_GMREDI */

      RETURN
      END