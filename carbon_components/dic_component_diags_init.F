#include "GCHEM_OPTIONS.h"
#include "DIC_OPTIONS.h"
C !INTERFACE: ==========================================================
      SUBROUTINE DIC_COMPONENT_DIAGS_INIT(myThid )

C !DESCRIPTION:
C define diagnostics for dic package
C experiment

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS

C     !LOCAL VARIABLES:
      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C     Define diagnostics Names :
#if ((defined DIC_CARBON_COMPONENTS) && (defined DIC_COMPONENT_DIAGS))
      diagName  = 'C_satINI'
      diagTitle = 'Saturated DIC at INItial atm pCO2'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM P    MR    '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)

      diagName  = 'C_resINI'
      diagTitle = 'Residual DIC at INItial atm pCO2'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)

      diagName  = 'C_satLOC'
      diagTitle = 'Saturated DIC at LOCal pCO2 (last surface contact)'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM P    MR    '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)

      diagName  = 'C_resLOC'
      diagTitle = 'Residual DIC at LOCal pCO2 (last surface contact)'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)

      diagName  = 'C_satATM'
      diagTitle = 'Saturated DIC at current ATM pCO2'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM P    MR    '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)
     
      diagName  = 'C_resATM'
      diagTitle = 'Residual DIC at current ATM pCO2'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, diagMate, myThid)
     
      diagName  = 'C_soft  '
      diagTitle = 'Soft-tissue carbon concentration'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'C_carb  '
      diagTitle = 'Carbonate concentration'
      diagUnits = 'mol C/m3      '
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'P_reg   '
      diagTitle = 'Regenerated Phosphate concentration'
      diagUnits = 'mol P/m3      '
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'A_reg   '
      diagTitle = 'Regenerated Alkalinity concentration'
      diagUnits = 'mol eq/m3      '
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'Si_pre  '
      diagTitle = 'Preformed Silica concentration'
      diagUnits = 'mol Si/m3'
      diagCode  = 'SM      MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'Silica  '
      diagTitle = 'Non-interactive 3d silica field'
      diagUnits = 'mol Si/m3 '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
       diagName  ='pH_glob '
      diagTitle = 'Global contemporary pH'
      diagUnits = ' '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'pH_ini   '
      diagTitle = 'Global pH at INItial pCO2'
      diagUnits = ' '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
     
      diagName  = 'pH_loc  '
      diagTitle = 'Global pH at LOCal pCO2 (last surface contact)'
      diagUnits = ' '
      diagCode  = 'SM P    MR    '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )     
#endif
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
