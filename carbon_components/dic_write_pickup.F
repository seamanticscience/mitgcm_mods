#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: DIC_WRITE_PICKUP

C !INTERFACE: ==========================================================
      SUBROUTINE DIC_WRITE_PICKUP( permPickup,
     I                             suff, myTime, myIter, myThid )

C !DESCRIPTION:
C     Writes DIC arrays (needed for a restart) to a pickup file

C !USES: ===============================================================
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIC_VARS.h"
#include "DIC_ATMOS.h"

C !INPUT PARAMETERS: ===================================================
C     permPickup :: write a permanent pickup
C     suff       :: suffix for pickup file (eg. ckptA or 0000000010)
C     myTime     :: Current time in simulation
C     myIter     :: Current iteration number in simulation
C     myThid     :: My Thread Id number
      LOGICAL permPickup
      CHARACTER*(*) suff
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_DIC

C     !LOCAL VARIABLES:
C     == Local variables ==
      CHARACTER*(MAX_LEN_FNAM) fn
      INTEGER prec
      INTEGER ioUnit
      _RL tmpFld(2)
      _RS dummyRS(1)
#ifdef DIC_BIOTIC
      LOGICAL glf
      _RL     timList(1)
      INTEGER j, nj
      INTEGER listDim, nWrFlds
      PARAMETER( listDim = 2 )
      CHARACTER*(8) wrFldList(listDim)
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#endif

#ifdef DIC_CARBON_COMPONENTS
      INTEGER pHlistDim, pHnWrFlds
#ifdef DIC_INSITU_COMPONENTS     
      PARAMETER( pHlistDim = 3 )
#else
      PARAMETER( pHlistDim = 2 )
#endif
      CHARACTER*(8) pHwrFldList(pHlistDim)
#endif

c     IF ( DIC_pickup_write_mdsio ) THEN
        prec = precFloat64

        IF ( dic_int1.GE.3 ) THEN
          WRITE(fn,'(A,A)') 'pickup_dic_co2atm.',suff
          ioUnit = 0
#ifdef ALLOW_OPENAD
          tmpFld(1) = total_atmos_carbon%v
          tmpFld(2) = atpco2%v
#else /* ALLOW_OPENAD */
          tmpFld(1) = total_atmos_carbon
          tmpFld(2) = atpco2
#endif /* ALLOW_OPENAD */
#ifdef ALLOW_MDSIO
          CALL MDS_WRITEVEC_LOC(
     I                         fn, prec, ioUnit,
     I                         'RL', 2, tmpFld, dummyRS,
     I                         0, 0, 1, myIter, myThid )
#endif
        ENDIF

#ifdef DIC_BIOTIC
        WRITE(fn,'(A,A)') 'pickup_dic.',suff
        j = 0

C     Firstly, write 3-D fields as consecutive records,

C-    switch to 2-D fields:
        nj = -j*Nr

C       record number < 0 : a hack not to write meta files now:
        j = j + 1
        nj = nj-1
        CALL WRITE_REC_3D_RL( fn, prec, 1, pH, nj, myIter, myThid )
        IF (j.LE.listDim) wrFldList(j) = 'DIC_pH2d'

C--------------------------
        nWrFlds = j
        IF ( nWrFlds.GT.listDim ) THEN
          WRITE(msgBuf,'(2A,I5,A)') 'DIC_WRITE_PICKUP: ',
     &     'trying to write ',nWrFlds,' fields'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A,I5,A)') 'DIC_WRITE_PICKUP: ',
     &     'field-list dimension (listDim=',listDim,') too small'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DIC_WRITE_PICKUP (list-size Pb)'
        ENDIF

#ifdef ALLOW_MDSIO
C     uses this specific S/R to write (with more informations) only meta
C     files
        j  = 1
        nj = ABS(nj)
        IF ( nWrFlds*Nr .EQ. nj ) THEN
          j  = Nr
          nj = nWrFlds
        ENDIF
        glf  = globalFiles
        timList(1) = myTime
        CALL MDS_WR_METAFILES( fn, prec, glf, .FALSE.,
     &                         0, 0, j, ' ',
     &                         nWrFlds, wrFldList,
     &                         1, timList, oneRL,
     &                         nj, myIter, myThid )
#endif /* ALLOW_MDSIO */
C--------------------------

#endif /* DIC_BIOTIC  */

#ifdef DIC_CARBON_COMPONENTS
          WRITE(fn,'(A,A)') 'pickup_dic_components.',suff

          j = 0

C     Firstly, write 3-D fields as consecutive records,

          j = j + 1
          CALL WRITE_REC_3D_RL( fn, prec, Nr,
     &         pH_ini, -j, myIter, myThid )
          IF (j.LE.pHlistDim)
     &      pHwrFldList(j) = 'pH_ini    '
    
          j = j + 1
          CALL WRITE_REC_3D_RL( fn, prec, Nr,
     &         pH_glob, -j, myIter, myThid )
          IF (j.LE.pHlistDim)
     &      pHwrFldList(j) = 'pH_glob   '

#ifdef DIC_INSITU_COMPONENTS        
          j = j + 1
          CALL WRITE_REC_3D_RL( fn, prec, Nr,
     &         pH_loc, -j, myIter, myThid )
          IF (j.LE.pHlistDim)
     &      pHwrFldList(j) = 'pH_loc  '
#endif      
C--------------------------
        pHnWrFlds = j
        
        IF ( pHnWrFlds.GT.pHlistDim ) THEN
          WRITE(msgBuf,'(2A,I5,A)') 'DIC_WRITE_PICKUP: ',
     &     'trying to write ',pHnWrFlds,' pH fields'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A,I5,A)') 'DIC_WRITE_PICKUP: ',
     &     'field-list dimension (pHlistDim=',pHlistDim,') too small'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R DIC_WRITE_PICKUP (list-size pH)'
        ENDIF

C     uses this specific S/R to write (with more informations) only meta
C     files
        
        glf  = globalFiles
        timList(1) = myTime
        CALL MDS_WR_METAFILES( fn, prec, glf, .FALSE.,
     &                         0, 0, Nr, ' ',
     &                         pHnWrFlds, pHwrFldList,
     &                         1, timList, oneRL,
     &                         j, myIter, myThid )
C--------------------------
#endif /* DIC_CARBON_COMPONENTS */

c     ENDIF /* DIC_pickup_write_mdsio */

#endif /* ALLOW_DIC  */

      RETURN
      END
