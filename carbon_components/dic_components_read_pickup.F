C $Header: /u/gcmpack/MITgcm/pkg/dic/dic_read_co2_pickup.F,v 1.2 2011/03/14 01:38:04 jmc Exp $
C $Name:  $

#include "DIC_OPTIONS.h"

      SUBROUTINE DIC_COMPONENTS_READ_PICKUP( pH_3disLoaded,
     I                                myIter, myThid )

      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DIC_VARS.h"

C     == Routine arguments ==
C     myThid  :: my Thread Id number
      INTEGER myIter
      INTEGER myThid

#ifdef ALLOW_DIC

C     fn          :: character buffer for creating filename
C     prec        :: precision of pickup files
C     filePrec    :: pickup-file precision (read from meta file)
C     nbFields    :: number of fields in pickup file (read from meta file)
C     fldName     :: Name of the field to read

C     !FUNCTIONS
      INTEGER  MDS_RECLEN
      EXTERNAL MDS_RECLEN
      
C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER i,j,k,bi,bj
      INTEGER filePrec, ioUnit, prec, nbFields, nj
      CHARACTER*(MAX_LEN_FNAM) fn, filNam
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(8) fldName
      LOGICAL useCurrentDir, fileExist, pH_3disLoaded
      INTEGER length_of_rec
      _RL vec(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy)
CEOP
#ifdef DIC_CARBON_COMPONENTS

      pH_3disLoaded =.FALSE.
       
        IF (pickupSuff.EQ.' ') THEN
          WRITE(fn,'(A,I10.10)') 'pickup_dic_components.', myIter
        ELSE
          WRITE(fn,'(A,A10)') 'pickup_dic_components.', pickupSuff
        ENDIF
        prec = precFloat64

C--   First check if pickup file exist
#ifdef ALLOW_MDSIO
        useCurrentDir = .FALSE.
        CALL MDS_CHECK4FILE(
     I                       fn, '.data', 'DIC_COMPONENTS_READ_PICKUP',
     O                       filNam, fileExist,
     I                       useCurrentDir, myThid )
#else
        STOP 'ABNORMAL END: S/R DIC_COMPONENTS_READ_PICKUP: Needs MDSIO pkg'
#endif
        IF ( fileExist ) THEN
C--   Read pickup file
    
         CALL READ_MFLDS_SET(
     I                      fn,
     O                      nbFields, filePrec,
     I                      Nr, myIter, myThid )
        
       _BEGIN_MASTER(myThid) 
       IF ( nbFields.GE.0 .AND. filePrec.NE.prec ) THEN
         WRITE(msgBuf,'(2A,I4)') 'DIC_COMPONENTS_READ_PICKUP:',
     &    'pickup-file binary precision do not match !'
         CALL PRINT_ERROR( msgBuf, myThid )
         WRITE(msgBuf,'(A,2(A,I4))') 'DIC_COMPONENTS_READ_PICKUP:',
     &    'file prec.=', filePrec, ' but expecting prec.=', prec
         CALL PRINT_ERROR( msgBuf, myThid )
         STOP 'ABNORMAL END: S/R DIC_COMPONENTS_READ_PICKUP'
       ENDIF
       _END_MASTER( myThid )
        
      IF ( nbFields.LE.0 ) THEN
C-      No meta-file or old meta-file without List of Fields
        ioUnit = errorMessageUnit
        IF ( pickupStrictlyMatch ) THEN
          WRITE(msgBuf,'(4A)') 'DIC_COMPONENTS_READ_PICKUP:',
     &      'no field-list found in meta-file',
     &      ' => cannot check for strick-matching'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(4A)') 'DIC_COMPONENTS_READ_PICKUP:',
     &      'try with " pickupStrictlyMatch=.FALSE.,"',
     &      ' in file: "data", NameList: "PARM03"'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
          STOP 'ABNORMAL END: S/R DIC_COMPONENTS_READ_PICKUP'
        ELSE
          WRITE(msgBuf,'(4A)') 'WARNING >> DIC_COMPONENTS_READ_PICKUP:',
     &      ' no field-list found'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         IF ( nbFields.EQ.-1 ) THEN
C-      No meta-file
          WRITE(msgBuf,'(4A)') 'WARNING >> ',
     &      ' try to read pickup as currently written'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         ELSE
C-      Old meta-file without List of Fields
          WRITE(msgBuf,'(4A)') 'WARNING >> ',
     &      ' try to read pickup as it used to be written'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(4A)') 'WARNING >> ',
     &      ' until checkpoint59l (2007 Dec 17)'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
         ENDIF
        ENDIF
       ENDIF
        
C---   New way to read ptracer pickup:
        nj = 0
C---    read pH 3-D fields for restart
          fldName = 'pH_INIT   '
          CALL READ_MFLDS_3D_RL( fldName,
     O                     vec,
     &                     nj, prec, Nr, myIter, myThid )
          CALL EXCH_3D_RL( vec, Nr, myThid )
     
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k=1,Nr
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
             pH_INIT(i,j,k,bi,bj) = vec(i,j,k,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        
          fldName = 'pH_GLOB   '
          CALL READ_MFLDS_3D_RL( fldName,
     O                     vec,
     &                     nj, prec, Nr, myIter, myThid )
          CALL EXCH_3D_RL( vec, Nr, myThid )

       DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k=1,Nr
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
             pH_GLOB(i,j,k,bi,bj) = vec(i,j,k,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        
#ifdef DIC_INSITU_COMPONENTS         
          fldName = 'pH_INSIT  '
          CALL READ_MFLDS_3D_RL( fldName,
     O                     vec,
     &                     nj, prec, Nr, myIter, myThid )
          CALL EXCH_3D_RL( vec, Nr, myThid )
    
           DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k=1,Nr
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
             pH_INSIT(i,j,k,bi,bj) = vec(i,j,k,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
#endif
        pH_3disLoaded =.TRUE.
      
      ELSE
        pH_3disLoaded = .FALSE.
        IF ( pickupStrictlyMatch ) THEN
          WRITE(msgBuf,'(4A)') 'DIC_COMPONENTS_READ_PICKUP:',
     &      'try with " pickupStrictlyMatch=.FALSE.,"',
     &      ' in file: "data", NameList: "PARM03"'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )
          STOP 'ABNORMAL END: S/R DIC_COMPONENTS_READ_PICKUP'
        ELSE
          WRITE(msgBuf,'(2A)') 'WARNING >> DIC_COMPONENTS_READ_PICKUP:',
     &      'will restart from approximated pH'
          CALL PRINT_MESSAGE( msgBuf, ioUnit, SQUEEZE_RIGHT, myThid )        
        ENDIF
      ENDIF
#endif /*  DIC_CARBON_COMPONENTS  */
#endif /*  ALLOW_DIC  */

      RETURN
      END
