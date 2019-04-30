C $Header: /u/gcmpack/MITgcm/pkg/dic/DIC_LOAD.h,v 1.8 2011/04/17 21:08:40 jmc Exp $
C $Name:  $

C--   COMMON /DIC_LOAD/
C     DIC_ldRec     :: time-record currently loaded (in temp arrays *[1])

      COMMON /DIC_LOAD_I/ DIC_ldRec
      COMMON /DIC_LOAD_RS/
     &    dicwind0, dicwind1, atmosp0,
     &    atmosp1, silica0, silica1
#ifdef READ_PAR
     &    ,par0, par1
#endif
#if !(defined ALLOW_SALT_PLUME) || !(defined SALT_PLUME_DIC)   
     &    , ice0, ice1
#endif
#ifdef ALLOW_FE
     &    ,feinput0,feinput1
#ifdef HYDROTHERMAL_FE
     &    ,htvinput0,htvinput1
#endif
#endif

      INTEGER DIC_ldRec(nSx,nSy)
      _RS dicwind0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS dicwind1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#if !(defined ALLOW_SALT_PLUME) || !(defined SALT_PLUME_DIC)   
      _RS ice0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS ice1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
      _RS atmosp0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS atmosp1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS silica0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS silica1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#ifdef READ_PAR
      _RS par0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS par1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
#ifdef ALLOW_FE
      _RS feinput0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS feinput1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      
#ifdef HYDROTHERMAL_FE
C JML include iron source due to hydrothermal input into bottom layer            
      _RS htvinput0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS htvinput1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif /* HYDROTHERMAL_FE */
#endif /* ALLOW_FE */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
