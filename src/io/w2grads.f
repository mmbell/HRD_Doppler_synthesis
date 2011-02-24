      PROGRAM CEDWRITE 
C
C     THIS ROUTINE IS A DRIVER FOR WRITING CEDRIC FORMAT FILES.
C     YOU WILL HAVE TO REPLACE THIS PROGRAM AND THE SUBROUTINE FETCHPLANE
C     IN THIS FILE WITH YOUR OWN. THESE ARE JUST EXAMPLES.
C
      PARAMETER (MAXFLD=25, NID=510)
      CHARACTER*4 PROJECT
      CHARACTER*6 SCINAME,NAMLND(7),RADSTN,TAPE
      CHARACTER*8 VOLNAM,FLDNAM(MAXFLD),SOURCE
      DIMENSION ISCLFLD(MAXFLD),XLND(7),YLND(7),ZLND(7)
      DIMENSION IDAT(512,512)
      INTEGER*2 IDATA(5,2000000)
      REAL U(2000000),V(2000000),W(2000000),DBZ(2000000),DIV(2000000)
      CHARACTER KEYWORD*4,FLTNAME*8,STMNAME*12,RADAR*4,EXPERIMENT*32
      CHARACTER CREATIME*32,EXTRA1*28,FILNAM*56
      INTEGER IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3
      REAL STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7
      CHARACTER NAME*56
      WRITE(6,*)'ENTER WIND3D FILE NAME EXCLUDING .W'
      READ(5,'(A56)')NAME
10    CALL READHEADER(99,NAME,KEYWORD,FLTNAME,
     + STMNAME,RADAR,EXPERIMENT,
     + CREATIME,EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     + SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     + POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)
      LUOUT=6
      WRITE(LUOUT,*)CREATIME
      WRITE(LUOUT,*)FLTNAME
      WRITE(LUOUT,*)STMNAME
      WRITE(LUOUT,*)'KEYWORD,RADAR,EXPERIMENT'
      WRITE(LUOUT,'(A4,2X,A4,2X,A32)')KEYWORD,RADAR,EXPERIMENT
      WRITE(LUOUT,*)'EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,IUNFLD,IATTEN'
      WRITE(LUOUT,*)EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,IUNFLD,IATTEN
      WRITE(LUOUT,*)'FLAG,EXTRA2'
      WRITE(LUOUT,*)FLAG,EXTRA2
      WRITE(LUOUT,*)'EXTRA3,STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ'
      WRITE(LUOUT,*)EXTRA3,STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ
      WRITE(LUOUT,*)'ROT,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL'
      WRITE(LUOUT,*)ROT,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL
      WRITE(LUOUT,*)'AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7'
      WRITE(LUOUT,*)AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7
      IF(KEYWORD.EQ.'RUV '.OR.KEYWORD.EQ.'ruv ')THEN
       IUV=1
      ELSE
       IUV=0
      ENDIF
      CALL LDPLN(U,V,W,DBZ,DIV,IMAX,JMAX,KMAX,IUV)
      CLOSE(99)
c      CALL WRITEVS(U,V,W,DBZ,DIV,IMAX,JMAX,KMAX)

C      CALL WRITEIDATA(IDATA,IMAX,JMAX,KMAX)
      SCINAME = 'GAMACH'
      BASANG  = 90.0+ROT
C      WRITE(6,*)'ENTER VOLUME NAME UP TO 8 CHARACTERS'
C      READ(5,'(A8)')VOLNAM
C      WRITE(6,*)'ENTER FLIGHT ID'
C      READ(5,'(A8)')SOURCE
C      WRITE(6,*)'ENTER YEAR (2-DIGIT)'
C      READ(5,*)IBEGYR
C      WRITE(6,*)'ENTER MONTH NUMBER'
C      READ(5,*)IBEGMNT
C      WRITE(6,*)'ENTER DAY OF MONTH'
C      READ(5,*)IBEGDAY
C      WRITE(6,*)'ENTER 4BYTE PROJECT NAME'
C      READ(5,'(A4)')PROJECT
      TAPE    = '000000'
      RADCON  = -99.9
      VNYQ    = -99.9
      IBEGHR=STIME/3600
      IREM=STIME-3600*IBEGHR
      IBEGMIN=IREM/60
      IBEGSEC=IREM-IBEGMIN*60
      IENDHR=STIME/3600
      IREM=STIME-3600*IENDHR
      IENDMIN=IREM/60
      IENDSEC=IREM-IENDMIN*60
      IENDYR  = IBEGYR
      IENDMNT = IBEGMNT
      IENDDAY = IBEGDAY

      LATDEG = INT(OLAT)
      REMLAT = OLAT-LATDEG
      LATMIN = INT(REMLAT*60.)
      REMLAT = REMLAT-LATMIN/60.
      LATSEC = NINT(REMLAT*3600.)
      LONDEG = INT(OLON)
      REMLON = OLON-LONDEG
      LONMIN = INT(REMLON*60.)
      REMLON = REMLON-LONMIN/60.
      LONSEC = NINT(REMLON*3600.)
      
c      XMIN  = -XZ+.5*SX
      XMIN  = -XZ
      XMAX  = XMIN+(IMAX-1)*SX
      NUMX  = IMAX
      ISPCX = NINT(SX*1000.)

c      YMIN  = -YZ+.5*SY
      YMIN  = -YZ
      YMAX  = YMIN+(JMAX-1)*SY
      NUMY  = JMAX
      ISPCY = NINT(SY*1000.)

      ZMIN  = ZZ*1000.
      ZMAX  = 1000*(ZZ+(KMAX-1)*SZ)
      NUMZ  = KMAX
      ISPCZ = NINT(SZ*1000.)


      CALL CHANGEFORM(U,V,W,DBZ,DIV,IMAX,JMAX,KMAX,IDATA,NAME,
     + OLAT,OLON,SX,SY,SZ,XMIN,YMIN)

      NFLD=5
      FLDNAM(1)  = 'U-COMPON'
      FLDNAM(2)  = 'V-COMPON'
      FLDNAM(3)  = 'W-COMPON'
      FLDNAM(4)  = 'DBZ     '
      FLDNAM(5)  = 'DIVERGEN'
      ISCLFLD(1) = 100
      ISCLFLD(2) = 100
      ISCLFLD(3) = 100
      ISCLFLD(4) = 100
      ISCLFLD(5) = 100
      RADSTN     = 'NOAA42'

      NUMLND = 1
      NUMRAD = 0
      NAMLND(1) = 'ANCHOR'
      XLND(1) = 0.0
      YLND(1) = 0.0
      ZLND(1) = 0.0


C      write(6,*)'nfld = ',nfld

C      print*, maxval(IDATA) 
C      CALL WRITCED(SCINAME,BASANG,VOLNAM,IBEGYR,
C     X     IBEGMNT,IBEGDAY,IBEGHR,IBEGMIN,IBEGSEC,IENDYR,
C     X     IENDMNT,IENDDAY,IENDHR,IENDMIN,IENDSEC,XMIN,
C     X     XMAX,NUMX,ISPCX,YMIN,YMAX,NUMY,ISPCY,ZMIN,
C     X     ZMAX,NUMZ,ISPCZ,NFLD,FLDNAM,ISCLFLD,NUMLND,
C     X     NUMRAD,NAMLND,XLND,YLND,ZLND,IDAT,RADSTN,SOURCE,
C     X     PROJECT,TAPE,RADCON,VNYQ,LATDEG,LATMIN,LATSEC,
C     X     LONDEG,LONMIN,LONSEC,IDATA)
      
      END


      SUBROUTINE CHANGEFORM(U,V,W,DBZ,DIV,IMAX,JMAX,KMAX,IDATA,NAME,
     + OLAT,OLON,SX,SY,SZ,XMIN,YMIN)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),W(IMAX,JMAX,KMAX)
      REAL DBZ(IMAX,JMAX,KMAX),DIV(IMAX,JMAX,KMAX)
      REAL FDATA(4,301,301,50)
      INTEGER*2 IDATA(5,IMAX,JMAX,KMAX)
      CHARACTER NAME*56, DATFILE*56, CTLFILE*56
      REAL latrad, fac_lat, fac_lon
      INTEGER LFN
      FLAG=-1.0E+10
      DO K=1,KMAX
       DO J=1,JMAX
        DO I=1,IMAX
         IF(U(I,J,K).GT.FLAG)THEN
          FDATA(1,I,J,K)=U(I,J,K)
         ELSE 
          FDATA(1,I,J,K)=-32768
         ENDIF
         IF(V(I,J,K).GT.FLAG)THEN
          FDATA(2,I,J,K)=V(I,J,K)
         ELSE
          FDATA(2,I,J,K)=-32768
         ENDIF
         IF(W(I,J,K).GT.FLAG)THEN
          FDATA(3,I,J,K)=W(I,J,K)
         ELSE
          FDATA(3,I,J,K)=-32768
         ENDIF
         IF(DBZ(I,J,K).GT.FLAG)THEN
          FDATA(4,I,J,K)=DBZ(I,J,K)
         ELSE
          FDATA(4,I,J,K)=-32768
         ENDIF
        ENDDO
       ENDDO
      ENDDO
c WRITE BINARY FORMAT FOR GRADS : modified by ms cccccccccccccccc
      LFN = index(NAME,' ')-1
      datfile = name(1:LFN)//'.dat'
      OPEN(40, file=datfile, FORM = 'unformatted',
     + ACCESS = 'direct', STATUS = 'unknown', recl=IMAX*JMAX*4)
      irec = 0 

      DO IVAR = 1, 4
      DO K=1,KMAX
      irec = irec + 1 
        WRITE(40,rec=irec)((FDATA(IVAR,I,J,K),I=1,IMAX)
     +     ,J=1,JMAX)
      ENDDO
      ENDDO 

C WRITE CTL FILE FOR GRADS 
      ctlfile = name(1:LFN)//'.ctl'      
      OPEN(50, file=ctlfile, STATUS = 'unknown')   
      write(50,'(a6,a56)') 'dset ^', DATFILE 
      write(50,'(a19)') 'title eldora 3d var' 
      write(50,'(a13)') 'undef -32768'   
      latrad = OLAT * 1.745329251994e-02
      fac_lat = 111.13209 - 0.56605 * cos(2.0 * latrad) 
     +     + 0.00012 * cos(4.0 * latrad) - 0.000002 * cos(6.0 * latrad)
      fac_lon = 111.41513 * cos(latrad)
     +     - 0.09455 * cos(3.0 * latrad) + 0.00012 * cos(5.0 * latrad)
      write(50,'(a5,I3,a19)'), 'xdef ', IMAX,' levels'
      DO I=1,IMAX
          write(50,'(f10.4)'),(XMIN + (I-1)*SX)/fac_lon + OLON
      enddo
      write(50,'(a5,I3,a19)'), 'ydef ', JMAX,' levels'
       DO J=1,JMAX
          write(50,'(f10.4)'),(YMIN + (J-1)*SY)/fac_lat + OLAT
      enddo
      write(50,'(a5,I3,a14)'), 'zdef ', KMAX,' linear 0. 0.5' 
      write(50,'(a33)'),'tdef 1 linear 23:50Z26SEP2008 1MN' 
      WRITE(50,'(a6)'),'vars 4' 
      WRITE(50,'(a4,I3,a12)'),'u   ', KMAX,' 99 variable'  
      WRITE(50,'(a4,I3,a12)'),'v   ', KMAX,' 99 variable'  
      WRITE(50,'(a4,I3,a12)'),'w   ', KMAX,' 99 variable'   
      WRITE(50,'(a4,I3,a12)'),'dbz ', KMAX,' 99 variable'  
      write(50, '(a7)') 'endvars' 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      RETURN
      END

      SUBROUTINE LDPLN(U,V,W,DBZ,DIV,IMAX,JMAX,KMAX,IUV)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),W(IMAX,JMAX,KMAX)
      REAL DBZ(IMAX,JMAX,KMAX),DIV(IMAX,JMAX,KMAX)
      REAL RWD(512),RWS(512),RWW(512),RDV(512)
      INTEGER*2 WD(512),WS(512),DB(512),WW(512),DV(512)         
      FLAG=-1.0E+10
      DO 1 K=1,KMAX                                                             
       WRITE(6,'("loading plane #",i3)') K                                    
       DO 2 J=1,JMAX                                                          
        IF(IUV.EQ.1)THEN
         READ(99)(RWD(I),RWS(I),RWW(I),DB(I),RDV(I),I=1,IMAX)
        ELSE
         READ(99) (WD(I),WS(I),WW(I),DB(I),DV(I),I=1,IMAX)                   
        ENDIF
        DO 3 I=1,IMAX                                                       
c         WRITE(6,*)'I,J,K,RWD,RWS,RWW,DB,RDV = ',I,J,K,
c     1              RWD(I),RWS(I),RWW(I),
c     1              DB(I),RDV(I)
         IF(IUV.EQ.0)THEN
          WDR=WD(I)*.1                                                     
          WSP=WS(I)*.1                                                     
          IF (WDR.LT.0.) THEN                                              
           U(I,J,K)=FLAG
           V(I,J,K)=FLAG
          ELSE                                                             
c           CALL COMP(XU,XV,WDR,WSP)                                       
           XU=-SIN(WDR*3.14159/180.)*WSP
           XV=-COS(WDR*3.14159/180.)*WSP
           U(I,J,K)=XU
           V(I,J,K)=XV
          ENDIF                                                            
          IF(WW(I).GT.-9000)THEN
           W(I,J,K)=WW(I)*.01
          ELSE
           W(I,J,K)=FLAG
          ENDIF
          IF(DV(I).LT.32767)THEN
           DIV(I,J,K)=DV(I)/100000.
          ELSE
           DIV(I,J,K)=FLAG
          ENDIF
         ELSE
          IF(RWW(I).LE.-100000.)RWD(I)=FLAG
          IF(RWS(I).LE.-100000.)RWS(I)=FLAG
          IF(RWW(I).LE.-100000.)RWW(I)=FLAG
          IF(RDV(I).LE.-100000.)RDV(I)=FLAG
          U(I,J,K)=RWD(I)
          V(I,J,K)=RWS(I)
          W(I,J,K)=RWW(I)
          DIV(I,J,K)=RDV(I)
         ENDIF
         IF(DB(I).GT.-9000)THEN
          DBZ(I,J,K)=DB(I)*.1
         ELSE
          DBZ(I,J,K)=FLAG
         ENDIF
C         WRITE(6,*)'I,J,K,U,V,W,DBZ,DIV = ',I,J,K,U(I,J,K),V(I,J,K),
C     1              W(I,J,K),DBZ(I,J,K),DIV(I,J,K)
3       CONTINUE                                                         
2      CONTINUE                                                               
1     CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       


      SUBROUTINE FETCHPLANE(IDAT,IMAX,JMAX,KMAX,IFIELD,K,IDATA)
C
C     SAMPLE ROUTINE FOR FETCHING Z PLANES OF CEDRIC DATA.
C     FILLS IN ARRAY WITH BOGUS DATA
C
      DIMENSION IDAT(IMAX,JMAX)
      INTEGER*2 IDATA(5,IMAX,JMAX,KMAX)
      DO 100 J=1,JMAX
         DO 50 I=1,IMAX
            IDAT(I,J)=IDATA(IFIELD,I,J,K)
C         write(6,*)'j,k,ifield,idat,idata = ',j,k,
C     1    ifield,(idat(i,j),idata(ifield,i,j,k),i=1,imax)
 50      CONTINUE
 100  CONTINUE

      RETURN

      END

         
      SUBROUTINE READHEADER(LU,FILNAM,KEYWORD,FLTNAME,
     + STMNAME,RADAR,EXPERIMENT,
     + CREATIME,EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     + SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     + POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)
      CHARACTER KEYWORD*4,FLTNAME*8,STMNAME*12,RADAR*4,EXPERIMENT*32
      CHARACTER CREATIME*32,EXTRA1*28,FILNAM*56,WFILE*58
      INTEGER IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,LFN
      REAL STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7
     
      LFN = index(FILNAM,' ')-1
      wfile = FILNAM(1:LFN)//'.w'
      OPEN(LU,file=wfile,IOSTAT=ISTAT,STATUS='OLD',FORM='UNFORMATTED')
      READ(LU,IOSTAT=ISTAT)KEYWORD,FLTNAME,STMNAME,RADAR,EXPERIMENT,
     + CREATIME,EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,IUNFLD,IATTEN,FLAG,
     + EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7
      RETURN
      END
      SUBROUTINE WRITEVS(U,V,W,DBZ,DIV,IMAX,JMAX,KMAX)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),W(IMAX,JMAX,KMAX)
      REAL DBZ(IMAX,JMAX,KMAX),DIV(IMAX,JMAX,KMAX)
      DO K=1,KMAX
       DO J=1,JMAX
C        WRITE(6,*)'J,K,U = ',(U(I,J,K),I=1,IMAX)
C        WRITE(6,*)'J,K,V = ',(V(I,J,K),I=1,IMAX)
C        WRITE(6,*)'J,K,W = ',(W(I,J,K),I=1,IMAX)
C        WRITE(6,*)'J,K,DBZ = ',(DBZ(I,J,K),I=1,IMAX)
C        WRITE(6,*)'J,K,DIV = ',(DIV(I,J,K),I=1,IMAX)
       ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE WRITEIDATA(IDATA,IMAX,JMAX,KMAX)
      INTEGER*2 IDATA(5,IMAX,JMAX,KMAX)
      DO K=1,KMAX
       DO J=1,JMAX
C        WRITE(6,*)'J,K,IDATU = ',(IDATA(1,I,J,K),I=1,IMAX)
C        WRITE(6,*)'J,K,IDATV = ',(IDATA(2,I,J,K),I=1,IMAX)
C        WRITE(6,*)'J,K,IDATW = ',(IDATA(3,I,J,K),I=1,IMAX)
C        WRITE(6,*)'J,K,IDATDBZ = ',(IDATA(4,I,J,K),I=1,IMAX)
C        WRITE(6,*)'J,K,IDATDIV = ',(IDATA(5,I,J,K),I=1,IMAX)
       ENDDO
      ENDDO
      RETURN
      END

      
