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
      INTEGER*2 IDATA(5,20000000)
      REAL U(20000000),V(20000000),W(20000000)
      REAL DBZ(20000000),DIV(20000000)
      CHARACTER KEYWORD*4,FLTNAME*8,STMNAME*12,RADAR*4,EXPERIMENT*32
      CHARACTER CREATIME*32,EXTRA1*28,FILNAM*56
      INTEGER IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,LFN
      REAL STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7
      CHARACTER NAME*56,FORT10*56,CEDFILE*56
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
      CALL CHANGEFORM(U,V,W,DBZ,DIV,IMAX,JMAX,KMAX,IDATA)
c      CALL WRITEIDATA(IDATA,IMAX,JMAX,KMAX)
      SCINAME = 'GAMACH'
      BASANG  = 90.0+ROT
      WRITE(6,*)'ENTER VOLUME NAME UP TO 8 CHARACTERS'
      READ(5,'(A8)')VOLNAM
      WRITE(6,*)'ENTER FLIGHT ID'
      READ(5,'(A8)')SOURCE
      WRITE(6,*)'ENTER YEAR (2-DIGIT)'
      READ(5,*)IBEGYR
      WRITE(6,*)'ENTER MONTH NUMBER'
      READ(5,*)IBEGMNT
      WRITE(6,*)'ENTER DAY OF MONTH'
      READ(5,*)IBEGDAY
      WRITE(6,*)'ENTER 4BYTE PROJECT NAME'
      READ(5,'(A4)')PROJECT
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
      
      XMIN  = -XZ+.5*SX
      XMAX  = XMIN+(IMAX-1)*SX
      NUMX  = IMAX
      ISPCX = NINT(SX*1000.)

      YMIN  = -YZ+.5*SY
      YMAX  = YMIN+(JMAX-1)*SY
      NUMY  = JMAX
      ISPCY = NINT(SY*1000.)

      ZMIN  = ZZ*1000.
      ZMAX  = 1000*(ZZ+(KMAX-1)*SZ)
      NUMZ  = KMAX
      ISPCZ = NINT(SZ*1000.)

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


      write(6,*)'nfld = ',nfld
      CALL WRITCED(SCINAME,BASANG,VOLNAM,IBEGYR,
     X     IBEGMNT,IBEGDAY,IBEGHR,IBEGMIN,IBEGSEC,IENDYR,
     X     IENDMNT,IENDDAY,IENDHR,IENDMIN,IENDSEC,XMIN,
     X     XMAX,NUMX,ISPCX,YMIN,YMAX,NUMY,ISPCY,ZMIN,
     X     ZMAX,NUMZ,ISPCZ,NFLD,FLDNAM,ISCLFLD,NUMLND,
     X     NUMRAD,NAMLND,XLND,YLND,ZLND,IDAT,RADSTN,SOURCE,
     X     PROJECT,TAPE,RADCON,VNYQ,LATDEG,LATMIN,LATSEC,
     X     LONDEG,LONMIN,LONSEC,IDATA)

      LFN = index(NAME,' ')-1
      FORT10= 'fort.10'
      cedfile = name(1:LFN)//'.ced'
      CALL RENAME(FORT10,cedfile)
      
      END


      SUBROUTINE CHANGEFORM(U,V,W,DBZ,DIV,IMAX,JMAX,KMAX,IDATA)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),W(IMAX,JMAX,KMAX)
      REAL DBZ(IMAX,JMAX,KMAX),DIV(IMAX,JMAX,KMAX)
      INTEGER*2 IDATA(5,IMAX,JMAX,KMAX)
* added 2/19/10 (ers)
      character*158 line
      character*6 cfltlvllat  
      character*7 cfltlvllon
*
      FLAG=-1.0E+10
      DO K=1,KMAX
       DO J=1,JMAX
        DO I=1,IMAX
         IF(U(I,J,K).GT.FLAG)THEN
          IDATA(1,I,J,K)=NINT(U(I,J,K)*100.)
         ELSE 
          IDATA(1,I,J,K)=-32768
         ENDIF
         IF(V(I,J,K).GT.FLAG)THEN
          IDATA(2,I,J,K)=NINT(V(I,J,K)*100.)
         ELSE
          IDATA(2,I,J,K)=-32768
         ENDIF
         IF(W(I,J,K).GT.FLAG)THEN
          IDATA(3,I,J,K)=NINT(W(I,J,K)*100.)
         ELSE
          IDATA(3,I,J,K)=-32768
         ENDIF
         IF(DBZ(I,J,K).GT.FLAG)THEN
          IDATA(4,I,J,K)=NINT(DBZ(I,J,K)*100.)
         ELSE
          IDATA(4,I,J,K)=-32768
         ENDIF
         IF(DIV(I,J,K).GT.FLAG)THEN
          IDATA(5,I,J,K)=NINT(DIV(I,J,K)*100000.)
         ELSE
          IDATA(5,I,J,K)=-32768
         ENDIF
        ENDDO
       ENDDO
      ENDDO
c 
c	go to 91919
c
c     Set w to missing .15 deg west and .10 deg east of the flight
c     longitude (beth sanabia 2/17/10)
c************
c  Open flight level data file to read from:
	open (22,file='rf13_fltlvl_0135_0205',status='old')
c************
c  Pertinent definitions (generic):  
c  XMIN  = -XZ+.5*SX
c  XMAX  = XMIN+(IMAX-1)*SX
c  YMIN  = -YZ+.5*SY
c  YMAX  = YMIN+(JMAX-1)*SY
c  OLAT = the lat of the origin
c  latrad = OLAT * 1.745329251994e-02
c  fac_lat = 111.13209 - 0.56605 * cos(2.0 * latrad)
c     +     + 0.00012 * cos(4.0 * latrad) - 0.000002 * cos(6.0 * latrad)
c  fac_lon = 111.41513 * cos(latrad)
c     +     - 0.09455 * cos(3.0 * latrad) + 0.00012 * cos(5.0 * latrad)
c  eldoralon=(YMIN + (J-1)*SY)/fac_lat + OLAT	
c  longitude that corresponds with "I"=(XMIN + (I-1)*SX)/fac_lon + OLON
c  latitude that corresponds with "J"=(YMIN + (J-1)*SY)/fac_lat + OLAT
c***************
c  Pertinent definitions (specific to analysis of rf13 from 0135-0205):
	rmyxmin=0.25
	rmyxmax=49.875
       	rmyymin=0.25
	rmyymax=49.875
	rmysx=0.5
	rmysy=0.5
	rmyolat=27.100000
	rmyolon=126.10000
        rmylatrad = myolat * 1.745329251994e-02
        rmyfac_lat = 111.13209 - 0.56605*cos(2.0*rmylatrad)
     +     + 0.00012*cos(4.0*rmylatrad)-0.000002*cos(6.0*rmylatrad)
	rmyfac_lon = 111.41513*cos(rmylatrad)
     +     - 0.09455*cos(3.0*rmylatrad)+0.00012*cos(5.0*rmylatrad)
c Loop
      DO K=1,KMAX
	write(6,*)'K=',K
       DO J=1,JMAX
c solve this to see latitude corresponding to this j (lat in eldora)
	 rmyeldlat=(rmyymin + (J-1)*rmysy)/rmyfac_lat + rmyolat

c read through the flt lvl lat's within the correct time frame to locate
c the flt lvl lat closest to the eldora longitude.
c
        rmymaxneg=-99999.9
	rmyminpos=99999.9
	rewind (22)
880	read (22,'(a158)',end=642) line
	  cfltlvllat = line(23:28) 
  	  read (cfltlvllat,'(F6.3)') fltlvllat
	  cfltlvllon = line(32:38) 
  	  read (cfltlvllon,'(F7.3)') fltlvllon
c  find min delta here
	  delta=rmyeldlat-fltlvllat
	  if (delta.gt.rmymaxneg) then 
	     rmymaxneg=delta 
	     rmaxneglat=fltlvllat
	     rmaxneglon=fltlvllon
	  endif 
	  if (delta.lt.rmyminpos) then
	     rmyminpos=delta
	     rminposlat=fltlvllat
	     rminposlon=fltlvllon
	  endif
	goto 880
642 	continue
	if (abs(rmymaxneg).ge.abs(rmyminpos)) then
		rmyfltlvllat=rminposlat
		rmyfltlvllon=rminposlon
	endif
	if (abs(rmyminpos).gt.abs(rmymaxneg)) then 
		rmyfltlvllat=rmaxneglat
		rmyfltlvllon=rmaxneglon
	endif
c note the flt lvl lon associated with the identified flt lvl lat	
         
c solve this for I (& then do an "NINT" on the result) to get the i
c closest to the flight level longitude
c	longitude that corresponds with "I"=(XMIN + (I-1)*SX)/fac_lon + OLON
c
	rifindi=(rmyfac_lon*(rmyfltlvllon-rmyolon)-rmyxmin)/rmysx + 1
	myi = NINT(rifindi)

	DO I=1,IMAX
	 if ((I.le.(myi-30)).or.(I.eq.myi).or.(I.ge.myi)) then
	     IDATA(3,I,J,K)=-32768
	 endif
	    
        ENDDO
       ENDDO
      ENDDO
      close(22) 
c91919 continue
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
c         WRITE(6,*)'I,J,K,U,V,W,DBZ,DIV = ',I,J,K,U(I,J,K),V(I,J,K),
c     1              W(I,J,K),DBZ(I,J,K),DIV(I,J,K)
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
c         write(6,*)'j,k,ifield,idat,idata = ',j,k,
c     1    ifield,(idat(i,j),idata(ifield,i,j,k),i=1,imax)
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
      CHARACTER CREATIME*32,EXTRA1*28,FILNAM*56,wfile*56
      INTEGER IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,LFN,ISTAT
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
        WRITE(6,*)'J,K,U = ',(U(I,J,K),I=1,IMAX)
        WRITE(6,*)'J,K,V = ',(V(I,J,K),I=1,IMAX)
        WRITE(6,*)'J,K,W = ',(W(I,J,K),I=1,IMAX)
        WRITE(6,*)'J,K,DBZ = ',(DBZ(I,J,K),I=1,IMAX)
        WRITE(6,*)'J,K,DIV = ',(DIV(I,J,K),I=1,IMAX)
       ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE WRITEIDATA(IDATA,IMAX,JMAX,KMAX)
      INTEGER*2 IDATA(5,IMAX,JMAX,KMAX)
      DO K=1,KMAX
       DO J=1,JMAX
        WRITE(6,*)'J,K,IDATU = ',(IDATA(1,I,J,K),I=1,IMAX)
        WRITE(6,*)'J,K,IDATV = ',(IDATA(2,I,J,K),I=1,IMAX)
        WRITE(6,*)'J,K,IDATW = ',(IDATA(3,I,J,K),I=1,IMAX)
        WRITE(6,*)'J,K,IDATDBZ = ',(IDATA(4,I,J,K),I=1,IMAX)
        WRITE(6,*)'J,K,IDATDIV = ',(IDATA(5,I,J,K),I=1,IMAX)
       ENDDO
      ENDDO
      RETURN
      END

      
