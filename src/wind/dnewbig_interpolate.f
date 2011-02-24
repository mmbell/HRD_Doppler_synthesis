      PROGRAM LEASTWIND
      CHARACTER KEYWORD*4,FLTNAME*8,STMNAME*12,RADAR*4,EXPERIMENT*32
      CHARACTER CREATIME*32,EXTRA1*28,FILNAM(40)*80,OFILNAM*80
      CHARACTER TESTFILE*80,KEYWORDOUT*4
      INTEGER IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3
      REAL STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7
      REAL TSTARTARRAY(40),TENDARRAY(40)
      REAL SXB,SYB,SZB,XZB,YZB,ZZB,ROTB
      CHARACTER KEYWORDB*4,FLNAM*80,TIMEMOVE*1
      CHARACTER COMYES*1,WINDFILE*80,CHARV*2,CHARDBZ*2
      CHARACTER SUMFILENAME*80,OSUMFILENAME(80)*80,IYES*1
      INTEGER IMAXB,JMAXB,KMAXB,linuxyes,nargs
      DOUBLE PRECISION DOTMIN
      REAL XSHIFT(80),YSHIFT(80),DBZFACT(80)
      REAL RHMAX,RZMAX,EHMAX,EZMAX
      REAL XACC,YACC,ZACC
c
c The next three commented lines have parameters set for a 48 Meg 
c maximum program.  If you wish to increase this you must have the 
c memory available.  If you do decide to increase the available points
c from 24,000, my first suggestion is that you
c multiply each of the parameters (N1,N2,NARRAY,NWORK,NPOINTS)
C by the same factor.  N1 will tell you how many grid points you are
c allowed in your analysis domain.  You can fiddle
c a little if you want to see if you can decrease NARRAY or NWORK a
c a little relative to the others, but you won't be able to do it much.
c I wish you success....John Gamache
c
c Added support for dynamic memory allocation in F90 -- Michael Bell
c
      INTEGER N1, N2, N8, NARRAY, NWORK, NPOINTS
c      PARAMETER(N1=24000,N2=48000,N8=192000)
c      PARAMETER(NARRAY=2300000,NWORK=550000)
c      PARAMETER (NPOINTS=72100)
      PARAMETER(IW0=0)

c      PARAMETER(LINUXYES=1)
c      REAL SUM1(N1),SUM2(N1),SUM3(N1)
c      REAL SUM4(N1),SUM5(N1),SUM6(N1)
c      REAL SUM7(N1),SUM8(N1),SUM9(N1)
c      REAL SUM10(N1),SUM11(N1),SUM12(N1)
c      REAL SUM13(N1),SUMWTSAVE(N1)
c      REAL DISTXYZ(N8)
c      REAL SUMDBZ(N1),SUMDB(N1)
c      REAL UGUESS(N1),VGUESS(N1),WGUESS(N1)
c      REAL XWEIGHT(N1),YWEIGHT(N1),ZWEIGHT(N1)
c      REAL U(N1),V(N1),W(N1),BAND1(NPOINTS)
c      REAL DIV(N1)
      CHARACTER RAMSFILES(80)*80,DIVERFILE*80
      DOUBLE PRECISION SXBD,SYBD,SZBD
c      INTEGER IA(NPOINTS)
c      INTEGER INDEX1(NPOINTS)
c      INTEGER*2 ICAST(N1)
      real, dimension (:), allocatable :: SUM1, SUM2, SUM3, SUM4,
     +   SUM5, SUM6, SUM7, SUM8, SUM9, SUM10, SUM11, SUM12, SUM13,
     +   SUMWTSAVE, SUMDBZ, SUMDB, UGUESS, VGUESS, WGUESS,
     +   XWEIGHT, YWEIGHT, ZWEIGHT, U, V, W, BAND1, DIV, DISTXYZ
      integer, dimension (:),  allocatable :: IA, INDEX1
      integer*2, dimension (:), allocatable :: ICAST
      integer :: err

      RFLAG=-1.0E+10
      OPEN(1,FILE='/dev/tty')
      CTR=3.14159/180.
      LUOUT=6
      LU=99
      LUSCAN=96
      KEYWORDB='WIND'
      RA=-999.
      CO1=-999.
      CO2=-999.
      AZMCOR=-999.
      ELCOR=-999.
      THRESH=-999.
      POWERT=-999.
      BIEL=-999.
      AZBIEL=-999.
      ETIME1=-999.
      STIME2=-999.
      EXTRA6=-999.
      EXTRA7=-999.
      IUNFLD=1
      IATTEN=0
      FLAG=0
      EXTRA2=-999
      EXTRA3=-999
C      WRITE(1,*)'WARNING--With this version you should not change'
C      WRITE(1,*)'your weighting parameters'
C      WRITE(1,*)'DBZ ONLY Y(1)/N(0)?'
C      READ(1,*)IDBZONLY
      IDBZONLY=0
C      WRITE(1,*)'Would you like to use a command file?'
C      READ(1,'(A)')COMYES
      COMYES='Y'
      IF(COMYES.EQ.'Y'.OR.COMYES.EQ.'y')THEN
       LUIN=15
       nargs = iargc()
       IF ( nargs.GT.0 ) THEN
          Call getarg(1,FLNAM)
          OPEN(15,FILE=FLNAM,STATUS='OLD')
       ELSE
          write(6,*)'Enter command file name'
          READ(*,'(A)')FLNAM
      ENDIF
        OPEN(LUIN,FILE=FLNAM,STATUS='OLD')
      ELSE
       LUIN=1
      ENDIF
c      WRITE(1,*)'Enter flight name 8 characters or less'
      READ(LUIN,'(A8)')FLTNAME
      WRITE(6,1501)FLTNAME
1501  FORMAT('FLIGHT NAME IS ',A8)
C      WRITE(1,*)'Enter Experiment Name 32 characters or less'
      READ(LUIN,'(A32)')EXPERIMENT
      WRITE(6,1502)EXPERIMENT
1502  FORMAT('EXPERIMENT NAME IS ',A32)
C      WRITE(1,*)'How many DOP3D files will be used in this program?'
      READ(LUIN,*)NFILE
      WRITE(6,1503)NFILE
1503  FORMAT('NUMBER OF INPUT UF FILES IS ',I3)
      IF(NFILE.GT.40)THEN
       WRITE(1,*)'Too many files'
       STOP
      ENDIF
      DO NF=1,NFILE
C       write(1,*)'Enter Doppler scanfile name for file # ',NF
       READ(LUIN,'(A80)')FILNAM(NF)
       WRITE(6,1504)NF,FILNAM(NF)
1504   FORMAT('FILE NO. ',I3,' IS ',A80)
C       WRITE(1,*)'How much would you like to multiply dbz by?'
c       READ(LUIN,*)DBZFACT(NF),XSHIFT(NF),YSHIFT(NF)
       DBZFACT(NF)=1.
       XSHIFT(NF)=0.
       YSHIFT(NF)=0.
       WRITE(6,1505)NF,DBZFACT(NF)
1505   FORMAT('DBZ FOR FILE NO. ',I3,' WILL BE MULTIPLIED BY ',
     +  F6.2)
       WRITE(6,1506)NF,XSHIFT(NF),YSHIFT(NF)
1506   FORMAT('OBSERVATIONS FOR FILE NO. ',I3,' WILL BE SHIFTED'/
     + F8.2,' AND ',F8.2,' KM IN X AND Y DIRECTION')
c       WRITE(1,*)'Enter beginning and ending times to use data in',
c     +  ' HHMMSS HHMMSS'
       READ(LUIN,'(3I2,X,3I2)')IB,JB,KB,IE,JE,KE
       WRITE(6,1507)NF,IB,JB,KB
1507   FORMAT('BEGIN TIME FOR FILE NO. ',I3,' IS ',3I2)
       WRITE(6,1508)NF,IE,JE,KE
1508   FORMAT('END TIME FOR FILE NO. ',I3,' IS ',3I2)
       TSTART=IB*3600.+JB*60+KB
       STIME=TSTART
       ETIME=IE*3600.+JE*60.+KE
       IF(ETIME.LT.STIME)THEN
        TEND=ETIME+86400.
       ELSE
        TEND=ETIME
       ENDIF
       WRITE(6,*)'STIME,ETIME,TSTART,TEND = ',
     +            STIME,ETIME,TSTART,TEND              
       TSTARTARRAY(NF)=TSTART
       TENDARRAY(NF)=TEND
C       WRITE(1,*)'Enter the SUMFILE for this flight leg'
       READ(LUIN,'(A80)')OSUMFILENAME(NF)
       WRITE(6,1509)NF,OSUMFILENAME(NF)
1509   FORMAT('OUTPUT SUMFILE FILE NAME FOR UF FILE NO. ',I3,' IS '/
     + A80)
      ENDDO
C      WRITE(1,*)'Mean storm velocity for this flight leg'
      READ(LUIN,*)SMOTIONU,SMOTIONV
      WRITE(6,1510)SMOTIONU,SMOTIONV
1510  FORMAT('STORM VELOCITY (U,V) IN M/S TO BE USED AS COMPOSITE '/
     + 'ADVECTION VELOCITY IS ',2F8.2)
C      WRITE(6,*)'Enter center time HHMMSS'
      READ(LUIN,'(3I2)')IH,IM,IS
      WRITE(6,1511)IH,IM,IS
1511  FORMAT('ADVECTION CENTER TIME FOR INTERPOLATION COMPOSITE IS '/
     + 3I2)
      CENTIME=IH*3600+IM*60+IS
      IF(CENTIME.LT.STIME) THEN
         PRINT *, 'Center time is prior to start time.'
         PRINT *, 'Do you want to move it forward one day? (Y/N)'
         READ (*,*) TIMEMOVE
         IF (TIMEMOVE.eq.'Y') THEN
            PRINT *, 'Adding 86400 seconds...'
            CENTIME=CENTIME+86400.
         ELSE
            PRINT *, 'Leaving it alone'
         ENDIF
      ENDIF
C      WRITE(1,*)'Enter latitude and longitude of anchor point'
      READ(LUIN,*)OLAT,OLON
      WRITE(6,1512)OLAT,OLON
1512  FORMAT('LATITUDE AND LONGITUDE FOR ORIGIN WILL BE ',
     + F8.3,' AND ',F8.3)
C      WRITE(6,*)'OLAT = ',OLAT,'   OLON = ',OLON
      I1REMOVE=0
C      WRITE(1,*)'Enter SX,SY,SZ OF OUTPUT FILE'
      READ(LUIN,*)SXBD,SYBD,SZBD
      WRITE(6,1513)SXBD,SYBD,SZBD
1513  FORMAT('X-, Y-, AND Z-RESOLUTIONS OF INTERPOLATION GRID WILL BE '/
     + 3F8.2)
      SXB=SXBD
      SYB=SYBD
      SZB=SZBD
C      WRITE(1,*)'Enter XZ,YZ,ZZ,ROT OF OUTPUT FILE'
      READ(LUIN,*)XZB,YZB,ZZB,ROTB
      WRITE(6,1514)XZB,YZB
1514  FORMAT('POSITION OF ORIGIN RELATIVE TO LOWER-LEFT CORNER OF '/
     + 'COMPOSITE WILL BE ',F8.2,', ',F8.2,' KM')
      WRITE(6,1515)ZZB
1515  FORMAT('HEIGHT OF FIRST VERTICAL LEVEL WILL BE ',F8.2,' KM')
      WRITE(6,1516)ROTB
1516  FORMAT('Y-AXIS OF INTERPOLATION GRID WILL BE ROTATED '/
     + F8.2,' DEGREES CLOCKWISE FROM DUE NORTH')
c      WRITE(1,*)'Enter IMAX,JMAX,KMAX of output file'
      READ(LUIN,*)IMAXB,JMAXB,KMAXB
      WRITE(6,1530)IMAXB,JMAXB,KMAXB

      N1=IMAXB*JMAXB*KMAXB
      N2=2*N1
      N8=8*N1
      allocate ( SUM1(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM2(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM3(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM4(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM5(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM6(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM7(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM8(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM9(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM10(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM11(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM12(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUM13(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUMWTSAVE(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUMDBZ(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( SUMDB(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( UGUESS(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( VGUESS(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( WGUESS(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( XWEIGHT(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( YWEIGHT(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( ZWEIGHT(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( U(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( V(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( W(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( DIV(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( DISTXYZ(N8), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( ICAST(N1), stat=err)
      if (err /=0) print *, "allocate failed."

      NPOINTS = N1*3+1
      allocate ( BAND1(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( IA(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( INDEX1(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."

      DO N=1,N1
       U(N)=RFLAG
       V(N)=RFLAG
       DIV(N)=0.     
       W(N)=RFLAG
       SUMDBZ(N)=RFLAG
      ENDDO

1530  FORMAT('INTEGER DOMAIN OF INTERPOLATION ARRAY IN X,Y,Z '/
     + 'WILL BE ',3I4)
c      write(1,*)'Enter horizontal, vertical maximum radii of ',
c     1 ' influence in km'       
      READ(LUIN,*)RHMAX,RZMAX
      WRITE(6,1517)RHMAX,RZMAX
1517  FORMAT('HORIZONTAL AND VERTICAL MAXIMUM RADII OF INFLUENCE '/
     + 'WILL BE ',2F8.2)
C      write(1,*)'Enter horizontal, vertical e-folding radius'
      READ(LUIN,*)EHMAX,EZMAX
      WRITE(6,1518)EHMAX,EZMAX
1518  FORMAT('HORIZONTAL AND VERTICAL E-FOLDING RADII FOR '/
     + 'INTERPOLATION WILL BE ',2F8.2)
C      write(1,*)'Enter acceptable x,y,z displacement'
C      READ(LUIN,*)XACC,YACC,ZACC
      XACC=50.
      YACC=50.
      ZACC=50.
      WRITE(6,1519)XACC,YACC,ZACC
1519  FORMAT('MAXIMUM X-, Y-, AND Z-DISPLACEMENT OF WEIGHTED MEAN '/
     + 'POSITION OF INTERPOLATED OBSERVATIONS FROM GRID POINT '/
     + 'LOCATION WILL BE ',3F8.2)
c      write(1,*)'Enter name of output file'
c      READ(LUIN,'(A80)')OFILNAM
C      write(1,*)'Enter ruv or wind file'
C      READ(1,*)KEYWORDOUT
      KEYWORDOUT='WIND'
      IF(KEYWORDB.EQ.'ruv '.OR.KEYWORDB.EQ.'RUV ')THEN
       IUVB=1
      ELSE
       IUVB=0
      ENDIF
C      WRITE(1,*)'XACC = ',XACC
C      WRITE(1,*)'YACC = ',YACC
C      WRITE(1,*)'ZACC = ',ZACC
C      WRITE(1,*)'Enter kmin and kmax for computation'
C      READ(LUIN,*)KMIN2,KMAX2
      KMIN2=1
      KMAX2=KMAXB
      WRITE(6,1520)KMIN2,KMAX2
1520  FORMAT('ENTERED MINIMUM AND MAXIMUM LEVEL TO ACTUALLY DO '/
     + 'INTERPOLATION WERE ',2I4)
      IF(KMAX2.GT.KMAXB)KMAX2=KMAXB
      IF(KMAX2.LT.KMAXB)THEN
       WRITE(6,*)'KMAX2 .LT. KMAXB'
       WRITE(6,*)'DID YOU REALLY MEAN TO DO THIS'
       READ(5,'(A)')IYES
       IF(IYES.NE.'Y'.OR.IYES.NE.'y')STOP
      ENDIF
      WRITE(6,1521)KMIN2,KMAX2
1521  FORMAT('MINIMUM AND MAXIMUM LEVEL TO ACTUALLY DO '/
     + 'INTERPOLATION WILL BE ',2I4)
C      WRITE(1,*)'Enter maximum acceptable elevation angle'
      SINTOL=90.
C      READ(LUIN,*)SINTOL
      WRITE(6,1522)SINTOL
1522  FORMAT('MAXIMUM ANGLE RELATIVE TO HORIZONTAL TO INCLUDE '/
     + 'OBSERVATIONS IN INTERPOLATION WILL BE ',F8.2,' DEGREES')
      IF(SINTOL.GE.90.)THEN
       SINTOL=1.
      ELSE
       SINTOL=SINTOL*3.14159/180.
       SINTOL=SIN(SINTOL)
      ENDIF
      WRITE(6,*)'Max abs(sin(elev)) = ',SINTOL      
C      WRITE(1,*)'Enter maximum acceptable condition number'
C      READ(LUIN,*)TOLERANCE
      TOLERANCE=1000.
C      WRITE(1,*)'Would you like to enter a fallspeed correction ',
C     + '(constant) y(1)/n(0)'
C      READ(LUIN,*)IFALLCOR
      IFALLCOR=0
      IF(IFALLCOR.EQ.1)THEN
       WRITE(6,*)'Enter correction amount'
       READ(LUIN,*)FALLCOR
      ELSE
       FALLCOR=0.
      ENDIF
C      WRITE(1,*)'Would you like to enter angle corrections other ',
C     +  'than those used in DOP3D or FAST3D y(1)/n(0)?'
C      READ(LUIN,*)IELCOR
      IELCOR=0
      IF(IELCOR.EQ.1)THEN
       WRITE(6,*)'Enter elevation correction in deg, approx ground ',
     +  'speed in m/s'
       READ(LUIN,*)ELCORR,GS
       WRITE(6,*)'Enter pitch and drift corrections in degrees'
       READ(LUIN,*)PCOR,DCOR
       PCOR=SIN(3.14159/180.*PCOR)*GS
       DCOR=SIN(3.14159/180.*DCOR)*GS
       EVRCOR=SIN(3.14159*ELCORR/180.)*GS
      ELSE
       GS=0.
       PCOR=0.
       DCOR=0.
       EVRCOR=0.
      ENDIF
C Test to make sure all files actually exist
C      WRITE(1,*)'Do you want to compute vertical wind y(1)/n(0)?'
C      READ(LUIN,*)IVWIND
      IVWIND=0
c      write(1,*)'i got here'
C      write(1,*)'Enter horizontal, vertical smoothing factors'
C      READ(LUIN,*)FACTRH,FACTRZ
C      WRITE(1,*)'Enter continuity factor,FACTWP,FACTWBOT'
C      READ(LUIN,*)FACTZ,FACTWP,FACTWBOT,FACTVR,FACTVR2
      FACTM=1.
C      WRITE(1,*)'Use the continuity equation (1) or solve ',
C     + 'directly (0, you need 3 Vrs'
C      READ(LUIN,*)ICONTINUITY
C      WRITE(1,*)'Is there a starting-guess wind3d file? y(1)/n(0)'
C      READ(1,*)IGUESS
      IGUESS=0
C      WRITE(1,*)'Use starting guess only as a mask? y(1)/n(0)'
C      READ(1,*)IGUESSER
      IGUESSER=0
C      WRITE(1,*)'Reverse the Vrs (very old file) y(1)/n(0)'
C      READ(1,*)IREVERSE
      IREVERSE=1
C      WRITE(1,*)'Do you want to read a SUMFILE y(1)/n(0)'
C      READ(LUIN,*)IREADABCD
C      WRITE(1,*)'Do you want to write a SUMFILE y(1)/n(0)'
C      READ(LUIN,*)IWRITEABCD
C      WRITE(1,*)'Do you want to use SUMFILE values? (1) or '
C      WRITE(1,*)'just pass them on (0)?'
C      READ(LUIN,*)IREADABCD
      IREADABCD=0
      IWRITEABCD=1
C      WRITE(1,*)'Joss and Waldvogel (0) or Willis gamma (1)'
C      READ(LUIN,*)IRSW
C      WRITE(1,*)'Enter melting band height, depth in original run'
C      READ(LUIN,*)HB1,DPB1
C      WRITE(1,*)'Enter new melting band height, depth'
C      READ(LUIN,*)HB2,DPB2
C      WRITE(1,*)'Enter dbz slope correction factor'
C      READ(LUIN,*)SLOPE1
C      WRITE(1,*)'Enter Dotmin'
C      READ(LUIN,*)DOTMIN
C      WRITE(1,*)'Enter weight for u-v and w boundary condition normal'
C      READ(LUIN,*)BOUND,WBOUND,WTOP
C      WRITE(1,*)'Enter range delay'
C      READ(LUIN,*)RDEL
C      WRITE(1,*)'Enter minimum dBZ and range to include data'
C      READ(1,*)DBZTEST,RRRR
      DBZTEST=-2000.
      RRRR=10.
C      WRITE(1,*)'Enter zlow and zhigh above melting level'
C      READ(LUIN,*)ZLOW,ZHIGH
C      WRITE(1,*)'Set a bottom boundary condition? y(1)/n(o)'
C      READ(1,*)IBOTTOP
      IBOTTOP=1
C      WRITE(1,*)'Enter maximum number of iterations and number of'
C      WRITE(1,*)'iterations between writing interim solution'
C      READ(1,*)ITMAX,ITMAXDIV
      ITMAX=50000 
      ITMAXDIV=50000
C      WRITE(1,*)'How many independent VRs are needed to find a '
C      WRITE(1,*)'wind solution (should be greater than or equal to 2)'
C      READ(1,*)IVRTEST
      IVRTEST=3
C      WRITE(1,*)'Use rams files for position--yes(1) or no (0)?'
C      READ(1,*)IRAMSFILES
c      WRITE(1,*)'ENTER 2 CHARACTER VOLUME NAME FOR VELOCITY'
      READ(LUIN,'(A2)')CHARV
C      WRITE(1,*)'ENTER 2 CHARACTER VOLUME NAME FOR DBZ'
      READ(LUIN,'(A2)')CHARDBZ
C      WRITE(1,*)'Subtract winds as echo motion yes(1) or no (0)?'
C      READ(LUIN,*)IECHO
      IECHO=0
      READ(LUIN,*)IBIL
      IF(IBIL.EQ.1)THEN
       WRITE(6,*)'Pseudo-bilinear interpolation will be used'
      ELSE
       WRITE(6,*)
     +  'Gaussian distance weighted interpolation will be used'    
      ENDIF
      IRAMSFILES=0
c      IF(IRAMSFILES.EQ.1)THEN
c       DO NF=1,NFILE
c        READ(LUIN,'(A80)')RAMSFILES(NF)
c       ENDDO
c      ENDIF
      IF(IGUESS.EQ.1)THEN
C
C If desired, enter first-guess wind field
C
       WRITE(6,*)'Enter WIND3D file name'
       READ(6,'(A80)')WINDFILE
       LULDPLN=99
       CALL READHEADER(LULDPLN,WINDFILE,KEYWORD,FLTNAME,
     +   STMNAME,RADAR,EXPERIMENT,
     +   CREATIME,EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,
     +   IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     +   SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     +   POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)
       WRITE(6,*)CREATIME
       WRITE(6,*)FLTNAME
       WRITE(6,*)STMNAME
       WRITE(6,*)'KEYWORD,RADAR,EXPERIMENT'
       WRITE(6,'(A4,2X,A4,2X,A32)')KEYWORD,RADAR,EXPERIMENT
       WRITE(6,*)'EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,IUNFLD,IATTEN'
       WRITE(6,*)EXTRA1,IMAXB,JMAXB,KMAXB,KOUNT,NMOSM,IUNFLD,IATTEN
       WRITE(6,*)'FLAG,EXTRA2'
       WRITE(6,*)FLAG,EXTRA2
       WRITE(6,*)'EXTRA3,STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ'
       WRITE(6,*)EXTRA3,STIME,ETIME,OLAT,OLON,SXB,SYB,SZB,
     +           XZB,YZB,ZZB
       WRITE(6,*)'ROT,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL'
       WRITE(6,*)ROTB,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL
       WRITE(6,*)'AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7'
       WRITE(6,*)AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7
       IF(KEYWORD.EQ.'ruv '.OR.KEYWORD.EQ.'RUV ')THEN
        IUV=1
       ELSE
        IUV=0
       ENDIF
       CALL LDPLNW(LULDPLN,UGUESS,VGUESS,WGUESS,IMAXB,JMAXB,KMAXB,IUV,
     1     ZZB,SZB)
       CLOSE(LULDPLN)
       CALL DIVER3(UGUESS,VGUESS,WGUESS,
     1             IMAXB,JMAXB,KMAXB,SXBD,SYBD,SZBD,DIV)
      ELSE
       DO IJK=1,N1
        UGUESS(IJK)=0.
        VGUESS(IJK)=0.
        WGUESS(IJK)=0.
       ENDDO
       IGUESS=1
      ENDIF

      DO NF=1,NFILE
       write(6,*)'at the beginning of nfile do loop nf = ',nf
C
C Read Scan File No. NF 
C
C
C sum effects from present wind file, and compute winds if last file
C in iteration
C
       XZD=XZ-XSHIFT(NF)
       YZD=YZ-YSHIFT(NF)
       IMAX2=IMAXB+2
       JMAX2=JMAXB+2
       XZD=XZ-XSHIFT(NF)
       YZD=YZ-YSHIFT(NF)
       WRITE(6,*)'CALLING ABCD'
       IF(IRAMSFILES.EQ.1)THEN
        CALL OPENRAMFILE(RAMSFILES(NF),LUOUT,STIME,ETIME,DATARATE,IERR)
       ENDIF
       write(6,'("about to open file ",A80)')filnam(nf)
       linuxyes=1
       IF(LINUXYES.EQ.1)THEN
        OPEN(80,FILE=FILNAM(NF),FORM='UNFORMATTED',RECL=2,
     1         ACCESS='DIRECT')
       ELSE
        OPEN(80,FILE=FILNAM(NF),FORM='UNFORMATTED')
       ENDIF
       write(6,*)'other side of open'
       TSTART=TSTARTARRAY(NF)
       TEND=TENDARRAY(NF)
       SLOPE1=DBZFACT(NF)
       CALL ABCD(U,V,W,
     1   RHMAX,RZMAX,EHMAX,EZMAX,XACC,YACC,ZACC,NF,NFILE,
     2   SUM1,SUM2,SUM3,SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,
     2   SUM10,SUM11,SUM12,SUM13,EVRCOR,
     3   SUMWEIGHT,TOLERANCE,KMIN2,KMAX2,
     4   SINTOL,SXB,SYB,SZB,XZB,YZB,ZZB,ROTB,IMAXB,JMAXB,KMAXB,SUMDB,
     6   BETATEST,SUMWTSAVE,XWEIGHT,YWEIGHT,ZWEIGHT,
     6   UNEW,IGUESS,FACTRH,
     7   ICONTINUITY,IVWIND,FACTZ,IREADABCD,OLAT,OLON,LUSCAN,
     9   IRSW,HB1,DPB1,HB2,DPB2,SLOPE1,IFFILTER,CENTIME,TSTART,TEND,
     1   SMOTIONU,SMOTIONV,FLTNAME,RDEL,DBZTEST,RRRR,ZLOW,ZHIGH,INDEX1,
     1   UGUESS,VGUESS,WGUESS,IREVERSE,IRAMSFILES,XSHIFT,YSHIFT,IECHO,
     1   DISTXYZ,CHARV,CHARDBZ,LINUXYES,IW0,IBIL)
c       WRITE(1,*)'AFTER LEAVING ABCD, FACTZ,FACTR = ',FACTZ,FACTRH
       CLOSE(80)
       IF(IWRITEABCD.EQ.1)THEN
        LUW=95
         OPEN(LUW,FILE=OSUMFILENAME(NF),FORM='UNFORMATTED')
         CALL WRITEABCD(LUW,IMAXB,JMAXB,KMAXB,SUM1,SUM2,SUM3,
     +   SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,
     +   SUM10,SUM11,SUM12,SUM13,SUMDB,      
     +   XWEIGHT,YWEIGHT,ZWEIGHT)
       write(6,*)'ended writeabcd for file ',nf
       ENDIF
      ENDDO
      STOP
      END
      SUBROUTINE WZERO(IMAXB,JMAXB,KMAXB,KMAX2,SUM1,SUM2,SUM3,
     + SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,SUM10,SUM11,SUM12,SUM13,INDEX1,
     + OSUMFILENAME,NFILE,SUMWTSAVE,SUMDBZ,SUMDB,
     + XWSAVE,YWSAVE,ZWSAVE,XACC,YACC,ZACC,IVRTEST,U,V)
      REAL SUM1(IMAXB,JMAXB,KMAXB),SUM2(IMAXB,JMAXB,KMAXB)
      REAL SUM3(IMAXB,JMAXB,KMAXB),SUM4(IMAXB,JMAXB,KMAXB)
      REAL SUM5(IMAXB,JMAXB,KMAXB),SUM6(IMAXB,JMAXB,KMAXB)
      REAL SUM7(IMAXB,JMAXB,KMAXB),SUM8(IMAXB,JMAXB,KMAXB)
      REAL SUM9(IMAXB,JMAXB,KMAXB)
      REAL SUMWTSAVE(IMAXB,JMAXB,KMAXB)
      REAL SUM10(IMAXB,JMAXB,KMAXB),SUM11(IMAXB,JMAXB,KMAXB)
      REAL SUM12(IMAXB,JMAXB,KMAXB),SUM13(IMAXB,JMAXB,KMAXB)
      REAL U(IMAXB,JMAXB,KMAXB),V(IMAXB,JMAXB,KMAXB)
      REAL SUMDBZ(IMAXB,JMAXB,KMAXB),SUMDB(IMAXB,JMAXB,KMAXB)
      REAL XWSAVE(IMAXB,JMAXB,KMAXB),YWSAVE(IMAXB,JMAXB,KMAXB)
      REAL ZWSAVE(IMAXB,JMAXB,KMAXB)
      INTEGER INDEX1(IMAXB,JMAXB,KMAXB)
      CHARACTER*80 OSUMFILENAME(80)
      FLAG=-1.0E+10
      DO K=1,KMAXB
       DO J=1,JMAXB
        DO I=1,IMAXB
         INDEX1(I,J,K)=0
         SUMWTSAVE(I,J,K)=0.
         SUMDBZ(I,J,K)=FLAG
        ENDDO
       ENDDO
      ENDDO
      LUW=95
      DO NF=1,NFILE
       OPEN(LUW,FILE=OSUMFILENAME(NF),FORM='UNFORMATTED')
       CALL READABCD(LUW,IMAXB,JMAXB,KMAXB,SUM1,SUM2,SUM3,SUM4,
     +               SUM5,SUM6,SUM7,SUM8,SUM9,
     +               SUM10,SUM11,SUM12,SUM13,
     +    SUMDB,KMAX2,XWSAVE,YWSAVE,ZWSAVE)
       DO K=1,KMAX2
        DO J=1,JMAXB
         DO I=1,IMAXB
          IF(NF.EQ.1)THEN
           U(I,J,K)=0.
           V(I,J,K)=0.
          ENDIF

          IF((SUM3(I,J,K).EQ.0.AND.SUM5(I,J,K).EQ.0.AND.
     +     SUM8(I,J,K).EQ.0.).OR.
     +      ABS(XWSAVE(I,J,K)).GT.XACC.OR.
     +      ABS(YWSAVE(I,J,K)).GT.YACC.OR.
     +      ABS(ZWSAVE(I,J,K)).GT.ZACC)THEN
           CONTINUE
          ELSE
           SUMWTSAVE(I,J,K)=SUMWTSAVE(I,J,K)+
     +                      SUM3(I,J,K)+SUM5(I,J,K)+SUM8(I,J,K)
           XVEC=SUM4(I,J,K)
           YVEC=SUM2(I,J,K)
           SIZE=SQRT(XVEC*XVEC+YVEC*YVEC)
           XVEC=XVEC/SIZE
           YVEC=YVEC/SIZE
           IF(NF.EQ.1)THEN
            U(I,J,K)=XVEC
            V(I,J,K)=YVEC
            INDEX1(I,J,K)=INDEX1(I,J,K)+1
           ELSE
            DOT=U(I,J,K)*XVEC+V(I,J,K)*YVEC
            IF(ABS(DOT).LT..95)THEN
             INDEX1(I,J,K)=INDEX1(I,J,K)+1
            ENDIF
           ENDIF
          ENDIF
          IF(ABS(XWSAVE(I,J,K)).LT.XACC.AND.
     +       ABS(YWSAVE(I,J,K)).LT.YACC.AND.
     +       ABS(ZWSAVE(I,J,K)).LT.ZACC)THEN
           IF(SUMDB(I,J,K).GT.SUMDBZ(I,J,K))SUMDBZ(I,J,K)=SUMDB(I,J,K)
           IF(SUMDBZ(I,J,K).GT.75.)SUMDBZ(I,J,K)=75.
          ENDIF
c         write(6,*)'i,j,k,sum3,sum5,sum8,indexbefore,indexafter = ',
c     1       i,j,k,sum3(i,j,k),sum5(i,j,k),sum8(i,j,k),
c     1       indexbefore,index1(i,j,k)
         ENDDO
        ENDDO
       ENDDO
       CLOSE(LUW)
      ENDDO
      DO K=1,KMAXB
       DO J=1,JMAXB
        DO I=1,IMAXB
         U(I,J,K)=0.
         V(I,J,K)=0.  
         indexbefore=index1(i,j,k)
         IF(INDEX1(I,J,K).GE.IVRTEST)THEN
          INDEX1(I,J,K)=1
         ELSEIF(INDEX1(I,J,K).GT.0)THEN
          INDEX1(I,J,K)=-1
         ELSE
          INDEX1(I,J,K)=0
         ENDIF
c         write(6,*)'i,j,k,indexbefore,indexafter,ivrtest = ',
c     1              i,j,k,indexbefore,index1(i,j,k),ivrtest
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END

