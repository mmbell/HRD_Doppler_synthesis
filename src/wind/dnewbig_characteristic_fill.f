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
      CHARACTER KEYWORDB*4,FLNAM*80
      CHARACTER COMYES*1,WINDFILE*80
      CHARACTER SUMFILENAME(2)*80,OSUMFILENAME(80)*80
      INTEGER KT(100000),ITOP(100000)
      INTEGER IMAXB,JMAXB,KMAXB
      INTEGER IVTSUB(80)
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
       INTEGER N1,N2,NARRAY,NWORK,NPOINTS
c      PARAMETER(N1=24000,N2=48000)
c      PARAMETER(NARRAY=2350000,NWORK=550000)
c      PARAMETER (NPOINTS=72100)

c      REAL SUM1(N1),SUM2(N1),SUM3(N1)
c      REAL SUM4(N1),SUM5(N1),SUM6(N1)
c      REAL SUM7(N1),SUM8(N1),SUM9(N1)
c      REAL SUM10(N1),SUM11(N1),SUM12(N1)
c      REAL SUM13(N1),SUMWTSAVE(N1)
c      REAL SUMDBZ(N1),SUMDB(N1)
c      REAL UGUESS(N1),VGUESS(N1),WGUESS(N1)
c      REAL XWEIGHT(N1),YWEIGHT(N1),ZWEIGHT(N1)
c      DOUBLE PRECISION RHS(NPOINTS)
c      REAL U(N1),V(N1),W(N1),BAND1(NPOINTS)
c      REAL DIV(N1)
      CHARACTER RAMSFILES(80)*80,DIVERFILE*80
c      DOUBLE PRECISION X12(NPOINTS)
c      DOUBLE PRECISION M(NARRAY)
c      DOUBLE PRECISION WORK(NWORK),
      DOUBLE PRECISION SXBD,SYBD,SZBD
c     INTEGER IWORK(NARRAY)
c      INTEGER IA(NPOINTS)
c      INTEGER JA(NARRAY)
c      REAL VELS(NPOINTS)
c      INTEGER INDEX44(N1),INDEX5(N1),INDX5(N1)
c      INTEGER INDEX1(NPOINTS),INDEX2(NPOINTS),INDEX4(N1)
c      INTEGER*2 ICAST(N1)

      real, dimension (:), allocatable :: SUM1, SUM2, SUM3, SUM4,
     +   SUM5, SUM6, SUM7, SUM8, SUM9, SUM10, SUM11, SUM12, SUM13,
     +   SUMWTSAVE, SUMDBZ, SUMDB, UGUESS, VGUESS, WGUESS,
     +   XWEIGHT, YWEIGHT, ZWEIGHT, U, V, W, BAND1, DIV, VELS
      double precision, dimension (:), allocatable :: M, X12, WORK, RHS
      integer, dimension (:),  allocatable :: IA, JA, IWORK,
     +   INDEX44, INDEX5, INDX5, INDEX1, INDEX2, INDEX4
      integer*2, dimension (:), allocatable :: ICAST
      integer :: err

      RFLAG=-1.0E+10
      OPEN(1,FILE='/dev/tty')
      CTR=3.14159/180.
      LUOUT=6
      LU=90
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
c      WRITE(1,*)'WARNING--With this version you should not change'
c      WRITE(1,*)'your weighting parameters'
c      WRITE(1,*)'DBZ ONLY (no winds) Y(1)/N(0)?'
c      READ(1,*)IDBZONLY
      IDBZONLY=0
c      WRITE(1,*)'Would you like to use a command file?'
C      READ(1,'(A)')COMYES
      COMYES='Y'
      IF(COMYES.EQ.'Y'.OR.COMYES.EQ.'y')THEN
       LUIN=16
       nargs = iargc()
       IF ( nargs.GT.0 ) THEN
          Call getarg(1,FLNAM)
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
1501  FORMAT('Flight name entered is: ',A8)
c      WRITE(1,*)'Enter Experiment Name 32 characters or less'
      READ(LUIN,'(A32)')EXPERIMENT
      WRITE(6,1502)EXPERIMENT
1502  FORMAT('Experiment name entered is: ',A32)
c      WRITE(1,*)'How many DOP3D files will be used in this program?'
      READ(LUIN,*)NFILE
      WRITE(6,1503)NFILE
1503  FORMAT('Number of scan files to be synthesized is ',I3)
      IF(NFILE.GT.40)THEN
       WRITE(1,*)'Too many files'
       STOP
      ENDIF
      ETIME=0
      STIME=200000
      DO NF=1,NFILE
c       write(1,*)'Enter Doppler scanfile name for file # ',NF
c       READ(LUIN,'(A80)')FILNAM(NF)
c       WRITE(1,*)'How much would you like to multiply dbz by?'
c       READ(LUIN,*)DBZFACT(NF),XSHIFT(NF),YSHIFT(NF)
c       WRITE(1,*)'Enter beginning and ending times to use data in',
c     +  ' HHMMSS HHMMSS'
       READ(LUIN,'(3I2,X,3I2)')IB,JB,KB,IE,JE,KE
       WRITE(6,1506)NF,IB,JB,KB
1506   FORMAT('Start time for file ',I3,' is ',3I2)
       WRITE(6,1507)NF,IE,JE,KE
1507   FORMAT('End time for file ',I3,' is ',3I2)
       TSTART=IB*3600.+JB*60+KB
       STIMENEW=TSTART
       TEND=IE*3600.+JE*60.+KE
       ETIMENEW=TEND
       IF(ETIMENEW.LT.STIMENEW)THEN
        TEND=ETIMENEW+86400.
        ETIMENEW=TEND
       ELSE
        TEND=ETIME
       ENDIF
       IF(STIMENEW.LT.STIME)STIME=STIMENEW
       IF(ETIMENEW.GT.ETIME)ETIME=ETIMENEW
       WRITE(6,*)'STIME,ETIME,TSTART,TEND = ',
     +            STIME,ETIME,TSTART,TEND 
       WRITE(6,*)'STIMENEW,STIME = ',STIMENEW,STIME
       WRITE(6,*)'ETIMENEW,ETIME = ',ETIMENEW,ETIME             
       TSTARTARRAY(NF)=TSTART
       TENDARRAY(NF)=TEND
c       WRITE(1,*)'Enter the SUMFILE for this flight leg'
c       READ(LUIN,*)IVTSUB(NF)
       IVTSUB(NF)=1
       IF(IVTSUB(NF).NE.0)THEN
        WRITE(6,*)'Fallspeed projection will be subtracted ',
     +   'from this file'
       ELSE
        WRITE(6,*)'Fallspeed projection will not be subtracted ',
     +   'from this file'
       ENDIF
       READ(LUIN,'(A80)')OSUMFILENAME(NF)
       WRITE(6,1508)NF,OSUMFILENAME(NF)
1508   FORMAT('Input sumfile name for file ',I3,' is ',A56)
      ENDDO
c      WRITE(1,*)'Mean storm velocity for this flight leg'
c      READ(LUIN,*)SMOTIONU,SMOTIONV
c      WRITE(6,*)'Enter center time HHMMSS'
c      READ(LUIN,'(3I2)')IH,IM,IS
c      CENTIME=IH*3600+IM*60+IS
c      IF(CENTIME.LT.STIME)CENTIME=CENTIME+86400.
c      WRITE(1,*)'Enter latitude and longitude of anchor point'
      READ(LUIN,*)OLAT,OLON
      WRITE(6,1511)OLAT,OLON
1511  FORMAT('Latitude and longitude of anchor point are ',2F9.3)
      WRITE(6,*)'OLAT = ',OLAT,'   OLON = ',OLON
      I1REMOVE=0
c      WRITE(1,*)'Enter SX,SY,SZ OF OUTPUT FILE'
      READ(LUIN,*)SXBD,SYBD,SZBD
      WRITE(6,1512)SXBD,SYBD,SZBD
1512  FORMAT('x-, y-, and z-resolution of grid (km) are ',3F8.2)
      SXB=SXBD
      SYB=SYBD
      SZB=SZBD
c      WRITE(1,*)'Enter XZ,YZ,ZZ,ROT OF OUTPUT FILE'
      READ(LUIN,*)XZB,YZB,ZZB,ROTB
      WRITE(6,1513)XZB,YZB
1513  FORMAT('Anchor position (x,y) relative to lower-left corner ',
     + 'in synthesis will be ',2F8.2)
      WRITE(6,1514)ZZB
1514  FORMAT('Vertical level 1 will be at height ',F8.2,' km')
      WRITE(6,1515)ROTB
1515  FORMAT('Clockwise rotation of analysis y-axis from north ',
     + 'will be ',F8.2,' degrees')
c      ZZB=0.
c      WRITE(1,*)'Enter IMAX,JMAX,KMAX of output file'
      READ(LUIN,*)IMAXB,JMAXB,KMAXB
      WRITE(6,1516)IMAXB,JMAXB,KMAXB
1516  FORMAT('Integral x-, y-, and z-domains will be ',3I4)
c      write(1,*)'Enter horizontal, vertical maximum radii of ',
c     1 ' influence in km'       
      NCHECK=IMAXB*JMAXB*KMAXB
      N1=NCHECK
      N2=2*N1
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
      allocate ( INDEX44(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( INDEX5(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( INDX5(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( INDEX4(N1), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( ICAST(N1), stat=err)
      if (err /=0) print *, "allocate failed."

      NPOINTS = N1*3+1
      allocate ( RHS(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( BAND1(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( X12(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( IA(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( VELS(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( INDEX1(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."
      allocate ( INDEX2(NPOINTS), stat=err)
      if (err /=0) print *, "allocate failed."

      NWORK = N1*18 + 100000
c NWORK based on 50000 ITMAX defined below -MMB
      allocate ( WORK(NWORK), stat=err)
      if (err /=0) print *, "WORK allocate failed."

      NARRAY = N1*135
      allocate ( M(NARRAY), stat=err)
      if (err /=0) print *, "M allocate failed."
      allocate ( IWORK(NARRAY), stat=err)
      if (err /=0) print *, "IWORK allocate failed."
      allocate ( JA(NARRAY), stat=err)
      if (err /=0) print *, "JA allocate failed."

      IF(NCHECK.GT.N1)THEN
       WRITE(6,*)'Total number of grid points is ',NCHECK
       WRITE(6,*)'Which is greater than the allowed ',N1
      ELSE
       WRITE(6,*)'N1,NCHECK = ',N1,NCHECK
      ENDIF

      DO N=1,100000
       ITOP(N)=0
      ENDDO
      DO N=1,NPOINTS
       VELS(N)=0.
      ENDDO
      DO N=1,N1
       U(N)=RFLAG
       V(N)=RFLAG
       DIV(N)=0.     
       W(N)=RFLAG
       SUMDBZ(N)=RFLAG
      ENDDO

c      READ(LUIN,*)RHMAX,RZMAX
c      write(1,*)'Enter horizontal, vertical e-folding radius'
c      READ(LUIN,*)EHMAX,EZMAX
c      write(1,*)'Enter acceptable x,y,z displacement'
c      READ(LUIN,*)XACC,YACC,ZACC
      XACC=50.
      YACC=50.
      ZACC=50.
c      write(1,*)'Enter name of output file'
c      READ(LUIN,'(A80)')OFILNAM
c      write(1,*)'Enter ruv or wind file'
C      READ(1,*)KEYWORDOUT
      KEYWORDOUT='WIND'
      IF(KEYWORDB.EQ.'ruv '.OR.KEYWORDB.EQ.'RUV ')THEN
       IUVB=1
      ELSE
       IUVB=0
      ENDIF
      WRITE(6,*)'XACC = ',XACC
      WRITE(6,*)'YACC = ',YACC
      WRITE(6,*)'ZACC = ',ZACC
c      WRITE(1,*)'Enter kmin and kmax for computation'
C      READ(LUIN,*)KMIN2,KMAX2
      KMIN2=1
      KMAX2=KMAXB
      WRITE(6,1519)KMIN2,KMAX2
1519  FORMAT('Entered kmin2, kmax2 are ',2I4)
      IF(KMAX2.GT.KMAXB)KMAX2=KMAXB
      WRITE(6,1520)KMIN2,KMAX2
1520  FORMAT('KMIN2,KMAX2 are ',2I4)
c      WRITE(1,*)'Enter maximum acceptable elevation angle'
c      READ(LUIN,*)SINTOL
c      IF(SINTOL.GE.90.)THEN
c       SINTOL=1.
c      ELSE
c       SINTOL=SINTOL*3.14159/180.
c       SINTOL=SIN(SINTOL)
c      ENDIF
c      WRITE(1,*)'Max abs(sin(elev)) = ',SINTOL      
c      WRITE(1,*)'Enter maximum acceptable condition number'
c      READ(LUIN,*)TOLERANCE
      TOLERANCE=1000.
C      WRITE(1,*)'Would you like to enter a fallspeed correction ',
C     + '(constant) y(1)/n(0)'
C      READ(LUIN,*)IFALLCOR
C      IF(IFALLCOR.EQ.1)THEN
C       WRITE(1,*)'Enter correction amount'
C       READ(LUIN,*)FALLCOR
C      ELSE
       FALLCOR=0.
C      ENDIF
C      WRITE(1,*)'Would you like to enter angle corrections other ',
C     +  'than those used in DOP3D or FAST3D y(1)/n(0)?'
C      READ(LUIN,*)IELCOR
C      IF(IELCOR.EQ.1)THEN
C       WRITE(1,*)'Enter elevation correction in deg, approx ground ',
C     +  'speed in m/s'
C       READ(LUIN,*)ELCORR,GS
C       WRITE(1,*)'Enter pitch and drift corrections in degrees'
C       READ(LUIN,*)PCOR,DCOR
C       PCOR=SIN(3.14159/180.*PCOR)*GS
C       DCOR=SIN(3.14159/180.*DCOR)*GS
C       EVRCOR=SIN(3.14159*ELCORR/180.)*GS
C      ELSE
       GS=0.
       PCOR=0.
       DCOR=0.
       EVRCOR=0.
C      ENDIF
C Test to make sure all files actually exist
C      WRITE(1,*)'Do you want to compute vertical wind y(1)/n(0)?'
C      READ(LUIN,*)IVWIND
      IVWIND=1
c      write(1,*)'i got here'
C      write(1,*)'Enter horizontal, vertical smoothing factors'
      FACTRH=.01
      FACTRZ=100.
c      READ(LUIN,*)FACTRH,FACTRZ
      WRITE(6,1540)FACTRH,FACTRZ
1540  FORMAT('Horizontal/vertical filter factors are ',F8.3)
c      WRITE(1,*)'Enter continuity factor,FACTWP,FACTWBOT'
C      READ(LUIN,*)FACTZ,FACTWP,FACTWBOT,FACTVR,FACTVR2,FACTWTOP
      READ(LUIN,*)FACTZ
C      FACTZ=1.
      FACTWP=0.
      FACTWBOT=1.
      FACTWTOP=1.
      FACTVR2=1.
      FACTVR=(SZB*1000.)+1.
      WRITE(6,*)'4th order continuity factor is ',FACTZ
      WRITE(6,*)'2nd order continuity factor is ',FACTWP
      WRITE(6,*)'bottom,top boundary factors are ',FACTWBOT,FACTWTOP
      WRITE(6,*)'Doppler radial weight is ',FACTVR2
      WRITE(6,*)'Thickness of top in meters is ',FACTVR
      FACTM=1.
c      WRITE(1,*)'Use the continuity equation (1) or solve ',
c     + 'directly (0, you need 3 Vrs'
C      READ(LUIN,*)ICONTINUITY
      ICONTINUITY=1
      IF(ICONTINUITY.EQ.1)THEN
       WRITE(6,*)'Synthesis will use continuity equation'
      ELSEIF(ICONTINUITY.EQ.0)THEN
       WRITE(6,*)'Synthesis will be determined from radials only'
      ELSE
       WRITE(6,*)'Problem with ICONTINUITY'
       STOP
      ENDIF
c      WRITE(1,*)'Is there a starting-guess wind3d file? y(1)/n(0)'
c      READ(1,*)IGUESS
c      IGUESS=0
c      WRITE(1,*)'Use starting guess only as a mask? y(1)/n(0)'
C      READ(1,*)IGUESSER
c      IGUESSER=0
c      WRITE(1,*)'Reverse the Vrs (very old file) y(1)/n(0)'
C      READ(1,*)IREVERSE
      IGUESS=0
      IGUESSER=0
      IREVERSE=1
C      WRITE(1,*)'Do you want to read a SUMFILE y(1)/n(0)'
C      READ(LUIN,*)IREADABCD
C      WRITE(1,*)'Do you want to write a SUMFILE y(1)/n(0)'
C      READ(LUIN,*)IWRITEABCD
C      WRITE(1,*)'Do you want to use SUMFILE values? (1) or '
C      WRITE(1,*)'just pass them on (0)?'
C      READ(LUIN,*)IREADABCD
      IREADABCD=1
      IWRITEABCD=0
C      WRITE(1,*)'Joss and Waldvogel (0) or Willis gamma (1)'
C      READ(LUIN,*)IRSW
      IRSW=0
      IF(IRSW.EQ.0)THEN
       WRITE(6,*)'Joss and Waldvogel will be used'
      ELSE
       WRITE(6,*)'Fall speeds from Willis gamma function will be used'  
      ENDIF
c      WRITE(1,*)'Enter melting band height, depth in original run'
c      READ(LUIN,*)HB1,DPB1
c      WRITE(1,*)'Enter new melting band height, depth'
      READ(LUIN,*)HB2,DPB2
      HB1=HB2
      DPB1=DPB2
      WRITE(6,1570)HB2,DPB2
1570  FORMAT('Melting-band height, depth are ',2F8.3,' km')
C      WRITE(1,*)'Enter dbz slope correction factor'
C      READ(LUIN,*)SLOPE1
C      WRITE(1,*)'Enter Dotmin'
C      READ(LUIN,*)DOTMIN
      DOTMIN=.9999999
C      WRITE(1,*)'Enter weight for u-v and w boundary condition normal'
C      READ(LUIN,*)BOUND,WBOUND,WTOP
C      WRITE(1,*)'Enter range delay'
       RDEL=0.
       BOUND=0.
       WBOUND=0.
       WTOP=0.
C      READ(LUIN,*)RDEL
C      WRITE(1,*)'Enter minimum dBZ and range to include data'
C      READ(1,*)DBZTEST,RRRR
C      DBZTEST=-2000.
C      RRRR=10.
C      WRITE(1,*)'Enter zlow and zhigh above melting level'
      READ(LUIN,*)ZLOW,ZHIGH
      WRITE(6,1560)ZLOW
1560  FORMAT('Reflecitivity below which all ice at T<0 is ',F8.3)
      WRITE(6,1561)ZHIGH
1561  FORMAT('Reflecitivity below which all water at T<0 is ',F8.3)
C      WRITE(1,*)'Set a bottom boundary condition? y(1)/n(o)'
C      READ(1,*)IBOTTOP
      IBOTTOP=1
C      WRITE(1,*)'Enter maximum number of iterations and number of'
C      WRITE(1,*)'iterations between writing interim solution'
C      READ(1,*)ITMAX,ITMAXDIV
      ITMAX=50000 
      ITMAXDIV=50000
c      WRITE(1,*)'How many independent VRs are needed to find a '
c      WRITE(1,*)'wind solution (should be greater than or equal to 2)'
c      READ(1,*)IVRTEST
      IVRTEST=2
C      WRITE(1,*)'Use rams files for position--yes(1) or no (0)?'
C      READ(1,*)IRAMSFILES
      IRAMSFILES=0
C      WRITE(1,*)'SUBTRACT VT YES (1) OR NO (0)?'
C      WRITE(1,*)'Subtract winds as echo motion yes(1) or no (0)?'
C      READ(LUIN,*)IECHO
      IECHO=0
c      WRITE(1,*)'BOTTOM (1), TOP (10), OR BOTH (11) '
c      WRITE(1,*)'OR EITHER (20) BOUND(S) NEEDED?'
c      READ(1,*)IBOUNDTEST
       IBOUNDTEST=20
c      WRITE(1,*)'ENTER HIGHEST LEVEL TO EXTEND BOTTOM B.C.'
c      READ(1,*)ZBOTMAX
      ZBOTMAX=40.
C      WRITE(1,*)'ENTER LOWEST LEVEL TO EXTEND TOP B.C.'
C      READ(1,*)ZTOPMIN
      ZTOPMIN=-1.
C      WRITE(1,*)'ENTER LOWEST LEVEL TO SET AS A TOP FOR B.C.s'
C      READ(1,*)ZTOPPER
      ZTOPPER=-1.
      IRAMSFILES=0
      IF(IRAMSFILES.EQ.1)THEN
       DO NF=1,NFILE
        READ(LUIN,'(A80)')RAMSFILES(NF)
       ENDDO
      ENDIF
      READ(LUIN,'(A80)')OFILNAM
C      WRITE(1,*)'Enter minimum angle between radials from ',
C     + 'different files'
      SINMIN=15.
c      SINMIN=20.
C      READ(1,*)SINMIN
      SINMIN=SIN(SINMIN*3.14159/180.)
C      WRITE(1,*)'ENTER WEIGHT FOR BACKGROUND FOR U,V,W'
C      READ(1,*)WTBACKU,WTBACKV,WTBACKW
      WTBACKU=0.
      WTBACKV=0.
      WTBACKW=.001
      IF(WTBACKU.GT.0..OR.WTBACKV.GT.0..OR.WTBACKW.GT.0.)THEN
C
C If desired, enter first-guess wind field
C
C       WRITE(1,*)'Set background wind to zero (0) or '
C       WRITE(1,*)'use guess field (1)?'
C       READ(1,*)ISETZERO
       ISETZERO=0
       IF(ISETZERO.NE.0)THEN
        WRITE(1,*)'Enter WIND3D file name'
        READ(1,'(A80)')WINDFILE
        LULDPLN=99

        CALL READHEADER(LULDPLN,WINDFILE,KEYWORD,FLTNAME,
     +   STMNAME,RADAR,EXPERIMENT,
     +   CREATIME,EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,
     +   IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME3,ETIME3,OLAT,OLON,
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
        WRITE(6,*)EXTRA3,STIME3,ETIME3,OLAT,OLON,SXB,SYB,SZB,
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
       ENDIF
      ELSE
       DO IJK=1,N1
        UGUESS(IJK)=0.
        VGUESS(IJK)=0.
        WGUESS(IJK)=0.
       ENDDO
       IGUESS=1
      ENDIF
      IF(IREADABCD.EQ.1)GO TO 92532
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
       OPEN(220,FILE=FILNAM(NF),FORM='UNFORMATTED')
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
     1   UGUESS,VGUESS,WGUESS,IREVERSE,IRAMSFILES,XSHIFT,YSHIFT,IECHO)
c       WRITE(1,*)'AFTER LEAVING ABCD, FACTZ,FACTR = ',FACTZ,FACTRH
       CLOSE(220)
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
92532 LUW=95
      NNN1=N1
      N4=NNN1
      NN1=NPOINTS
      NNNN1=NPOINTS
      NNNNN=NARRAY
      LUOUT=6
C
C Determine horizontal divergence of wind field
C
      CALL DIVER(U,V,DIV,IMAXB,JMAXB,KMAX2,SXB,SYB)
       LUO=91
       KEYWORDB='WIND'
      IF(KEYWORDB.EQ.'ruv '.OR.KEYWORDB.EQ.'RUV ')THEN
       IUVB=1
      ELSE
       IUVB=0
      ENDIF
       TESTFILE='TESTFILE'
       CALL WRITEHEADER(LUO,TESTFILE,KEYWORDB,FLTNAME,
     +  STMNAME,RADAR,EXPERIMENT,
     +  CREATIME,EXTRA1,IMAXB,JMAXB,KMAXB,KOUNT,NMOSM,
     +  IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     +  SXB,SYB,SZB,XZB,YZB,ZZB,ROTB,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     +  POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)
C
C Write result to output file
C
c       WRITE(1,*)'ENTERING WRPLN'
       CALL WRPLN(LUO,IMAXB,JMAXB,KMAXB,U,V,W,DIV,SUMDB,
     +  IUVB,ZZB,SZB)
c       WRITE(1,*)'LEAVING WRPLN'
      CLOSE(LUO)
c      WRITE(1,*)'BEFORE ENTERING SPARSE, FACTZ,FACTR = ',FACTZ,FACTR
      DO I=1,N1
       INDEX5(I)=0
       SUMDBZ(I)=RFLAG
      ENDDO
C      WRITE(1,*)'MULTIPLE VRS (0) OR MULTIPLE PAIRS (1)?'
C      READ(1,*)IMULTI
      IMULTI=1 
      IF(IMULTI.EQ.0)THEN
       CALL WCHECK(IMAXB,JMAXB,KMAXB,KMAX2,SUM1,SUM2,SUM3,SUM4,
     +   SUM5,SUM6,SUM7,SUM8,SUM9,SUM10,SUM11,SUM12,SUM13,INDEX4,
     +   OSUMFILENAME,NFILE,SUMWTSAVE,SUMDBZ,SUMDB,XWEIGHT,YWEIGHT,
     +  ZWEIGHT,XACC,YACC,ZACC,IVRTEST,U,V,W,SZB,ZZB)
       DO NTOP=1,100000
        ITOP(NTOP)=1
       ENDDO
      ELSE
       DO L=1,N1
        INDEX4(L)=0
       ENDDO
       DO L=NFILE,2,-1
        SUMFILENAME(2)=OSUMFILENAME(L)
        DO LL=1,L-1
         write(6,*)'calling wzero'
         SUMFILENAME(1)=OSUMFILENAME(LL)
         write(6,'(A80)')SUMFILENAME(1)
         write(6,'(A80)')SUMFILENAME(2)
         CALL WZERO(IMAXB,JMAXB,KMAXB,KMAX2,SUM1,SUM2,SUM3,SUM4,
     +   SUM5,SUM6,SUM7,SUM8,SUM9,SUM10,SUM11,SUM12,SUM13,INDEX4,
     +   SUMFILENAME,NFILE,SUMWTSAVE,SUMDBZ,SUMDB,XWEIGHT,YWEIGHT,
     +   ZWEIGHT,XACC,YACC,ZACC,IVRTEST,U,V,W,SZB,ZZB,SINMIN,
     +   VELS(1),VELS(N1+1),VELS(N1+N1+1))
        ENDDO
       ENDDO
      ENDIF
      NZ=NARRAY
      NW=NWORK
      NINSIDE=ITMAX/ITMAXDIV
       IF(IWRITEABCD.EQ.1)THEN
        DO IJK=1,N1
         UGUESS(IJK)=0.
         VGUESS(IJK)=0.
         WGUESS(IJK)=0.
        ENDDO
       ENDIF
      IF(IDBZONLY.NE.1)THEN
      write(6,*)'call sparse kmaxb,nnn1,nn1 = ',kmaxb,nnn1,nn1
       CALL SPARSE(M,RHS,U,V,IMAXB,JMAXB,KMAXB,SXBD,SYBD,SZBD,NNN1,NN1,
     + INDEX1,INDEX2,IPOS,SUM1,SUM2,SUM3,SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,
     + SUM10,SUM11,SUM12,SUM13,
     + XWEIGHT,YWEIGHT,ZWEIGHT,XACC,YACC,ZACC,
     + OSUMFILENAME,NZ,NFILE,IWORK,NW,INDEX4,
     + FACTZ,FACTRH,FACTRZ,ICAST,W,BOUND,WBOUND,ZZB,NNNN1,
     + SUMDBZ,SUMWTSAVE,DBZTEST,NORM,NNNNN,INDXN,WTOP,IBOTTOP,
     + IA,JA,IGUESSER,UGUESS,VGUESS,WGUESS,
     + BAND1,HB2,DPB2,IRSW,ZLOW,ZHIGH,
     1 FACTWP,DUMMY1,FACTVR,FACTVR2,FACTWTOP,
     1 IGUESS,FACTWBOT,DIV,KT,KMAX2,ITOP,IVTSUB,INDX5,ZTOPPER,
     1 WTBACKU,WTBACKV,WTBACKW,VELS(1),VELS(1+N1))    
c ,INDEX6)
       DO JIJ=1,NPOINTS
        VELS(JIJ)=0.
       ENDDO
       write(6,*)'index1(1) = ',index1(1)
        KEYWORDB=KEYWORDOUT
        IF(KEYWORDB.EQ.'ruv '.OR.KEYWORDB.EQ.'RUV ')THEN
         IUVB=1
        ELSE
         IUVB=0
        ENDIF
        DIVERFILE='diverfile'

        CALL WRITEHEADER(LUO,DIVERFILE,KEYWORDB,FLTNAME,
     +  STMNAME,RADAR,EXPERIMENT,
     +  CREATIME,EXTRA1,IMAXB,JMAXB,KMAXB,KOUNT,NMOSM,
     +  IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     +  SXB,SYB,SZB,XZB,YZB,ZZB,ROTB,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     +  POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)

C Write result to output file

        CALL WRPLN(LUO,IMAXB,JMAXB,KMAXB,UGUESS,VGUESS,WGUESS,
     +  DIV,SUMDBZ,
     +  IUVB,ZZB,SZB)
        CLOSE(LUO)
       IF(IGUESSER.EQ.1)THEN
        DO IJK=1,N1
         UGUESS(IJK)=0.
         VGUESS(IJK)=0.
         WGUESS(IJK)=0.
        ENDDO
       ENDIF
       DO JKI=1,NINSIDE
        WRITE(6,*)'BEGINNING SOLUTION LOOP NUMBER ',JKI
        DO IJK=1,N1
         JJK=IJK+N1
         KJK=IJK+N2
         IF((JKI.EQ.1.AND.IGUESSER.EQ.1).OR.JKI.GT.1)THEN
          VELS(IJK)=UGUESS(IJK)
          VELS(JJK)=VGUESS(IJK)
          VELS(KJK)=WGUESS(IJK)
         ENDIF
        ENDDO
        KEYWORDB=KEYWORDOUT
        IF(KEYWORDB.EQ.'ruv '.OR.KEYWORDB.EQ.'RUV ')THEN
         IUVB=1
        ELSE
         IUVB=0
        ENDIF
        WRITE(6,*)'SOLVESPARSE,STIME,ETIME = ',STIME,ETIME
        CALL SOLVESPARSE(IMAXB,JMAXB,KMAXB,NNN1,NN1,NNNN1,IPOS,SXB,
     + SYB,SZB,VELS,LUOUT,M,RHS,C,R,P,AP,INDEX1,INDEX2,DOTMIN,
     + LUO,OFILNAM,KEYWORDB,FLTNAME,
     +  STMNAME,RADAR,EXPERIMENT,
     +  CREATIME,EXTRA1,KOUNT,NMOSM,
     +  IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     +  XZB,YZB,ZZB,ROTB,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     +  POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7,
     +  DIV,SUMDBZ,SUMWTSAVE,IGUESS,X12,C12,U,V,W,IUVB,NORM,NNNNN,
     +  INDXN,UGUESS,VGUESS,WGUESS,
     +  NZ,NW,IA,JA,IWORK,WORK,ITMAXDIV,JKI,INDEX4,N4)
        DO IJK=1,N1
         JJK=IJK+N1
         KJK=IJK+N2
         U(IJK)=VELS(IJK)
         V(IJK)=VELS(JJK)
         W(IJK)=VELS(KJK)
        ENDDO
        CALL DIVER3(U,V,W,IMAXB,JMAXB,KMAXB,SXBD,SYBD,SZBD,DIV)
        LUO=91
        KEYWORDB=KEYWORDOUT
        WRITE(6,*)'FINAL FILE STIME,ETIME = ',STIME,ETIME
        CALL WRITEHEADER(LUO,OFILNAM,KEYWORDB,FLTNAME,
     +  STMNAME,RADAR,EXPERIMENT,
     +  CREATIME,EXTRA1,IMAXB,JMAXB,KMAXB,KOUNT,NMOSM,
     +  IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     +  SXB,SYB,SZB,XZB,YZB,ZZB,ROTB,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     +  POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)
C
C Write result to output file
C

        CALL WRPLN4(LUO,IMAXB,JMAXB,KMAXB,U,V,W,DIV,SUMDBZ,
     +  IUVB,ZZB,SZB,INDEX4)
        CLOSE(LUO)
        DO IJK=1,N1
         UGUESS(IJK)=U(IJK)
         VGUESS(IJK)=V(IJK)
         WGUESS(IJK)=W(IJK)
        ENDDO
       ENDDO
      ELSE
       DO IJK=1,N1
        U(IJK)=1.
        V(IJK)=1.
        W(IJK)=1.
        DIV(IJK)=0.
       ENDDO
      ENDIF
      IF(IDBZONLY.EQ.1)THEN
       OFILNAM='DBZ'
       LUO=91
       KEYWORDB=KEYWORDOUT
       CALL WRITEHEADER(LUO,OFILNAM,KEYWORDB,FLTNAME,
     +  STMNAME,RADAR,EXPERIMENT,
     +  CREATIME,EXTRA1,IMAXB,JMAXB,KMAXB,KOUNT,NMOSM,
     +  IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     +  SXB,SYB,SZB,XZB,YZB,ZZB,ROTB,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     +  POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)
C
C Write result to output file
C
       CALL MINMAX(U,V,W,IMAXB,JMAXB,KMAXB)
       CALL WRPLN(LUO,IMAXB,JMAXB,KMAXB,U,V,W,DIV,SUMDBZ,
     +  IUVB,ZZB,SZB)
      CLOSE(LUO)
      ENDIF

      STOP
      END

      SUBROUTINE WZERO(IMAXB,JMAXB,KMAXB,KMAX2,SUM1,SUM2,SUM3,
     + SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,SUM10,SUM11,SUM12,SUM13,INDEX1,
     + OSUMFILENAME,NFILE,SUMWTSAVE,SUMDBZ,SUMDB,
     + XWSAVE,YWSAVE,ZWSAVE,XACC,YACC,ZACC,IVRTEST,U,V,W,SZB,ZZB,
     + SINMIN,UTEST,VTEST,WTEST)
      REAL UTEST(IMAXB,JMAXB,KMAXB),VTEST(IMAXB,JMAXB,KMAXB)
      REAL WTEST(IMAXB,JMAXB,KMAXB)
      REAL SUM1(IMAXB,JMAXB,KMAXB),SUM2(IMAXB,JMAXB,KMAXB)
      REAL SUM3(IMAXB,JMAXB,KMAXB),SUM4(IMAXB,JMAXB,KMAXB)
      REAL SUM5(IMAXB,JMAXB,KMAXB),SUM6(IMAXB,JMAXB,KMAXB)
      REAL SUM7(IMAXB,JMAXB,KMAXB),SUM8(IMAXB,JMAXB,KMAXB)
      REAL SUM9(IMAXB,JMAXB,KMAXB)
      REAL SUMWTSAVE(IMAXB,JMAXB,KMAXB)
      REAL SUM10(IMAXB,JMAXB,KMAXB),SUM11(IMAXB,JMAXB,KMAXB)
      REAL SUM12(IMAXB,JMAXB,KMAXB),SUM13(IMAXB,JMAXB,KMAXB)
      REAL U(IMAXB,JMAXB,KMAXB),V(IMAXB,JMAXB,KMAXB)
      REAL W(IMAXB,JMAXB,KMAXB)
      REAL SUMDBZ(IMAXB,JMAXB,KMAXB),SUMDB(IMAXB,JMAXB,KMAXB)
      REAL XWSAVE(IMAXB,JMAXB,KMAXB),YWSAVE(IMAXB,JMAXB,KMAXB)
      REAL ZWSAVE(IMAXB,JMAXB,KMAXB)
      REAL RHOK(100)
      INTEGER INDEX1(IMAXB,JMAXB,KMAXB)
      CHARACTER*80 OSUMFILENAME(80)

      FLAG=-1.0E+10
      write(6,*)'entered wzero'
      HSINMIN=.5*SINMIN
      DO K=1,KMAXB
       HEIGHT=ZZB+(K-1)*SZB
       RHOK(K)=EXP(-HEIGHT*.10437052)
       DO J=1,JMAXB
        DO I=1,IMAXB
         SUMWTSAVE(I,J,K)=0.
        ENDDO
       ENDDO
      ENDDO
      LUW=95
      DO NF=1,2
       OPEN(LUW,FILE=OSUMFILENAME(NF),FORM='UNFORMATTED')
       CALL READABCD(LUW,IMAXB,JMAXB,KMAXB,SUM1,SUM2,SUM3,SUM4,
     +               SUM5,SUM6,SUM7,SUM8,SUM9,
     +               SUM10,SUM11,SUM12,SUM13,
     +    SUMDB,KMAX2,XWSAVE,YWSAVE,ZWSAVE)
       DO K=1,KMAX2
        DO J=1,JMAXB
         DO I=1,IMAXB
          IF(NF.EQ.1)THEN
           UTEST(I,J,K)=0.
           VTEST(I,J,K)=0.
           WTEST(I,J,K)=0.
           U(I,J,K)=0.
           V(I,J,K)=0.
           W(I,J,K)=0.
          ENDIF
          IF(XWSAVE(I,J,K).LT.XACC.AND.YWSAVE(I,J,K).LT.YACC.AND.
     +       ZWSAVE(I,J,K).LT.ZACC)THEN
           IF(NF.EQ.2)THEN
            RLENGTH12=UTEST(I,J,K)*UTEST(I,J,K)
     +      +VTEST(I,J,K)*VTEST(I,J,K)+WTEST(I,J,K)*WTEST(I,J,K)
            RLENGTH22=SUM10(I,J,K)*SUM10(I,J,K)
     +      +SUM11(I,J,K)*SUM11(I,J,K)+SUM12(I,J,K)*SUM12(I,J,K)
            IF(RLENGTH12.GT.0..AND.RLENGTH22.GT.0.)THEN
             RLENGTH1=SQRT(RLENGTH12)
             RLENGTH2=SQRT(RLENGTH22)
             U1=UTEST(I,J,K)/RLENGTH1
             V1=VTEST(I,J,K)/RLENGTH1
             W1=WTEST(I,J,K)/RLENGTH1
             U2=SUM10(I,J,K)/RLENGTH2
             V2=SUM11(I,J,K)/RLENGTH2
             W2=SUM12(I,J,K)/RLENGTH2
             UX=V1*W2-W1*V2
             VX=W1*U2-U1*W2
             WX=U1*V2-V1*U2
             UVW2=UX*UX+VX*VX+WX*WX
             RLENGTH=SQRT(UVW2)
             SINTEST=RLENGTH
            ELSE
             SINTEST=0.
            ENDIF
           ENDIF
           U(I,J,K)=U(I,J,K)+SUM5(I,J,K)
           V(I,J,K)=V(I,J,K)+SUM3(I,J,K)
           W(I,J,K)=W(I,J,K)+SUM8(I,J,K)
           UTEST(I,J,K)=UTEST(I,J,K)+SUM10(I,J,K)
           VTEST(I,J,K)=VTEST(I,J,K)+SUM11(I,J,K)
           WTEST(I,J,K)=WTEST(I,J,K)+SUM12(I,J,K)
          ELSE
           SINTEST=0.
          ENDIF
          IF(NF.EQ.2)THEN
           INDEXX=0
           RLENGTH2=U(I,J,K)+V(I,J,K)+W(I,J,K)
           TESTMIN=1.
           IF(RLENGTH2.GT.0.)THEN
            U1=SQRT(U(I,J,K)/RLENGTH2-(UTEST(I,J,K)/RLENGTH2)**2.)
            IF(U1.GT.HSINMIN)INDEXX=1
            V1=SQRT(V(I,J,K)/RLENGTH2-(VTEST(I,J,K)/RLENGTH2)**2.)
            IF(V1.GT.HSINMIN)INDEXX=INDEXX+1
            W1=SQRT(W(I,J,K)/RLENGTH2-(WTEST(I,J,K)/RLENGTH2)**2.)
            IF(W1.GT.HSINMIN)INDEXX=INDEXX+1
            UTEST(I,J,K)=W1
            VTEST(I,J,K)=SUM12(I,J,K)/RLENGTH2
           ELSE
            UTEST(I,J,K)=0.
            VTEST(I,J,K)=0.
           ENDIF
C           IF(J.GT.35.AND.I.LT.5)THEN
C            WRITE(6,*)'I,J,K,U1,V1,W1,UX,VX,WX,RLENGTH2,SINTEST = ',
C     +      I,J,K,U1,V1,W1,UX,VX,WX,RLENGTH2,SINTEST
C           ENDIF
           IF(INDEXX.LT.2.AND.SINTEST.GT.SINMIN)INDEXX=2
           IF(INDEXX.GE.IVRTEST)THEN
            INDEX1(I,J,K)=1
           ENDIF
c           if(k.eq.1)write(6,*)'i,j,sintest,indexx,index1,ivrtest = ',
c     1     i,j,sintest,indexx,index1(i,j,k),ivrtest
          ENDIF
          IF(SUMDB(I,J,K).GT.-100.)THEN
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
      write(6,*)'leaving wzero'

      RETURN
      END

      SUBROUTINE WCHECK(IMAXB,JMAXB,KMAXB,KMAX2,SUM1,SUM2,SUM3,
     + SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,SUM10,SUM11,SUM12,SUM13,INDEX1,
     + OSUMFILENAME,NFILE,SUMWTSAVE,SUMDBZ,SUMDB,
     + XWSAVE,YWSAVE,ZWSAVE,XACC,YACC,ZACC,IVRTEST,U,V,W,SZB,ZZB)
      REAL SUM1(IMAXB,JMAXB,KMAXB),SUM2(IMAXB,JMAXB,KMAXB)
      REAL SUM3(IMAXB,JMAXB,KMAXB),SUM4(IMAXB,JMAXB,KMAXB)
      REAL SUM5(IMAXB,JMAXB,KMAXB),SUM6(IMAXB,JMAXB,KMAXB)
      REAL SUM7(IMAXB,JMAXB,KMAXB),SUM8(IMAXB,JMAXB,KMAXB)
      REAL SUM9(IMAXB,JMAXB,KMAXB)
      REAL SUMWTSAVE(IMAXB,JMAXB,KMAXB)
      REAL SUM10(IMAXB,JMAXB,KMAXB),SUM11(IMAXB,JMAXB,KMAXB)
      REAL SUM12(IMAXB,JMAXB,KMAXB),SUM13(IMAXB,JMAXB,KMAXB)
      REAL U(IMAXB,JMAXB,KMAXB),V(IMAXB,JMAXB,KMAXB)
      REAL W(IMAXB,JMAXB,KMAXB)
      REAL SUMDBZ(IMAXB,JMAXB,KMAXB),SUMDB(IMAXB,JMAXB,KMAXB)
      REAL XWSAVE(IMAXB,JMAXB,KMAXB),YWSAVE(IMAXB,JMAXB,KMAXB)
      REAL ZWSAVE(IMAXB,JMAXB,KMAXB)
      REAL RHOK(100)
      INTEGER INDEX1(IMAXB,JMAXB,KMAXB)
      CHARACTER*80 OSUMFILENAME(80)

      FLAG=-1.0E+10
      write(6,*)'entered wcheck'
      DO K=1,KMAXB
       DO J=1,JMAXB
        DO I=1,IMAXB
         INDEX1(I,J,K)=0
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
          IF(SUMDB(I,J,K).GT.SUMDBZ(I,J,K))SUMDBZ(I,J,K)=SUMDB(I,J,K)  
          IF((SUM3(I,J,K).EQ.0.AND.SUM5(I,J,K).EQ.0.AND.
     +     SUM8(I,J,K).EQ.0.).OR.
     +      XWSAVE(I,J,K).GT.XACC.OR.
     +      YWSAVE(I,J,K).GT.YACC.OR.ZWSAVE(I,J,K).GT.ZACC)THEN
            U(I,J,K)=FLAG
            V(I,J,K)=FLAG
            W(I,J,K)=FLAG
          ELSE 
            U(I,J,K)=0
            V(I,J,K)=0
            W(I,J,K)=0
            INDEX1(I,J,K)=INDEX1(I,J,K)+1
c            WRITE(6,*)'I,J,K,INDEX1 = ',I,J,K,INDEX1(I,J,K)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
       CLOSE(LUW)
      ENDDO
      write(6,*)'leaving wcheck'
      DO K=1,KMAXB
       DO J=1,JMAXB
        DO I=1,IMAXB
         IF(INDEX1(I,J,K).GE.IVRTEST)THEN
          INDEX1(I,J,K)=1
c          WRITE(6,*)'PASSED TEST I,J,K,INDEX1 = ',
c     1     I,J,K,INDEX1(I,J,K)
         ELSE
          INDEX1(I,J,K)=0
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END



