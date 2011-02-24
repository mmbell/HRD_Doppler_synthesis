      SUBROUTINE ABCD(U,V,W,
     1   RHMAX,RZMAX,EHMAX,EZMAX,XACC,YACC,ZACC,NF,NFILE,
     2   SUM1,SUM2,SUM3,SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,
     2   SUM10,SUM11,SUM12,SUM13,EVRCOR,
     3   SUMWEIGHT,TOLERANCE,KMIN2,KMAX2,
     4   SINTOL,SXB,SYB,SZB,XZB,YZB,ZZZB,
     5   ROTB,IMAXB,JMAXB,KMAXB,SUMDBZ,
     6   BETATEST,SUMWTSAVE,XWSAVE,YWSAVE,ZWSAVE,UNEW,IGUESS,FACTR,
     7   ICONTINUITY,IVWIND,FACTZ,IREADABCD,OLAT,OLON,LUSCAN,
     9   IRSW,HB1,DPB1,HB2,DPB2,SLOPE1,IFFILTER,CENTIME,TSTART,TEND,
     1   SMOTIONU,SMOTIONV,FLIGHTID,RDEL,DBZTEST,RRRR,ZLOW,ZHIGH,INDEX1,
     1   UGUESS,VGUESS,WGUESS,IREVERSE,IRAMSFILES,XSHIFT,YSHIFT,IECHO,
     1   DISTXYZ,ALPHADIST,BETADIST,GAMMADIST,VRDIST,DBZDIST,
     1   CHARV,CHARDBZ)

      CHARACTER FLIGHTID*8, CHARV*2, CHARDBZ*2
      REAL SUMWTSAVE(IMAXB,JMAXB,KMAXB)
      REAL XWSAVE(IMAXB,JMAXB,KMAXB)
      REAL YWSAVE(IMAXB,JMAXB,KMAXB),ZWSAVE(IMAXB,JMAXB,KMAXB)
      REAL U(IMAXB,JMAXB,KMAXB),V(IMAXB,JMAXB,KMAXB)
      REAL W(IMAXB,JMAXB,KMAXB),SUMDBZ(IMAXB,JMAXB,KMAXB)
      REAL VR(1024),RANGE(1024),DBZ(1024),XR(1024),YR(1024),ZR(1024)
      REAL SUM1(IMAXB,JMAXB,KMAXB),SUM2(IMAXB,JMAXB,KMAXB)
      REAL SUM3(IMAXB,JMAXB,KMAXB),SUM4(IMAXB,JMAXB,KMAXB)
      REAL SUM5(IMAXB,JMAXB,KMAXB),SUM6(IMAXB,JMAXB,KMAXB)
      REAL SUM7(IMAXB,JMAXB,KMAXB),SUM8(IMAXB,JMAXB,KMAXB)
      REAL SUM9(IMAXB,JMAXB,KMAXB),SUM10(IMAXB,JMAXB,KMAXB)
      REAL SUM11(IMAXB,JMAXB,KMAXB),SUM12(IMAXB,JMAXB,KMAXB)
      REAL SUM13(IMAXB,JMAXB,KMAXB)
      REAL UGUESS(IMAXB,JMAXB,KMAXB),VGUESS(IMAXB,JMAXB,KMAXB)
      REAL WGUESS(IMAXB,JMAXB,KMAXB),RHOK(100),WGT(200)
      REAL DISTXYZ(IMAXB,JMAXB,KMAXB,2,2,2)
      INTEGER INDEX1(IMAXB,JMAXB,KMAXB)
      INTEGER NS(3)
      REAL XSHIFT(80),YSHIFT(80)
      FLAG=-1.0E+10
      MDIM=IMAXB*JMAXB*KMAXB
      WRITE(6,'(A2)')CHARV
      WRITE(6,'(A2)')CHARDBZ
      WRITE(6,*)'MDIM = IMAXB*JMAXB*KMAXB = ',MDIM
      PI=2.*ASIN(1.)
      PI2=2.*PI
      DTR=PI/180.
      ISCANDIM=1024
      IBINS=7
      SXB2=.5*SXB*SXB
      SYB2=.5*SYB*SYB
      SZB2=.5*SZB*SZB
      JBINS=3*IBINS+1
      ZZB=0.
      RDEL=0.
      DO K=1,KMAXB
       DO J=1,JMAXB
        DO I=1,IMAXB
         DO KK=1,2
          DO JJ=1,2
           DO II=1,2
            DISTXYZ(I,J,K,II,JJ,KK)=1000.
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      DO I=1,JBINS
       RATIOI=(FLOAT(I)/FLOAT(IBINS))**2.
       WGT(I)=EXP(-RATIOI)
      ENDDO
      DO K=1,KMAXB
       HEIGHT=ZZB+(K-1)*SZB
       RHOK(K)=EXP(-HEIGHT*.10437052)
      ENDDO
      ZTOP=(KMAXB-1)*SZB
      RLZ=ZTOP
      RLH=25.*SXB
      UBAR=10.
      DO I=1,1024
       A=I
       RANGE(I)=A
       VR(I)=-A
       DBZ(I)=A*2.
      ENDDO
      write(6,*)'entered abcd'
      NS(1)=0
      NS(2)=0
      NS(3)=1
      NSTEP=2
      NSTEP2=1
      MDIM2=IMAXB*JMAXB
      WRITE(6,*)'IMAXB,JMAXB,KMAXB,SXB,SYB,SZB,ROTB'
      WRITE(6,*)IMAXB,JMAXB,KMAXB,SXB,SYB,SZB,ROTB
      WRITE(6,*)'Entered ABCD'
      RHWEIGHT=EHMAX*EHMAX
      RZWEIGHT=EZMAX*EZMAX
      RXMAX=RHMAX
      RYMAX=RHMAX
      WRITE(6,*)'EHMAX,EZMAX,RXMAX,RYMAX,RZMAX'
      WRITE(6,*)EHMAX,EZMAX,RXMAX,RYMAX,RZMAX    
      WRITE(6,*)'XACC = ',XACC
      WRITE(6,*)'YACC = ',YACC
      WRITE(6,*)'ZACC = ',ZACC
      TIMER=0.
c      WRITE(1,*)'AFTER ENTERING ABCD, FACTZ,FACTR = ',FACTZ,FACTR

      IF(IREADABCD.NE.1)THEN
       DO KK=KMIN2,KMAX2
        DO JJ=1,JMAXB
         DO II=1,IMAXB
C
C First call to ABCD:  zero all angle-projection summations
C
           SUMDBZ(II,JJ,KK)=0
           SUMWTSAVE(II,JJ,KK)=0
           SUM1(II,JJ,KK)=0
           SUM2(II,JJ,KK)=0
           SUM3(II,JJ,KK)=0
           SUM4(II,JJ,KK)=0
           SUM5(II,JJ,KK)=0
           SUM6(II,JJ,KK)=0
           SUM7(II,JJ,KK)=0
           SUM8(II,JJ,KK)=0
           SUM9(II,JJ,KK)=0
           SUM10(II,JJ,KK)=0
           SUM11(II,JJ,KK)=0
           SUM12(II,JJ,KK)=0
           SUM13(II,JJ,KK)=0
           XWSAVE(II,JJ,KK)=0
           YWSAVE(II,JJ,KK)=0
           ZWSAVE(II,JJ,KK)=0
           INDEX1(II,JJ,KK)=0
C
C If w found during first iteration, set to zero to begin second iteration
C
         ENDDO
        ENDDO
       ENDDO
      ENDIF

      RLATBEFORE=0.
      IBEGIN=0
      IF(IREADABCD.EQ.1)GO TO 72735
c      write(6,*)'calling readuf_disc'
92923 RLATBEFORE=RLAT
      CALL READUF_DISC(TIMES,NUMSWEEP,FLIGHTID,NRECORD,RLAT,RLON,
     + RALT,AZIM,ELEV,RANGE,DBZ,VR,ISCANDIM,IERR,CHARV,CHARDBZ)
c      WRITE(6,*)'VE = ',(VR(I),I=1,ISCANDIM)
c      WRITE(6,*)'DZ = ',(DBZ(I),I=1,ISCANDIM)
      IF(RLAT.EQ.0.)RLAT=RLATBEFORE
       write(6,*)'times,rlat,rlon from reader = ',times,rlat,rlon
       IF(IRAMSFILES.EQ.1)THEN
        CALL READRAMFILE(TIMES,TIME_RAM,RLAT,RLON,DM1,DM2,DM3,
     *                      DM4,DM5,DM6,DM7,IOS)
       ENDIF
      IF(IBEGIN.EQ.1.AND.TIMES.LT.TSTART)TIMES=TIMES+86400.
      IF(TIMES.GT.TEND.OR.IERR.EQ.-1)GO TO 92929
      IF(TIMES.GT.TIMEOLD)WRITE(6,*)TIMES,RALT
      TIMEOLD=TIMES
      IF(IREADABCD.NE.1.AND.TIMES.GE.TSTART.AND.TIMES.LE.TEND)THEN
       JJJSTART=1
       IF(RDEL.LT.0)THEN
        DO JJJ=JJJSTART,ISCANDIM
         RANGE(JJJ)=RANGE(JJJ)+RDEL
        ENDDO
       ENDIF
       RHIGH=300.
       RLOW=.1
       WRITE(6,*)'RANGE(1) = ',RANGE(1)
c      WRITE(6,*)TIMES,RLAT,RLON,ELEV,AZIM,JJJSTART,RANGE(JJJSTART),
c     +RLOW,RHIGH
C
C Angle projection summing follows.  
C
       IBEGIN=1
       XWEIGHT=0
       YWEIGHT=0
       ZWEIGHT=0
       SUMWEIGHT=0
c       ISCANDIM=512
c       DO I=1,ISCANDIM
c        DBZ(I)=0.
c        VR(I)=0.
c        RANGE(I)=I*.150
c       ENDDO
       CALL VRLOCATE(OLAT,OLON,TIMES,
     + RLAT,RLON,RALT,ELEV,AZIM,RANGE,XR,YR,ZR,
     + SMOTIONU,SMOTIONV,ISCANDIM,CENTIME,IBINS,WGT,VR)
       AZIMM=AZIM-ROTB
       GAMMA=SIN(ELEV*DTR)
       ALPHA=SIN(AZIMM*DTR)*COS(ELEV*DTR)
       BETA=COS(AZIMM*DTR)*COS(ELEV*DTR)
c       WRITE(6,*)'RLAT,RLON = ',RLAT,RLON
c       WRITE(6,*)'ZZZB = ',ZZZB
c       WRITE(6,*)'TIMES,RANGE(1),XR(1),YR(1),ZR(1) = ',
c     +  TIMES,RANGE(1),XR(1),YR(1),ZR(1)
       DO I=JJJSTART,ISCANDIM
        ZR(I)=ZR(I)-ZZZB
       ENDDO
c       WRITE(6,*)'RLAT,RLON = ',RLAT,RLON
c       WRITE(6,*)'TIMES,XR(1),YR(1),ZR(1) = ',TIMES,XR(1),YR(1),ZR(1)
c       do i=1,iscandim
c        write(6,*)i,xr(i),yr(i),zr(i)
c       enddo
       DO I=JJJSTART,ISCANDIM
        DBZCHK=DBZTEST+20.*ALOG10(RANGE(I)/RRRR)
c        WRITE(1,*)'RANGE(I),DBZTEST,RRRR,DBZCHK,DBZ(I) = '
c        WRITE(1,*)RANGE(I),DBZTEST,RRRR,DBZCHK,DBZ(I)
c        write(6,*)'i,range,dbz,vr = ',range(i),dbz(i),vr(i)
        RHO=EXP(-ZR(I)*.10437052)
c        RHO=EXP(-ZR(I)*.11)
c        ZRCHECK=RANGE(I)/57.*.95*COS(ELEV*DTR)
        ZRCHECK=0.
c        DISTCHECK=RANGE(I)*.08
C        RLOW=0.1
c        write(6,*)'rlow,rhigh,ralt,i,range(i) = ',
c     +             rlow,rhigh,ralt,i,range(i)
        XR(I)=XR(I)+XSHIFT(NF)
        YR(I)=YR(I)+YSHIFT(NF)
c        WRITE(6,*)'VR,DBZ,ZR,RANGE,RLOW,RHIGH = ',
c     + VR(I),DBZ(I),ZR(I),RANGE(I),RLOW,RHIGH
c        X=XR(I)
c        Y=YR(I)
c        Z=ZR(I)
c        UAN=UBAR*SIN(PI2*X/RLH)*COS(PI2*Y/RLH)*COS(PI2*Z/RLZ)
c        VAN=UBAR*COS(PI2*X/RLH)*SIN(PI2*Y/RLH)*COS(PI2*Z/RLZ)
c        WAN=-2.*UBAR*RLZ/RLH*COS(PI2*X/RLH)*
c     +                     COS(PI2*Y/RLH)*SIN(PI2*Z/RLZ)
c        WRITE(6,*)'I,RLZ,RLH,PI2,UBAR = ',I,RLZ,RLH,PI2,UBAR
c        WRITE(6,*)'I,X,Y,Z,UAN,VAN,WAN = ',I,X,Y,Z,UAN,VAN,WAN
c        UAN=UAN/RHO
c        VAN=VAN/RHO
c        WAN=WAN/RHO
c        UAN=UAN
c        VAN=VAN
c        WAN=WAN
c        VRVR=UAN*ALPHA+VAN*BETA+WAN*GAMMA
        VRVR=VR(I)
C        DBZ(I)=0.
        IF(VR(I).GT.-900..OR.DBZ(I).GT.DBZCHK.AND.
     +     ZR(I).GT.ZRCHECK.AND.ZR(I).GT..075.
C     +      ZR(I).GT.-4.
     +     AND.RANGE(I).GT.RLOW.AND.RANGE(I).LT.RHIGH)THEN
c         write(6,*)'I GOT HERE'
         IF(GAMMA.LE.SINTOL)THEN
C
C Determine actual distance from location of datum in data file
C to the point being evaluated in output file 
C
          CALL DIFFOBLONG(XR(I),YR(I),ZR(I),
     +     RXMAX,RYMAX,RZMAX,SXB,SYB,SZB,
     +     XZB,YZB,ZZB,IMAXB,JMAXB,KMAXB,
     +     IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX,ROTB)
          DO KK=IZMIN,IZMAX
           RHOVRVR=VRVR*RHOK(KK)
           DO JJ=IYMIN,IYMAX
            DO II=IXMIN,IXMAX  
             CALL DIFFDIST(XR(I),YR(I),ZR(I),
     +       II,JJ,KK,INEAR,JNEAR,KNEAR,
     +       SXB,SYB,SZB,XZB,YZB,ZZB,RX,RY,RZ,ROTB)
c             WRITE(6,*)'II,JJ,KK = ',II,JJ,KK
c             WRITE(6,*)'RX,RY,RZ = ',RX,RY,RZ
             XX=(FLOAT(II-1)+.5)*SXB-XZB
             YY=(FLOAT(JJ-1)+.5)*SYB-YZB
             ZZ=(KK-1)*SZB+ZZB
             IF(RX.LE.0)IA=1
             IF(RX.GT.0)IA=2
             IF(RY.LE.0)JA=1
             IF(RY.GT.0)JA=2
             IF(RZ.LE.0)KA=1
             IF(RZ.GT.0)KA=2
             RR=SQRT(RX*RX+RY*RY+RZ*RZ)
             IF(RR.LT.DISTXYZ(II,JJ,KK,IA,JA,KA))THEN
              DISTXYZ(II,JJ,KK,IA,JA,KA)=RR
             ENDIF
             ABC=(RX*RX+RY*RY)/RHWEIGHT+RZ*RZ/RZWEIGHT
             WEIGHT=EXP(-ABC)
             P=WEIGHT
C
C Accumulate weighted sums of x, y, and z distances of input data
C
             SUMWEIGHT=P
             IF(DBZ(I).GT.-100)THEN
              DBZZ=DBZ(I)*SLOPE1
              ELSE
              DBZZ=-200.
             ENDIF
C
C Perform DBZ weighted sum
C
             ZE=10.**(DBZZ/10.)
             IF(ABS(VR(I)).LT.900.)THEN
              VTEST=1.
             ELSE
              VTEST=0.
             ENDIF
             IF(DBZ(I).GT.DBZCHK)THEN
              SUMDBZ(II,JJ,KK)=
     +          SUMDBZ(II,JJ,KK)+P*ZE
              SUMWTSAVE(II,JJ,KK)=SUMWTSAVE(II,JJ,KK)+SUMWEIGHT
c              IF(KK.EQ.5)THEN
c               WRITE(6,*)'II,JJ,KK,I,RANGE(I),ELEV,XR,YR,ZR,DBZ = ',
c     + II,JJ,KK,I,RANGE(I),ELEV,XR(I),YR(I),ZR(I),DBZ(I)
c              ENDIF
             ENDIF
C
C Perform weighted sums of 9 direction cosine/Vr terms
C
c             IF(VTEST.GT.0.)THEN
c              IF(KK.EQ.1)THEN
c               WRITE(6,*)'II,JJ,KK,I,RANGE(I),ELEV,XR,YR,ZR,VR = ',
c     + II,JJ,KK,I,RANGE(I),ELEV,XR(I),YR(I),ZR(I),VR(I)
c              ENDIF
c             ENDIF
             SUM1(II,JJ,KK)=
     +         SUM1(II,JJ,KK)+P*ALPHA*BETA*VTEST
             SUM2(II,JJ,KK)=
     +         SUM2(II,JJ,KK)+P*BETA*RHOVRVR*VTEST
             SUM3(II,JJ,KK)=
     +         SUM3(II,JJ,KK)+P*BETA*BETA*VTEST
             SUM4(II,JJ,KK)=
     +         SUM4(II,JJ,KK)+P*ALPHA*RHOVRVR*VTEST
             SUM5(II,JJ,KK)=
     +         SUM5(II,JJ,KK)+P*ALPHA*ALPHA*VTEST
             SUM6(II,JJ,KK)=
     +         SUM6(II,JJ,KK)+P*ALPHA*GAMMA*VTEST
             SUM7(II,JJ,KK)=
     +         SUM7(II,JJ,KK)+P*BETA*GAMMA*VTEST
             SUM8(II,JJ,KK)=
     +         SUM8(II,JJ,KK)+P*GAMMA*GAMMA*VTEST
             SUM9(II,JJ,KK)=
     +         SUM9(II,JJ,KK)+P*GAMMA*RHOVRVR*VTEST
             SUM10(II,JJ,KK)=
     +         SUM10(II,JJ,KK)+P*ALPHA*VTEST
             SUM11(II,JJ,KK)=
     +         SUM11(II,JJ,KK)+P*BETA*VTEST
             SUM12(II,JJ,KK)=
     +         SUM12(II,JJ,KK)+P*GAMMA*VTEST
             SUM13(II,JJ,KK)=
     +         SUM13(II,JJ,KK)+P*RHOVRVR*VTEST
             XWSAVE(II,JJ,KK)=XWSAVE(II,JJ,KK)+P*RX*VTEST
             YWSAVE(II,JJ,KK)=YWSAVE(II,JJ,KK)+P*RY*VTEST
             ZWSAVE(II,JJ,KK)=ZWSAVE(II,JJ,KK)+P*RZ*VTEST
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDIF
       ENDDO
      ENDIF
      GO TO 92923
92929 WRITE(6,*)'JUST FINISHED VR SUMMATION'
c5554  OPEN(103,FILE='/hrd/dat/gustav/tracking/insideabcd')
c      WRITE(103,*)'AT STATEMENT 5554--NF=',NF
c      WRITE(103,*)'JUST FINISHED VR SUMMATION'
c      CLOSE(103)
      DO KK=1,KMAXB
       DO JJ=1,JMAXB
        DO II=1,IMAXB
         ITEST=0
         IF(KK.LE.3)THEN
          KABOT=2
         ELSE
          KABOT=1
         ENDIF
         DO KA=KABOT,2
          DO JA=1,2
           DO IA=1,2
C            IF(KK.LE.2.AND.KA.EQ.1)THEN
C             DIST=DISTXYZ(II,JJ,KK,IA,JA,2)
C             IF(DIST.GT.EZMAX)DIST=1000.
C            ELSE
             DIST=DISTXYZ(II,JJ,KK,IA,JA,KA)
C            ENDIF
            IF(DIST.GE.999.)THEN
             SUM1(II,JJ,KK)=0
             SUM2(II,JJ,KK)=0
             SUM3(II,JJ,KK)=0
             SUM4(II,JJ,KK)=0
             SUM5(II,JJ,KK)=0
             SUM6(II,JJ,KK)=0
             SUM7(II,JJ,KK)=0
             SUM8(II,JJ,KK)=0
             SUM9(II,JJ,KK)=0
             SUM10(II,JJ,KK)=0
             SUM11(II,JJ,KK)=0
             SUM12(II,JJ,KK)=0
             SUM13(II,JJ,KK)=0
             XWSAVE(II,JJ,KK)=1000.*SUMWTSAVE(II,JJ,KK)
             YWSAVE(II,JJ,KK)=1000.*SUMWTSAVE(II,JJ,KK)
             ZWSAVE(II,JJ,KK)=1000.*SUMWTSAVE(II,JJ,KK)
             GO TO 10097
            ENDIF
           ENDDO
          ENDDO
         ENDDO
10097   ENDDO
       ENDDO
      ENDDO
      IOKAY=1
72735 IF(IOKAY.EQ.1)THEN
C
C At end of each complete iteration
C
       WRITE(6,*)'IN ABCD'
       WRITE(6,*)'SUM1,SUM2,SUM3 = ',
     + SUM1(20,20,4),SUM2(20,20,4),SUM3(20,20,4),
     + SUMDBZ(20,20,4),SUMWTSAVE(20,20,4)
       DO KK=KMIN2,KMAX2
c        OPEN(102,FILE='/hrd/dat/gustav/tracking/abcd')
c        WRITE(102,*)'IN ABCD--K,NF = ',KK,NF
c        CLOSE(102)
        DO JJ=1,JMAXB
         DO II=1,IMAXB
          SUMWEIGHT=SUM3(II,JJ,KK)+SUM5(II,JJ,KK)+SUM8(II,JJ,KK)
          IF(SUMWEIGHT.GT.0..AND.SUMWTSAVE(II,JJ,KK).GT.0)THEN
           SUM1(II,JJ,KK)=SUM1(II,JJ,KK)/SUMWEIGHT
           SUM2(II,JJ,KK)=SUM2(II,JJ,KK)/SUMWEIGHT
           SUM3(II,JJ,KK)=SUM3(II,JJ,KK)/SUMWEIGHT
           SUM4(II,JJ,KK)=SUM4(II,JJ,KK)/SUMWEIGHT
           SUM5(II,JJ,KK)=SUM5(II,JJ,KK)/SUMWEIGHT
           SUM6(II,JJ,KK)=SUM6(II,JJ,KK)/SUMWEIGHT
           SUM7(II,JJ,KK)=SUM7(II,JJ,KK)/SUMWEIGHT
           SUM8(II,JJ,KK)=SUM8(II,JJ,KK)/SUMWEIGHT
           SUM9(II,JJ,KK)=SUM9(II,JJ,KK)/SUMWEIGHT
           SUM10(II,JJ,KK)=SUM10(II,JJ,KK)/SUMWEIGHT
           SUM11(II,JJ,KK)=SUM11(II,JJ,KK)/SUMWEIGHT
           SUM12(II,JJ,KK)=SUM12(II,JJ,KK)/SUMWEIGHT
           SUM13(II,JJ,KK)=SUM13(II,JJ,KK)/SUMWEIGHT
           XWSAVE(II,JJ,KK)=XWSAVE(II,JJ,KK)/SUMWEIGHT
           YWSAVE(II,JJ,KK)=YWSAVE(II,JJ,KK)/SUMWEIGHT
           ZWSAVE(II,JJ,KK)=ZWSAVE(II,JJ,KK)/SUMWEIGHT
           ZE=SUMDBZ(II,JJ,KK)/SUMWTSAVE(II,JJ,KK)
           SUMDBZ(II,JJ,KK)=10.*ALOG10(ZE)
          ELSE
           SUMDBZ(II,JJ,KK)=FLAG
          ENDIF
         ENDDO
        ENDDO
       ENDDO      
       WRITE(6,*)'SUM1,SUM2,SUM3,SUMDBZ,SUMWTSAVE = ',
     + SUM1(20,20,4),SUM2(20,20,4),SUM3(20,20,4),
     + SUMDBZ(20,20,4),SUMWTSAVE(20,20,4)
      ENDIF
      RETURN
      END
      SUBROUTINE VRLOCATE(OLAT,OLON,TIMES,
     + RLAT,RLON,ALT,ELEV,AZIM,RANGE,XR,YR,ZR,
     + SMOTIONU,SMOTIONV,ISCANDIM,CENTIME,IBINS,WEIGHT,VR)
      REAL RANGE(ISCANDIM),XR(ISCANDIM),YR(ISCANDIM),ZR(ISCANDIM)
      REAL SUM(1024),WSUM(1024),WEIGHT(200),VR(1024)
      RALT=ALT/1000.
      PI=ASIN(1.0)*2.
      DTR=PI/180.
      REARTH=6366.
      TDIFF=TIMES-CENTIME
      JBINS=IBINS*3
c      DO I=1,ISCANDIM
c       SUM(I)=0
c       WSUM(I)=0
c       ISTART=I-JBINS
c       IEND=I+JBINS
c       IF(ISTART.LT.1)ISTART=1
c       IF(IEND.GT.ISCANDIM)IEND=ISCANDIM
c       DO J=ISTART,IEND
c        IF(VR(J).GT.-900.)THEN
c         IDIFF=J-I
c         IDIFF=IABS(IDIFF)
c         WSUM(I)=WSUM(I)+WEIGHT(IDIFF)
c         SUM(I)=SUM(I)+WEIGHT(IDIFF)*VR(J)
c        ENDIF
c       ENDDO
c      ENDDO
c      DO I=1,ISCANDIM
c       IF(WSUM(I).LE.0..AND.VR(I).LE.-900.)THEN
c        VR(I)=-999.9
c       ELSE
c        VR(I)=SUM(I)/WSUM(I)
c       ENDIF
c      ENDDO
      XSHIFT=TDIFF*SMOTIONU/1000.
      YSHIFT=TDIFF*SMOTIONV/1000.     
      RAZIM=AZIM*DTR
      RELEV=ELEV*DTR
      DO I=1,ISCANDIM
c xradar=horizontal distance from aircraft
c yradar=vertical distance from center of earth
       XRADAR=RANGE(I)*COS(RELEV)
       YRADAR=REARTH+RALT+RANGE(I)*SIN(RELEV)
c radnew=distance of datum from center of earth
       RADNEW=SQRT(XRADAR*XRADAR+YRADAR*YRADAR)
c delangle=angle between vertical at aircraft and vertical at datum
       DELANGLE=ATAN(XRADAR/YRADAR)
c distance along great circle from radar to datum
       DISTANCE=.5*(RADNEW+REARTH+RALT)*DELANGLE
c height=altitude of datum
       HEIGHT=RADNEW-REARTH
c distance2=xradar
       DISTANCE2=RANGE(I)*COS(RELEV)
c height2=altitude of datum not accounting for earth's curvature
       HEIGHT2=RALT+RANGE(I)*SIN(RELEV)
c delx=distance (aircraft/datum) along great circle due to east west
       DELX=DISTANCE*SIN(RAZIM)
c dely=distance(aircraft/datum)along great circle due to north-south
       DELY=DISTANCE*COS(RAZIM)
c dx=east west distance from origin to aircraft
       XLAT=.5*(RLAT+OLAT)*DTR
       DX=(RLON-OLON)*DTR*REARTH*COS(XLAT)
c dy=north south distance from origin to aircraft
       DY=(RLAT-OLAT)*DTR*REARTH
c xshift=east-west distance from lower left corner to origin
c yshift=north-south distance from lower left corner to origin
c XR=distance of datum from lower left corner
c YR=distance of datum from lower left corner
       XR(I)=DELX+DX-XSHIFT
c yshift=distance from lower left corner to origin
       YR(I)=DELY+DY-YSHIFT
       ZR(I)=HEIGHT
c       WRITE(1,*)'I,RANGE,DISTANCE,HEIGHT,DELX,DELY,DX,DY,XSHIFT,'
c     +  'YSHIFT = ',
c     +  I,RANGE(I),DISTANCE,HEIGHT,DELX,DELY,DX,DY,XSHIFT,YSHIFT
      ENDDO
      CONTINUE
      RETURN
      END

