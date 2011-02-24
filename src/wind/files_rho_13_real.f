      SUBROUTINE WRPLN(LU,IMAX,JMAX,KMAX,U,V,W,DIV,SUMDBZ,
     + IUV,ZZ,SZ)
      INTEGER*2 ID(5,512)
      REAL RD(5,512)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),W(IMAX,JMAX,KMAX)
      REAL DIV(IMAX,JMAX,KMAX),SUMDBZ(IMAX,JMAX,KMAX)

      FLAG=-1.0E+10
c      WRITE(1,*)'ENTERED WRPLN'
      ITOTALTEST=21*IMAX*JMAX+48*IMAX
      DO K=1,KMAX
       H=ZZ+(K-1)*SZ
       RHO=EXP(-H*.10437052)
c       write(1,*)'k,zz,sz,h,rho = ',k,zz,sz,h,rho
       WRITE(6,*)'Writing plane ',K
       DO J=1,JMAX
c        WRITE(6,*)'K = ',K,'J = ',J
        DO I=1,IMAX
         IF(U(I,J,K).GT.FLAG.AND.V(I,J,K).GT.FLAG.AND.
     +     (U(I,J,K).NE.0..OR.V(I,J,K).NE.0.))THEN
          IF(IUV.EQ.0)THEN
           UU=U(I,J,K)/RHO
           VV=V(I,J,K)/RHO
           WDR=AMOD(630.0-57.29578*ATAN2(VV,UU),360.0)
           WSP=SQRT(UU*UU+VV*VV)
           ID(1,I)=WDR*10.
           ID(2,I)=WSP*10.
          ELSE
           RD(1,I)=U(I,J,K)/RHO
           RD(2,I)=V(I,J,K)/RHO
          ENDIF
         ELSE
          ID(1,I)=-32000
          ID(2,I)=-32000
          RD(1,I)=FLAG
          RD(2,I)=FLAG
         ENDIF
         IF(IUV.EQ.1)THEN
          IF(W(I,J,K).GT.FLAG)THEN
           RD(3,I)=W(I,J,K)/RHO
          ELSE
           RD(3,I)=FLAG
          ENDIF
         ELSE
          IF(W(I,J,K).GT.FLAG)THEN
           ID(3,I)=W(I,J,K)*100./RHO
          ELSE
           ID(3,I)=-32000
          ENDIF
         ENDIF
         IF(IUV.EQ.1)THEN
          IF(DIV(I,J,K).GT.FLAG)THEN
           RD(5,I)=DIV(I,J,K)/RHO
          ELSE
           RD(5,I)=FLAG
          ENDIF
         ELSE
          IF(ABS(DIV(I,J,K)).LT..32766)THEN
           ID(5,I)=DIV(I,J,K)*1.0E+05/RHO
          ELSE
           ID(5,I)=32767
          ENDIF
         ENDIF
         IF(SUMDBZ(I,J,K).LT.-100.)THEN
           ID(4,I)=-32000
         ELSE
           DBZ=SUMDBZ(I,J,K)
           ID(4,I)=NINT(DBZ*10.)
         ENDIF
        ENDDO
        IF(IUV.NE.1)THEN
         WRITE(LU)((ID(N,I),N=1,5),I=1,IMAX)
        ELSE
         WRITE(LU)((RD(N,I),N=1,3),ID(4,I),RD(5,I),I=1,IMAX)
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END



      SUBROUTINE WRPLN4(LU,IMAX,JMAX,KMAX,U,V,W,DIV,SUMDBZ,
     + IUV,ZZ,SZ,INDEX4)
      INTEGER*2 ID(5,512)
      REAL RD(5,512)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),W(IMAX,JMAX,KMAX)
      REAL DIV(IMAX,JMAX,KMAX),SUMDBZ(IMAX,JMAX,KMAX)
      INTEGER INDEX4(IMAX,JMAX,KMAX)

      FLAG=-1.0E+10
c      WRITE(1,*)'ENTERED WRPLN'
      ITOTALTEST=21*IMAX*JMAX+48*IMAX
      DO K=1,KMAX
       H=ZZ+(K-1)*SZ
       RHO=EXP(-H*.10437052)
       WRITE(6,*)'Writing plane ',K
       DO J=1,JMAX
C        WRITE(6,*)'K = ',K,'J = ',J
        DO I=1,IMAX
         IF(U(I,J,K).GT.FLAG.AND.V(I,J,K).GT.FLAG.AND.
     +     (U(I,J,K).NE.0..OR.V(I,J,K).NE.0.).AND.
     +      INDEX4(I,J,K).EQ.1)THEN
          IF(IUV.EQ.0)THEN
           UU=U(I,J,K)/RHO
           VV=V(I,J,K)/RHO
           WDR=AMOD(630.0-57.29578*ATAN2(VV,UU),360.0)
           WSP=SQRT(UU*UU+VV*VV)
           ID(1,I)=WDR*10.
           ID(2,I)=WSP*10.
          ELSE
           RD(1,I)=U(I,J,K)/RHO
           RD(2,I)=V(I,J,K)/RHO
          ENDIF
         ELSE
          ID(1,I)=-32000
          ID(2,I)=-32000
          RD(1,I)=FLAG
          RD(2,I)=FLAG
         ENDIF
         IF(IUV.EQ.1)THEN
          IF(W(I,J,K).GT.FLAG.AND.INDEX4(I,J,K).EQ.1)THEN
           RD(3,I)=W(I,J,K)/RHO
          ELSE
           RD(3,I)=FLAG
          ENDIF
         ELSE
          IF(W(I,J,K).GT.FLAG.AND.INDEX4(I,J,K).EQ.1)THEN
           ID(3,I)=W(I,J,K)*100./RHO
          ELSE
           ID(3,I)=-32000
          ENDIF
         ENDIF
         IF(IUV.EQ.1)THEN
          IF(DIV(I,J,K).GT.FLAG.AND.INDEX4(I,J,K).EQ.1)THEN
           RD(5,I)=DIV(I,J,K)/RHO
          ELSE
           RD(5,I)=FLAG
          ENDIF
         ELSE
          IF(ABS(DIV(I,J,K)).LT..32766.AND.INDEX4(I,J,K).EQ.1)THEN
           ID(5,I)=DIV(I,J,K)*1.0E+05/RHO
          ELSE
           ID(5,I)=32767
          ENDIF
         ENDIF
         IF(SUMDBZ(I,J,K).LT.-100.)THEN
           ID(4,I)=-32000
         ELSE
           DBZ=SUMDBZ(I,J,K)
           ID(4,I)=NINT(DBZ*10.)
         ENDIF
        ENDDO
        IF(IUV.NE.1)THEN
         WRITE(LU)((ID(N,I),N=1,5),I=1,IMAX)
        ELSE
         WRITE(LU)((RD(N,I),N=1,3),ID(4,I),RD(5,I),I=1,IMAX)
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE READHEADER(LU,FILNAM,KEYWORD,FLTNAME,
     + STMNAME,RADAR,EXPERIMENT,
     + CREATIME,EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     + SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     + POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)
      CHARACTER KEYWORD*4,FLTNAME*8,STMNAME*12,RADAR*4,EXPERIMENT*32
      CHARACTER CREATIME*32,EXTRA1*28,FILNAM*80
      INTEGER IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3
      REAL STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7

      OPEN(LU,FILE=FILNAM,STATUS='OLD',FORM='UNFORMATTED')
      READ(LU)KEYWORD,FLTNAME,STMNAME,
     + RADAR,EXPERIMENT,CREATIME,
     + EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,
     + EXTRA3,STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7

      RETURN
      END


      SUBROUTINE WRITEHEADER(LU,FILNAM,KEYWORD,FLTNAME,
     + STMNAME,RADAR,EXPERIMENT,
     + CREATIME,EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3,STIME,ETIME,OLAT,OLON,
     + SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,CO2,AZMCOR,ELCOR,THRESH,
     + POWERT,BIEL,AZBIEL,ETIME1,STIME2,EXTRA6,EXTRA7)
      CHARACTER KEYWORD*4,FLTNAME*8,STMNAME*12,RADAR*4,EXPERIMENT*32
      CHARACTER CREATIME*32,EXTRA1*28,FILNAM*80
      INTEGER IMAX,JMAX,KMAX,KOUNT,NMOSM,
     + IUNFLD,IATTEN,FLAG,EXTRA2,EXTRA3
      REAL STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7

      write(6,*)'WRITEHEADER ENTERED'
      OPEN(LU,FILE=FILNAM,FORM='UNFORMATTED')
      write(6,*)'OPENED LU ',LU,' AS '
      write(6,'(A)')FILNAM
      KEYWORD='WIND'
      CREATIME='                                '
      EXTRA1='                            '
      WRITE(LU)KEYWORD,FLTNAME,STMNAME,RADAR,EXPERIMENT,CREATIME,
     + EXTRA1,IMAX,JMAX,KMAX,KOUNT,NMOSM,IUNFLD,IATTEN,FLAG,EXTRA2,
     + EXTRA3,STIME,ETIME,OLAT,OLON,SX,SY,SZ,XZ,YZ,ZZ,ROT,RA,CO1,
     + CO2,AZMCOR,ELCOR,THRESH,POWERT,BIEL,AZBIEL,ETIME1,STIME2,
     + EXTRA6,EXTRA7

      RETURN
      END



      SUBROUTINE WRITEABCD(LUW,IMAX,JMAX,KMAX,SUM1,SUM2,SUM3,
     +   SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,SUM10,SUM11,SUM12,SUM13,SUMDBZ,
     +   XWEIGHT,YWEIGHT,ZWEIGHT)
      REAL SUM1(IMAX,JMAX,KMAX),SUM2(IMAX,JMAX,KMAX)
      REAL SUM3(IMAX,JMAX,KMAX),SUM4(IMAX,JMAX,KMAX)
      REAL SUM5(IMAX,JMAX,KMAX),SUM6(IMAX,JMAX,KMAX)
      REAL SUM7(IMAX,JMAX,KMAX),SUM8(IMAX,JMAX,KMAX)
      REAL SUM10(IMAX,JMAX,KMAX),SUM11(IMAX,JMAX,KMAX)
      REAL SUM12(IMAX,JMAX,KMAX),SUM13(IMAX,JMAX,KMAX)
      REAL SUM9(IMAX,JMAX,KMAX)
      REAL XWEIGHT(IMAX,JMAX,KMAX)
      REAL YWEIGHT(IMAX,JMAX,KMAX),ZWEIGHT(IMAX,JMAX,KMAX)
      REAL SUMDBZ(IMAX,JMAX,KMAX)

      WRITE (LUW)IMAX,JMAX,KMAX
      DO K=1,KMAX
       DO J=1,JMAX
        WRITE(LUW)(SUM1(I,J,K),SUM2(I,J,K),SUM3(I,J,K),
     +             SUM4(I,J,K),SUM5(I,J,K),SUM6(I,J,K),
     +             SUM7(I,J,K),SUM8(I,J,K),SUM9(I,J,K),
     +             SUM10(I,J,K),SUM11(I,J,K),SUM12(I,J,K),SUM13(I,J,K),
     +              SUMDBZ(I,J,K),XWEIGHT(I,J,K),YWEIGHT(I,J,K),
     +              ZWEIGHT(I,J,K),I=1,IMAX)
       ENDDO
      ENDDO
      write(6,*)'about to close sumfile'
      CLOSE(LUW)
      write(6,*)'closing sumfile and returning'

      RETURN
      END


      SUBROUTINE READABCD(LUW,IMAX,JMAX,KMAX,SUM1,SUM2,SUM3,SUM4,SUM5,
     +   SUM6,SUM7,SUM8,SUM9,SUM10,SUM11,SUM12,SUM13,
     +   SUMDBZ,KMAX2,XWEIGHT,YWEIGHT,ZWEIGHT)
      REAL SUM1(IMAX,JMAX,KMAX),SUM2(IMAX,JMAX,KMAX)
      REAL SUM3(IMAX,JMAX,KMAX),SUM4(IMAX,JMAX,KMAX)
      REAL SUM5(IMAX,JMAX,KMAX),SUM6(IMAX,JMAX,KMAX)
      REAL SUM7(IMAX,JMAX,KMAX),SUM8(IMAX,JMAX,KMAX)
      REAL SUM10(IMAX,JMAX,KMAX),SUM11(IMAX,JMAX,KMAX)
      REAL SUM12(IMAX,JMAX,KMAX),SUM13(IMAX,JMAX,KMAX)
      REAL SUM9(IMAX,JMAX,KMAX)
      REAL XWEIGHT(IMAX,JMAX,KMAX)
      REAL YWEIGHT(IMAX,JMAX,KMAX),ZWEIGHT(IMAX,JMAX,KMAX)
      REAL SUMDBZ(IMAX,JMAX,KMAX)
      REAL SUMM1(512),SUMM2(512),SUMM3(512)
      REAL SUMM4(512),SUMM5(512),SUMM6(512)
      REAL SUMM7(512),SUMM8(512),SUMM9(512)
      REAL SUMM10(512),SUMM11(512),SUMM12(512),SUMM13(512)
      REAL SUMDBZZ(512),XW(512),YW(512),ZW(512)
      INTEGER*2 ITEST(512,512,80)

      READ(LUW)IMAXTEST,JMAXTEST,KMAXTEST
      IF(IMAXTEST.GT.512.OR.JMAXTEST.GT.512.OR.KMAXTEST.GT.80)THEN
       WRITE(6,*)'ITEST FAILED',IMAXTEST,JMAXTEST,KMAXTEST
       STOP
      ENDIF
      write(6,*)'entered readabcd'
      IF(IMAX.NE.IMAXTEST.OR.JMAX.NE.JMAXTEST.OR.KMAX.NE.KMAXTEST)THEN
       WRITE(6,*)'SUMFILE ERROR'
       WRITE(6,*)'IMAX,JMAX,KMAX = ',IMAX,JMAX,KMAX
       WRITE(6,*)'IMAXTEST,JMAXTEST,KMAXTEST ',
     +            IMAXTEST,JMAXTEST,KMAXTEST
       IF(KMAX.NE.KMAXTEST.AND.IMAX.EQ.IMAXTEST.AND.JMAX.EQ.JMAXTEST)
     + THEN
        WRITE(6,*)'SUMFILE ERROR'
        WRITE(6,*)'IMAX,JMAX,KMAX = ',IMAX,JMAX,KMAX
        WRITE(6,*)'IMAXTEST,JMAXTEST,KMAXTEST ',
     +            IMAXTEST,JMAXTEST,KMAXTEST
        ICONTINUE=1
        IF(ICONTINUE.EQ.0)STOP
       ELSE
        STOP
       ENDIF
      ENDIF
      DO K=1,KMAX2
       DO J=1,JMAX
c        write(6,*)'j,k,imax,jmax,kmax2 = ',j,k,imax,jmax,kmax2
        READ(LUW)(SUMM1(I),SUMM2(I),SUMM3(I),
     +   SUMM4(I),SUMM5(I),SUMM6(I),SUMM7(I),SUMM8(I),SUMM9(I),
     +   SUMM10(I),SUMM11(I),SUMM12(I),SUMM13(I),
     +              SUMDBZZ(I),XW(I),YW(I),ZW(I),I=1,IMAX)
        DO I=1,IMAX
         SUM1(I,J,K)=SUMM1(I)
         SUM2(I,J,K)=SUMM2(I)
         SUM3(I,J,K)=SUMM3(I)
         SUM4(I,J,K)=SUMM4(I)
         SUM5(I,J,K)=SUMM5(I)
         SUM6(I,J,K)=SUMM6(I)
         SUM7(I,J,K)=SUMM7(I)
         SUM8(I,J,K)=SUMM8(I)
         SUM9(I,J,K)=SUMM9(I)
         SUM10(I,J,K)=SUMM10(I)
         SUM11(I,J,K)=SUMM11(I)
         SUM12(I,J,K)=SUMM12(I)
         SUM13(I,J,K)=SUMM13(I)
         SUMDBZ(I,J,K)=SUMDBZZ(I)
         XWEIGHT(I,J,K)=XW(I)
         YWEIGHT(I,J,K)=YW(I)
         ZWEIGHT(I,J,K)=ZW(I)
c         write(6,*)'i,j,k,sum1,sum2,sum3,sumdbz = '
c       write(6,*)i,j,k,sum1(i,j,k),sum2(i,j,k),sum3(i,j,k),sumdbz(i,j,k)
c         write(6,*)'sum4,sum5,sum6,sum7,sum8,sum9 = '
c        write(6,*)sum4(i,j,k),sum5(i,j,k),sum6(i,j,k),sum7(i,j,k),
c     +            sum8(i,j,k),sum9(i,j,k)
c        write(6,*)'xw,yw,zw = ',xweight(i,j,k),yweight(i,j,k),
c     + zweight(i,j,k)
        ENDDO
       ENDDO
      ENDDO
      CLOSE(LUW)
      DO K=1,KMAX
       DO J=1,JMAX
        DO I=1,IMAX
         SUMTEST=SUM3(I,J,K)+SUM5(I,J,K)+SUM8(I,J,K)
         IF(SUMTEST.GT.0.)THEN
          ITEST(I,J,K)=1
         ELSE
          ITEST(I,J,K)=0
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      NFILLED=0
      DO K=1,KMAX
       K1=K-1
       IF(K1.LE.0)K1=1
       K3=K+1
       IF(K3.GT.KMAX)K3=KMAX
       DO J=1,JMAX
        J1=J-1
        IF(J1.LE.0)J1=1
        J3=J+1
        IF(J3.GT.JMAX)J3=JMAX
        DO I=1,IMAX
         IF(ITEST(I,J,K).EQ.0)THEN
          DIVIDE=0.
          SUMM1(1)=0.
          SUMM2(1)=0.
          SUMM3(1)=0.
          SUMM4(1)=0.
          SUMM5(1)=0.
          SUMM6(1)=0.
          SUMM7(1)=0.
          SUMM8(1)=0.
          SUMM9(1)=0.
          SUMM10(1)=0.
          SUMM11(1)=0.
          SUMM12(1)=0.
          SUMM13(1)=0.
c          SUMDBZZ(1)=0.
          XW(1)=0.
          YW(1)=0.
          ZW(1)=0.
          I1=I-1
          IF(I1.LE.0)I1=1
          I3=I+1
          IF(I3.GT.IMAX)I3=IMAX
          IF(K1.NE.K.AND.K3.NE.K)THEN
           IF(ITEST(I,J,K1).EQ.1.AND.ITEST(I,J,K3).EQ.1)THEN
            DIVIDE=DIVIDE+2.
            SUMM1(1)=SUM1(I,J,K1)+SUM1(I,J,K3)
            SUMM2(1)=SUM2(I,J,K1)+SUM2(I,J,K3)
            SUMM3(1)=SUM3(I,J,K1)+SUM3(I,J,K3)
            SUMM4(1)=SUM4(I,J,K1)+SUM4(I,J,K3)
            SUMM5(1)=SUM5(I,J,K1)+SUM5(I,J,K3)
            SUMM6(1)=SUM6(I,J,K1)+SUM6(I,J,K3)
            SUMM7(1)=SUM7(I,J,K1)+SUM7(I,J,K3)
            SUMM8(1)=SUM8(I,J,K1)+SUM8(I,J,K3)
            SUMM9(1)=SUM9(I,J,K1)+SUM9(I,J,K3)
            SUMM10(1)=SUM10(I,J,K1)+SUM10(I,J,K3)
            SUMM11(1)=SUM11(I,J,K1)+SUM11(I,J,K3)
            SUMM12(1)=SUM12(I,J,K1)+SUM12(I,J,K3)
            SUMM13(1)=SUM13(I,J,K1)+SUM13(I,J,K3)
            XW(1)=XWEIGHT(I,J,K1)+XWEIGHT(I,J,K3)
            YW(1)=YWEIGHT(I,J,K1)+YWEIGHT(I,J,K3)
            ZW(1)=ZWEIGHT(I,J,K1)+ZWEIGHT(I,J,K3)
           ENDIF
          ENDIF
          IF(J1.NE.J.AND.J3.NE.J)THEN
           IF(ITEST(I,J1,K).EQ.1.AND.ITEST(I,J3,K).EQ.1)THEN
            DIVIDE=DIVIDE+2.
            SUMM1(1)=SUMM1(1)+SUM1(I,J1,K)+SUM1(I,J3,K)
            SUMM2(1)=SUMM2(1)+SUM2(I,J1,K)+SUM2(I,J3,K)
            SUMM3(1)=SUMM3(1)+SUM3(I,J1,K)+SUM3(I,J3,K)
            SUMM4(1)=SUMM4(1)+SUM4(I,J1,K)+SUM4(I,J3,K)
            SUMM5(1)=SUMM5(1)+SUM5(I,J1,K)+SUM5(I,J3,K)
            SUMM6(1)=SUMM6(1)+SUM6(I,J1,K)+SUM6(I,J3,K)
            SUMM7(1)=SUMM7(1)+SUM7(I,J1,K)+SUM7(I,J3,K)
            SUMM8(1)=SUMM8(1)+SUM8(I,J1,K)+SUM8(I,J3,K)
            SUMM9(1)=SUMM9(1)+SUM9(I,J1,K)+SUM9(I,J3,K)
            SUMM10(1)=SUMM10(1)+SUM10(I,J1,K)+SUM10(I,J3,K)
            SUMM11(1)=SUMM11(1)+SUM11(I,J1,K)+SUM11(I,J3,K)
            SUMM12(1)=SUMM12(1)+SUM12(I,J1,K)+SUM12(I,J3,K)
            SUMM13(1)=SUMM13(1)+SUM13(I,J1,K)+SUM13(I,J3,K)
            XW(1)=XW(1)+XWEIGHT(I,J1,K)+XWEIGHT(I,J3,K)
            YW(1)=YW(1)+YWEIGHT(I,J1,K)+YWEIGHT(I,J3,K)
            ZW(1)=ZW(1)+ZWEIGHT(I,J1,K)+ZWEIGHT(I,J3,K)
           ENDIF
          ENDIF
          IF(I1.NE.I.AND.I3.NE.I)THEN
           IF(ITEST(I1,J,K).EQ.1.AND.ITEST(I3,J,K).EQ.1)THEN
            DIVIDE=DIVIDE+2.
            SUMM1(1)=SUMM1(1)+SUM1(I1,J,K)+SUM1(I3,J,K)
            SUMM2(1)=SUMM2(1)+SUM2(I1,J,K)+SUM2(I3,J,K)
            SUMM3(1)=SUMM3(1)+SUM3(I1,J,K)+SUM3(I3,J,K)
            SUMM4(1)=SUMM4(1)+SUM4(I1,J,K)+SUM4(I3,J,K)
            SUMM5(1)=SUMM5(1)+SUM5(I1,J,K)+SUM5(I3,J,K)
            SUMM6(1)=SUMM6(1)+SUM6(I1,J,K)+SUM6(I3,J,K)
            SUMM7(1)=SUMM7(1)+SUM7(I1,J,K)+SUM7(I3,J,K)
            SUMM8(1)=SUMM8(1)+SUM8(I1,J,K)+SUM8(I3,J,K)
            SUMM9(1)=SUMM9(1)+SUM9(I1,J,K)+SUM9(I3,J,K)
            SUMM10(1)=SUMM10(1)+SUM10(I1,J,K)+SUM10(I3,J,K)
            SUMM11(1)=SUMM11(1)+SUM11(I1,J,K)+SUM11(I3,J,K)
            SUMM12(1)=SUMM12(1)+SUM12(I1,J,K)+SUM12(I3,J,K)
            SUMM13(1)=SUMM13(1)+SUM13(I1,J,K)+SUM13(I3,J,K)
            XW(1)=XWEIGHT(I1,J,K)+XWEIGHT(I3,J,K)
            YW(1)=YWEIGHT(I1,J,K)+YWEIGHT(I3,J,K)
            ZW(1)=ZWEIGHT(I1,J,K)+ZWEIGHT(I3,J,K)
           ENDIF
          ENDIF
          IF(DIVIDE.GE.3.999)THEN
           SUM1(I,J,K)=SUMM1(1)/DIVIDE
           SUM2(I,J,K)=SUMM2(1)/DIVIDE
           SUM3(I,J,K)=SUMM3(1)/DIVIDE
           SUM4(I,J,K)=SUMM4(1)/DIVIDE
           SUM5(I,J,K)=SUMM5(1)/DIVIDE
           SUM6(I,J,K)=SUMM6(1)/DIVIDE
           SUM7(I,J,K)=SUMM7(1)/DIVIDE
           SUM8(I,J,K)=SUMM8(1)/DIVIDE
           SUM9(I,J,K)=SUMM9(1)/DIVIDE
           SUM10(I,J,K)=SUMM10(1)/DIVIDE
           SUM11(I,J,K)=SUMM11(1)/DIVIDE
           SUM12(I,J,K)=SUMM12(1)/DIVIDE
           SUM13(I,J,K)=SUMM13(1)/DIVIDE
c           SUMDBZ(I,J,K)=SUMDBZZ(1)/DIVIDE
           XWEIGHT(I,J,K)=XW(1)/DIVIDE
           YWEIGHT(I,J,K)=YW(1)/DIVIDE
           ZWEIGHT(I,J,K)=ZW(1)/DIVIDE
           NFILLED=NFILLED+1
          ENDIF
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      WRITE(6,*)'number of points filled is ',nfilled
      write(6,*)'leaving readabcd'

      RETURN
      END

      SUBROUTINE LDPLNW(LU,U,V,W,IMAX,JMAX,KMAX,IUV,ZZ,SZ)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),W(IMAX,JMAX,KMAX)
      INTEGER*2 IU(512),IV(512),IW(512),IDB(512),IDV(512)
      REAL RU(512),RV(512),RW(512),RDV(512)

      FLAG=-1.0E+10
      WRITE(6,*)'INSIDE LDPLNW'
      WRITE(6,*)'IMAX,JMAX,KMAX,IUV = ',IMAX,JMAX,KMAX,IUV
      DO K=1,KMAX
       H=ZZ+(K-1)*SZ
       RHO=EXP(-H*.10437052)
       DO J=1,JMAX
        IF(IUV.EQ.0)THEN
         READ(LU)(IU(I),IV(I),IW(I),IDB(I),IDV(I),I=1,IMAX)
        ELSE
         READ(LU)(RU(I),RV(I),RW(I),IDB(I),RDV(I),I=1,IMAX)
        ENDIF
        DO I=1,IMAX
         IF(IUV.EQ.0)THEN
          IF(IU(I).LT.0.)THEN
           U(I,J,K)=FLAG
           V(I,J,K)=FLAG
          ELSE
           WDR=IU(I)*.1
           WSP=IV(I)*.1
           XU=-SIN(WDR*3.14159/180.)*WSP
           XV=-COS(WDR*3.14159/180.)*WSP
           U(I,J,K)=XU*RHO
           V(I,J,K)=XV*RHO
          ENDIF
          IF(IW(I).GT.-31999)THEN
           W(I,J,K)=IW(I)*.01*RHO
          ELSE
           W(I,J,K)=FLAG
          ENDIF
         ELSE
          IF(RU(I).GT.FLAG.AND.RV(I).GT.FLAG.AND.RW(I).GT.FLAG)THEN
           U(I,J,K)=RU(I)*RHO
           V(I,J,K)=RV(I)*RHO
           W(I,J,K)=RW(I)*RHO
          ELSE
           U(I,J,K)=FLAG
           V(I,J,K)=FLAG
           W(I,J,K)=FLAG
          ENDIF
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END

