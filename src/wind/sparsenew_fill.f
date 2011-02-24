       SUBROUTINE SPARSE(M,RHS,U,V,IMAX,JMAX,KMAX,SX,SY,SZ,N,NN,INDEX1,
     + INDEX2,IPOS,SUM1,SUM2,SUM3,SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,
     + SUM10,SUM11,SUM12,SUM13,
     + XWEIGHT,YWEIGHT,ZWEIGHT,XACC,YACC,ZACC,
     + OSUMFILENAME,NZ,NFILE,IWORK,NW,INDEX4,
     + FACTZZ,FACTRH,FACTRZ,ICAST,W,BOUNDUV,WBOUND,Z0,NNNN,
     + SUMDBZ,SUMWTSAVE,DBZTEST,NORM,NNNNN,INDXN,WTOP,IBOTTOP,IA,JA,
     + IGUESSER,UGUESS,VGUESS,WGUESS,BAND1,HB,DPB,IRSW,ZLOW,ZHIGH,
     + FACTWP,WPRECIP,FACTVR,FACTVR2,FACTWTOP,IGUESS,
     + FACTWBOT,DIV,KT,KMAX2,ITOP,IVTSUB,INDEX5,ZTOPPER,
     + WTBACKU,WTBACKV,WTBACKW,TESTRMS,TESTGAMMA)
      INTEGER*2 ICAST(IMAX,JMAX,KMAX)
      INTEGER IVTSUB(80)
      INTEGER KT(IMAX,JMAX)
      INTEGER KTTOP(IMAX,JMAX),ITOP(IMAX,JMAX)
      REAL TESTRMS(IMAX,JMAX,KMAX),TESTGAMMA(IMAX,JMAX,KMAX)
      REAL SUM1(IMAX,JMAX,KMAX),SUM2(IMAX,JMAX,KMAX)
      REAL SUM3(IMAX,JMAX,KMAX),SUM4(IMAX,JMAX,KMAX)
      REAL SUM5(IMAX,JMAX,KMAX),SUM6(IMAX,JMAX,KMAX)
      REAL SUM7(IMAX,JMAX,KMAX),SUM8(IMAX,JMAX,KMAX)
      REAL SUM10(IMAX,JMAX,KMAX),SUM11(IMAX,JMAX,KMAX)
      REAL SUM12(IMAX,JMAX,KMAX),SUM13(IMAX,JMAX,KMAX)
      REAL SUM9(IMAX,JMAX,KMAX)
      REAL XWEIGHT(IMAX,JMAX,KMAX)
      REAL RHOK(80),VGUESS(IMAX,JMAX,KMAX),WGUESS(IMAX,JMAX,KMAX)
      REAL YWEIGHT(IMAX,JMAX,KMAX),ZWEIGHT(IMAX,JMAX,KMAX)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),UGUESS(IMAX,JMAX,KMAX)
      REAL SUMDBZ(IMAX,JMAX,KMAX)
      REAL SUMWTSAVE(IMAX,JMAX,KMAX)
      REAL DIV(IMAX,JMAX,KMAX)
      INTEGER INDEX4(IMAX,JMAX,KMAX),INDEX5(IMAX,JMAX,KMAX)
      DOUBLE PRECISION M(NZ),SM,SX,SY,SZ,DIVSUM,DIVSUM2
      INTEGER IA(NNNN),JA(NZ),IWORK(NZ)
      REAL BAND1(NNNN)
      CHARACTER*80 OSUMFILENAME(80)
      DOUBLE PRECISION DIVDZ1,DIVERGENCE,DZ1,DY1,DX1
      DOUBLE PRECISION VIWMINUS2,VIWMINUS1,VIWMIDDLE,VIWPLUS1,VIWPLUS2
      DOUBLE PRECISION VIVMINUS2,VIVMINUS1,VIVMIDDLE,VIVPLUS1,VIVPLUS2
      DOUBLE PRECISION VIUMINUS2,VIUMINUS1,VIUMIDDLE,VIUPLUS1,VIUPLUS2
      DOUBLE PRECISION V2IWMINUS2,V2IWMINUS1
      DOUBLE PRECISION V2IWMIDDLE,V2IWPLUS1,V2IWPLUS2
      DOUBLE PRECISION V2IVMINUS2,V2IVMINUS1
      DOUBLE PRECISION V2IVMIDDLE,V2IVPLUS1,V2IVPLUS2
      DOUBLE PRECISION V2IUMINUS2,V2IUMINUS1
      DOUBLE PRECISION V2IUMIDDLE,V2IUPLUS1,V2IUPLUS2
      DOUBLE PRECISION UTEST(9),FILTU,FILTV,FILTW
      DOUBLE PRECISION RHS(NNNN)
      REAL W(IMAX,JMAX,KMAX)
      INTEGER INDEX1(NN),INDEX2(NN),IDIF(115)
      INTEGER INDEXNN(100),INDEXSAVE(20)      
      DOUBLE PRECISION VALUES(100),VALUESAVE(20)
      write(6,*)'entered sparse imax,jmax,kmax,n,nn,nnnn = ',
     1 imax,jmax,kmax,n,nn,nnnn
      FLAG=-1.0E+10
      ALFA=-.0001043705
      BOUND=0.
      FACTORWP=10.
      FACTRH10=0.
      FACTRH10=FACTRH
      FACTRH=0.
      HIGH=HB+.5*DPB
      HLOW=HB-.5*DPB
      IJMAX=IMAX*JMAX
      WRITE(6,*)'POSITION 1'
      GAMMAMAX=0.0
      write(6,*)'nz = ',nz
      DO K=1,100
       INDEXNN(K)=0
      ENDDO
      write(6,*)'entered sparse imax,jmax,kmax,n,nn,nnnn = ',
     1 imax,jmax,kmax,n,nn,nnnn
      write(6,*)'entered sparse imax,jmax,kmax,n,nn,nnnn = ',
     1 imax,jmax,kmax,n,nn,nnnn
      DO K=1,KMAX
       HEIGHT=Z0+(K-1)*SZ
       RHO=EXP(-HEIGHT*.10437052)
       RHOK(K)=RHO
      ENDDO
      write(6,*)'entered sparse imax,jmax,kmax,n,nn,nnnn = ',
     1 imax,jmax,kmax,n,nn,nnnn
      DO I=1,NNNN
       BAND1(I)=0.
       RHS(I)=0.
      ENDDO
      write(6,*)'entered sparse imax,jmax,kmax,n,nn,nnnn = ',
     1 imax,jmax,kmax,n,nn,nnnn
      DO I=1,NN
       INDEX1(I)=0
       INDEX2(I)=0
      ENDDO
      write(6,*)'entered sparse imax,jmax,kmax,n,nn,nnnn = ',
     1 imax,jmax,kmax,n,nn,nnnn
      WRITE(6,*)'AFTER ENTERING SPARSE, FACTZ,FACTRH = ',FACTZZ,FACTRH
      WRITE(6,*)'SUM1,SUM2,SUM3 = ',
     + SUM1(15,15,2),SUM2(15,15,2),SUM3(15,15,2)
      JPOSITIONU=15+14*IMAX+IMAX*JMAX
      JPOSITIONV=JPOSITIONU+N
      JPOSITIONW=JPOSITIONV+N
      WRITE(6,*)'n,imax,jmax,POSITIONS FOR TEST POINT ARE ',
     + n,imax,jmax,JPOSITIONU,JPOSITIONV,JPOSITIONW
      DX1=SX*1000.
      DY1=SY*1000.
      DZ1=SZ*1000.
      DX2=DX1*DX1
      DY2=DY1*DY1
      DZ2=DZ1*DZ1
      DX4=DX2*DX2
      DY4=DY2*DY2
      DZ4=DZ2*DZ2
      XFACTOR=1.
      YFACTOR=1.
c      WRITE(1,*)'DX1,DY1,DZ1,DX2,DY2,DZ2,DX4,DY4,DZ4'
c      WRITE(1,*)DX1,DY1,DZ1,DX2,DY2,DZ2,DX4,DY4,DZ4
      CENTERPROD=2.*(DX2+DY2)/(DX2*DY2)
      IDIF(1)=0
      IDIF(2)=1
      IDIF(3)=2
      IDIF(4)=IMAX-1
      IDIF(5)=IMAX
      IDIF(6)=IMAX+1
      IDIF(7)=2*IMAX
      IDIF(8)=IMAX*JMAX-2*IMAX
      IDIF(9)=IMAX*JMAX-IMAX
      IDIF(10)=IMAX*JMAX-2
      IDIF(11)=IMAX*JMAX-1
      IDIF(12)=IMAX*JMAX
      IDIF(13)=IMAX*JMAX+1
      IDIF(14)=IMAX*JMAX+2
      IDIF(15)=IMAX*JMAX+IMAX
      IDIF(16)=IMAX*JMAX+2*IMAX
      IDIF(17)=2*IMAX*JMAX
      IDIF(18)=N-IMAX*JMAX-2*IMAX-2
      IDIF(19)=N-IMAX*JMAX-2*IMAX-1
      IDIF(20)=N-IMAX*JMAX-2*IMAX
      IDIF(21)=N-IMAX*JMAX-2*IMAX+1
      IDIF(22)=N-IMAX*JMAX-2*IMAX+2
      IDIF(23)=N-IMAX*JMAX-IMAX-2
      IDIF(24)=N-IMAX*JMAX-IMAX-1
      IDIF(25)=N-IMAX*JMAX-IMAX
      IDIF(26)=N-IMAX*JMAX-IMAX+1
      IDIF(27)=N-IMAX*JMAX-IMAX+2
      IDIF(28)=N-IMAX*JMAX-2
      IDIF(29)=N-IMAX*JMAX-1
      IDIF(30)=N-IMAX*JMAX
      IDIF(31)=N-IMAX*JMAX+1
      IDIF(32)=N-IMAX*JMAX+2
      IDIF(33)=N-IMAX*JMAX+IMAX-2
      IDIF(34)=N-IMAX*JMAX+IMAX-1
      IDIF(35)=N-IMAX*JMAX+IMAX
      IDIF(36)=N-IMAX*JMAX+IMAX+1
      IDIF(37)=N-IMAX*JMAX+IMAX+2
      IDIF(38)=N-IMAX*JMAX+2*IMAX-2
      IDIF(39)=N-IMAX*JMAX+2*IMAX-1
      IDIF(40)=N-IMAX*JMAX+2*IMAX
      IDIF(41)=N-IMAX*JMAX+2*IMAX+1
      IDIF(42)=N-IMAX*JMAX+2*IMAX+2
      IDIF(43)=N-2*IMAX-2
      IDIF(44)=N-2*IMAX-1
      IDIF(45)=N-2*IMAX
      IDIF(46)=N-2*IMAX+1
      IDIF(47)=N-2*IMAX+2
      IDIF(48)=N-IMAX-2
      IDIF(49)=N-IMAX-1
      IDIF(50)=N-IMAX
      IDIF(51)=N-IMAX+1
      IDIF(52)=N-IMAX+2
      IDIF(53)=N-1
      IDIF(54)=N
      IDIF(55)=N+1
      IDIF(56)=N+IMAX-2
      IDIF(57)=N+IMAX-1
      IDIF(58)=N+IMAX
      IDIF(59)=N+IMAX+1
      IDIF(60)=N+IMAX+2
      IDIF(61)=N+2*IMAX-2
      IDIF(62)=N+2*IMAX-1
      IDIF(63)=N+2*IMAX
      IDIF(64)=N+2*IMAX+1
      IDIF(65)=N+2*IMAX+2
      IDIF(66)=N+IMAX*JMAX-2*IMAX-2
      IDIF(67)=N+IMAX*JMAX-2*IMAX-1
      IDIF(68)=N+IMAX*JMAX-2*IMAX
      IDIF(69)=N+IMAX*JMAX-2*IMAX+1
      IDIF(70)=N+IMAX*JMAX-2*IMAX+2
      IDIF(71)=N+IMAX*JMAX-IMAX-2
      IDIF(72)=N+IMAX*JMAX-IMAX-1
      IDIF(73)=N+IMAX*JMAX-IMAX
      IDIF(74)=N+IMAX*JMAX-IMAX+1
      IDIF(75)=N+IMAX*JMAX-IMAX+2
      IDIF(76)=N+IMAX*JMAX-2
      IDIF(77)=N+IMAX*JMAX-1
      IDIF(78)=N+IMAX*JMAX
      IDIF(79)=N+IMAX*JMAX+1
      IDIF(80)=N+IMAX*JMAX+2
      IDIF(81)=N+IMAX*JMAX+IMAX-2
      IDIF(82)=N+IMAX*JMAX+IMAX-1
      IDIF(83)=N+IMAX*JMAX+IMAX
      IDIF(84)=N+IMAX*JMAX+IMAX+1
      IDIF(85)=N+IMAX*JMAX+IMAX+2
      IDIF(86)=N+IMAX*JMAX+2*IMAX-2
      IDIF(87)=N+IMAX*JMAX+2*IMAX-1
      IDIF(88)=N+IMAX*JMAX+2*IMAX
      IDIF(89)=N+IMAX*JMAX+2*IMAX+1
      IDIF(90)=N+IMAX*JMAX+2*IMAX+2
      IDIF(91)=2*N-IMAX*JMAX-2
      IDIF(92)=2*N-IMAX*JMAX-1
      IDIF(93)=2*N-IMAX*JMAX
      IDIF(94)=2*N-IMAX*JMAX+1
      IDIF(95)=2*N-IMAX*JMAX+2
      IDIF(96)=2*N-2
      IDIF(97)=2*N-1
      IDIF(98)=2*N
      IDIF(99)=2*N+1
      IDIF(100)=2*N+2
      IDIF(101)=2*N+IMAX*JMAX-2
      IDIF(102)=2*N+IMAX*JMAX-1
      IDIF(103)=2*N+IMAX*JMAX
      IDIF(104)=2*N+IMAX*JMAX+1
      IDIF(105)=2*N+IMAX*JMAX+2
      IDIF(106)=N-2
      IDIF(107)=N+2
      IDIF(108)=IMAX*JMAX-4*IMAX
      IDIF(109)=IMAX*JMAX-3*IMAX
      IDIF(110)=IMAX*JMAX-4
      IDIF(111)=IMAX*JMAX-3
      IDIF(112)=IMAX*JMAX+3
      IDIF(113)=IMAX*JMAX+4
      IDIF(114)=IMAX*JMAX+3*IMAX
      IDIF(115)=IMAX*JMAX+4*IMAX
      WRITE(6,*)'POSITION 2'
      DO K=1,KMAX
       DO J=1,JMAX
c        WRITE(6,*)'K,J = ',K,J,'  INDEX4 FOLLOWS'
c        WRITE(6,17248)(INDEX4(I,J,K),I=1,IMAX)
c17248   FORMAT(35I2)
       ENDDO
      ENDDO
      DO J=1,JMAX
       DO I=1,IMAX
        KTTOP(I,J)=-1
       ENDDO
      ENDDO
      DO J=1,JMAX
       DO I=1,IMAX
        DO K=1,KMAX
         KK=K
         IF(TESTRMS(I,J,K).GT.0.)GO TO 7653
        ENDDO
7653    TESTRMS(I,J,1)=TESTRMS(I,J,KK)
        TESTGAMMA(I,J,1)=TESTGAMMA(I,J,KK)
        DO K=KMAX,1,-1
         KK=K
         IF(TESTRMS(I,J,K).GT.0.)GO TO 7654
        ENDDO
7654    TESTRMS(I,J,KMAX)=TESTRMS(I,J,KK)
        TESTGAMMA(I,J,KMAX)=TESTGAMMA(I,J,KK)
       ENDDO
      ENDDO
      DO K=1,KMAX
       DO J=1,JMAX
        DO I=1,IMAX
         IF(INDEX4(I,J,K).EQ.1)THEN
           U(I,J,K)=0.
           V(I,J,K)=0.
           W(I,J,K)=0.
           INDEX5(I,J,K)=1
c           write(6,*)'i,j,k,index4 = ',i,j,k,index4(i,j,k)
         ELSE
           U(I,J,K)=FLAG
           INDEX5(I,J,K)=0
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      DO K=1,KMAX
       DO J=1,JMAX
        DO I=1,IMAX
         KT(I,J)=0
         KTTOP(I,J)=-1
        ENDDO
       ENDDO
      ENDDO
      DO J=1,JMAX
       DO I=1,IMAX
        DO K=1,KMAX
         IF(INDEX4(I,J,K).GT.0)KT(I,J)=K
        ENDDO
       ENDDO
      ENDDO
      KTADD=FACTVR/DZ1
      WRITE(6,*)'KTADD = ',KTADD
      DO J=1,JMAX
       DO I=1,IMAX
        IF(KT(I,J).GT.0)KTTOP(I,J)=KT(I,J)+KTADD
        IF(KTTOP(I,J).GT.KMAX)KTTOP(I,J)=KMAX
c        write(6,*)'i,j,kt,kttop = ',i,j,kt(i,j),kttop(i,j)
       ENDDO
      ENDDO
8787  DO K=1,KMAX
       DO J=1,JMAX
        DO I=1,IMAX
c         write(6,*)'i,j,k,u(i,j,k) = ',i,j,k,u(i,j,k)
         IPOSITION=I+(J-1)*IMAX+(K-1)*IJMAX
         IPOSITIONU=IPOSITION
         IPOSITIONV=IPOSITION+N
         IPOSITIONW=IPOSITIONV+N
         U(I,J,K)=0
         IF(U(I,J,K).GT.FLAG)THEN
          INDEX2(IPOSITIONU)=1
          INDEX2(IPOSITIONV)=1
          INDEX2(IPOSITIONW)=1
         ELSE
          INDEX2(IPOSITIONU)=0
          INDEX2(IPOSITIONV)=0
          INDEX2(IPOSITIONW)=0
         ENDIF
        ENDDO
       ENDDO
      ENDDO       
      IPOS=0
      DO I=1,NN
c       write(6,*)'i,index2(i) = ',i,index2(i)
       IF(INDEX2(I).NE.0)THEN
        IPOS=IPOS+1 
c        write(6,*)'ipos = ',ipos
        IF(IPOS.GT.NNNN)THEN
         WRITE(6,*)'I= ',I,'  IPOS = ',IPOS
         STOP
        ENDIF
        INDEX1(IPOS)=I
        INDEX2(I)=IPOS
       ENDIF
      ENDDO 
      ZEROWRITE=0.
      MODE=0
      LEVEL=0
      NOUT=6
      write(6,*)'calling sbini with '
      write(6,*)'ipos,nz = ',ipos,nz
      write(6,*)'index1(1)= ',index1(1)
      WRITE(6,*)'SBINI--IPOS,NZ = ',
     + IPOS,NZ
      CALL SBINI(IPOS,NZ,IA,JA,M,IWORK)
1995  LEVEL=0
      MODE=1
      FACTOR=1.
c      WRITE(6,*)'U = ',U(15,15,3),U(15,14,4),U(14,15,4),
c     +                 U(16,15,4),U(15,16,4),U(15,15,5)
      MODE=1
      IF(FACTZZ.LE.0..AND.FACTWP.LE.0)GO TO 77778
      DIVSUM=0.
      DIVSUM2=0.
      NDIVSUM=0.
      DO K=1,KMAX
       WRITE(6,*)'SUMMING FOR DIV AND VR AT LEVEL ',K
       HEIGHT=(Z0+(K-1)*SZ)*1000.
       RHO=RHOK(K)
c       WRITE(6,*)'HEIGHT,Z0,SZ,K = ',HEIGHT,Z0,SZ,K
       write(6,*)'height, rho = ',height,rho
c       WRITE(6,*)'RHO = ',RHO
       FACTPOINT=1.
       DO J=1,JMAX
        DO I=1,IMAX 
C
C Determine data positions for 6 divergence data points
C
         DIV(I,J,K)=FLAG
         IF(U(I,J,K).LE.FLAG)GO TO 9998
         CALL BOUND5(I,J,K,IMAX,JMAX,KMAX,U,IBOUNDX,IBOUNDY,IBOUNDZ)
         IPOSITION=I+(J-1)*IMAX+(K-1)*IJMAX
         IPOSITIONU=IPOSITION
         IPOSITIONV=IPOSITION+N
         IPOSITIONW=IPOSITIONV+N
         INPOSITIONU=INDEX2(IPOSITIONU)
         INPOSITIONV=INDEX2(IPOSITIONV)
         INPOSITIONW=INDEX2(IPOSITIONW)
         BAND1(INPOSITIONU)=1
         BAND1(INPOSITIONV)=1
         BAND1(INPOSITIONW)=1
c         FACTZ=FACTZZ/RHO/RHO
         FACTZ=FACTZZ
         IF(IBOUNDX.GT.-99.AND.IBOUNDY.GT.-99.AND.
     +      IBOUNDZ.GT.-99.)THEN
          ICAST(I,J,K)=0
         ELSE
          ICAST(I,J,K)=1
         ENDIF
         CALL VALUESET(IPOSITION,IPOSITIONU,IPOSITIONV,IPOSITIONW,
     +    IMAX,JMAX,KMAX,IBOUNDX,IBOUNDY,IBOUNDZ,DX1,DY1,DZ1,
     +    IUPLUS2,IUPLUS1,IUMIDDLE,IUMINUS1,IUMINUS2,
     +    VIUPLUS2,VIUPLUS1,VIUMIDDLE,VIUMINUS1,VIUMINUS2,
     +    IVPLUS2,IVPLUS1,IVMIDDLE,IVMINUS1,IVMINUS2,
     +    VIVPLUS2,VIVPLUS1,VIVMIDDLE,VIVMINUS1,VIVMINUS2,
     +    IWPLUS2,IWPLUS1,IWMIDDLE,IWMINUS1,IWMINUS2,
     +    VIWPLUS2,VIWPLUS1,VIWMIDDLE,VIWMINUS1,VIWMINUS2,
     +    I,J,K,IP2,IP1,IP0,IM1,IM2,JP2,JP1,JP0,JM1,JM2,
     +    KP2,KP1,KP0,KM1,KM2)
C         write(6,*)'i,j,k,iboundx,iboundy,iboundz = ',
C     +    i,j,k,iboundx,iboundy,iboundz
C         write(6,*)'iuplus2,iuplus1,iumiddle,iuminus1,iuminus2 = ',
C     +    iuplus2,iuplus1,iumiddle,iuminus1,iuminus2
C         write(6,*)'ivplus2,ivplus1,ivmiddle,ivminus1,ivminus2 = ',
C     +    iuplus2,ivplus1,ivmiddle,ivminus1,ivminus2
C         write(6,*)'iwplus2,iwplus1,iwmiddle,iwminus1,iwminus2 = ',
C     +    iwplus2,iwplus1,iwmiddle,iwminus1,iwminus2
C         write(6,*)
C     +  'viuplus2,viuplus1,viumiddle,viuminus1,viuminus2 = ',
C     +    viuplus2,viuplus1,viumiddle,viuminus1,viuminus2
C         write(6,*)
C     +   'vivplus2,vivplus1,vivmiddle,vivminus1,vivminus2 = ',
C     +    vivplus2,vivplus1,vivmiddle,vivminus1,vivminus2
C         write(6,*)
C     +   'viwplus2,viwplus1,viwmiddle,viwminus1,viwminus2 = ',
C     +    viwplus2,viwplus1,viwmiddle,viwminus1,viwminus2
         NINDEX=15
         DIVERGENCE=0.
         DIVERGENCEH=0.
         IF(KM2.GT.0)THEN
          DIVERGENCE=DIVERGENCE+(WGUESS(I,J,KP2)*VIWPLUS2
     +                         +WGUESS(I,J,KP1)*VIWPLUS1
     +                         +WGUESS(I,J,KP0)*VIWMIDDLE
     +                         +WGUESS(I,J,KM1)*VIWMINUS1
     +                         +WGUESS(I,J,KM2)*VIWMINUS2)/DZ1
         ELSE
          DIVERGENCE=DIVERGENCE+(WGUESS(I,J,KP2)*VIWPLUS2
     +                         +WGUESS(I,J,KP1)*VIWPLUS1
     +                         +WGUESS(I,J,KP0)*VIWMIDDLE
     +                         +WGUESS(I,J,KM1)*VIWMINUS1)/DZ1
         ENDIF
         DIVERGENCEH=DIVERGENCEH+(VGUESS(I,JP2,K)*VIVPLUS2
     +                         +VGUESS(I,JP1,K)*VIVPLUS1
     +                         +VGUESS(I,JP0,K)*VIVMIDDLE
     +                         +VGUESS(I,JM1,K)*VIVMINUS1
     +                         +VGUESS(I,JM2,K)*VIVMINUS2)/DZ1
         DIVERGENCEH=DIVERGENCEH+(UGUESS(IP2,J,K)*VIUPLUS2
     +                         +UGUESS(IP1,J,K)*VIUPLUS1
     +                         +UGUESS(IP0,J,K)*VIUMIDDLE
     +                         +UGUESS(IM1,J,K)*VIUMINUS1
     +                         +UGUESS(IM2,J,K)*VIUMINUS2)/DZ1
         DIVERGENCE=DIVERGENCE+DIVERGENCEH
         DIVDZ1=DIVERGENCE*DZ1
         DIVDZ1=0.
c         WRITE(6,*)'I,J,K,DIVERGENCE = ',I,J,K,DIVERGENCE
         if(k.gt.100)then
          write(6,*)'v = ',vguess(i,jm2,k),vguess(i,jm1,k),
     + vguess(i,jp0,k),vguess(i,jp1,k),vguess(i,jp2,k)
          write(6,*)'vval = ',vivminus2,vivminus1,vivmiddle,
     +    vivplus1,vivplus2
          write(6,*)'u = ',uguess(im2,j,k),uguess(im1,j,k),
     + uguess(ip0,j,k)+uguess(ip1,j,k)+uguess(ip2,j,k)
          write(6,*)'uval = ',viuminus2,viuminus1,viumiddle,
     +    viuplus1,viuplus2
         endif
         DIV(I,J,K)=DIVERGENCEH
         DIVSUM=DIVSUM+DIVERGENCEH
         DIVSUM2=DIVSUM2+DIVERGENCEH*DIVERGENCEH
         NDIVSUM=NDIVSUM+1
         INDEXNN(1)=IWPLUS2
         INDEXNN(2)=IWPLUS1
         INDEXNN(3)=IWMIDDLE
         INDEXNN(4)=IWMINUS1
         INDEXNN(5)=IWMINUS2
         INDEXNN(6)=IVPLUS2
         INDEXNN(7)=IVPLUS1
         INDEXNN(8)=IVMIDDLE
         INDEXNN(9)=IVMINUS1
         INDEXNN(10)=IVMINUS2
         INDEXNN(11)=IUPLUS2
         INDEXNN(12)=IUPLUS1
         INDEXNN(13)=IUMIDDLE
         INDEXNN(14)=IUMINUS1
         INDEXNN(15)=IUMINUS2
         VALUES(1)=VIWPLUS2
         VALUES(2)=VIWPLUS1
         VALUES(3)=VIWMIDDLE
         VALUES(4)=VIWMINUS1
         VALUES(5)=VIWMINUS2
         VALUES(6)=VIVPLUS2
         VALUES(7)=VIVPLUS1
         VALUES(8)=VIVMIDDLE
         VALUES(9)=VIVMINUS1
         VALUES(10)=VIVMINUS2
         VALUES(11)=VIUPLUS2
         VALUES(12)=VIUPLUS1
         VALUES(13)=VIUMIDDLE
         VALUES(14)=VIUMINUS1
         VALUES(15)=VIUMINUS2
         FACTOR=FACTZ
         DO NNN1=1,15
          IF(INDEXNN(NNN1).GT.0)THEN
           INDNNN=INDEX2(INDEXNN(NNN1))
           RHS(INDNNN)=RHS(INDNNN)-DIVDZ1*VALUES(NNN1)*FACTOR
          ENDIF
         ENDDO         
c         IF(IBOUNDTOP.EQ.3)THEN
C          WRITE(6,*)'IBOUNDTOP=3,I,J,K = ',I,J,K
C          WRITE(6,*)'INDEXNN = ',(INDEXNN(ILIL),ILIL=1,NINDEX)
C          WRITE(6,*)'VALUES = ',(VALUES(ILIL),ILIL=1,NINDEX)
C         ENDIF       
c      write(6,*)'i,j,k,indexnn = ',i,j,k,(indexnn(llll),
c     1  llll=1,nindex)
c      write(6,*)'i,j,k,values = ',i,j,k,(values(llll),llll=1,nindex)
c         WRITE(6,*)'MATRIX I,J,K,DX,DY,DZ = ',I,J,K,DX,DY,DZ
c         WRITE(6,*)'DXI,DYI = ',DXI,DYI
         CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTOR,IPOS,NZ,IA,JA,M,IWORK,NW)
         FACTRH10R=FACTRH10/RHOK(K)/RHOK(K)
         IF(INDEX4(I,J,K).EQ.1)THEN
          RF=1.
         ELSE
          RF=100.
         ENDIF
         IF(IABS(IBOUNDX).EQ.2)RF=4.
         IF(IABS(IBOUNDY).EQ.2)RF=RF*4.
         IF(IABS(IBOUNDZ).EQ.2)RF=RF*4.
         FACTRH10R=FACTRH10R*RF
         FACTCORNER=FACTRH10R*2.
         NINDEX=9
         DO JI=1,3
          NADD=(JI-1)*N
          IF(IBOUNDX.EQ.2.AND.IBOUNDY.EQ.2)THEN
           INDEXNN(1)=IPOSITIONU+2+2*IMAX+NADD
           INDEXNN(2)=IPOSITIONU+1+2*IMAX+NADD
           INDEXNN(3)=IPOSITIONU+2*IMAX+NADD
           INDEXNN(4)=IPOSITIONU+2+IMAX+NADD
           INDEXNN(5)=IPOSITIONU+1+IMAX+NADD
           INDEXNN(6)=IPOSITIONU+IMAX+NADD
           INDEXNN(7)=IPOSITIONU+2+NADD
           INDEXNN(8)=IPOSITIONU+1+NADD
           INDEXNN(9)=IPOSITIONU+NADD
           VALUES(1)=.25
           VALUES(2)=-1.0
           VALUES(3)=.75
           VALUES(4)=-1.0
           VALUES(5)=4.0
           VALUES(6)=-3.0
           VALUES(7)=.75
           VALUES(8)=-3.0
           VALUES(9)=2.25
          ELSEIF(IBOUNDX.EQ.2.AND.IBOUNDY.EQ.-2)THEN
           INDEXNN(1)=IPOSITIONU+2+NADD
           INDEXNN(2)=IPOSITIONU+1+NADD
           INDEXNN(3)=IPOSITIONU+NADD
           INDEXNN(4)=IPOSITIONU+2-IMAX+NADD
           INDEXNN(5)=IPOSITIONU+1-IMAX+NADD
           INDEXNN(6)=IPOSITIONU-IMAX+NADD
           INDEXNN(7)=IPOSITIONU+2-2*IMAX+NADD
           INDEXNN(8)=IPOSITIONU+1-2*IMAX+NADD
           INDEXNN(9)=IPOSITIONU-2*IMAX+NADD
           VALUES(1)=-.75
           VALUES(2)=3.0
           VALUES(3)=-2.25
           VALUES(4)=1.0
           VALUES(5)=-4.0
           VALUES(6)=3.0
           VALUES(7)=-.25
           VALUES(8)=1.0
           VALUES(9)=-.75
          ELSEIF(IBOUNDX.EQ.-2.AND.IBOUNDY.EQ.-2)THEN
           INDEXNN(1)=IPOSITIONU+NADD
           INDEXNN(2)=IPOSITIONU-1+NADD
           INDEXNN(3)=IPOSITIONU-2+NADD
           INDEXNN(4)=IPOSITIONU-IMAX+NADD
           INDEXNN(5)=IPOSITIONU-1-IMAX+NADD
           INDEXNN(6)=IPOSITIONU-2-IMAX+NADD
           INDEXNN(7)=IPOSITIONU-2*IMAX+NADD
           INDEXNN(8)=IPOSITIONU-1-2*IMAX+NADD
           INDEXNN(9)=IPOSITIONU-2-2*IMAX+NADD
           VALUES(1)=2.25
           VALUES(2)=-3.0
           VALUES(3)=.75
           VALUES(4)=-3.0
           VALUES(5)=4.0
           VALUES(6)=-1.0
           VALUES(7)=.75
           VALUES(8)=-1.0
           VALUES(9)=.25
          ELSEIF(IBOUNDX.EQ.-2.AND.IBOUNDY.EQ.2.)THEN
           INDEXNN(1)=IPOSITIONU+2*IMAX+NADD
           INDEXNN(2)=IPOSITIONU-1+2*IMAX+NADD
           INDEXNN(3)=IPOSITIONU-2+2*IMAX+NADD
           INDEXNN(4)=IPOSITIONU+IMAX+NADD
           INDEXNN(5)=IPOSITIONU-1+IMAX+NADD
           INDEXNN(6)=IPOSITIONU-2+IMAX+NADD
           INDEXNN(7)=IPOSITIONU+NADD
           INDEXNN(8)=IPOSITIONU-1+NADD
           INDEXNN(9)=IPOSITIONU-2+NADD
           VALUES(1)=-.75
           VALUES(2)=1.0
           VALUES(3)=-.25
           VALUES(4)=3.0
           VALUES(5)=-4.0
           VALUES(6)=1.0
           VALUES(7)=-2.25
           VALUES(8)=3.0
           VALUES(9)=-.75
          ELSEIF(IBOUNDX.EQ.2.AND.IABS(IBOUNDY).LE.1)THEN
           INDEXNN(1)=IPOSITIONU+2+IMAX+NADD
           INDEXNN(2)=IPOSITIONU+1+IMAX+NADD
           INDEXNN(3)=IPOSITIONU+IMAX+NADD
           INDEXNN(4)=0
           INDEXNN(5)=0
           INDEXNN(6)=0
           INDEXNN(7)=IPOSITIONU+2-IMAX+NADD
           INDEXNN(8)=IPOSITIONU+1-IMAX+NADD
           INDEXNN(9)=IPOSITIONU-IMAX+NADD
           VALUES(1)=-.25
           VALUES(2)=1.00
           VALUES(3)=-.75
           VALUES(4)=0.
           VALUES(5)=0.
           VALUES(6)=0.
           VALUES(7)=0.25
           VALUES(8)=-1.0
           VALUES(9)=0.75
          ELSEIF(IBOUNDX.EQ.-2.AND.IABS(IBOUNDY).LE.1)THEN
           INDEXNN(1)=IPOSITIONU+IMAX+NADD
           INDEXNN(2)=IPOSITIONU-1+IMAX+NADD
           INDEXNN(3)=IPOSITIONU-2+IMAX+NADD
           INDEXNN(4)=0
           INDEXNN(5)=0
           INDEXNN(6)=0
           INDEXNN(7)=IPOSITIONU-IMAX+NADD
           INDEXNN(8)=IPOSITIONU-1-IMAX+NADD
           INDEXNN(9)=IPOSITIONU-2-IMAX+NADD
           VALUES(1)=0.75
           VALUES(2)=-1.0
           VALUES(3)=0.25
           VALUES(4)=0.
           VALUES(5)=0.
           VALUES(6)=0.
           VALUES(7)=-.75
           VALUES(8)=1.00
           VALUES(9)=-.25
          ELSEIF(IABS(IBOUNDX).LE.1.AND.IBOUNDY.EQ.-2)THEN
           INDEXNN(1)=IPOSITIONU+1+NADD
           INDEXNN(2)=0
           INDEXNN(3)=IPOSITIONU-1+NADD
           INDEXNN(4)=IPOSITIONU+1-IMAX+NADD
           INDEXNN(5)=0
           INDEXNN(6)=IPOSITIONU-1-IMAX+NADD
           INDEXNN(7)=IPOSITIONU+1-2*IMAX+NADD
           INDEXNN(8)=0
           INDEXNN(9)=IPOSITIONU-1-2*IMAX+NADD
           VALUES(1)=.75
           VALUES(2)=0.
           VALUES(3)=-.75
           VALUES(4)=-1.
           VALUES(5)=0.
           VALUES(6)=1.0
           VALUES(7)=.25
           VALUES(8)=0.
           VALUES(9)=-.25
          ELSEIF(IABS(IBOUNDX).LE.1.AND.IBOUNDY.EQ.2.)THEN
           INDEXNN(1)=IPOSITIONU+1+2*IMAX+NADD
           INDEXNN(2)=0
           INDEXNN(3)=IPOSITIONU-1+2*IMAX+NADD
           INDEXNN(4)=IPOSITIONU+1+IMAX+NADD
           INDEXNN(5)=0
           INDEXNN(6)=IPOSITIONU-1+IMAX+NADD
           INDEXNN(7)=IPOSITIONU+1+NADD
           INDEXNN(8)=0
           INDEXNN(9)=IPOSITIONU-1+NADD
           VALUES(1)=-.25
           VALUES(2)=0.
           VALUES(3)=0.25
           VALUES(4)=1.00
           VALUES(5)=0.
           VALUES(6)=-1.0
           VALUES(7)=-.75
           VALUES(8)=0.
           VALUES(9)=0.75
          ELSE
           GO TO 234
          ENDIF
C          write(6,*)i,j,k,iboundx,iboundy,factcorner
C          write(6,*)nindex,(indexnn(kkj),values(kkj),kkj=1,nindex)
          CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTCORNER,IPOS,NZ,IA,JA,M,IWORK,NW)
         ENDDO
234      IF(IABS(IBOUNDX).LE.1.AND.IABS(IBOUNDY).LE.1)THEN
          NINDEX=4
          DO JI=1,3
           NADD=(JI-1)*N
           INDEXNN(1)=IPOSITIONU+1+IMAX+NADD
           INDEXNN(2)=IPOSITIONU-1+IMAX+NADD
           INDEXNN(3)=IPOSITIONU+1-IMAX+NADD
           INDEXNN(4)=IPOSITIONU-1-IMAX+NADD
           VALUES(1)=.25
           VALUES(2)=-.25
           VALUES(3)=-.25
           VALUES(4)=.25
            FACTXY=2.*FACTRH10R
            CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTXY,IPOS,NZ,IA,JA,M,IWORK,NW)
          ENDDO
         ENDIF
         IF(IABS(IBOUNDX).LE.1.AND.IABS(IBOUNDZ).LE.1)THEN
          NINDEX=4
          DO JI=1,3
           NADD=(JI-1)*N
           INDEXNN(1)=IPOSITIONU+1+IMAX*JMAX+NADD
           INDEXNN(2)=IPOSITIONU-1+IMAX*JMAX+NADD
           INDEXNN(3)=IPOSITIONU+1-IMAX*JMAX+NADD
           INDEXNN(4)=IPOSITIONU-1-IMAX*JMAX+NADD
           VALUES(1)=.25
           VALUES(2)=-.25
           VALUES(3)=-.25
           VALUES(4)=.25
            FACTXY=2.*FACTRH10R
            CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTXY,IPOS,NZ,IA,JA,M,IWORK,NW)
          ENDDO
         ENDIF
         IF(IABS(IBOUNDY).LE.1.AND.IABS(IBOUNDZ).LE.1)THEN
          NINDEX=4
          DO JI=1,3
           NADD=(JI-1)*N
           INDEXNN(1)=IPOSITIONU+IMAX+IMAX*JMAX+NADD
           INDEXNN(2)=IPOSITIONU-IMAX+IMAX*JMAX+NADD
           INDEXNN(3)=IPOSITIONU+IMAX-IMAX*JMAX+NADD
           INDEXNN(4)=IPOSITIONU-IMAX-IMAX*JMAX+NADD
           VALUES(1)=.25
           VALUES(2)=-.25
           VALUES(3)=-.25
           VALUES(4)=.25
            FACTXY=2.*FACTRH10R
            CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTXY,IPOS,NZ,IA,JA,M,IWORK,NW)
          ENDDO
         ENDIF
         NINDEX=5
         DO IJI=1,3
          NADD=(IJI-1)*N
          IF(IBOUNDX.EQ.-2)THEN
           VALUES(1)=35./12.*XFACTOR
           VALUES(2)=-26./3.*XFACTOR
           VALUES(3)=19./2.*XFACTOR
           VALUES(4)=-14./3.*XFACTOR
           VALUES(5)=11./12.*XFACTOR
          ELSEIF(IBOUNDX.EQ.2)THEN
           VALUES(5)=35./12.*XFACTOR
           VALUES(4)=-26./3.*XFACTOR
           VALUES(3)=19./2.*XFACTOR
           VALUES(2)=-14./3.*XFACTOR
           VALUES(1)=11./12.*XFACTOR
          ELSEIF(IBOUNDX.EQ.-1)THEN
           VALUES(1)=1.*XFACTOR
           VALUES(2)=-2.*XFACTOR
           VALUES(3)=1.*XFACTOR
           VALUES(4)=0.
           VALUES(5)=0.
          ELSEIF(IBOUNDX.EQ.1)THEN
           VALUES(5)=1.*XFACTOR
           VALUES(4)=-2.*XFACTOR
           VALUES(3)=1.*XFACTOR
           VALUES(2)=0.
           VALUES(1)=0.
          ELSE
           VALUES(1)=0.
           VALUES(2)=1.*XFACTOR
           VALUES(3)=-2.*XFACTOR
           VALUES(4)=1.*XFACTOR
           VALUES(5)=0.
          ENDIF
          INDEXNN(1)=IUPLUS2+NADD
          INDEXNN(2)=IUPLUS1+NADD
          INDEXNN(3)=IUMIDDLE+NADD
          INDEXNN(4)=IUMINUS1+NADD
          INDEXNN(5)=IUMINUS2+NADD
           CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTRH10R,IPOS,NZ,IA,JA,M,IWORK,NW)
          IF(IBOUNDY.EQ.-2)THEN
           VALUES(1)=35./12.*YFACTOR
           VALUES(2)=-26./3.*YFACTOR
           VALUES(3)=19./2.*YFACTOR
           VALUES(4)=-14./3.*YFACTOR
           VALUES(5)=11./12.*YFACTOR
          ELSEIF(IBOUNDY.EQ.2)THEN
           VALUES(5)=35./12.*YFACTOR
           VALUES(4)=-26./3.*YFACTOR
           VALUES(3)=19./2.*YFACTOR
           VALUES(2)=-14./3.*YFACTOR
           VALUES(1)=11./12.*YFACTOR
          ELSEIF(IBOUNDY.EQ.-1)THEN
           VALUES(1)=1.*YFACTOR
           VALUES(2)=-2.*YFACTOR
           VALUES(3)=1.*YFACTOR
           VALUES(4)=0.
           VALUES(5)=0.
          ELSEIF(IBOUNDY.EQ.1)THEN
           VALUES(5)=1.*YFACTOR
           VALUES(4)=-2.*YFACTOR
           VALUES(3)=1.*YFACTOR
           VALUES(2)=0.
           VALUES(1)=0.
          ELSE
           VALUES(1)=0.
           VALUES(2)=1.*YFACTOR
           VALUES(3)=-2.*YFACTOR
           VALUES(4)=1.*YFACTOR
           VALUES(5)=0.
          ENDIF
          INDEXNN(1)=IVPLUS2-N+NADD
          INDEXNN(2)=IVPLUS1-N+NADD
          INDEXNN(3)=IVMIDDLE-N+NADD
          INDEXNN(4)=IVMINUS1-N+NADD
          INDEXNN(5)=IVMINUS2-N+NADD
           CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTRH10R,IPOS,NZ,IA,JA,M,IWORK,NW)
          IF(IBOUNDZ.EQ.-2)THEN
           VALUES(1)=35./12.
           VALUES(2)=-26./3.
           VALUES(3)=19./2.
           VALUES(4)=-14./3.
           VALUES(5)=11./12.
          ELSEIF(IBOUNDZ.EQ.2)THEN
           VALUES(5)=35./12.
           VALUES(4)=-26./3.
           VALUES(3)=19./2.
           VALUES(2)=-14./3.
           VALUES(1)=11./12.
          ELSEIF(IBOUNDZ.EQ.-1)THEN
           VALUES(1)=1.
           VALUES(2)=-2.
           VALUES(3)=1.
           VALUES(4)=0.
           VALUES(5)=0.
          ELSEIF(IBOUNDZ.EQ.1)THEN
           VALUES(5)=1.
           VALUES(4)=-2.
           VALUES(3)=1.
           VALUES(2)=0.
           VALUES(1)=0.
          ELSE
           VALUES(1)=0.
           VALUES(2)=1.
           VALUES(3)=-2.
           VALUES(4)=1.
           VALUES(5)=0.
          ENDIF
          INDEXNN(1)=IWPLUS2-2*N+NADD
          INDEXNN(2)=IWPLUS1-2*N+NADD
          INDEXNN(3)=IWMIDDLE-2*N+NADD
          INDEXNN(4)=IWMINUS1-2*N+NADD
          INDEXNN(5)=IWMINUS2-2*N+NADD
           CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTRH10R,IPOS,NZ,IA,JA,M,IWORK,NW)
         ENDDO
         IF(IBOUNDZ.EQ.-2000)THEN
          INDEXNN(1)=IWPLUS2
          VALUES(1)=1.
          FACTOR=FACTZ
          NINDEX=1
          DO NNN1=1,NINDEX
           IF(INDEXNN(NNN1).GT.0)THEN
            INDNNN=INDEX2(INDEXNN(NNN1))
            RHS(INDNNN)=RHS(INDNNN)-DIVERGENCEH*FACTVR*
     1       VALUES(NNN1)*FACTOR
           ENDIF
          ENDDO         
          CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTOR,IPOS,NZ,IA,JA,M,IWORK,NW)
         ENDIF
         IF(IBOUNDX.EQ.0)THEN
          IUPLUS2=0
          IUPLUS1=IPOSITIONU+1
          IUMIDDLE=0
          IUMINUS1=IPOSITIONU-1
          IUMINUS2=0
          VIUPLUS2=0.
          VIUPLUS1=1./2.*DZ1/DX1
          VIUMIDDLE=0.
          VIUMINUS1=-1./2.*DZ1/DX1
          VIUMINUS2=0.
          IP2=I+2
          IP1=I+1
          IP0=I
          IM1=I-1
          IM2=I-2         
         ELSEIF(IBOUNDX.EQ.1)THEN
          IUPLUS2=0
          IUPLUS1=0
          IUMIDDLE=IPOSITIONU+1
          IUMINUS1=0
          IUMINUS2=IPOSITIONU-1
          VIUPLUS2=0.
          VIUPLUS1=0.
          VIUMIDDLE=1./2.*DZ1/DX1
          VIUMINUS1=0.
          VIUMINUS2=-1./2.*DZ1/DX1
          IP2=I+3
          IP1=I+2
          IP0=I+1
          IM1=I
          IM2=I-1
         ELSEIF(IBOUNDX.EQ.-1)THEN
          IUPLUS2=IPOSITIONU+1
          IUPLUS1=0
          IUMIDDLE=IPOSITIONU-1
          IUMINUS1=0
          IUMINUS2=0
          VIUPLUS2=1./2.*DZ1/DX1
          VIUPLUS1=0.
          VIUMIDDLE=-1./2.*DZ1/DX1
          VIUMINUS1=0.
          VIUMINUS2=0.
          IP2=I+1
          IP1=I
          IP0=I-1
          IM1=I-2
          IM2=I-3
         ELSEIF(IBOUNDX.EQ.2)THEN
          IUPLUS2=0
          IUPLUS1=0
          IUMIDDLE=IPOSITIONU+2
          IUMINUS1=IPOSITIONU+1
          IUMINUS2=IPOSITIONU
          VIUPLUS2=0.
          VIUPLUS1=0.
          VIUMIDDLE=-.5*DZ1/DX1
          VIUMINUS1=2.*DZ1/DX1
          VIUMINUS2=-1.5*DZ1/DX1
          IP2=I+4
          IP1=I+3
          IP0=I+2
          IM1=I+1
          IM2=I
         ELSEIF(IBOUNDX.EQ.-2)THEN
          IUPLUS2=IPOSITIONU
          IUPLUS1=IPOSITIONU-1
          IUMIDDLE=IPOSITIONU-2
          IUMINUS1=0
          IUMINUS2=0
          VIUPLUS2=1.5*DZ1/DX1
          VIUPLUS1=-2.*DZ1/DX1
          VIUMIDDLE=.5*DZ1/DX1
          VIUMINUS1=0.
          VIUMINUS2=0.
          IP2=I
          IP1=I-1
          IP0=I-2
          IM1=I-3
          IM2=I-4
         ENDIF
         IF(IBOUNDY.EQ.0)THEN
          IVPLUS2=0
          IVPLUS1=IPOSITIONV+IMAX
          IVMIDDLE=0
          IVMINUS1=IPOSITIONV-IMAX
          IVMINUS2=0
          VIVPLUS2=0
          VIVPLUS1=1./2.*DZ1/DY1
          VIVMIDDLE=0.
          VIVMINUS1=-1./2.*DZ1/DY1
          VIVMINUS2=0.
          JP2=J+2
          JP1=J+1
          JP0=J
          JM1=J-1
          JM2=J-2
         ELSEIF(IBOUNDY.EQ.1)THEN
          IVPLUS2=0
          IVPLUS1=0
          IVMIDDLE=IPOSITIONV+IMAX
          IVMINUS1=0
          IVMINUS2=IPOSITIONV-IMAX
          VIVPLUS2=0.
          VIVPLUS1=0.
          VIVMIDDLE=1./2.*DZ1/DY1
          VIVMINUS1=0.
          VIVMINUS2=-1./2.*DZ1/DY1
          JP2=J+3
          JP1=J+2
          JP0=J+1
          JM1=J
          JM2=J-1
         ELSEIF(IBOUNDY.EQ.-1)THEN
          IVPLUS2=IPOSITIONV+IMAX
          IVPLUS1=0
          IVMIDDLE=IPOSITIONV-IMAX
          IVMINUS1=0
          IVMINUS2=0
          VIVPLUS2=1./2.*DZ1/DY1
          VIVPLUS1=0.
          VIVMIDDLE=-1./2.*DZ1/DY1
          VIVMINUS1=0.
          VIVMINUS2=0.
          JP2=J+1
          JP1=J
          JP0=J-1
          JM1=J-2
          JM2=J-3
         ELSEIF(IBOUNDY.EQ.2)THEN
          IVPLUS2=0
          IVPLUS1=0
          IVMIDDLE=IPOSITIONV+2*IMAX
          IVMINUS1=IPOSITIONV+IMAX
          IVMINUS2=IPOSITIONV
          VIVPLUS2=0.
          VIVPLUS1=0.
          VIVMIDDLE=-.5*DZ1/DY1
          VIVMINUS1=2.*DZ1/DY1
          VIVMINUS2=-1.5*DZ1/DY1
          JP2=J+4
          JP1=J+3
          JP0=J+2
          JM1=J+1
          JM2=J
         ELSEIF(IBOUNDY.EQ.-2)THEN
          IVPLUS2=IPOSITIONV
          IVPLUS1=IPOSITIONV-IMAX
          IVMIDDLE=IPOSITIONV-2*IMAX
          IVMINUS1=0
          IVMINUS2=0
          VIVPLUS2=1.5*DZ1/DY1
          VIVPLUS1=-2.*DZ1/DY1
          VIVMIDDLE=.5*DZ1/DY1
          VIVMINUS1=0.
          VIVMINUS2=0.
          JP2=J
          JP1=J-1
          JP0=J-2
          JM1=J-3
          JM2=J-4
         ENDIF
         IF(IBOUNDZ.EQ.0)THEN
          IWPLUS2=0
          IWPLUS1=IPOSITIONW+IMAX*JMAX
          IWMIDDLE=0
          IWMINUS1=IPOSITIONW-IMAX*JMAX
          IWMINUS2=0
          VIWPLUS2=0.
          VIWPLUS1=1./2.
          VIWMIDDLE=0.
          VIWMINUS1=-1./2.
          VIWMINUS2=0.
          KP2=K+2
          KP1=K+1
          KP0=K
          KM1=K-1
          KM2=K-2
         ELSEIF(IBOUNDZ.EQ.1)THEN
          IWPLUS2=0
          IWPLUS1=0
          IWMIDDLE=IPOSITIONW+IMAX*JMAX
          IWMINUS1=0
          IWMINUS2=IPOSITIONW-IMAX*JMAX
          VIWPLUS2=0.
          VIWPLUS1=0.
          VIWMIDDLE=1./2.
          VIWMINUS1=0.
          VIWMINUS2=-1./2.
          KP2=K+3
          KP1=K+2
          KP0=K+1
          KM1=K
          KM2=K-1
         ELSEIF(IBOUNDZ.EQ.-1)THEN
          IWPLUS2=IPOSITIONW+IMAX*JMAX
          IWPLUS1=0
          IWMIDDLE=IPOSITIONW-IMAX*JMAX
          IWMINUS1=0
          IWMINUS2=0
          VIWPLUS2=1./2.
          VIWPLUS1=0.
          VIWMIDDLE=-1./2.
          VIWMINUS1=0.
          VIWMINUS2=0.
          KP2=K+1
          KP1=K
          KP0=K-1
          KM1=K-2
          KM2=K-3
         ELSEIF(IBOUNDZ.EQ.2)THEN
          IWPLUS2=0
          IWPLUS1=0
          IWMIDDLE=IPOSITIONW+2*IMAX*JMAX
          IWMINUS1=IPOSITIONW+IMAX*JMAX
          IWMINUS2=IPOSITIONW
          VIWPLUS2=0.
          VIWPLUS1=0.
          VIWMIDDLE=-.5
          VIWMINUS1=2.
          VIWMINUS2=-1.5
          KP2=K+4
          KP1=K+3
          KP0=K+2
          KM1=K+1
          KM2=K
         ELSEIF(IBOUNDZ.EQ.-2)THEN
          IWPLUS2=IPOSITIONW
          IWPLUS1=IPOSITIONW-IMAX*JMAX
          IWMIDDLE=IPOSITIONW-2*IMAX*JMAX
          IWMINUS1=0
          IWMINUS2=0
          VIWPLUS2=1.5
          VIWPLUS1=-2.
          VIWMIDDLE=.5
          VIWMINUS1=0.
          VIWMINUS2=0.
          KP2=K
          KP1=K-1
          KP0=K-2
          KM1=K-3
          KM2=K-4
         ENDIF
         NINDEX=15
         DIVERGENCE=0
         DIVDZ1=0.
         INDEXNN(1)=IWPLUS2
         INDEXNN(2)=IWPLUS1
         INDEXNN(3)=IWMIDDLE
         INDEXNN(4)=IWMINUS1
         INDEXNN(5)=IWMINUS2
         INDEXNN(6)=IVPLUS2
         INDEXNN(7)=IVPLUS1
         INDEXNN(8)=IVMIDDLE
         INDEXNN(9)=IVMINUS1
         INDEXNN(10)=IVMINUS2
         INDEXNN(11)=IUPLUS2
         INDEXNN(12)=IUPLUS1
         INDEXNN(13)=IUMIDDLE
         INDEXNN(14)=IUMINUS1
         INDEXNN(15)=IUMINUS2
         VALUES(1)=VIWPLUS2
         VALUES(2)=VIWPLUS1
         VALUES(3)=VIWMIDDLE
         VALUES(4)=VIWMINUS1
         VALUES(5)=VIWMINUS2
         VALUES(6)=VIVPLUS2
         VALUES(7)=VIVPLUS1
         VALUES(8)=VIVMIDDLE
         VALUES(9)=VIVMINUS1
         VALUES(10)=VIVMINUS2
         VALUES(11)=VIUPLUS2
         VALUES(12)=VIUPLUS1
         VALUES(13)=VIUMIDDLE
         VALUES(14)=VIUMINUS1
         VALUES(15)=VIUMINUS2
         FACTOR=FACTWP

         DO NNN1=1,15
          IF(INDEXNN(NNN1).GT.0)THEN
           INDNNN=INDEX2(INDEXNN(NNN1))
           RHS(INDNNN)=RHS(INDNNN)-DIVDZ1*VALUES(NNN1)*FACTOR
          ENDIF
         ENDDO         
c         IF(IBOUNDTOP.EQ.3)THEN
C          WRITE(6,*)'IBOUNDTOP=3,I,J,K = ',I,J,K
C          WRITE(6,*)'INDEXNN = ',(INDEXNN(ILIL),ILIL=1,NINDEX)
C          WRITE(6,*)'VALUES = ',(VALUES(ILIL),ILIL=1,NINDEX)
C         ENDIF       
c      write(6,*)'i,j,k,indexnn = ',i,j,k,(indexnn(llll),
c     1  llll=1,nindex)
c      write(6,*)'i,j,k,values = ',i,j,k,(values(llll),llll=1,nindex)
c         WRITE(6,*)'MATRIX I,J,K,DX,DY,DZ = ',I,J,K,DX,DY,DZ
c         WRITE(6,*)'DXI,DYI = ',DXI,DYI
         CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1                  NNNN,FACTOR,IPOS,NZ,IA,JA,M,IWORK,NW)
         NINDEX=3
         IF(IBOUNDX.EQ.1.OR.IBOUNDX.EQ.-1)THEN
          VALUES(1)=1.
          VALUES(2)=-2.
          VALUES(3)=1.
          INDEXNN(1)=IPOSITIONU-1
          INDEXNN(2)=IPOSITIONU
          INDEXNN(3)=IPOSITIONU+1
         ENDIF
         IF(IBOUNDY.EQ.1.OR.IBOUNDY.EQ.-1)THEN
          IF(INDEX4(I,J,K).NE.0)THEN
           VALUES(1)=1.
           VALUES(2)=-2.
           VALUES(3)=1.
           INDEXNN(1)=IPOSITIONV-IMAX
           INDEXNN(2)=IPOSITIONV
           INDEXNN(3)=IPOSITIONV+IMAX
          ELSE
           VALUES(1)=VIVPLUS2
           VALUES(2)=VIVPLUS1
           VALUES(3)=VIVMIDDLE
           VALUES(4)=VIVMINUS1
           VALUES(5)=VIVMINUS2
           INDEXNN(1)=IVPLUS2
           INDEXNN(2)=IVPLUS1
           INDEXNN(3)=IVMIDDLE
           INDEXNN(4)=IVMINUS1
           INDEXNN(5)=IVMINUS2
          ENDIF
         ENDIF
         IF(IBOUNDZ.EQ.1.OR.IBOUNDZ.EQ.-1)THEN
          IF(K.NE.1)THEN
           VALUES(1)=1.
           INDEXNN(1)=IPOSITIONW-IMAX*JMAX
          ELSE
           VALUES(1)=1.
           INDEXNN(1)=0
          ENDIF
          VALUES(2)=-2.
          VALUES(3)=1.
          INDEXNN(2)=IPOSITIONW
          INDEXNN(3)=IPOSITIONW+IMAX*JMAX
         ENDIF
9998    ENDDO
       ENDDO
      ENDDO
C Make sure all points used are anchored
77778 LUW=95

      FACTOR=1.
      MODE=1
      DIVSUM=DIVSUM/NDIVSUM
      DIVSUM2=DSQRT(DIVSUM2/NDIVSUM)
      WRITE(6,*)'DIVERGENCE ERROR, RMS ERROR = ',DIVSUM,DIVSUM2
      DO K=1,KMAX
       DO J=1,JMAX
        DO I=1,IMAX
         V(I,J,K)=0.
         W(I,J,K)=0.
        ENDDO
       ENDDO
      ENDDO
      DO NF=1,NFILE
       IF(NF.GT.2)THEN
        FACTORPOINT=.01
       ELSE
        FACTORPOINT=1.
       ENDIF
       OPEN(LUW,FILE=OSUMFILENAME(NF),FORM='UNFORMATTED')
       CALL READABCD(LUW,IMAX,JMAX,KMAX,SUM1,SUM2,SUM3,
     +    SUM4,SUM5,SUM6,SUM7,SUM8,SUM9,SUM10,SUM11,SUM12,SUM13,
     +    V,KMAX2,XWEIGHT,YWEIGHT,ZWEIGHT)
       DO K=1,KMAX
        write(6,*)'anchoring points at level ',K
        HEIGHT=Z0+(K-1)*SZ
c        WRITE(6,*)'HEIGHT,Z0,SZ,K',HEIGHT,Z0,SZ,K
        RHO=RHOK(K)
        IF(K.GT.3)THEN
         FACTPOINT=FACTVR2/RHO/RHO*FACTORPOINT
        ELSEIF(K.EQ.1)THEN
         FACTPOINT=1.*FACTVR2/RHO/RHO*FACTORPOINT
        ELSEIF(K.EQ.2)THEN
         FACTPOINT=1.*FACTVR2/RHO/RHO*FACTORPOINT
        ELSEIF(K.EQ.3)THEN
         FACTPOINT=1.*FACTVR2/RHO/RHO*FACTORPOINT
        ENDIF
        write(6,*)'imax,jmax = ',imax,jmax
        DO J=1,JMAX
         DO I=1,IMAX
          if(i.eq.30.and.j.eq.35)then
c           write(6,*)'index4 at level ',k,' is ',index4(i,j,k)
          endif
c          RHOZ=EXP(ZWEIGHT(I,J,K)*.10437052)
          RHOZ=1.0
          IPOSITION=I+(J-1)*IMAX+(K-1)*IJMAX
          IPOSITIONU=IPOSITION
          IPOSITIONV=IPOSITION+N
          IPOSITIONW=IPOSITION+2*N
          INPOSITION=INDEX2(IPOSITION)
          INPOSITIONU=INDEX2(IPOSITIONU)
          INPOSITIONV=INDEX2(IPOSITIONV)
          INPOSITIONW=INDEX2(IPOSITIONW)
          SUMSTEST=SUM3(I,J,K)+SUM5(I,J,K)+SUM8(I,J,K)
          IF(SUMSTEST.GT.0.AND.FACTZZ.LE.0..AND.FACTWP.LE.0.)THEN
           BAND1(INPOSITIONU)=1.
           BAND1(INPOSITIONV)=1.
           BAND1(INPOSITIONW)=1.
          ENDIF
          IF(U(I,J,K).GT.FLAG)THEN
           NINDEX=1
           IF(INDEX4(I,J,K).GT.0.OR.SUMSTEST.GT.0.OR.
     1       K.LE.KTTOP(I,J))THEN
            VALUES(1)=1.
            WTBACKUU=WTBACKU
            WTBACKVV=WTBACKV
            WTBACKWW=WTBACKW
           ELSE
            VALUES(1)=1.
            WTBACKUU=WTBACKU*10.
            WTBACKVV=WTBACKV*10.
            WTBACKWW=WTBACKW*10.
           ENDIF
           INDEXNN(1)=IPOSITIONU
           CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1      NNNN,WTBACKUU,IPOS,NZ,IA,JA,M,IWORK,NW)
c           write(6,*)'before,i,j,k,rhsu,rhsv,rhsw = ',
c     1      i,j,k,rhs(inpositionu),rhs(inpositionv),rhs(inpositionw)
           RHS(INPOSITIONU)=RHS(INPOSITIONU)
     1      +WTBACKUU*UGUESS(I,J,K)
           INDEXNN(1)=IPOSITIONV
           CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1      NNNN,WTBACKVV,IPOS,NZ,IA,JA,M,IWORK,NW)
           RHS(INPOSITIONV)=RHS(INPOSITIONV)
     1      +WTBACKVV*VGUESS(I,J,K)
           INDEXNN(1)=IPOSITIONW
           CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1      NNNN,WTBACKWW,IPOS,NZ,IA,JA,M,IWORK,NW)
           RHS(INPOSITIONW)=RHS(INPOSITIONW)
     1      +WTBACKWW*WGUESS(I,J,K)
c           write(6,*)'i,j,k,values,u,v,w = ',i,j,k,values(1),
c     1      uguess(i,j,k),vguess(i,j,k),wguess(i,j,k)
c           write(6,*)'after,i,j,k,rhsu,rhsv,rhsw = ',
c     1      i,j,k,rhs(inpositionu),rhs(inpositionv),rhs(inpositionw)
           VALUES(1)=1.
           IF(K.EQ.1.OR.(FACTZZ.LE.0..AND.FACTWP.LE.0.))THEN
            NINDEX=1
            INDEXNN(1)=IPOSITIONW
            VALUES(1)=1.
            FACTZ=FACTWBOT
            IF(TESTRMS(I,J,K).GT..25)THEN
             FACTRMS=0.
            ELSEIF(TESTRMS(I,J,K).GT..15)THEN
             FACTRMS=100*(.25-TESTRMS(I,J,K))**2.
            ELSE
             FACTRMS=1.
            ENDIF
            FACTGAMMA=(1.-TESTGAMMA(I,J,K)*TESTGAMMA(I,J,K))
            FACTZ=FACTZ*FACTRMS*FACTGAMMA
c            if(factz.lt.0..or.factz.gt.1.)then
c             write(6,*)'i,j,k,rms,gamma,factgamma,factrms,factz = ',
c     1  i,j,k,testrms(i,j,k),testgamma(i,j,k),factgamma,factrms,factz
c            endif
            IF(FACTZ.GT.0.)THEN
             CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1      NNNN,FACTZ,IPOS,NZ,IA,JA,M,IWORK,NW)
            ENDIF
           ENDIF
          ENDIF
          IF(K.EQ.KMAX.OR.K.EQ.KTTOP(I,J))THEN
           KTT=KT(I,J)
           NINDEX=1
           INDEXNN(1)=IPOSITIONW
           VALUES(1)=1.
           FACTZ=FACTWTOP
           IF(FACTZ.LT.0.)FACTZ=0
c           write(6,*)'applying top at ',i,j,k
C CHANGED BY MBELL 12/28/05 IN CASE KTT < 1
           IF(KTT.LE.0) KTT = 0

            IF(TESTRMS(I,J,KTT).GT..25)THEN
             FACTRMS=0.
            ELSEIF(TESTRMS(I,J,KTT).GT..15)THEN
             FACTRMS=100*(.25-TESTRMS(I,J,KTT))**2.
            ELSE
             FACTRMS=1.
            ENDIF
            FACTGAMMA=(1.-TESTGAMMA(I,J,KTT)*TESTGAMMA(I,J,KTT))
            FACTZ=FACTZ*FACTRMS*FACTGAMMA/RHO/RHO
c            if(factz.lt.0..or.factz.gt.1.)then
c             write(6,*)'i,j,k,rms,gamma,factgamma,factrms,factz = ',
c     1  i,j,k,testrms(i,j,k),testgamma(i,j,k),factgamma,factrms,factz
c            endif
           IF(FACTZ.GT.0.)THEN
            CALL MATRIXSUM(INDEXNN,INDEX2,VALUES,NINDEX,NN,NORM,IDIF,
     1      NNNN,FACTZ,IPOS,NZ,IA,JA,M,IWORK,NW)
           ENDIF
          ENDIF
          IF(K.GT.KTTOP(I,J).AND.SUMSTEST.LE.0)GO TO 5434
          K1=K
          IF(K.GT.KT(I,J).AND.SUMSTEST.LE.0.AND.FACTWTOP.GT.0.)THEN
           K1=KT(I,J)
          ELSE
           K1=K
          ENDIF
c          write(6,*)'i,j,k,k1 = ',i,j,k,k1
          IF(INPOSITIONU.GT.0.AND.INPOSITIONV.GT.0.AND.
     +      INPOSITIONW.GT.0.AND.
     +      ABS(XWEIGHT(I,J,K1)).LT.XACC.AND.
     +      ABS(YWEIGHT(I,J,K1)).LT.YACC.AND.
     +      ABS(ZWEIGHT(I,J,K1)).LT.ZACC)THEN

           IF(BAND1(INPOSITIONU).GT.0..AND.BAND1(INPOSITIONV).GT.0.
     +      .AND.BAND1(INPOSITIONW).GT.0.)THEN
             IF(K1.EQ.K)THEN
              SUM2(I,J,K1)=-SUM2(I,J,K1)
              SUM4(I,J,K1)=-SUM4(I,J,K1)
              SUM9(I,J,K1)=-SUM9(I,J,K1)
             ENDIF
            SM=FACTPOINT*SUM5(I,J,K1)
            CALL SBSIJ(IPOS,NZ,IA,JA,M,IWORK,INPOSITIONU,
     1      INPOSITIONU,SM,
     1      MODE,LEVEL,NOUT,IER)
            SM=FACTPOINT*SUM3(I,J,K1)
            CALL SBSIJ(IPOS,NZ,IA,JA,M,IWORK,INPOSITIONV,
     1      INPOSITIONV,SM,
     1      MODE,LEVEL,NOUT,IER)
            SM=FACTPOINT*SUM8(I,J,K1)
            CALL SBSIJ(IPOS,NZ,IA,JA,M,IWORK,INPOSITIONW,
     1      INPOSITIONW,SM,
     1      MODE,LEVEL,NOUT,IER)
            SM=SUM1(I,J,K1)*FACTPOINT
            CALL SBSIJ(IPOS,NZ,IA,JA,M,IWORK,INPOSITIONU,
     1      INPOSITIONV,SM,
     1      MODE,LEVEL,NOUT,IER)
            SM=SUM6(I,J,K1)*FACTPOINT
            CALL SBSIJ(IPOS,NZ,IA,JA,M,IWORK,INPOSITIONU,
     1      INPOSITIONW,SM,
     1      MODE,LEVEL,NOUT,IER)
            SM=SUM7(I,J,K1)*FACTPOINT
            CALL SBSIJ(IPOS,NZ,IA,JA,M,IWORK,INPOSITIONV,
     1      INPOSITIONW,SM,
     1      MODE,LEVEL,NOUT,IER)
            IF(IVTSUB(NF).NE.0)THEN
             VT=VTERM_NEW(SUMDBZ(I,J,K1),HEIGHT,
     1         hb,dpb,IRSW,ZLOW,ZHIGH)*RHOZ
            ELSE
             VT=0.
            ENDIF
            RHS(INPOSITIONU)=RHS(INPOSITIONU)-(SUM4(I,J,K1)*RHOZ
     +       +SUM6(I,J,K1)*VT*RHOK(K1))*FACTPOINT
c            WRITE(6,*)'I,J,K,RHS = ',I,J,K,RHS(INPOSITIONV)
            RHS(INPOSITIONV)=RHS(INPOSITIONV)-(SUM2(I,J,K1)*RHOZ
     +       +SUM7(I,J,K1)*VT*RHOK(K1))*FACTPOINT
c        WRITE(6,*)'I,J,K,U,V,W,SUM1,SUM7,SUM3,SUM2 = ',I,J,K,
c     +  UGUESS(I,J,K1),VGUESS(I,J,K1),WGUESS(I,J,K1),SUM1(I,J,K1),
c     +  SUM3(I,J,K1),SUM7(I,J,K1),SUM2(I,J,K1)
            RHS(INPOSITIONW)=RHS(INPOSITIONW)-(SUM9(I,J,K1)*RHOZ
     +       +SUM8(I,J,K1)*VT*RHOK(K1))*FACTPOINT
           ENDIF
          ENDIF
5434     ENDDO
        ENDDO
       ENDDO
      ENDDO
C
      WRITE(6,*)'GAMMAMAX = ',GAMMAMAX
      WRITE(6,*)'AT I,J,K,NF = ',IGAMMA,JGAMMA,KGAMMA,NFGAMMA
      FACTOR=1.
1996  WRITE(6,*)'Reforming w vector'
      CALL SBEND(IPOS,NZ,IA,JA,M,IWORK)
      write(6,*)'ia(i),i=1,100 = ',(ia(i),i=1,100)
      write(6,*)'ja(i),i=1,100 = ',(ja(i),i=1,100)
      write(6,*)'m(i),i=1,100 = ',(m(i),i=1,100)
      DO K=1,KMAX
       DO J=1,JMAX
        DO I=1,IMAX
         IF(DIV(I,J,K).GT.FLAG)THEN
          WGUESS(I,J,K)=DIV(I,J,K)*1.0E+06
         ELSE
          WGUESS(I,J,K)=FLAG
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END 
      SUBROUTINE SORTER(INDX,VALUES,N)
      INTEGER INDX(N)
      DOUBLE PRECISION VALUES(N)
      DO I=1,N
       INTEST=I
       DO L=I,N
        IF(INDX(L).GT.INDX(INTEST))INTEST=L
       ENDDO
       IA=INDX(INTEST)
       INDX(INTEST)=INDX(I)
       INDX(I)=IA
       A=VALUES(INTEST)
       VALUES(INTEST)=VALUES(I)
       VALUES(I)=A
      ENDDO

      RETURN
      END



      SUBROUTINE MATRIXSUM(INDX,INDEX2,VALUES,N,NN,NORM,IDIF,NNNN,
     1                     FACTOR,IPOS,NZ,IA,JA,A,IWORK,NW)

      INTEGER IA(NNNN),JA(NZ),IWORK(NZ)
      DOUBLE PRECISION A(NZ),SM,VALUES(N)
      INTEGER INDEX2(NN),INDX(N),IDIF(115)
      CALL SORTER(INDX,VALUES,N)
c        WRITE(6,*)'INDX = ',INDX
c        WRITE(6,*)'VALUES = ',VALUES
      NOUT=6
      LEVEL=0
      MODE=1
      IF(FACTOR.EQ.0)RETURN
      DO I=1,N
       DO L=1,N
        IF(INDX(I).GE.INDX(L))THEN
         IF(INDX(I).GT.0.AND.INDX(L).GT.0)THEN
          NDIFF=INDX(I)-INDX(L)
          NPOS=100*L+I
          IN=INDEX2(INDX(I))
          INN=INDEX2(INDX(L))
          IF(IN.GT.0.AND.INN.GT.0)THEN
           SM=VALUES(I)*VALUES(L)*
     1                     FACTOR
         INI=INDEX2(INDX(I))
         INL=INDEX2(INDX(L))
         CALL SBSIJ(IPOS,NZ,IA,JA,A,IWORK,INL,INI,SM,
     1    MODE,LEVEL,NOUT,IER)
          ENDIF
         ENDIF
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE BOUND5(I,J,K,IMAX,JMAX,KMAX,U,
     + IBOUNDX,IBOUNDY,IBOUNDZ)
      REAL UTEST(9),U(IMAX,JMAX,KMAX)
      FLAG=-1.0E+10
      IBOUNDX=-99
      IBOUNDY=-99
      IBOUNDZ=-99
      DO LKLK=1,9
       UTEST(LKLK)=FLAG
      ENDDO
      IF(I.GE.5)UTEST(1)=U(I-4,J,K)
      IF(I.GE.4)UTEST(2)=U(I-3,J,K)
      IF(I.GE.3)UTEST(3)=U(I-2,J,K)
      IF(I.GE.2)UTEST(4)=U(I-1,J,K)
      UTEST(5)=U(I,J,K)
      IF(I.LE.IMAX-1)UTEST(6)=U(I+1,J,K)
      IF(I.LE.IMAX-2)UTEST(7)=U(I+2,J,K)
      IF(I.LE.IMAX-3)UTEST(8)=U(I+3,J,K)
      IF(I.LE.IMAX-4)UTEST(9)=U(I+4,J,K)
      IF(UTEST(3).GT.FLAG.AND.UTEST(4).GT.FLAG.AND.
     +      UTEST(5).GT.FLAG.AND.UTEST(6).GT.FLAG.AND.
     +      UTEST(7).GT.FLAG)THEN
       IBOUNDX=0
       ELSEIF(UTEST(2).GT.FLAG.AND.UTEST(3).GT.FLAG.AND.
     +          UTEST(4).GT.FLAG.AND.UTEST(5).GT.FLAG.AND.
     +          UTEST(6).GT.FLAG)THEN
          IBOUNDX=-1
       ELSEIF(UTEST(4).GT.FLAG.AND.UTEST(5).GT.FLAG.AND.
     +          UTEST(6).GT.FLAG.AND.UTEST(7).GT.FLAG.AND.
     +          UTEST(8).GT.FLAG)THEN
        IBOUNDX=1
       ELSEIF(UTEST(1).GT.FLAG.AND.UTEST(2).GT.FLAG.AND.
     +          UTEST(3).GT.FLAG.AND.UTEST(4).GT.FLAG.AND.
     +          UTEST(5).GT.FLAG)THEN
        IBOUNDX=-2
       ELSEIF(UTEST(5).GT.FLAG.AND.UTEST(6).GT.FLAG.AND.
     +          UTEST(7).GT.FLAG.AND.UTEST(8).GT.FLAG.AND.
     +          UTEST(9).GT.FLAG)THEN
        IBOUNDX=2
      ENDIF
      DO LKLK=1,9
       UTEST(LKLK)=FLAG
      ENDDO
      IF(J.GE.5)UTEST(1)=U(I,J-4,K)
      IF(J.GE.4)UTEST(2)=U(I,J-3,K)
      IF(J.GE.3)UTEST(3)=U(I,J-2,K)
      IF(J.GE.2)UTEST(4)=U(I,J-1,K)
      UTEST(5)=U(I,J,K)
      IF(J.LE.JMAX-1)UTEST(6)=U(I,J+1,K)
      IF(J.LE.JMAX-2)UTEST(7)=U(I,J+2,K)
      IF(J.LE.JMAX-3)UTEST(8)=U(I,J+3,K)
      IF(J.LE.JMAX-4)UTEST(9)=U(I,J+4,K)
      IF(UTEST(3).GT.FLAG.AND.UTEST(4).GT.FLAG.AND.
     +      UTEST(5).GT.FLAG.AND.UTEST(6).GT.FLAG.AND.
     +      UTEST(7).GT.FLAG)THEN
       IBOUNDY=0
      ELSEIF(UTEST(2).GT.FLAG.AND.UTEST(3).GT.FLAG.AND.
     +          UTEST(4).GT.FLAG.AND.UTEST(5).GT.FLAG.AND.
     +          UTEST(6).GT.FLAG)THEN
       IBOUNDY=-1
      ELSEIF(UTEST(4).GT.FLAG.AND.UTEST(5).GT.FLAG.AND.
     +          UTEST(6).GT.FLAG.AND.UTEST(7).GT.FLAG.AND.
     +          UTEST(8).GT.FLAG)THEN
       IBOUNDY=1
      ELSEIF(UTEST(1).GT.FLAG.AND.UTEST(2).GT.FLAG.AND.
     +          UTEST(3).GT.FLAG.AND.UTEST(4).GT.FLAG.AND.
     +          UTEST(5).GT.FLAG)THEN
       IBOUNDY=-2
      ELSEIF(UTEST(5).GT.FLAG.AND.UTEST(6).GT.FLAG.AND.
     +          UTEST(7).GT.FLAG.AND.UTEST(8).GT.FLAG.AND.
     +          UTEST(9).GT.FLAG)THEN
       IBOUNDY=2
      ENDIF
      DO LKLK=1,9
       UTEST(LKLK)=FLAG
      ENDDO
      IF(K.GE.5)UTEST(1)=U(I,J,K-4)
      IF(K.GE.4)UTEST(2)=U(I,J,K-3)
      IF(K.GE.3)UTEST(3)=U(I,J,K-2)
      IF(K.GE.2)UTEST(4)=U(I,J,K-1)
      UTEST(5)=U(I,J,K)
      IF(K.LE.KMAX-1)UTEST(6)=U(I,J,K+1)
      IF(K.LE.KMAX-2)UTEST(7)=U(I,J,K+2)
      IF(K.LE.KMAX-3)UTEST(8)=U(I,J,K+3)
      IF(K.LE.KMAX-4)UTEST(9)=U(I,J,K+4)
      IF(UTEST(3).GT.FLAG.AND.UTEST(4).GT.FLAG.AND.
     +      UTEST(5).GT.FLAG.AND.UTEST(6).GT.FLAG.AND.
     +      UTEST(7).GT.FLAG)THEN
       IBOUNDZ=0
      ELSEIF(UTEST(2).GT.FLAG.AND.UTEST(3).GT.FLAG.AND.
     +          UTEST(4).GT.FLAG.AND.UTEST(5).GT.FLAG.AND.
     +          UTEST(6).GT.FLAG)THEN
       IBOUNDZ=-1
      ELSEIF(UTEST(4).GT.FLAG.AND.UTEST(5).GT.FLAG.AND.
     +          UTEST(6).GT.FLAG.AND.UTEST(7).GT.FLAG.AND.
     +          UTEST(8).GT.FLAG)THEN
       IBOUNDZ=1
      ELSEIF(UTEST(1).GT.FLAG.AND.UTEST(2).GT.FLAG.AND.
     +          UTEST(3).GT.FLAG.AND.UTEST(4).GT.FLAG.AND.
     +          UTEST(5).GT.FLAG)THEN
       IBOUNDZ=-2
      ELSEIF(UTEST(5).GT.FLAG.AND.UTEST(6).GT.FLAG.AND.
     +          UTEST(7).GT.FLAG.AND.UTEST(8).GT.FLAG.AND.
     +          UTEST(9).GT.FLAG)THEN
       IBOUNDZ=2
      ENDIF

      RETURN
      END


      SUBROUTINE VALUESET(IPOSITION,IPOSITIONU,IPOSITIONV,IPOSITIONW,
     +    IMAX,JMAX,KMAX,IBOUNDX,IBOUNDY,IBOUNDZ,DX1,DY1,DZ1,
     +    IUPLUS2,IUPLUS1,IUMIDDLE,IUMINUS1,IUMINUS2,
     +    VIUPLUS2,VIUPLUS1,VIUMIDDLE,VIUMINUS1,VIUMINUS2,
     +    IVPLUS2,IVPLUS1,IVMIDDLE,IVMINUS1,IVMINUS2,
     +    VIVPLUS2,VIVPLUS1,VIVMIDDLE,VIVMINUS1,VIVMINUS2,
     +    IWPLUS2,IWPLUS1,IWMIDDLE,IWMINUS1,IWMINUS2,
     +    VIWPLUS2,VIWPLUS1,VIWMIDDLE,VIWMINUS1,VIWMINUS2,
     +    I,J,K,IP2,IP1,IP0,IM1,IM2,JP2,JP1,JP0,JM1,JM2,
     +    KP2,KP1,KP0,KM1,KM2)
      DOUBLE PRECISION DZ1,DY1,DX1
      DOUBLE PRECISION VIWMINUS2,VIWMINUS1,VIWMIDDLE,VIWPLUS1,VIWPLUS2
      DOUBLE PRECISION VIVMINUS2,VIVMINUS1,VIVMIDDLE,VIVPLUS1,VIVPLUS2
      DOUBLE PRECISION VIUMINUS2,VIUMINUS1,VIUMIDDLE,VIUPLUS1,VIUPLUS2

      IF(IBOUNDX.EQ.0)THEN
       IUPLUS2=IPOSITIONU+2
       IUPLUS1=IPOSITIONU+1
       IUMIDDLE=IPOSITIONU
       IUMINUS1=IPOSITIONU-1
       IUMINUS2=IPOSITIONU-2
       VIUPLUS2=-1./12.*DZ1/DX1
       VIUPLUS1=2./3.*DZ1/DX1
       VIUMIDDLE=0.
       VIUMINUS1=-2./3.*DZ1/DX1
       VIUMINUS2=1./12.*DZ1/DX1
       IP2=I+2
       IP1=I+1
       IP0=I
       IM1=I-1
       IM2=I-2         
      ELSEIF(IBOUNDX.EQ.1)THEN
       IUPLUS2=IPOSITIONU+3
       IUPLUS1=IPOSITIONU+2
       IUMIDDLE=IPOSITIONU+1
       IUMINUS1=IPOSITIONU
       IUMINUS2=IPOSITIONU-1
       VIUPLUS2=1./12.*DZ1/DX1
       VIUPLUS1=-.5*DZ1/DX1
       VIUMIDDLE=3./2.*DZ1/DX1
       VIUMINUS1=-5./6.*DZ1/DX1
       VIUMINUS2=-.25*DZ1/DX1
       IP2=I+3
       IP1=I+2
       IP0=I+1
       IM1=I
       IM2=I-1
      ELSEIF(IBOUNDX.EQ.-1)THEN
       IUPLUS2=IPOSITIONU+1
       IUPLUS1=IPOSITIONU
       IUMIDDLE=IPOSITIONU-1
       IUMINUS1=IPOSITIONU-2
       IUMINUS2=IPOSITIONU-3
       VIUPLUS2=.25*DZ1/DX1
       VIUPLUS1=5./6.*DZ1/DX1
       VIUMIDDLE=-3./2.*DZ1/DX1
       VIUMINUS1=.5*DZ1/DX1
       VIUMINUS2=-1./12.*DZ1/DX1
       IP2=I+1
       IP1=I
       IP0=I-1
       IM1=I-2
       IM2=I-3
      ELSEIF(IBOUNDX.EQ.2)THEN
       IUPLUS2=IPOSITIONU+4
       IUPLUS1=IPOSITIONU+3
       IUMIDDLE=IPOSITIONU+2
       IUMINUS1=IPOSITIONU+1
       IUMINUS2=IPOSITIONU
       VIUPLUS2=-.25*DZ1/DX1
       VIUPLUS1=4./3.*DZ1/DX1
       VIUMIDDLE=-3.*DZ1/DX1
       VIUMINUS1=4.*DZ1/DX1
       VIUMINUS2=-25./12.*DZ1/DX1
       IP2=I+4
       IP1=I+3
       IP0=I+2
       IM1=I+1
       IM2=I
      ELSEIF(IBOUNDX.EQ.-2)THEN
       IUPLUS2=IPOSITIONU
       IUPLUS1=IPOSITIONU-1
       IUMIDDLE=IPOSITIONU-2
       IUMINUS1=IPOSITIONU-3
       IUMINUS2=IPOSITIONU-4
       VIUPLUS2=25./12.*DZ1/DX1
       VIUPLUS1=-4.*DZ1/DX1
       VIUMIDDLE=3.*DZ1/DX1
       VIUMINUS1=-4./3.*DZ1/DX1
       VIUMINUS2=.25*DZ1/DX1
       IP2=I
       IP1=I-1
       IP0=I-2
       IM1=I-3
       IM2=I-4
      ENDIF
      IF(IBOUNDY.EQ.0)THEN
       IVPLUS2=IPOSITIONV+2*IMAX
       IVPLUS1=IPOSITIONV+IMAX
       IVMIDDLE=IPOSITIONV
       IVMINUS1=IPOSITIONV-IMAX
       IVMINUS2=IPOSITIONV-2*IMAX
       VIVPLUS2=-1./12.*DZ1/DY1
       VIVPLUS1=2./3.*DZ1/DY1
       VIVMIDDLE=0.
       VIVMINUS1=-2./3.*DZ1/DY1
       VIVMINUS2=1./12.*DZ1/DY1
       JP2=J+2
       JP1=J+1
       JP0=J
       JM1=J-1
       JM2=J-2
      ELSEIF(IBOUNDY.EQ.1)THEN
       IVPLUS2=IPOSITIONV+3*IMAX
       IVPLUS1=IPOSITIONV+2*IMAX
       IVMIDDLE=IPOSITIONV+IMAX
       IVMINUS1=IPOSITIONV
       IVMINUS2=IPOSITIONV-IMAX
       VIVPLUS2=1./12.*DZ1/DY1
       VIVPLUS1=-.5*DZ1/DY1
       VIVMIDDLE=3./2.*DZ1/DY1
       VIVMINUS1=-5./6.*DZ1/DY1
       VIVMINUS2=-.25*DZ1/DY1
       JP2=J+3
       JP1=J+2
       JP0=J+1
       JM1=J
       JM2=J-1
      ELSEIF(IBOUNDY.EQ.-1)THEN
       IVPLUS2=IPOSITIONV+IMAX
       IVPLUS1=IPOSITIONV
       IVMIDDLE=IPOSITIONV-IMAX
       IVMINUS1=IPOSITIONV-2*IMAX
       IVMINUS2=IPOSITIONV-3*IMAX
       VIVPLUS2=.25*DZ1/DY1
       VIVPLUS1=5./6.*DZ1/DY1
       VIVMIDDLE=-3./2.*DZ1/DY1
       VIVMINUS1=.5*DZ1/DY1
       VIVMINUS2=-1./12.*DZ1/DY1
       JP2=J+1
       JP1=J
       JP0=J-1
       JM1=J-2
       JM2=J-3
      ELSEIF(IBOUNDY.EQ.2)THEN
       IVPLUS2=IPOSITIONV+4*IMAX
       IVPLUS1=IPOSITIONV+3*IMAX
       IVMIDDLE=IPOSITIONV+2*IMAX
       IVMINUS1=IPOSITIONV+IMAX
       IVMINUS2=IPOSITIONV
       VIVPLUS2=-.25*DZ1/DY1
       VIVPLUS1=4./3.*DZ1/DY1
       VIVMIDDLE=-3.*DZ1/DY1
       VIVMINUS1=4.*DZ1/DY1
       VIVMINUS2=-25./12.*DZ1/DY1
       JP2=J+4
       JP1=J+3
       JP0=J+2
       JM1=J+1
       JM2=J
      ELSEIF(IBOUNDY.EQ.-2)THEN
       IVPLUS2=IPOSITIONV
       IVPLUS1=IPOSITIONV-IMAX
       IVMIDDLE=IPOSITIONV-2*IMAX
       IVMINUS1=IPOSITIONV-3*IMAX
       IVMINUS2=IPOSITIONV-4*IMAX
       VIVPLUS2=25./12.*DZ1/DY1
       VIVPLUS1=-4.*DZ1/DY1
       VIVMIDDLE=3.*DZ1/DY1
       VIVMINUS1=-4./3.*DZ1/DY1
       VIVMINUS2=.25*DZ1/DY1
       JP2=J
       JP1=J-1
       JP0=J-2
       JM1=J-3
       JM2=J-4
      ENDIF
      IF(IBOUNDZ.EQ.0)THEN
       IWPLUS2=IPOSITIONW+2*IMAX*JMAX
       IWPLUS1=IPOSITIONW+IMAX*JMAX
       IWMIDDLE=IPOSITIONW
       IWMINUS1=IPOSITIONW-IMAX*JMAX
       IWMINUS2=IPOSITIONW-2*IMAX*JMAX
       VIWPLUS2=-1./12.
       VIWPLUS1=2./3.
       VIWMIDDLE=0.
       VIWMINUS1=-2./3.
       VIWMINUS2=1./12.
       KP2=K+2
       KP1=K+1
       KP0=K
       KM1=K-1
       KM2=K-2
      ELSEIF(IBOUNDZ.EQ.1)THEN
       IWPLUS2=IPOSITIONW+3*IMAX*JMAX
       IWPLUS1=IPOSITIONW+2*IMAX*JMAX
       IWMIDDLE=IPOSITIONW+IMAX*JMAX
       IWMINUS1=IPOSITIONW
       IWMINUS2=IPOSITIONW-IMAX*JMAX
       VIWPLUS2=1./12.
       VIWPLUS1=-.5
       VIWMIDDLE=3./2.
       VIWMINUS1=-5./6.
       VIWMINUS2=-.25
       KP2=K+3
       KP1=K+2
       KP0=K+1
       KM1=K
       KM2=K-1
      ELSEIF(IBOUNDZ.EQ.-1)THEN
       IWPLUS2=IPOSITIONW+IMAX*JMAX
       IWPLUS1=IPOSITIONW
       IWMIDDLE=IPOSITIONW-IMAX*JMAX
       IWMINUS1=IPOSITIONW-2*IMAX*JMAX
       IWMINUS2=IPOSITIONW-3*IMAX*JMAX
       VIWPLUS2=.25
       VIWPLUS1=5./6.
       VIWMIDDLE=-3./2.
       VIWMINUS1=.5
       VIWMINUS2=-1./12.
       KP2=K+1
       KP1=K
       KP0=K-1
       KM1=K-2
       KM2=K-3
      ELSEIF(IBOUNDZ.EQ.2)THEN
       IWPLUS2=IPOSITIONW+4*IMAX*JMAX
       IWPLUS1=IPOSITIONW+3*IMAX*JMAX
       IWMIDDLE=IPOSITIONW+2*IMAX*JMAX
       IWMINUS1=IPOSITIONW+IMAX*JMAX
       IWMINUS2=IPOSITIONW
       VIWPLUS2=-.25
       VIWPLUS1=4./3.
       VIWMIDDLE=-3.
       VIWMINUS1=4.
       VIWMINUS2=-25./12.
       KP2=K+4
       KP1=K+3
       KP0=K+2
       KM1=K+1
       KM2=K
      ELSEIF(IBOUNDZ.EQ.-2)THEN
       IWPLUS2=IPOSITIONW
       IWPLUS1=IPOSITIONW-IMAX*JMAX
       IWMIDDLE=IPOSITIONW-2*IMAX*JMAX
       IWMINUS1=IPOSITIONW-3*IMAX*JMAX
       IWMINUS2=IPOSITIONW-4*IMAX*JMAX
       VIWPLUS2=25./12.
       VIWPLUS1=-4.
       VIWMIDDLE=3.
       VIWMINUS1=-4./3.
       VIWMINUS2=.25
       KP2=K
       KP1=K-1
       KP0=K-2
       KM1=K-3
       KM2=K-4
      ENDIF

      RETURN
      END
