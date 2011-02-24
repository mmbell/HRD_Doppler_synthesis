      SUBROUTINE DIVER3(U,V,W,IMAX,JMAX,KMAX,SX,SY,SZ,
     + DIV)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX)
      REAL DIV(IMAX,JMAX,KMAX)
      DOUBLE PRECISION SX,SY,SZ
      DOUBLE PRECISION DIVDZ1,DIVERGENCE,DZ1,DY1,DX1
      DOUBLE PRECISION VIW1,VIW2,VIW3,VIW4,VIW5
      DOUBLE PRECISION VIV1,VIV2,VIV3,VIV4,VIV5
      DOUBLE PRECISION VIU1,VIU2,VIU3,VIU4,VIU5
      DOUBLE PRECISION UTEST(9),SUMLEVEL,SUM,SUM2,SUMLEVEL2
      DOUBLE PRECISION BOUNDSUM,BOUNDSUM2
      DOUBLE PRECISION BOUNDSUMLEVEL,BOUNDSUMLEVEL2
      REAL W(IMAX,JMAX,KMAX)

      FLAG=-1.0E+10
      WRITE(6,*)'IMAX,JMAX,KMAX = ',IMAX,JMAX,KMAX
      WRITE(6,*)'SX,SY,SZ = ',SX,SY,SZ
      DX1=SX*1000.
      DY1=SY*1000.
      DZ1=SZ*1000.
      NPTS=0
      NPTSBOUND=0
      SUM=0
      SUM2=0
      BOUNDSUM=0
      BOUNDSUM2=0
      DO K=1,KMAX
       write(6,*)'divergence at level ',k
       NPTSLEVEL=0
       NPTSBOUNDLEVEL=0
       BOUNDSUMLEVEL=0.
       BOUNDSUMLEVEL2=0.
       SUMLEVEL=0
       SUMLEVEL2=0
       DO J=1,JMAX
        DO I=1,IMAX 
c         write(6,*)'i,j,k,u = ',i,j,k,u(i,j,k)
C
C Determine data positions for 6 divergence data points
C
         IF(U(I,J,K).LE.FLAG)GO TO 9999
         IBOUNDX=99
         IBOUNDY=99
         IBOUNDZ=99
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
         IF(K.GE.5)THEN
          UTEST(1)=U(I,J,K-4)
         ELSEIF(K.EQ.4)THEN
          UTEST(1)=0.
         ENDIF
         IF(K.GE.4)THEN
          UTEST(2)=U(I,J,K-3)
         ELSEIF(K.EQ.3)THEN
          UTEST(2)=0.
         ENDIF
         IF(K.GE.3)THEN
          UTEST(3)=U(I,J,K-2)
         ELSEIF(K.EQ.2)THEN
          UTEST(3)=0.
         ENDIF
         IF(K.GE.2)THEN
          UTEST(4)=U(I,J,K-1)
         ELSEIF(K.EQ.1)THEN
          UTEST(4)=0.
         ENDIF
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
         IF(IBOUNDX.EQ.99.OR.IBOUNDY.EQ.99.OR.
     +      IBOUNDZ.EQ.99)GO TO 9999
         IF(IBOUNDX.EQ.0)THEN
          VIU5=-1./12.*DZ1/DX1
          VIU4=2./3.*DZ1/DX1
          VIU3=0.
          VIU2=-2./3.*DZ1/DX1
          VIU1=1./12.*DZ1/DX1
          I5=I+2
          I4=I+1
          I3=I
          I2=I-1
          I1=I-2
         ELSEIF(IBOUNDX.EQ.1)THEN
          VIU5=1./12.*DZ1/DX1
          VIU4=-.5*DZ1/DX1
          VIU3=3./2.*DZ1/DX1
          VIU2=-5./6.*DZ1/DX1
          VIU1=-.25*DZ1/DX1
          I5=I+3
          I4=I+2
          I3=I+1
          I2=I
          I1=I-1
         ELSEIF(IBOUNDX.EQ.-1)THEN
          VIU5=.25*DZ1/DX1
          VIU4=5./6.*DZ1/DX1
          VIU3=-3./2.*DZ1/DX1
          VIU2=.5*DZ1/DX1
          VIU1=-1./12.*DZ1/DX1
          I5=I+1
          I4=I
          I3=I-1
          I2=I-2
          I1=I-3
         ELSEIF(IBOUNDX.EQ.2)THEN
          VIU5=-.25*DZ1/DX1
          VIU4=4./3.*DZ1/DX1
          VIU3=-3.*DZ1/DX1
          VIU2=4.*DZ1/DX1
          VIU1=-25./12.*DZ1/DX1
          I5=I+4
          I4=I+3
          I3=I+2
          I2=I+1
          I1=I
         ELSEIF(IBOUNDX.EQ.-2)THEN
          VIU5=25./12.*DZ1/DX1
          VIU4=-4.*DZ1/DX1
          VIU3=3.*DZ1/DX1
          VIU2=-4./3.*DZ1/DX1
          VIU1=.25*DZ1/DX1
          I5=I
          I4=I-1
          I3=I-2
          I2=I-3
          I1=I-4
         ENDIF
         IF(IBOUNDY.EQ.0)THEN
          VIV5=-1./12.*DZ1/DY1
          VIV4=2./3.*DZ1/DY1
          VIV3=0.
          VIV2=-2./3.*DZ1/DY1
          VIV1=1./12.*DZ1/DY1
          J5=J+2
          J4=J+1
          J3=J
          J2=J-1
          J1=J-2
         ELSEIF(IBOUNDY.EQ.1)THEN
          VIV5=1./12.*DZ1/DY1
          VIV4=-.5*DZ1/DY1
          VIV3=3./2.*DZ1/DY1
          VIV2=-5./6.*DZ1/DY1
          VIV1=-.25*DZ1/DY1
          J5=J+3
          J4=J+2
          J3=J+1
          J2=J
          J1=J-1
         ELSEIF(IBOUNDY.EQ.-1)THEN
          VIV5=.25*DZ1/DY1
          VIV4=5./6.*DZ1/DY1
          VIV3=-3./2.*DZ1/DY1
          VIV2=.5*DZ1/DY1
          VIV1=-1./12.*DZ1/DY1
          J5=J+1
          J4=J
          J3=J-1
          J2=J-2
          J1=J-3
         ELSEIF(IBOUNDY.EQ.2)THEN
          VIV5=-.25*DZ1/DY1
          VIV4=4./3.*DZ1/DY1
          VIV3=-3.*DZ1/DY1
          VIV2=4.*DZ1/DY1
          VIV1=-25./12.*DZ1/DY1
          J5=J+4
          J4=J+3
          J3=J+2
          J2=J+1
          J1=J
         ELSEIF(IBOUNDY.EQ.-2)THEN
          VIV5=25./12.*DZ1/DY1
          VIV4=-4.*DZ1/DY1
          VIV3=3.*DZ1/DY1
          VIV2=-4./3.*DZ1/DY1
          VIV1=.25*DZ1/DY1
          J5=J
          J4=J-1
          J3=J-2
          J2=J-3
          J1=J-4
         ENDIF
         IF(IBOUNDZ.EQ.0)THEN
          VIW5=-1./12.
          VIW4=2./3.
          VIW3=0.
          VIW2=-2./3.
          VIW1=1./12.
          K5=K+2
          K4=K+1
          K3=K
          K2=K-1
          K1=K-2
         ELSEIF(IBOUNDZ.EQ.1)THEN
          VIW5=1./12.
          VIW4=-.5
          VIW3=3./2.
          VIW2=-5./6.
          VIW1=-.25
          K5=K+3
          K4=K+2
          K3=K+1
          K2=K
          K1=K-1
         ELSEIF(IBOUNDZ.EQ.-1)THEN
          VIW5=.25
          VIW4=5./6.
          VIW3=-3./2.
          VIW2=.5
          VIW1=-1./12.
          K5=K+1
          K4=K
          K3=K-1
          K2=K-2
          K1=K-3
         ELSEIF(IBOUNDZ.EQ.2)THEN
          VIW5=-.25
          VIW4=4./3.
          VIW3=-3.
          VIW2=4.
          VIW1=-25./12.
          K5=K+4
          K4=K+3
          K3=K+2
          K2=K+1
          K1=K
         ELSEIF(IBOUNDZ.EQ.-2)THEN
          VIW5=25./12.
          VIW4=-4.
          VIW3=3.
          VIW2=-4./3.
          VIW1=.25
          K5=K
          K4=K-1
          K3=K-2
          K2=K-3
          K1=K-4
         ENDIF
         DIVERGENCE=0.
         DIVERGENCEW=0
         IF(K1.GT.0)THEN
          DIVERGENCEW=(VIW1*W(I,J,K1)+VIW2*W(I,J,K2)+VIW3*W(I,J,K3)+
     +              VIW4*W(I,J,K4)+VIW5*W(I,J,K5))/DZ1
         ELSE
          DIVERGENCEW=(VIW2*W(I,J,K2)+VIW3*W(I,J,K3)+
     +              VIW4*W(I,J,K4)+VIW5*W(I,J,K5))/DZ1
         ENDIF
         DIVERGENCEV=(VIV1*V(I,J1,K)+VIV2*V(I,J2,K)+VIV3*V(I,J3,K)+
     +              VIV4*V(I,J4,K)+VIV5*V(I,J5,K))/DZ1
         DIVERGENCEU=(VIU1*U(I1,J,K)+VIU2*U(I2,J,K)+VIU3*U(I3,J,K)+
     +               VIU4*U(I4,J,K)+VIU5*U(I5,J,K))/DZ1
c         WRITE(6,*)'I,J,K,DIVERGENCEU,DIVERGENCEV,DIVERGENCEW = ',
c     +    I,J,K,DIVERGENCEU,DIVERGENCEV,DIVERGENCEW
         DIVERGENCE=DIVERGENCEU+DIVERGENCEV+DIVERGENCEW
c         WRITE(6,*)'I,J,K,DIVERGENCE = ',I,J,K,DIVERGENCE
c         WRITE(6,*)'K1,K2,K3,K4,K5 = ',K1,K2,K3,K4,K5
c         WRITE(6,*)'W = ',W(I,J,K2),W(I,J,K3),W(I,J,K4),W(I,J,K5)
c         WRITE(6,*)'VIW = ',VIW2,VIW3,VIW4,VIW5
c         WRITE(6,*)'V = ',V(I,J1,K),V(I,J2,K),V(I,J3,K),V(I,J4,K),
c     +     V(I,J5,K)
         DIV(I,J,K)=DIVERGENCE
         SUM=SUM+DIVERGENCE
         SUM2=SUM2+DIVERGENCE*DIVERGENCE
         SUMLEVEL=SUMLEVEL+DIVERGENCE
         SUMLEVEL2=SUMLEVEL2+DIVERGENCE*DIVERGENCE
         NPTS=NPTS+1
         NPTSLEVEL=NPTSLEVEL+1
         IF(IBOUNDX.EQ.0.AND.IBOUNDY.EQ.0.AND.IBOUNDZ.EQ.0)THEN
          NPTSBOUND=NPTSBOUND+1
          NPTSBOUNDLEVEL=NPTSBOUNDLEVEL+1
          BOUNDSUM=BOUNDSUM+DIVERGENCE
          BOUNDSUM2=BOUNDSUM2+DIVERGENCE*DIVERGENCE
          BOUNDSUMLEVEL=BOUNDSUMLEVEL+DIVERGENCE
          BOUNDSUMLEVEL2=BOUNDSUMLEVEL2+DIVERGENCE*DIVERGENCE
         ENDIF
9999    ENDDO
       ENDDO
       IF(NPTSLEVEL.GT.0)THEN
        DIVERGENCE=SUMLEVEL/NPTSLEVEL
       ELSE
        DIVERGENCE=FLAG
       ENDIF
       WRITE(6,*)'MEAN DIVERGENCE AT LEVEL ',K,' IS ',DIVERGENCE
       IF(NPTSLEVEL.GT.0)THEN
        DIVERGENCE=DSQRT(SUMLEVEL2/NPTSLEVEL)
       ELSE
        DIVERGENCE=FLAG
       ENDIF
       WRITE(6,*)'RMS DIVERGENCE AT LEVEL ',K,' IS ',DIVERGENCE
       IF(NPTSBOUNDLEVEL.GT.0)THEN
        DIVERGENCE=BOUNDSUMLEVEL/NPTSBOUNDLEVEL
       ELSE
        DIVERGENCE=FLAG
       ENDIF
       WRITE(6,*)'MEAN INNER DIVERGENCE AT LEVEL ',K,' IS ',
     +   DIVERGENCE
       IF(NPTSBOUNDLEVEL.GT.0)THEN
        DIVERGENCE=DSQRT(BOUNDSUMLEVEL2/NPTSBOUNDLEVEL)
       ELSE
        DIVERGENCE=FLAG
       ENDIF
       WRITE(6,*)'RMS INNER DIVERGENCE AT LEVEL ',K,' IS ',
     +   DIVERGENCE
      ENDDO
      DIVERGENCE=SUM/NPTS
      WRITE(6,*)'MEAN VOLUME DIVERGENCE IS ',DIVERGENCE
      DIVERGENCE=DSQRT(SUM2/NPTS)
      WRITE(6,*)'RMS VOLUME DIVERGENCE IS ',DIVERGENCE
      DIVERGENCE=BOUNDSUM/NPTSBOUND
      WRITE(6,*)'MEAN INNER VOLUME DIVERGENCE IS ',DIVERGENCE
      DIVERGENCE=DSQRT(BOUNDSUM2/NPTSBOUND)
      WRITE(6,*)'RMS INNER VOLUME DIVERGENCE IS ',DIVERGENCE
      RETURN
      END

