      SUBROUTINE DIFFOBLONG(XIN,YIN,Z,XDIFFMAX,YDIFFMAX,ZDIFFMAX,
     + SXB,SYB,SZB,
     + XZB,YZB,ZZB,IMAXB,JMAXB,KMAXB,
     + IMN,IMX,JMN,JMX,KMN,KMX,ROTB)

      IMPLICIT NONE
      REAL X, Y, Z, XB, YB, XIN, YIN, ROTB
      REAL XMINREL, YMINREL, ZMINREL
      REAL XMAXREL, YMAXREL, ZMAXREL
      REAL XDIFFMAX, YDIFFMAX, ZDIFFMAX
      REAL SXB,SYB,SZB,XZB,YZB,ZZB
      INTEGER IMN, JMN, KMN, IMX, JMX, KMX, IMAXB, JMAXB, KMAXB
      
      X=0.
      Y=0.
      XB=0.
      YB=0.
      CALL ROTATE(XIN,YIN,XB,YB,-ROTB)
      X=XB
      Y=YB
      XMINREL=X-XDIFFMAX
      YMINREL=Y-YDIFFMAX
      ZMINREL=Z-ZDIFFMAX
      XMAXREL=X+XDIFFMAX
      YMAXREL=Y+YDIFFMAX
      ZMAXREL=Z+ZDIFFMAX
      IMN=INT((XMINREL+XZB)/SXB+.5)
      JMN=INT((YMINREL+YZB)/SYB+.5)
      KMN=1+INT((ZMINREL-ZZB)/SZB)
      IF(IMN.LT.1)IMN=1
      IF(JMN.LT.1)JMN=1
      IF(KMN.LT.1)KMN=1
      IF(IMN.GT.IMAXB)IMN=IMAXB+1
      IF(JMN.GT.JMAXB)JMN=JMAXB+1
      IF(KMN.GT.KMAXB)KMN=KMAXB+1
      IMX=INT((XMAXREL+XZB)/SXB+.5)
      JMX=INT((YMAXREL+YZB)/SYB+.5)
      KMX=1+INT((ZMAXREL-ZZB)/SZB)
      IF(IMX.LT.1)IMX=0
      IF(JMX.LT.1)JMX=0
      IF(KMX.LT.1)KMX=0
      IF(IMX.GT.IMAXB)IMX=IMAXB
      IF(JMX.GT.JMAXB)JMX=JMAXB
      IF(KMX.GT.KMAXB)KMX=KMAXB

      RETURN
      END


      SUBROUTINE ROTATE(X,Y,XROT,YROT,ROT)
      IMPLICIT NONE
      REAL ROT, XROT, YROT, X, Y, ROTD

      ROTD=3.14159/180.*ROT
      XROT=X*COS(ROTD)+Y*SIN(ROTD)
      YROT=-X*SIN(ROTD)+Y*COS(ROTD)

      RETURN
      END


      SUBROUTINE DIFFDIST(X,Y,Z,IB,JB,KB,INEAR,JNEAR,KNEAR,
     + SXB,SYB,SZB,
     + XZB,YZB,ZZB,DISTX,DISTY,DISTZ,ROTB)
      ROTBD=-3.14159/180.*ROTB
      A=COS(ROTBD)
      B=SIN(ROTBD)
C      WRITE(100,*)'ROTB,ROTBD,COS,SIN = ',ROTB,ROTBD,A,B
      XBROT=A*X+B*Y
      YBROT=-B*X+A*Y
      ZBROT=(KB-1)*SZB+ZZB
      XB=(FLOAT(IB)-.5)*SXB-XZB
      YB=(FLOAT(JB)-.5)*SYB-YZB
      ZB=(FLOAT(KB)-1)*SZB+ZZB
      INEAR=1+INT((XBROT+XZB)/SXB)
      JNEAR=1+INT((YBROT+YZB)/SYB)
      KNEAR=1+INT((Z-ZZB)/SZB)
      DISTX=XB-XBROT
      DISTY=YB-YBROT
      DISTZ=Z-ZBROT
C      WRITE(100,*)'IB,JB,KB,XB,YB,ZB = ',IB,JB,KB,XB,YB,ZB
C      WRITE(100,*)'XBROT,YBROT,ZBROT = ',XBROT,YBROT,ZBROT

      RETURN
      END


      SUBROUTINE DIVER(U,V,DIV,IMAX,JMAX,KMAX,SX,SY)
      REAL U(IMAX,JMAX,KMAX),V(IMAX,JMAX,KMAX),DIV(IMAX,JMAX,KMAX)
      FLAG=-1.0E+10
      DO 1 K=1,KMAX                                                             
c        WRITE(6,'("- computing divergence for plane #",i3)') K  
        DO 2 J=1,JMAX                                                           
          DO 3 I=1,IMAX                                                         
            DIV(I,J,K)=FLAG
            SCALX=SX*2000.0      !TWICE THE RES IN M                            
            SCALY=SY*2000.0      !TWICE THE RES IN M                            
            IF(I.EQ.1) THEN                                                     
              U1=U(I,J,K)
              SCALX=SCALX*.5                                                    
            ELSE                                                                
              U1=U(I-1,J,K)
            ENDIF                                                               
            IF(I.EQ.IMAX) THEN                                                  
              U2=U(I,J,K)
              SCALX=SCALX*.5                                                    
            ELSE                                                                
              U2=U(I+1,J,K)
            ENDIF                                                               
            IF(J.EQ.1) THEN                                                     
              V1=V(I,J,K)
              SCALY=SCALY*.5                                                    
            ELSE                                                                
              V1=V(I,J-1,K)
            ENDIF                                                               
            IF(J.EQ.JMAX) THEN                                                  
              V2=V(I,J,K)
              SCALY=SCALY*.5                                                    
            ELSE                                                                
              V2=V(I,J+1,K)
            ENDIF                                                               
            IF(U1.GT.FLAG.AND.U2.GT.FLAG.AND.                                 
     1       V1.GT.FLAG.AND.V2.GT.FLAG)THEN                                   
             DIV(I,J,K)=(U2-U1)/SCALX + (V2-V1)/SCALY
             IF (ABS(DIV(I,J,K)).GT.0.32767) 
     1         DIV(I,J,K)=SIGN(0.32767,DIV(I,J,K))  
            ELSE                                                                
             DIV(I,J,K)=FLAG
            ENDIF                                                               
C            IF(IDATA(5,I,J,K).GE.32767)
C     1       write(6,*)'I,J,K,V1,V2,U1,U2',I,J,K,V1,V2,U1,U2
3           CONTINUE                                                            
2         CONTINUE                                                              
1       CONTINUE                                                                
      RETURN                                                                    
      END                                                                       
      FUNCTION VTERM_NEW(Z,H,hb,dpb,IRSW,ZLOW,ZHIGH) 
C Terminal velocity from dBZ and height 
C This function computes mean terminal velocity from the reflectivity 
C according to Paul Willis' 2-parameter gamma distribution and
C a snow relationship developed by Atlas et al. (1973)
C Reflectivity, in terms of dBZ, must be passed to this routine.
C HEIGHT MUST BE IN KM
c irsw - 0 Joss and Waldvogel; >0 Willis Gamma
C hb is the height of the bright band, and dpb is the depth of the bright 
C   band both in km 
C 
      ZZ=10.0**(Z*0.1)
      hlow= hb - dpb * .5 
      hhi= hlow + dpb 
C density correction term (rhoo/rho)*0.45 [rho(Z)=rhoo exp-(z/H), where 
C  H is the scale height = 9.58125 from Gray's inner 2 deg composite] 
C 0.45 density correction from Beard (1985, JOAT pp 468-471)
      DCOR=EXP(0.45*H*.10437052)
C The snow relationship (Atlas et al., 1973) --- VT=0.817*Z**0.063  (m/s) 
      VTS=-DCOR * (0.817*ZZ**0.063) 
      if(irsw.gt.0) then
C The rain relationship --- from Willis analytical-gamma distribution 
         TERM1=7.331/ZZ**0.010022 
         TERM2=0.14034*ZZ**0.095238 
         VTR=-DCOR * (5.5011E+09/(TERM1+TERM2)**10.5) 
      else
C The rain relationship (Joss and Waldvogel,1971) --- VT=2.6*Z**.107 (m/s)
         VTR=-DCOR * (2.6*ZZ**.107) 
      endif 
C test if height is in the transition region between SNOW and RAIN
C  defined as hlow in km < H < hhi in km
C  if in the transition region do a linear weight of VTR and VTS
      IF(Z.GT.ZLOW.AND.Z.LE.ZHIGH)THEN
       WEIGHTR=(Z-ZLOW)/(ZHIGH-ZLOW)
       WEIGHTS=1.-WEIGHTR
       VTS=(VTR*WEIGHTR+VTS*WEIGHTS)/(WEIGHTR+WEIGHTS)
      ELSEIF(Z.GT.ZHIGH)THEN
       VTS=VTR
      ENDIF
      VTERM_NEW=VTR*(hhi-H)/dpb + VTS*(H-hlow)/dpb  
      IF(H.LT.hlow) VTERM_NEW=VTR 
      IF(H.GT.hhi) VTERM_NEW=VTS

      RETURN
      END 

