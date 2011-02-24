      SUBROUTINE DATEE(DATVAL)
C
C     RETURNS THE DATE IN ASCII FORM MM/DD/YY
C
      CHARACTER*8 DATVAL
      CHARACTER*15 DATIM

      CALL CDATE(DATIM)

      WRITE(DATVAL,20)DATIM(5:6),DATIM(7:8),DATIM(3:4)
 20   FORMAT(A2,'/',A2,'/',A2)

      RETURN

      END

      SUBROUTINE CLOCK(TIMVAL)
C
C     RETURNS THE CURRENT TIME IN ASCII FORM HH:MM:SS
C
      CHARACTER*8 TIMVAL
      CHARACTER*15 DATIM

      CALL CDATE(DATIM)


      WRITE(TIMVAL,20)DATIM(9:10),DATIM(11:12),DATIM(13:14)
 20   FORMAT(A2,':',A2,':',A2)

      RETURN

      END


      FUNCTION CVMGT(TRANS,FLANS,LTEST)
C
C     THIS FUNCTION RETURNS TRANS IF LTEST IS TRUE
C                   RETURNS FLANS IF LTEST IS FALSE
C
      LOGICAL LTEST
C      INTEGER TRANS,FLANS

      IF (LTEST) THEN
         CVMGT=TRANS
      ELSE
         CVMGT=FLANS
      END IF

      RETURN

      END


      INTEGER FUNCTION CVMGP(I1,I2,I3)
C
C     RETURNS I1 IF I3 >= 0
C     RETURNS I2 IF I3 <  0
C
      IF (I3.GE.0.0) THEN
         CVMGP=I1
      ELSE
         CVMGP=I2
      END IF

      RETURN

      END

      INTEGER FUNCTION CVMGZ(I1,I2,I3)
C
C     RETURNS I1 IF I3 = 0.0
C     RETURNS I2 IF I3 != 0.0
C
      IF (I3.EQ.0.0) THEN
         CVMGZ=I1
      ELSE
         CVMGZ=I2
      END IF

      RETURN

      END

      INTEGER FUNCTION CVMGN(I1,I2,I3)
C
C     RETURNS I1 IF I3 != 0.0
C     RETURNS I2 IF I3 =  0.0
C
      IF (I3.NE.0.0) THEN
         CVMGN=I1
      ELSE
         CVMGN=I2
      END IF

      RETURN

      END
