C **********
      subroutine readuf_disc(time,numsweep,flightid,nrecord,rlat,rlon,
     +            radalt,azmgr,elegr,radii,dbzbuf,velbuf,nbins,ierr,
     +            charv,chardbz,numrec,linuxyes)
c 
c This routine reads the Universal Format Tapes from disc and returns
c the following parameters:
c   time -    time in seconds. (real)
c   flightid - 8 character variable containing the project name. (character)
c   nrecord - physical record number relative to begining of file.
c             If the record number > 32767, then the following # is 
c             -32768, -32767, etc.(integer)
c   rlat -    radar latitude in degrees.(real)
c   rlon -    radar longitude in degrees.(real)
c   radalt -  radar altitude in meters.(real)
c   azmgr -   azimuth using the same orientation as ground based radars.(real)
c   elegr -   elev. using the same orientation as ground based radars.(real)
c   radii -   buffer containing the ranges in km for every data point.(real)
c   dbzbuf-   buffer containing dbz values, -1.0e+10 for no data.(real)
c   velbuf-   buffer containing velocities, -1.0e+10 for no data.(real)
c   nbins -   number of data points.(integer)
c   ierr -    error code > 0 record length in bytes (integer)
c                        -1 for end of data.
c                        -2 no optional header block.
c                        -3 no local use header.
c                        -4 no data header.
c
c variables use in the subroutine:
c   ibuf -    buffer read from disc
c   filename -Universal Format disc file name.
c   jfile -   Disc file name with /hrd/dat // filename
c   readlen - record length return from disc
c   fname-    field name: VE for velocities, DZ for dbz
c
      implicit none
      integer*2 ibuf(32768),ihs,iss,ims,idh,istart,iend,scale
      integer*2 nfields,rdel,rdel2,i1,UFASCII
      integer*2 readlen,rl_upper1,rl_upper2,rl_lower1,rl_lower2 
      real radii(1024), dbzbuf(1024), velbuf(1024)
      character flightid*8, fname*2, charv*2, chardbz*2
      character filename*80, jfile*60
      logical first, linux
      integer lp,lup,istat,numsweep,nrecord,nbins,ierr
      integer iremain,itimehh,itimemm,itimess,itimes,i,j,jk,k
      integer NBLK,linuxyes,rlu,rll,numrec
      real time,rlat,rlon,radalt,azmgr,elegr,flag,gtln
      save lp, ibuf
      data first/.true./
c
C Read data from disc:
c      write(6,*)'statement 50 in ufdisc_io.f'
c50    call disk_stream(220,first,ibuf,4,length)
c      NUMREC=1
      if(linuxyes.eq.1)then
       linux=.true.
      else
       linux=.false.
      endif
      FLAG=-1.0E+10
      DO JK=1,1024
       VELBUF(JK)=FLAG
       DBZBUF(JK)=FLAG
      ENDDO
c      write(6,*)'im inside reader'
c      WRITE(6,'(A2)')CHARV
c      WRITE(6,'(A2)')CHARDBZ
      IF(LINUXYES.EQ.1)THEN
c       write(6,*)'about to read length and numrec = ',numrec
       READ(80,REC=numrec,IOSTAT=ISTAT,ERR=99999)rl_upper1
       NUMREC=NUMREC+1
       READ(80,REC=numrec,IOSTAT=ISTAT,ERR=99999)rl_upper2
       CALL BYTE_SWAP_I2(RL_UPPER1,1)
       CALL BYTE_SWAP_I2(RL_UPPER2,1)
       RLU=RL_UPPER1+RL_UPPER2
       NBLK = RLU/2
c       write(6,*)'nblk = ',nblk
       NUMREC=NUMREC+1
c       write(6,*)'numrec1 = ',numrec
       DO I=1,NBLK
        READ(80,REC=NUMREC,ERR=99999) IBUF(I)
        NUMREC=NUMREC+1
c        write(6,*)'rec = ',numrec
       ENDDO
c       write(6,*)'numrec2 = ',numrec
       CALL BYTE_SWAP_I2(IBUF,NBLK)
       NUMREC=NUMREC+1
       READ(80,REC=NUMREC,ERR=99999)RL_LOWER1
       NUMREC=NUMREC+1
       READ(80,REC=NUMREC,ERR=99999)RL_LOWER2
       CALL BYTE_SWAP_I2(RL_LOWER1,1)
       CALL BYTE_SWAP_I2(RL_LOWER2,1)
       RLL=RL_LOWER1+RL_LOWER2
       IF (RLL.ne.RLU) THEN
         write(6,*)'Error in UF read: ',rlu,rll
       ENDIF
      ELSE
       READ(80,IOSTAT=ISTAT,END=99998,err=99999)IBUF(1),IBUF(2),
     +    (IBUF(I),I=3,IBUF(2))
      ENDIF
c      write(6,*)'numrec3 = ',numrec
      ierr=0
      readlen=ibuf(2)*2 - 4

c 
c Get the time, azimuth, elevation, record number, lat and lon from
c the Mandatory Header Block.
      nrecord=ibuf(6)                    !record number
      numsweep=ibuf(10)              
      rlat=real(ibuf(19)) + real(ibuf(20))/60.0 +
     +     (real(ibuf(21))/64.0)/3600.  !lat in degrees
      rlon=real(ibuf(22)) + real(ibuf(23))/60.0 +
     +     (real(ibuf(24))/64.0)/3600.  !lon in degrees
c      write(6,*)'rlat ibuf = ',rlat,ibuf(19),ibuf(20),ibuf(21)
c      write(6,*)'rlon ibuf = ',rlon,ibuf(22),ibuf(23),ibuf(24)
      ihs=ibuf(29)
      ims=ibuf(30)
      iss=ibuf(31)
      time=3600.*ihs+60.*ims+iss
      azmgr=real(ibuf(33))/64.0
      elegr=real(ibuf(34))/64.0
      radalt=real(ibuf(25)) 
c
c Get the Optional Header Block, ierr= -2 if not present
      if (ibuf(3) .eq. ibuf(4)) ierr= -2 !no optional header block. 
c
c Get the local use header, ierr= -3 if not present,
c returns flight id and radar altitude if present
      if (ibuf(4) .eq. ibuf(5)) then     !no local use header.
         ierr= -3
      else
         write(flightid,'(4a2)') (ibuf(i),i=ibuf(4)+5,ibuf(4)+8)
      endif
c
c Get the data headers and field headers 
      if (ibuf(2) .le. ibuf(5)) then     !no data headers
         ierr= -4
      else
         idh=ibuf(5)                      !1st word of data header block
         nfields=ibuf(idh)                !number of fields in the ray
         do i=1, nfields
            if(linuxyes.eq.1)then
             call char_swap(ibuf(1+i*2+idh),linuxyes,fname)
            else
c             write(fname,'(a2)') achar(ibuf(1+(i*2)+idh)) !field name.
              write(fname,'(a2)') ibuf(1+(i*2)+idh)
            endif
            i1=ibuf(2+(i*2)+idh)              !1st word of field header
            scale=ibuf(i1+1)                  !data divided by this factor
            rdel=ibuf(i1+2)                   !range to 1st gate(km)      
            rdel2=ibuf(i1+3)                  !adjustment to rdel (m)
            rdel=rdel+rdel2/1000.             !full range delay (km)
            gtln=real(ibuf(i1+4))/1000.      !sample volume spacing(km)
            nbins=ibuf(i1+5)                  !number of gates
            istart=ibuf(i1)                   !position of 1st data word
            iend=(istart + nbins) - 1         !position of last data word) 
c       write(6,'(a2,x,a2)')fname,charv
c       write(6,*)'i1,scale,rdel,gtln,nbins,istart,iend'
c       write(6,*)i1,scale,rdel,gtln,nbins,istart,iend
c       write(6,*)'ibuf = ',(ibuf(k),k=istart,iend)
c           get the velocities and dbz values
c           set no data to -1.0e+10
c            if (fname .eq. 'VU' .or. fname .eq. 'VE') then         !velocities
c            write(6,'(3A2,X)')fname,charv,chardbz
c               write(6,*)istart,iend
            if (fname .eq. charv) then         !velocities
               j=0
               if(istart.gt.32768.or.iend.gt.32768)then
                write(6,*)'uf buffer too large for velocity'
             write(6,*)'increase dimension of ibuf in ufdisc_io_new.f'     
               endif
               do k=istart, iend
                  j=j + 1
c                  write(6,*)'k,j = ',k,j
                  if (ibuf(k) .eq. -32768) then
                     velbuf(j)= -1.0e+10
                  else
                     velbuf(j)= real(ibuf(k))/scale
                  endif
c                  write(6,*)'k,ibuf(k),j,velbuf(j) = '
c                  write(6,*)k,ibuf(k),j,velbuf(j)
               end do
            else if (fname .eq. chardbz) then    !dbz values
               j=0
c               write(6,*)istart,iend
               if(istart.gt.32768.or.iend.gt.32768)then
                write(6,*)'uf buffer too large for dbz'
             write(6,*)'increase dimension of ibuf in ufdisc_io_new.f'     
               endif
               do k=istart, iend
                  j=j + 1
                  if (ibuf(k) .eq. -32768) then
                     dbzbuf(j)= -1.0e+10
                  else
                     dbzbuf(j)= real(ibuf(k))/scale
                  endif
c                  write(6,*)'k,ibuf(k),j,dbzbuf(j) = '
c                  write(6,*)k,ibuf(k),j,dbzbuf(j)
               end do
            endif
         end do
c        fill ranges for every data points
         radii(1)= rdel 
         if (radii(1).lt.0) then
c            write(6,*)'Bad range!', radii(1)
         endif
         do k=2, 1024
            radii(k)=radii(k-1) + gtln
         end do
      endif
      ITIMES=NINT(TIME)
      ITIMEHH=ITIMES/3600
      IREMAIN=ITIMES-3600*ITIMEHH
      ITIMEMM=IREMAIN/60
      ITIMESS=IREMAIN-ITIMEMM*60
c      if(itimes.gt.200000)then
c       WRITE(6,10103)ITIMEHH,ITIMEMM,ITIMESS
10103 FORMAT('TIME = ',3I2)
c       WRITE(6,*)'LAT,LON,ALT,AZM,ELE',RLAT,RLON,RADALT,AZMGR,ELEGR
c       WRITE(6,*)'RANGE = ',(RADII(ijk),ijk=1,512)
c       WRITE(6,*)'VR = ',(VELBUF(ijk),ijk=1,512)
c       WRITE(6,*)'DBZ = ',(DBZBUF(ijk),ijk=1,512)
c      endif
      return
c **********
c      subroutine openuf_disc(filename,lup,istat)
C This section opens the universal format file.  It searches for files
C on /hrd/dat/.
C filename - file name; format: sub-directory/name.UF no need to specify
C           /hrd/dat/ if the file resides there.  {input, character}
C lup - logical unit number to display messages   {input, integer}
C istat - error condition flag: 0 - no errors      {integer}
C                             >0 - error condition exists.
      lp=lup
      open(80,file=filename,iostat=istat,err=900,status='old')
10    write(lp,'(" UF file is opened.")')
      return
c ***** errors *****
900   jfile=filename
c      open(80,file=jfile,iostat=istat,err=999,status='old')
      go to 10
999   write(lp,'(" error opening file:",i6)') istat
      return
99998 ierr=-1
      return
99999 IF(LINUXYES.EQ.1)THEN
       IERR=-1
      ELSE
       IERR=-5
      ENDIF
c I'll assume that ierr = 5 is the beginning of a new volume
      return
      end

