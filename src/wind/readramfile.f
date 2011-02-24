      subroutine readramfile(timein,time_ram,rlat,rlon,ra,wd,ws,
     *                       ucom,vcom,vws,wg,ios)
C
C This subroutine reads the RAM files and returns the following parameters:
C      time_ram - time from file (sec)  {real}
C      rlat - radar latitude (deg)      {real}
C      rlon - radar longitude(deg)      {real}
C      ra - radar altitude (meters)     {real}
C      wd - wind direction (deg)        {real}
C      ws - wind speed (m/s)            {real}
C      ucom - E/W wind speed (m/s)      {real}
C      vcom - N/S wind speed (m/s)      {real}
C      vws - vertical wind speed (m/s)  {real}
C      wg - vertical ground speed (m/s) {real}
C      ios - error condition: 0 good read {integer}
C                            -1 time gap on ram file
C                            -2 end of file or timein is not found on rams
C                            >0 HP run time eror condition
C timein - time requested in seconds;send by the calling routine
C
C stime_ram - start time of the ram file
C etime_ram - end time of the ram file, 48 hour clock.
C rate - data rate in seconds
C reclen - record length of the file
C fileid - file identification name, 'ram1'
C maxrec - last record number
C irec - record number
C ngoodrec - next good record after a time jump
C lgoodrec - last good record before a time jump.
C             FILE FORMAT
C first record:
C             stime_ram,etime_ram,rate,reclen,'ram1',0,0
C good read:
C             time_ram,rlat,rlon,ra,u,v,w,wg,extra
C time jump:
C             time_ram,-999.0,ngoodrec,lgoodrec,-999,-999,-999,-999
C last record:
C             -999.,0.,0.,0,0,0,0,0,0
C
      integer*2 idat(6)
      character nameram*(*), jfile*60, fileid*4
      save stime_ram, etime_ram, rate, reclen, fileid, maxrec, irec, lp
      equivalence (dat3,next), (idat(1),last),
     *            (idat(1),reclen), (idat(3),fileid)    
      data itimein_old/-99/
C
      time_send=timein
      itime_send=timein
      irec= nint((time_send - stime_ram)/rate) + 2 !calculate record number.
      ios= 0
C add 24 hours to the time and try again:
      if (irec.lt.0) then
         time_send=time_send + 86400.0
         irec=nint((time_send - stime_ram)/rate) + 2
      endif
C return if the time does not match a record number:
      if (irec.lt.2 .or. irec.gt.maxrec .or. timein.lt.0.0) then
         call ctme(timein,ihs,imn,isc)
         write(lp,'(" time= ",3i2.2," not on rams")') ihs,imn,isc
         ios= -2
         return
      endif
C read file:
10    read(92,err=99,iostat=iose,rec=irec)
     +        time_ram,rlat,dat3,(idat(i),i=1,6)
      if (time_ram .le. -999.) go to 90
      if (rlat .le. -999.) then            !time gap.
         lgoodrec=last
         ngoodrec=next
         read(92,err=99,iostat=iose,rec=lgoodrec)
     +           time_ram,rlat,dat3,(idat(i),i=1,6)
         diflast=abs(time_send - time_ram)
         read(92,err=99,iostat=iose,rec=ngoodrec)
     +           time_ram,rlat,dat3,(idat(i),i=1,6)
         difnext=abs(time_send - time_ram)
         call ctme(time_ram,ihs,imn,isc)
         call ctme(time_send,ih,im,is)
         if (itimein_old.lt.itime_send .or. itimein_old.gt.itime_send)
     +   write(lp,'(1x,"time jump on ram:",3i2.2," tape time=",3i2.2)')
     +        ihs,imn,isc,ih,im,is
         ios= -1
         if (diflast .lt. difnext) irec=lgoodrec
         if (difnext .le. diflast) irec=ngoodrec
         go to 10
      endif
      itimein_old=itime_send
      rlon=dat3
      ra=idat(1)
      ucom=real(idat(2))/100.0
      vcom=real(idat(3))/100.0
      vws=real(idat(4))/100.0
      wg=real(idat(5))/100.0
      wd=amod(630.0-57.29578*atan2(vv,uu),360.0)
      ws=sqrt(uu*uu+vv*vv)
      return
C
C open file and read first record
      entry openramfile(nameram,lup,stime,etime,datarate,ierr)
C This section opens a ram file.  It searches for files on /hrd/dat/.
C nameram - ram file name; format: sub-directory/name.ram no need to specify
C           /hrd/dat/ if the file resides there.  {input, character}
C lup - logical unit number to display messages   {input, integer}
C stime - start time of the file                  {output, real}
C etime - end time of the file                    {output, real}
C datarate - file data rate                       {output, real}
C ierr - error condition flag: 0 - no errors      {integer}
C                             >0 - error condition exists.
C                             <0 - not a ram file
C logical unit 92 is connected to the file.
      lp=lup
      open(92,err=100,file=nameram,iostat=ierr,status='old',
     +     access='direct',recl=24)
100   if (ierr .ne. 0) then
         jfile='/hrd/dat/' // nameram
         open(92,err=200,file=jfile,iostat=ierr,status='old',
     +        access='direct',recl=24)
      endif
      read(92,err=200,iostat=ierr,rec=1)
     *     stime_ram,etime_ram,rate,reclen,fileid
      if (fileid .eq. 'ram1') then
         write(lp,'(" ram file opened: ",a)') nameram
         call ctme(stime_ram,ihs,imn,isc)
         call ctme(etime_ram,ih,im,is)
         write(lp,'(" start,end & rate: ",3i2.2,1x,3i2.2,1x,f5.0)')
     +         ihs,imn,isc,ih,im,is,rate
         stime=stime_ram
         etime=etime_ram
         datarate=rate
         maxrec= ((etime_ram - stime_ram)/rate) + 2
      else
         ierr= -1
      endif
      return
C
C close ram file.
      entry closeramfile
      close(92)
      return
c ***** error messages *****
90    write(lp,'(" end of rams file")')
      ios= -2
      return
99    ios=iose
      write(lp,'(" error on reading rams: ",i6)') ios
200   return
      end
      subroutine ctme(time,itimehh,itimemm,itimess)
      itimes=nint(time)
      itimehh=itimes/3600
      iremain=itimes-3600*itimehh
      itimemm=iremain/60
      itimess=iremain-itimemm*60
      return
      end

