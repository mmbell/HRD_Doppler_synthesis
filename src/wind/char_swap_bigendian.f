      subroutine char_swap(data,linux,text)
c  this is a routine that should not actually be called if LINUXYES=0
c  and should not be used if LINUXYES=1
c  LINUXYES is a parameter in main program indicating whether machine
c  being used is littleendian (LINUXYES=1) or bigendian (LINUXYES=0)
c
      integer*2 data
c      integer*1 i1(2)
      character*2 text
      data=0
      text='  '
      return
      end
