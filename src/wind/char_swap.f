      subroutine char_swap(data,linux,text)
      integer*2 data,i2
      integer*1 i1(2)
      character*2 text
      logical linux
      equivalence (i1(1),i2)
      i2=data
      if (linux) then
      text(1:1)=char(i1(2))
      text(2:2)=char(i1(1))
      else
      text(1:1)=char(i1(1))
      text(2:2)=char(i1(2))
      endif
      return
      end
