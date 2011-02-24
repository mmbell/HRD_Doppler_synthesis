      subroutine byte_swap_i2(data,max)
c********0*********0*********0*********0*********0*********0*********
c
c Byte swapping for an integer*2 data
c    
c     Allen Zhao                             May 2002
c
c********0*********0*********0*********0*********0*********0*********
c
      implicit none
c
      integer max
      integer n

      integer*2 data(max)
      integer*2 i2
      integer*1 i1(2)
      integer*1 byte
      equivalence (i1(1),i2)
c
c********0*********0*********0*********0*********0*********0*********
c
      do n=1,max
         i2=data(n)
         byte=i1(2)
         i1(2)=i1(1)
         i1(1)=byte
         data(n)=i2
      enddo
c
c********0*********0*********0*********0*********0*********0*********
c
      return
      end
      
