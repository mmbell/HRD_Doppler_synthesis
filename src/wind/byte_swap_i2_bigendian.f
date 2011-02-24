      subroutine byte_swap_i2(data,max)
c********0*********0*********0*********0*********0*********0*********
c
c Byte swapping for an integer*2 data
c    
c     Allen Zhao                             May 2002
c
c********0*********0*********0*********0*********0*********0*********
c
c This is a dummy unused routine for a unix compilation of interpolation
c programs
      integer*2 data(max)
      data(1)=1
      data(2)=data(1)
      return
      end
      
