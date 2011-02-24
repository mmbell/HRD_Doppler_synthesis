      subroutine byte_swap_i4(data,max)
c********0*********0*********0*********0*********0*********0*********
c
c Byte swapping for an integer*4 data
c    
c     Allen Zhao                             May 2002
c
c********0*********0*********0*********0*********0*********0*********
c  This is an unused routine for a unix compilation compatible 
c  with the programs and make files
      integer*4 data(max)
      data(2)=data(1)
      return
      end
      
