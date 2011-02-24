      subroutine byte_swap_i4(data,max)
c********0*********0*********0*********0*********0*********0*********
c
c Byte swapping for an integer*4 data
c    
c     Allen Zhao                             May 2002
c
c********0*********0*********0*********0*********0*********0*********
c
      implicit none
c
      integer max
      integer n

      integer*4 data(max)
      integer*4 i4
      integer*1 i1(4)
      integer*1 byte
      equivalence (i1(1),i4)
c
c********0*********0*********0*********0*********0*********0*********
c
c      write(6,*)'i1 = ',i1
      do n=1,max
         i4=data(n)
         byte=i1(2)
         i1(2)=i1(3)
         i1(3)=byte
         byte=i1(4)
         i1(4)=i1(1)
         i1(1)=byte
         data(n)=i4
      enddo
c      write(6,*)'i1 = ',i1
c
c********0*********0*********0*********0*********0*********0*********
c
      return
      end
      
