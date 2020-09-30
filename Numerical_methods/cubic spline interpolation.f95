      Program spline_interpolation
      implicit none

      real(kind=8),allocatable::finaldat(:,:),grid(:)
      real(kind=8)::xy(5,5),x(5),y(5),b(4),c(5),d(4),s(4,4),h(4),f(4),e(4,4)
      real(kind=8)::g(4)
      integer(kind=4)::i,j,steps,pivot(4),ok
        
      !defining x and f(x)
      do i=1,5
         x(i)=i
         y(i)=i**2
      enddo
 
      !calculating h
      do i=1,4
         h(i)=x(i+1)-x(i)
      enddo

      !calculating f
      do i=1,4
         f(i)=y(i+1)-y(i)
      enddo 

      !calculating e
      do i=1,4
         do j=1,4
            e(i,j)=0
            e(1,1)=1
            if(i .ge. 2) then
              e(i,i)=2*(h(i-1)+h(i))
              if(i-j .eq. 1) then
                e(i,j)=h(i-1)
              endif
              if(i-j .eq. -1) then
                e(i,j)=h(i)
              endif
            endif
         enddo
      enddo

      !calculate g
      do i=1,4
         g(1)=0
         if(i .ge. 2) then
           g(i)=3*(f(i)/h(i) - f(i-1)/h(i-1))
         endif
      enddo 
 

      ! find the solution using the LAPACK routine DGESV
      call DGESV(4,1,e,4,pivot,g,4,ok)

      !calculating c
      do i=1,4
         c(i)=g(i)
         c(5)=0
      enddo

      !calculating d
      do i=1,4
         d(i)=(c(i+1)-c(i))/3/h(i)
      enddo

      !calculating b
      do i=1,4
         b(i)=(y(i+1)-y(i))/h(i) - (2*c(i)+c(i+1))*h(i)/3
      enddo
      
      !writing the spline coeffecient into array
      do i=1,4
         s(i,1)=y(i)
         s(i,2)=b(i)
         s(i,3)=c(i)
         s(i,4)=d(i)
      enddo
      
      !setting grid for interpolation 
      steps=(x(5)-x(1))/0.1
      allocate(grid(steps+1))
      do i=1,steps+1
         grid(i)=x(1)+0.1*(i-1)
      enddo 
      
      !calculating interpolated values from grid from the splines  
      allocate(finaldat(steps+1,2))
      do i=1,steps+1
            if(x(1) .le. grid(i) .and. grid(i) .lt. x(2)) then
               finaldat(i,1)=grid(i)
               finaldat(i,2)=s(1,1)+s(1,2)*(grid(i)-x(1))+s(1,3)*(grid(i)-x(1))**2+s(1,4)*(grid(i)-x(1))**3 
            endif 
            if(x(2) .le. grid(i) .and. grid(i) .lt. x(3)) then
               finaldat(i,1)=grid(i)
               finaldat(i,2)=s(2,1)+s(2,2)*(grid(i)-x(2))+s(2,3)*(grid(i)-x(2))**2+s(2,4)*(grid(i)-x(2))**3 
            endif
            if(x(3) .le. grid(i) .and. grid(i) .lt. x(4)) then
               finaldat(i,1)=grid(i)
               finaldat(i,2)=s(3,1)+s(3,2)*(grid(i)-x(3))+s(3,3)*(grid(i)-x(3))**2+s(3,4)*(grid(i)-x(3))**3 
            endif
            if(x(4) .le. grid(i) .and. grid(i) .lt. x(5)) then
               finaldat(i,1)=grid(i)
               finaldat(i,2)=s(4,1)+s(4,2)*(grid(i)-x(4))+s(4,3)*(grid(i)-x(4))**2+s(4,4)*(grid(i)-x(4))**3 
            endif
      enddo
      open(30,file="spli.txt")
      do i=1,steps+1 
         write(30,*)(finaldat(i,j),j=1,2)
      enddo
      close(30)    
      end
