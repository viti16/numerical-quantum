      PROGRAM idsafj
       implicit none

       real(kind=8)::f,y,t,k1,k2,k3,k4
       integer(kind=4)::i,j,k
      
       f(t,y) = t/y

       y=1.0d0
       k1=0;k2=0;k3=0;k4=0

       open(10,file="rkdat.txt")
  
       do i=1,200
          t = 1.0d0 + (i-1)*0.1d0
          k1=f(t,y)
          k2=f(t+0.05d0,y+0.05d0*k1)
          k3=f(t+0.05d0,y+0.05d0*k2)
          k4=f(t+0.1d0,y+0.1d0*k3) 

          write(10,*)t,y
          y = y + (0.1d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
          
       enddo
       close(10)

       end

