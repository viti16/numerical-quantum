       PROGRAM integration
       implicit none

       real(kind=8),allocatable::p(:,:),eig(:,:),w(:),work(:),wt(:)
       integer(kind=8)::i,j,k,n,info,lwork,lda
       real(kind=8)::x,f,fint,a,b,s,t

       !points to calculate
       write(*,*)"enter points to compute"
       read(*,*)n
       
       
       !calculating legendre matrix
       allocate(p(n,n))
       do i=1,n
          do j=1,n
             if(i-j .eq. -1) then
               p(i,j)=i/((2.0d0*i-1.0d0)*(2.0d0*i+1.0d0))**0.50d0 
               p(j,i)=p(i,j)
             endif
             if(abs(i-j) .ne. 1) then
               p(i,j)=0
             endif
          enddo 
       enddo    
       
       !calulating eigen values and vectors
       allocate(eig(n,n))
       do i=1,n
          do j=1,n
             eig(i,j)=p(i,j)
          enddo
       enddo
       allocate(w(n))
       allocate(work(1000))
       lwork=-1
       lda=n
       call dsyev('Vectors','Upper',n,eig,lda,w,work,lwork,info)
       lwork=min(1000,int(work(1)))
        
       call dsyev('Vectors','Upper',n,eig,lda,w,work,lwork,info)
        
       
       !calculating weights and function value at nodes
       allocate(wt(n))
       wt=0 
       do i=1,n
          wt(i)=2*(eig(1,i)**2)
       enddo   


       open(10,file="intdat.txt")
       do i=1,200
          x=-40+i*0.4
          do j=1,n
             do k=1,n
                t=((3.14d0)**2)*w(k) + (3.14d0)**2
                s=(50*w(j) + 50) 
                f=wt(k)*(3.14)*wt(j)*(50)*1/((3.14d0)**0.5)*(s**(-0.5))*exp(-1*s*(x-cos(t))**2)
                fint=f+fint
             enddo
          enddo
          write(10,*)x,fint
       enddo


       close(10)

       end


