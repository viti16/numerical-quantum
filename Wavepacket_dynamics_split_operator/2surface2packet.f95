!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   To compute dynamics of 2 wavepackets in 2 Harmonic Potential  !!!!
!!                   Using split opeartor method                   !!!! 
!!                     command to  run the code                    !!!!   
!!        gfortran 2surface2packet.f95   fourn.f95                 !!!!
!!                  ./a.out                                        !!!!
!!    Author: Vishal Tiwari
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                 !-------------------------!
                                   PROGRAM FFT
                                 !-------------------------!
        use   dftfft
                                        
        implicit none
  
        REAL*8, allocatable :: ftdata1(:),ftdata2(:)
        COMPLEX*8,allocatable :: PSIO(:),TEMPVEC(:)
        COMPLEX*8,allocatable :: EVOLVE(:,:),EVOLTE(:,:)
        COMPLEX*8 :: IOTA=(0d0,1d0),XEXP,X2EXP,ENERGY
        COMPLEX*8 :: PE,KE,PEXP,OVERLAP
        REAL*8,allocatable :: PGRID(:),XGRID(:)
        REAL*8 :: HBAROMGA,HBAR,MASS,PI=3.1415926535d0
        REAL*8 :: DX,ALPHA,X1,X2,PO,AC,XC,DP,XP,NORM
        REAL*8 :: DELT,K12,LAMBDA,V1,V2,V12,DD
        INTEGER*4:: i,j,k,l,m,n,o,p
  
        n=int(2**10)
                
        allocate(PSIO(2*n))

        allocate(PGRID(n))
        allocate(XGRID(n))

        allocate(EVOLVE(2*n,2*n))
        allocate(EVOLTE(2*n,2*n))

        allocate(TEMPVEC(2*n))
        allocate(ftdata1(2*n))
        allocate(ftdata2(2*n))
  
        HBAR=0.6582d0; HBAROMGA=0.44d0; MASS=86.65d0
        K12=0.02d0; LAMBDA=0.19d0; V12=LAMBDA
     

!!!!!!   Defining box and starting point of potentials !!!!!!

        AC=-24.4d0; DELT=0.10d0; X1=-8.67d0; X2=5.62d0

!!!!!!   Grid spacing for space and momentum           !!!!!!

        DX=abs(35.6d0+24.4d0)/real(n-1); DP=2d0*pi*HBAR/DX/n

        open(20,file='calculatedparameter.txt')

!!!!!!   Grid for space and momentum            !!!!!!!!!!!!!

        do i=1,n
           XGRID(i)=AC+DX*real(i-1)
        enddo

        do i=1,n/2+1
           PGRID(i)=2.d0*pi*(i-1)/(DX*N)*HBAR
        enddo
        do i=n/2+1,n-1
           PGRID(i+1)=2.d0*pi*(i-n)/(DX*N)*HBAR
        enddo

!!!!!!!!  Wavepacket matrix                   !!!!!!!!!!!!!!!
            
        do i=1,n
              PSIO(i)=(1d0/pi)**(0.25d0) & 
                              *dexp(-0.5d0*(XGRID(i)-X1)**2) 
              PSIO(i+n)=(1d0/pi)**(0.25d0) & 
                              *dexp(-0.5d0*(XGRID(i)-X1)**2) 
        enddo

!!!!!!!   Kinetic energy evolution opeartor matrix    !!!!!!!

        EVOLTE=0.d0
        do i=1,n
           EVOLTE(i,i)=cdexp(-1.d0*IOTA*(PGRID(i)**2)*DELT/4.d0/MASS/HBAR)
           EVOLTE(i+n,i+n)=cdexp(-1.d0*IOTA*(PGRID(i)**2)*DELT/4.d0/MASS/HBAR)
        enddo

!!!!!!!   Potential energy evolution opeartor matrix    !!!!!!!
       
        EVOLVE=0.d0

        do i=1,n
             V1=0.5d0*K12*((XGRID(i)-X1)**2)
             V2=HBAROMGA+0.5d0*K12*((XGRID(i)-X2)**2)
             DD=4d0*(V12)**2 + (V1-V2)**2
             EVOLVE(i,i)=cdexp(-0.5d0*IOTA*(V1+V2)*DELT/HBAR) &
                   *( dcos(dsqrt(DD)*DELT/2d0/HBAR) + &
                   IOTA*dsin(dsqrt(DD)*DELT/2d0/HBAR)/dsqrt(DD) & 
                   *(V2-V1))
         EVOLVE(i+n,i+n)=cdexp(-0.5d0*IOTA*(V1+V2)*DELT/HBAR) &
                   *( dcos(dsqrt(DD)*DELT/2d0/HBAR) + &
                   IOTA*dsin(dsqrt(DD)*DELT/2d0/HBAR)/dsqrt(DD) & 
                   *(V1-V2))
        enddo

        do i=1,n
             V1=0.5d0*K12*((XGRID(i)-X1)**2)
             V2=HBAROMGA+0.5d0*K12*((XGRID(i)-X2)**2)
             DD=4d0*(V12)**2 + (V1-V2)**2
           EVOLVE(i,i+n)=cdexp(-0.5d0*IOTA*(V1+V2)*DELT/HBAR) &
                   *( IOTA*dsin(dsqrt(DD)*DELT/2d0/HBAR)/dsqrt(DD) & 
                   *(-2d0*V12))
           EVOLVE(i+n,i)=cdexp(-0.5d0*IOTA*(V1+V2)*DELT/HBAR) &
                   *( IOTA*dsin(dsqrt(DD)*DELT/2d0/HBAR)/dsqrt(DD) & 
                   *(-2d0*V12))
       enddo

       ftdata1=0.d0
       ftdata2=0.d0
       TEMPVEC=0.d0

!!!!!!!!!!!   Time loop starts    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do k=1,10000

                NORM=0.d0; PE=0.d0; KE=0.d0; 

!!!!!!!!!!!!!!  Normalizing at each timestep    !!!!!!!!!!!!!!!
 
                do i=1,2*n
                   NORM=NORM + real(conjg(PSIO(i))*PSIO(i)*DX)
                enddo
                PSIO=PSIO/dsqrt(NORM)

               do i=1,n
                  V1=0.5d0*K12*((XGRID(i)-X1)**2)
                  V2=0.5d0*K12*((XGRID(i)-X2)**2)
                  PE=PE + conjg(PSIO(i))*V1*PSIO(i)*DX
                  PE=PE + conjg(PSIO(i+n))*V2*PSIO(i+n)*DX
               enddo

               do i=1,n
                  ftdata1(2*i-1)=real(PSIO(i))
                  ftdata1(2*i)=aimag(PSIO(i))
                  ftdata2(2*i-1)=real(PSIO(i+n))
                  ftdata2(2*i)=aimag(PSIO(i+n))
               enddo

!!!!!!!!!!!!   Fourier transform subroutine     !!!!!!!!!!!!!!        

               call four1(ftdata1,n,1)
               call four1(ftdata2,n,1)

               do i=1,n
                  TEMPVEC(i)=ftdata1(2*i-1)+IOTA*ftdata1(2*i)
                  TEMPVEC(i+n)=ftdata2(2*i-1)+IOTA*ftdata2(2*i)
               enddo

               NORM=0.d0
               do i=1,2*n
                  NORM=NORM + real(conjg(TEMPVEC(i))*TEMPVEC(i)*DP)
               enddo
               TEMPVEC=TEMPVEC/dsqrt(NORM)/dsqrt(HBAR)

!!!!!!!!!!!!   Calculating Kinetic,potential energy     !!!!!!!!!!!!!!

               do i=1,n
                  KE=KE + conjg(TEMPVEC(i))*0.5d0*PGRID(i)**2* &
                          TEMPVEC(i)*DP/MASS
                  KE=KE + conjg(TEMPVEC(i+n))*0.5d0*PGRID(i)**2* & 
                          TEMPVEC(i+n)*DP/MASS
               enddo
              
               OVERLAP=0.d0
               do i=1,n
                  OVERLAP=OVERLAP+ real(abs(conjg(PSIO(i))*PSIO(i+n)*DX)**2) 
               enddo
       
                write(20,*)DELT*real(k),real(PE),real(KE) &
                         ,real(PE)+real(KE)


!!!!!!!!!!!!!!!!  Evolution of wavepacket by evolution operator !!!!!!!
 
                TEMPVEC=PSIO
              
                do i=1,n
                   ftdata1(2*i-1)=real(TEMPVEC(i))
                   ftdata1(2*i)=aimag(TEMPVEC(i))
                   ftdata2(2*i-1)=real(TEMPVEC(i+n))
                   ftdata2(2*i)=aimag(TEMPVEC(i+n))
                enddo
         
                call four1(ftdata1,n,1)
                call four1(ftdata2,n,1)

                  
                do i=1,n
                   TEMPVEC(i)=ftdata1(2*i-1)+IOTA*ftdata1(2*i)
                   TEMPVEC(i+n)=ftdata2(2*i-1)+IOTA*ftdata2(2*i)
                enddo


                TEMPVEC=matmul(EVOLTE,TEMPVEC)
        
                do i=1,n
                   ftdata1(2*i-1)=real(TEMPVEC(i))
                   ftdata1(2*i)=aimag(TEMPVEC(i))
                   ftdata2(2*i-1)=real(TEMPVEC(i+n))
                   ftdata2(2*i)=aimag(TEMPVEC(i+n))
                enddo
        
                call four1(ftdata1,n,-1)
                call four1(ftdata2,n,-1)

                do i=1,n
                   TEMPVEC(i)=ftdata1(2*i-1) + IOTA*ftdata1(2*i)
                   TEMPVEC(i+n)=ftdata2(2*i-1) + IOTA*ftdata2(2*i)
                enddo
                   TEMPVEC=TEMPVEC/real(n)

                TEMPVEC=matmul(EVOLVE,TEMPVEC)

                do i=1,n
                   ftdata1(2*i-1)=real(TEMPVEC(i))
                   ftdata1(2*i)=aimag(TEMPVEC(i))
                   ftdata2(2*i-1)=real(TEMPVEC(i+n))
                   ftdata2(2*i)=aimag(TEMPVEC(i+n))
                enddo
        
                call four1(ftdata1,n,1)
                call four1(ftdata2,n,1)

                do i=1,n
                   TEMPVEC(i)=ftdata1(2*i-1) + IOTA*ftdata1(2*i)
                   TEMPVEC(i+n)=ftdata2(2*i-1) + IOTA*ftdata2(2*i)
                enddo

                TEMPVEC=matmul(EVOLTE,TEMPVEC)

                do i=1,n
                   ftdata1(2*i-1)=real(TEMPVEC(i))
                   ftdata1(2*i)=aimag(TEMPVEC(i))
                   ftdata2(2*i-1)=real(TEMPVEC(i+n))
                   ftdata2(2*i)=aimag(TEMPVEC(i+n))
                enddo
        
                call four1(ftdata1,n,-1)
                call four1(ftdata2,n,-1)

                do i=1,n
                   TEMPVEC(i)=ftdata1(2*i-1) + IOTA*ftdata1(2*i)
                   TEMPVEC(i+n)=ftdata2(2*i-1) + IOTA*ftdata2(2*i)
                enddo
                TEMPVEC=TEMPVEC/real(n)

                PSIO=TEMPVEC

        enddo
        
        close(20)

        end


