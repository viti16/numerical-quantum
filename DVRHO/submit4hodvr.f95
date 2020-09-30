                                 !-------------------------!
                                   PROGRAM gridbasHO
                                 !-------------------------!
                                        
      implicit none

      REAL*8,allocatable :: HAM(:,:), EIGVAL(:,:), EIG(:), TE(:,:)
      REAL*8,allocatable :: VEX(:,:),VEY(:,:)
!      REAL*8,allocatable :: TEMPVEC1(:,:),TEMPVEC2(:,:)
      REAL*8 :: OMGA,HBAR,MASS,DX,TECONST,POT1,POT2,AC,pi=3.1415926535d0
      REAL*8 :: DY
      INTEGER*4:: i,j,k,l,n,o,p,a


      HBAR=1.d0;OMGA=1.d0;MASS=1.d0

      
!!!!!!!!!!!!!!!!!  grid is from -10 to 10 with grid spacing 0.01 unit~~~~~~~~~~~

   do o=1,1

      AC=-1.d0*10.d0
      DX=0.2d0*real(o)
      DY=DX
     
      n=int(abs(AC)*2.d0/DX)

      allocate(EIG(n*n))
      allocate(HAM(n*n,n*n))
      allocate(VEX(n,n))
      allocate(VEY(n,n))
      allocate(TE(n,n))



!!!!!!!!!!!!!!!!  setting up hamiltonian   ~~~~~~~~~~~~~~~~~

            do i=1,n
               do j=1,n
                  TE(i,j)=0.d0
                  VEX(i,j)=0.d0
                  VEY(i,j)=0.d0
               enddo
            enddo
            do i=1,n*n
               do j=1,n*n
                  HAM(i,j)=0.d0
               enddo
            enddo

!!!!!!!!!!!!!!!! kinetic energy ~~~~~~~~~~~~~~~~~~~~~~

            TECONST=0.5d0/(DX)**2
  
            do i=1,n
               TE(i,i)=TECONST*(pi**2)/3.d0
            enddo 
            do i=1,n
               do j=i+1,n
                  TE(i,j)=TECONST*(-1.d0)**(i-j)*2.d0/(i-j)**2
               enddo
            enddo 
            do i=1,n
               do j=i+1,n
                  TE(j,i)=TE(i,j)
               enddo
            enddo 
            
!!!!!!!!!!!!!!!!!!! potential energy ~~~~~~~~~~~~~~~~~~~

            do i=1,n
               VEX(i,i)=0.5d0*(-10.d0+real(i-1)*DX)**2
               VEY(i,i)=0.5d0*(-10.d0+real(i-1)*DY)**2
            enddo

!!!!!!!!!!!!!!!!!!! hamiltonian ~~~~~~~~~~~~~~~~~~~~~~~~
!!!!!!!!!!!!!!!!!!! diagonal term ~~~~~~~~~~~~~~~~~~~~~~

            do k=1,n
               do i=1,n
                  do j=1,n
                     HAM(i+(k-1)*n,j+(k-1)*n)=TE(i,j)+VEX(i,j)
                  enddo
                     HAM(i+(k-1)*n,i+(k-1)*n)=VEY(k,k)+HAM(i+(k-1)*n,i+(k-1)*n)+TE(k,k)
               enddo
            enddo

!!!!!!!!!!!!!!!!!!! hamiltonian ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!!!!!!!!!!!!!!!!! off- diagonal term ~~~~~~~~~~~~~~~~~~~~~~

            do k=1,n-1
               do l=k-1,0,-1
                  do i=1,n
                        HAM(i+l*n,i+k*n)=TE(l+1,k+1)
                        HAM(i+k*n,i+l*n)=TE(k+1,l+1)
                  enddo
               enddo
            enddo

!!!!!!!!!!!!!!!  diagonalization ~~~~~~~~~~~~~~~~~~~~~~~~~
      
       call CALL_DSYEV(HAM,n*n,EIG)

!!!!!!!!!!!!!!!  writing terms ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
         write(*,*)EIG(1),EIG(2),EIG(3),EIG(4)

          do i=1,n
             do j=1,n
                POT1=-10.d0+real(j-1)*DX
                POT2=-10.d0+real(i-1)*DY
                write(*,*)POT1,POT2,(HAM(j+(i-1)*n,k)**2,k=1,4)
             enddo
          enddo

   

    enddo

         

  end



     SUBROUTINE CALL_DSYEV(H, NDIM, eigvals)

       IMPLICIT NONE

       CHARACTER (LEN = 1), PARAMETER  :: JOBZ="V"
       CHARACTER (LEN = 1), PARAMETER  :: UPLO="L"

       REAL    (KIND = 8), ALLOCATABLE :: WORK(:)

       INTEGER (KIND = 4) :: NDIM,INFO
       INTEGER (KIND = 4) :: LDA, LDWORK

       REAL    (KIND = 8) :: H(NDIM,NDIM)
       REAL    (KIND = 8) :: EIGVALS(NDIM)

        LDA = NDIM
        LDWORK = 3*NDIM - 1

        ALLOCATE(WORK(LDWORK))
        CALL DSYEV (JOBZ, UPLO, NDIM, H, LDA, &
                          EIGVALS, WORK, LDWORK, INFO)
        DEALLOCATE(WORK)

     END




