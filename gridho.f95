!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   Program to compute eigenvalues and eigenvectors     !!!!!
!!!!!!!   for 1-D Harmonic oscillator using grid based method !!!!!
!!!!!!!   Author: Vishal Tiwari                               !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                 !-------------------------!
                                   PROGRAM gridbasHO
                                 !-------------------------!
                                        
      implicit none

      REAL*8,allocatable :: HAM(:,:), EIGVAL(:,:), EIG(:), TE(:,:),VE(:,:)
      REAL*8,allocatable :: TEMPVEC1(:,:),TEMPVEC2(:,:)
      REAL*8 :: OMGA,HBAR,MASS,DX,TECONST,POT,AC
      INTEGER*4:: i,j,k,n,o,p,a


      HBAR=1.d0;OMGA=1.d0;MASS=1.d0

      
!!!!!!!!!!!!!!!!!  grid is from -16 to 16 with grid spacing 0.1 unit~~~~~~~~~~~

   do o=1,1

      AC=-16d0
      n=300
      DX=2.d0*abs(AC)/(n-1)
    

!      n=int(abs(AC)*2.d0/DX)

      allocate(EIG(n))
      allocate(HAM(n,n))
      allocate(VE(n,n))
      allocate(TE(n,n))
      allocate(TEMPVEC1(n,n))



!!!!!!!!!!!!!!!!  setting up hamiltonian   ~~~~~~~~~~~~~~~~~

            do i=1,n
               do j=1,n
                  HAM(i,j)=0.d0
                  TE(i,j)=0.d0
                  VE(i,j)=0.d0
               enddo
                  EIG(i)=0.d0
            enddo
  
!!!!!!!!!!!!!!!!  kinetic energy matrix   ~~~~~~~~~~~~~~~~~

            TECONST=0.5d0/(DX)**2

            do i=1,n-1
               TE(i,i)=TECONST*2.d0
               TE(i,i+1)=TECONST*(-1.d0)
               TE(i+1,i)=TECONST*(-1.d0)
            enddo 
               TE(n,n)=TECONST*2.d0
      
!!!!!!!!!!!!!!!!  potential energy matrix   ~~~~~~~~~~~~~~~~~

                VE=0.d0
            do i=1,n
               VE(i,i)=0.5d0*(AC+real(i-1)*DX)**2
            enddo

            do i=1,n 
               do j=1,n
                  HAM(i,j)=TE(i,j)+VE(i,j)
               enddo
            enddo

!!!!!!!!!!!!!!!  diagonalization ~~~~~~~~~~~~~~~~~~~~~~~~~
      
       call CALL_DSYEV(HAM,n,EIG)

!!!!!!!!!!!!!!! writing the energy and eigenvector ~~~~~~~

            do i=1,n
               do j=1,n
                  TEMPVEC1(i,j)=EIG(j)+150.d0*HAM(i,j)**2
               enddo
            enddo
       
     
            do i=1,n 
!               write(*,*)DX,EIG(3),abs(EIG(3)-(2.5d0)),dlog(abs(EIG(3)-(2.5d0)))* &
!                         0.4342944819d0
                         
               POT=AC+real(i-1)*DX
               write(*,*)POT,(TEMPVEC1(i,j),j=1,100)
!               write(*,*)POT,VE(i,i),EIG(i),(i-0.5d0),(TEMPVEC1(i,j),j=1,n)
            enddo
    
   
      deallocate(EIG)
      deallocate(HAM)
      deallocate(VE)
      deallocate(TE)
      deallocate(TEMPVEC1)

    enddo

         

  end


!!!!!!!!!!!!!!!!!!! Diagonalizatin subroutine ~~~~~~~~~~~~~~~

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




