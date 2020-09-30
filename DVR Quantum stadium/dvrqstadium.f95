!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!   Program to compute eigenvalues and eigenvectors  !!!!!
!!!!!!!   for quantum stadium using dvr method             !!!!!
!!!!!!!   Author: Vishal Tiwari                            !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                 !-------------------------!
                                   PROGRAM ddvrHO
                                 !-------------------------!
                                        
      implicit none

      REAL*8,allocatable :: HAM(:,:), EIGVAL(:,:), EIG(:), TE(:,:)
      REAL*8,allocatable :: VEX(:,:),VEY(:,:)
!      REAL*8,allocatable :: TEMPVEC1(:,:),TEMPVEC2(:,:)
      REAL*8 :: OMGA,HBAR,MASS,TECONST,POT1,POT2,AC,pi=3.1415926535d0
      REAL*8 :: DX,DY,XC,YC
      INTEGER*4:: i,j,k,l,n,o,p,a


      HBAR=1.d0;OMGA=1.d0;MASS=1.d0
      open(10,file="recqstadeval2.txt")
      open(20,file="recqstadevec2.txt")

      
!!!!!!!!!!!!!!!!!  grid is from -2.5 to 2.5 with grid spacing 0.025 unit~~~~~~~~~~~

       write(*,*)"jobs starts"
      AC=-1.d0*2.5d0
      n=200
      DX=2*abs(AC)/(n-1)
      DY=DX
     
!!!!!!!!!!!!!!!!  setting up hamiltonian   ~~~~~~~~~~~~~~~~~

      allocate(EIG(n*n))
      allocate(HAM(n*n,n*n))
      allocate(TE(n,n))



            do i=1,n
               do j=1,n
                  TE(i,j)=0.d0
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
            

!!!!!!!!!!!!!!!!!!! hamiltonian    ~~~~~~~~~~~~~~~~~~~~~~~~
!!!!!!!!!!!!!!!!!!! diagonal block ~~~~~~~~~~~~~~~~~~~~~~~~
!!!!!!!!!!!!!!!!!!! with potential ~~~~~~~~~~~~~~~~~~~~~~~~

            do k=1,n
               do i=1,n
                  do j=1,n
                     HAM(i+(k-1)*n,j+(k-1)*n)=TE(i,j)
                  enddo
                     XC=AC+real(i-1)*DX
                     YC=AC+real(k-1)*DY
                     if ( ((XC-0.5d0)**2 + YC**2) .le. (0.5d0)**2 .or. &
                          ((XC+0.5d0)**2 + YC**2) .le. (0.5d0)**2 .or. &
                          (abs(XC)) .le. 0.5d0 .and. (abs(YC)) .le. 0.5d0 ) then
!!!!!!!!!!!!!!!!!!! inside potential zero ~~~~~~~~~~~~~~~~~~~
                               a=0.d0
                     else 
!!!!!!!!!!!!!!!!!!! outside potential infinite ~~~~~~~~~~~~~~
                               a=100000000d0
                     endif
                     HAM(i+(k-1)*n,i+(k-1)*n)=a+HAM(i+(k-1)*n,i+(k-1)*n)+TE(k,k)
               enddo
            enddo

!!!!!!!!!!!!!!!!!!! hamiltonian ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!!!!!!!!!!!!!!!!! off-diagonal block ~~~~~~~~~~~~~~~~~~~~~~

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

!!!!!!!!!!!!!!!  writing energy and eigenvectors ~~~~~~~~~~~~~~~~~~~~~~~

       do i=1,n 
         write(10,*)EIG(i)
       enddo

          do i=1,n
             do j=1,n
                POT1=AC + real(j-1)*DX
                POT2=AC + real(i-1)*DY
                write(20,*)POT1,POT2,(HAM(j+(i-1)*n,k),k=1,n,20)
             enddo
                write(20,*)
                write(20,*)
          enddo

   close(10)
   close(20)

         

  end

!!!!!!!!!!!!!!!!!!!!!!! subroutine for diagonalization ~~~~~~~~~~~~~~

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




