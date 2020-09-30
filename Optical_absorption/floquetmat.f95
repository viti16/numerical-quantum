    

      module floquetsub

         implicit none

         contains

      subroutine floquetmat(n,EO,OMGA,thet,ZEROHAM & 
                                ,MUMAT,m,HAM,EIG,EIGVEC)
       
         implicit none
       
             COMPLEX :: HAM((2*n+1)*m,(2*n+1)*m), EIG((2*n+1)*m)
             COMPLEX :: EIGVEC((2*n+1)*m,(2*n+1)*m)
             COMPLEX*8 :: phase
             COMPLEX*8 :: IOTA=(0d0,1d0)
             REAL*8 :: OMGA, EO,ZEROHAM(m,m),MUMAT(m,m),thet
             INTEGER*4:: i,j,k,l,m,n,di
       
             di=(2*n+1)*m
       
             HAM=(0d0,0d0)
             EIG=(0d0,0d0)
             EIGVEC=(0d0,0d0)

             phase=dcos(thet)-IOTA*dsin(thet)
       
       !!!!!!!!!!!!!!!  diagonal elements
       

             j=-n
             do k=1,2*n+1
                do i=1,m
                   HAM((k-1)*m+i,(k-1)*m+i)=ZEROHAM(i,i)+real(j)*OMGA
                enddo
             j=j+1
             enddo
        
       !!!!!!!!!!!!!!!  offdiagonal elements
             
          do k=1,2*n+1
             do i=1,m
                do j=1,m
                   if ( k .le. 2*n) then
                        HAM((k-1)*m+i,k*m+j)=MUMAT(i,j)*EO/2d0*phase
                   endif
                enddo
             enddo
          enddo
               
             do i=1,di
                do j=i,di
                   HAM(j,i)=conjg(HAM(i,j))
                enddo
             enddo
       
       !!!!!!!!!!!!!!!  diagonalization
       
               call CALL_CGEEV (HAM,di,EIG,EIGVEC)
                      
            return 

             end  subroutine


       !!!!!!!!!!!!!!!  diagonalization subroutine

                      SUBROUTINE CALL_CGEEV(H, NDIM,eigvals,VR)
                 
                        IMPLICIT NONE
                 
                        CHARACTER (LEN = 1), PARAMETER  :: JOBVL="N"
                        CHARACTER (LEN = 1), PARAMETER  :: JOBVR="V"
                 
                        COMPLEX    (KIND = 8), ALLOCATABLE :: WORK(:)
                 
                 
                        INTEGER (KIND = 4) :: NDIM,INFO
                        INTEGER (KIND = 4) :: LDA,LDVL,LDVR,LWORK
                 
                        REAL    (KIND = 8):: RWORK(2*NDIM)

                        COMPLEX  :: H(NDIM,NDIM)
                        COMPLEX  :: VL(NDIM,NDIM)
                        COMPLEX  :: VR(NDIM,NDIM)
                        COMPLEX  :: EIGVALS(NDIM)
                 
                         LDA = NDIM
                         LDVL = NDIM
                         LDVR = NDIM

                         LWORK=2*NDIM+10
                 
                         ALLOCATE(WORK(LWORK))

                         CALL CGEEV (JOBVL,JOBVR, NDIM, H, LDA, &
                                   EIGVALS,VL,LDVL,VR,LDVR,WORK,LWORK &
                                        ,RWORK,INFO)
                         DEALLOCATE(WORK)
                 
                        return

                      END



                      END MODULE
