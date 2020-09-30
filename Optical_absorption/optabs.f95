!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!     Optical absorption code for n level system         !!!!!!!!
!!!!!!!        monochromatic laser on single processor         !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!        run as         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!        gfortran optabs.f95 floquetmat.f95 -llapack     !!!!!!!!!
!!!!!!        ./a.out > <outputfilename>               !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       PROGRAM optabsnlvl

                use floquetsub

               implicit none
               
               real*8,parameter :: pi=3.1415962535
               complex*8,parameter :: IOTA=(0d0,1d0)
               
               complex,allocatable :: HAM(:,:),EIG(:),EIGVEC(:,:)
               real*8,allocatable :: ZEROHAM(:,:),MUMAT(:,:)
               real*8,allocatable ::absorb(:,:,:,:),finabs(:,:)
               complex*8,allocatable :: MUNLLP(:,:)
               complex*8,allocatable :: MUNSUM(:,:),MUNF(:,:),TEMPVEC(:)
               real*8 :: phi,pp,TIME,THET
               real*8 :: deltaf1,deltaf2
               real*8 :: OMEGA,ED
               complex*8:: popfact1,popfact2,popfact3,popfact4
               complex*8:: pop
               
               integer*4::i,j,k,l,m,n,o,p,t,di,mun
               integer*4::q,n1,n2,n3
               integer*4,allocatable::testmat(:),LGRID(:)

               n=21
               m=2
               mun=n*4+1
               di=(2*n+1)*m

               allocate(EIG(di))
               allocate(TEMPVEC(di))
               allocate(HAM(di,di))
               allocate(EIGVEC(di,di))
               allocate(MUNLLP(di*di,2))
               allocate(MUNSUM(di*di/m,2))
               allocate(MUNF(mun,2))
               allocate(LGRID(di))
               allocate(ZEROHAM(m,m))
               allocate(MUMAT(m,m))
               allocate(testmat(m))
              allocate(absorb(m,m,mun,2))
              allocate(finabs(mun,2))

               k=-n
               do i=1,2*n+1
                  do j=1,m
                       LGRID((i-1)*m+j)=k
                  enddo
                  k=k+1
               enddo

!!!!!!!!!!!!!!!!!!!!      laser parameters     !!!!!!!!!!!

               ED=0.1d0
               THET=0d0
               OMEGA=0.85d0

!!!!!!!!!!!!!!!!!!!!  Field free Hamiltonian   !!!!!!!!!!!

               ZEROHAM=0d0
               ZEROHAM(1,1)= 0.5d0
               ZEROHAM(2,2)=-0.5d0

!!!!!!!!!!!!!!!!!!!!!!!    Dipole matrix      !!!!!!!!!!!!!!!!!!!!!!!!

               MUMAT=0d0
               MUMAT(1,2)=1d0
               MUMAT(2,1)=1d0
            
               
           do p=0,0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Subroutine for Floqeut matrix      !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               call floquetmat(n,ED,OMEGA,THET,ZEROHAM, & 
                        MUMAT,m,HAM,EIG,EIGVEC)

                testmat=0
                k=1
                do i=1,di

                    if (abs(EIG(i)) .lt. OMEGA/2d0) then
                           testmat(k)=i
                            k=k+1
                    endif
                enddo
                if (k .gt. m+1) then
                        write(*,*)k,'Error FBZ'
                endif
                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Calculating Population Factor      !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                absorb=0d0
                HAM=EIGVEC

                do n2=1,m
                do n1=1,m

                   t=testmat(n2)
                   o=testmat(n1)

                popfact1=0d0
                popfact2=0d0
                popfact3=0d0

                do q=1,m-1
                  do i=m,di,m
                     do j=m,di,m
                        do k=q,di,m
                           do l=q,di,m
                              if (LGRID(i) .eq. LGRID(j) .and. & 
                               LGRID(k) .eq. LGRID(l) ) then
                            popfact1=popfact1+conjg(HAM(i,o))*HAM(j,o) &
                                    *conjg(HAM(k,t))*HAM(l,t)
                              else  if (LGRID(i) .eq. LGRID(l) .and. & 
                               LGRID(j) .eq. LGRID(k) ) then
                            popfact2=popfact2+conjg(HAM(i,o))*HAM(j,o) &
                                    *conjg(HAM(k,t))*HAM(l,t)
                              else  if (LGRID(i) .eq. -LGRID(k) .and. & 
                               LGRID(j) .eq. -LGRID(l) ) then
                            popfact3=popfact3+conjg(HAM(i,o))*HAM(j,o) &
                                    *conjg(HAM(k,t))*HAM(l,t)
                               endif
                            enddo
                         enddo
                    enddo
                  enddo
                enddo

                             
                                          
               pop=(popfact1+popfact2+popfact3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!  Calculating time dependent Dipole moments  !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                MUNLLP=0d0
                MUNF=0d0
                MUNSUM=0d0
               
               do i=1,di
                    do j=1,di
                         MUNLLP((i-1)*di+j,1)=LGRID(j)-LGRID(i)
                    enddo
               enddo

               do i=1,2*n+1
                  do k=1,m
                    do j=1,2*n+1
                           do l=1,m
                             MUNLLP((i-1)*di*m+(k-1)*di+(j-1)*m+l,2) = &
                                           MUMAT(k,l) &
                                    *conjg(HAM((i-1)*m+k,t)) & 
                                       *(HAM((j-1)*m+l,o))
                            enddo
                         enddo
                    enddo
               enddo


               j=1 
               do i=1,di*di/m
                  MUNSUM(j,1)=MUNLLP((i-1)*m+1,1)
                    do k=1,m
                       MUNSUM(j,2)=MUNSUM(j,2)+MUNLLP((i-1)*m+k,2)
                    enddo
                  j=j+1
               enddo


                do k=-2*n,2*n
                do i=1,di*di/m
                   if (abs(MUNSUM(i,1)-real(k)) .lt. 0.00001d0) then
                          MUNF(k+2*n+1,1)=real(k)
                          MUNF(k+2*n+1,2)=MUNF(k+2*n+1,2)+MUNSUM(i,2)
                   endif
                enddo
                enddo
                

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  calculating and storing absorptions !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            

         do i=1,mun

            absorb(n1,n2,i,1)=real(EIG(t)-EIG(o))+MUNF(i,1)*OMEGA
            absorb(n1,n2,i,2)=real(conjg(MUNF(i,2)*MUNF(i,2)))*pop

         enddo

         enddo
         enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!  diagonal sorting !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         finabs=0d0

         do i=1,m
            do j=1,mun
            finabs(j,1)=absorb(i,i,j,1)
            finabs(j,2)=finabs(j,2)+absorb(i,i,j,2)
            enddo
         enddo

         do i=1,(mun-1)/2
            finabs(mun+1-i,2)=finabs(mun+1-i,2)-finabs(i,2)
         enddo

         do i=(mun-1)/2+1,mun
            write(*,*)finabs(i,1),finabs(i,2)
         enddo
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!  offf diagonal sorting and printing !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          finabs=0d0

       do o=1,m-1
          do t=o+1,m

             finabs=0d0
             k=0
            do i=1,mun
               if (absorb(o,t,mun+1-i,1) .ge. 0d0) then
                   finabs(i,1)=absorb(o,t,mun+1-i,1)
                   finabs(i,2)=absorb(o,t,mun+1-i,2)-absorb(t,o,i,2)
                   k=k+1
                endif
            enddo

            do i=1,mun
               if (absorb(t,o,mun+1-i,1) .ge. 0d0) then
                   finabs(i+k,1)=absorb(t,o,mun+1-i,1)
                   finabs(i+k,2)=absorb(t,o,mun+1-i,2)-absorb(o,t,i,2)
                endif
            enddo

            do i=1,mun
               write(*,*)finabs(i,1),finabs(i,2)
            enddo

          enddo
       enddo


       write(*,*)
       write(*,*)

       !!!!!!!!!!! p loop end!!!!!!!!!!!!
       enddo
            
          
            






      end


