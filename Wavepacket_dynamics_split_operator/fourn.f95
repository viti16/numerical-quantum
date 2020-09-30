 
 
          module dftfft
 
          implicit none
 
          contains

!* Subroutine that handles the factor of n that appears in the FFT
!* routine of 
!* Numerical Recipes (see below).
!*
!        subroutine fft(data,n,sign)
!        implicit none
!        integer i,n,sign
!        double complex data(n)
!        double precision fact
!
!        call four1(data,n,sign)
!        fact = 1d0
!        if (sign.eq.-1) fact = 1d0/dble(n)
!        do i = 1,n
!           data(i) = data(i)*fact
!        enddo
!
!        return
!        end
!*
!*
!* The order of frequencies is extremely important. (See Numerical
!* Recipes 
!* for further information.) The routine uses the following map:
!*
!* Term #:               Frequency:
!* 1 through N/2+1       positive [from 0 to Nyquist]
!* N/2+1 through N       negative [-Nyquist to 0)
!*
!* P.S. We have made use of the mathematical meaning of the brackets.
!*
        subroutine four1(data,nn,isign)
        implicit none
        integer isign, nn
        double precision data(2*nn)
        integer i, istep, j, m, mmax, n
        double precision tempi, tempr
        double precision theta, wi, wpi, wpr, wr, wtemp
        n = 2*nn
        j = 1
        do i = 1,n,2
           if(j .gt. i)then
              tempr = data(j)
              tempi = data(j+1)
              data(j) = data(i)
              data(j+1) = data(i+1)
              data(i) = tempr
              data(i+1) = tempi
           endif
           m = n/2
1          if((m .ge. 2).and.(j .gt. m))then
              j = j - m
              m = m/2
              goto 1
           endif
           j = j + m
        enddo
        mmax = 2
2       if(n .gt. mmax)then
           istep = 2*mmax
           theta = 6.28318530717959d0/dble(isign*mmax)
           wpr = -2.d0*dsin(0.5d0*theta)**2
           wpi = dsin(theta)
           wr = 1.d0
           wi = 0.d0
           do m = 1,mmax,2
              do i = m,n,istep
                 j = i + mmax
                 tempr = sngl(wr)*data(j) - sngl(wi)*data(j+1)
                 tempi = sngl(wr)*data(j+1) + sngl(wi)*data(j)
                 data(j) = data(i) - tempr
                 data(j+1) = data(i+1) - tempi
                 data(i) = data(i) + tempr
                 data(i+1) = data(i+1) + tempi
              enddo
              wtemp = wr
              wr = wr*wpr - wi*wpi + wr
              wi = wi*wpr + wtemp*wpi + wi
           enddo
           mmax = istep
           goto 2
        endif
        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        end module
