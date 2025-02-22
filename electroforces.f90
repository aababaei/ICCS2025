!     Subroutines for computing electrostatic interaction between two spheres
!     Developed by Ahmad Ababaei (ahmad.ababaei@imgw.pl) and Antoine Michel (antoine.michel@cea.fr)

      PROGRAM ELECTROFORCES
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PARAMETER ( acu = 1d-10 )
      INTEGER sample

      pi = 4d0*DATAN(1d0)
       e = 1.602176634d-19 ! C

! ============ I N P U T S ==============

        alam = 1.00d-0   ! Radius ratio
          a1 = 1d+0      ! Larger drop radius [μm]
          a2 = alam * a1 ! Smaller drop radius [μm]
          E0 = 0d-1      ! Electric field intensity [V/cm]
         psi = 0d0*pi    ! Electric field angle
          qr = 1.0d0     ! Electric charge ratio
          q1 = 200*e     ! Electric charge [C]
          q2 = qr * q1   ! Electric charge [C]
        epsr = 80d0      ! Dielectric constant / relative permittivity (water = 80, perfect conductor: ∞)

! ============= U N I T S ===============

          a1 = a1 * 1d-4 ! [μm] to [cm]
          a2 = a2 * 1d-4 ! [μm] to [cm]
          q1 = q1 * 2997924580d0 ! [C] to [statC]
          q2 = q2 * 2997924580d0 ! [C] to [statC]

! ========= L O G   D I S T R. ========== 
! Logarithmic distribution of normalized
! gap size ξ = s — 2 in JO84 notation:
      xi_min = 1d-2
      xi_max = 1d+2
      sample = 10
      dlt_xi = DLOG ( xi_max / xi_min ) / DBLE(sample-1)
           s = 2d0 + xi_min
        DO i = 1, sample
           r = s * ( a1 + a2 ) / 2d0

! ============ M E T H O D ==============

!     CALL COULOMB(q1,q2,r,F12)
!     CALL D64(a1,a2,r,q1,q2,E0,psi,F1,F2,acu)
      CALL KCSB14_T(a1,a2,r,q1,q2,epsr,F12,acu)
!     CALL KCSB14_H(a1,a2,r,q1,q2,epsr,F12,acu)
!     CALL KCSB14_G(a1,a2,r,q1,q2,epsr,F12,acu)
!     CALL BPBS21(a1,a2,r,q1,q2,F12,acu)


! ============ O U T P U T ==============
!     WRITE(*,*) s-2d0, F1, F2, F12
      STOP

! ========= L O G   D I S T R. ==========
      s = DLOG ( s - 2d0 ) + dlt_xi
      s = DEXP ( s ) + 2d0
      ENDDO

      END PROGRAM ELECTROFORCES

! ======================= S U B R O U T I N E S ========================

! === Coulomb (1785) ===================================================

      SUBROUTINE COULOMB(q1,q2,r,F12)
      IMPLICIT DOUBLE PRECISION (A-Z)
      F12 = q1 * q2 / r**2
      END SUBROUTINE COULOMB

! ======================================================================

! === Davis (1964) =====================================================
!     Davis, M. H. (1964). Two charged spherical conductors in a uniform electric field: Forces and field strength. The Quarterly Journal of Mechanics and Applied Mathematics, 17(4), 499-511.

      SUBROUTINE D64(a1,a2,r,q1,q2,E0,psi,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,L-Z)

!     eps0 = 8.854187817d-12 * 2.99792458d9**2 * 1d-6 ! esu
      eps0 = 1d0

      IF ( DABS(E0) .LT. 1d-10 ) E0 = 1d-10
      Ex = E0*DSIN(psi)
      Ez = E0*DCOS(psi)

          h = r - ( a1 + a2 )
        eps =  h / a1
       alam = a2 / a1
     coshal = 1d0 + eps   * ( alam + eps / 2d0 ) / ( 1d0 + alam + eps )
     coshbe = 1d0 + eps/alam*( 1d0 + eps / 2d0 ) / ( 1d0 + alam + eps )
        mu1 = DACOSH(coshal)
        mu2 = DACOSH(coshbe)
        m12 = mu1 + mu2
          c = a1 * DSINH(mu1)

        Qs1 = 2d0 * eps0 * c**2 * E0 * DCOS(psi) * ( S(1,mu2,m12,acu) + S(1,0d0,m12,acu) )
        Qs2 =-2d0 * eps0 * c**2 * E0 * DCOS(psi) * ( S(1,mu1,m12,acu) + S(1,0d0,m12,acu) )
        C11 = 2d0 * eps0 * c * S(0,mu2,m12,acu)
        C12 =-2d0 * eps0 * c * S(0,0d0,m12,acu)
        C22 = 2d0 * eps0 * c * S(0,mu1,m12,acu)
        del = C11 * C22 - C12**2
        P11 = C22 / del
        P12 =-C12 / del
        P22 = C11 / del
         v1 =-( P11*Qs1 + P12*Qs2 ) / ( E0 * c * DCOS(psi) )
         v2 =-( P12*Qs1 + P22*Qs2 ) / ( E0 * c * DCOS(psi) )
         w1 = ( P11*(Q1-Qs1) + P12*(Q2-Qs2) ) / ( E0 * c * DCOS(psi) )
         w2 = ( P12*(Q1-Qs1) + P22*(Q2-Qs2) ) / ( E0 * c * DCOS(psi) )
        p11 = c * P11
        p12 = c * P12
        p22 = c * P22

       F2zo = 0d0
       F2xo = 0d0
        rel = 1d0
         n  = 0d0
      DO WHILE ( rel .GT. acu )
         Yn =-DSQRT(2d0)*(2d0*n+1d0)*DEXP((n+5d-1)*mu2)
         Yn = Yn * ( (2d0*n+1d0)*(DEXP((2d0*n+1d0)*mu1)+1d0) - w2*DEXP((2d0*n+1d0)*mu1) + w1 )
         Yn = Yn / ( DEXP((2d0*n+1d0)*m12) - 1d0 )
        Ynp =-DSQRT(2d0)*(2d0*n+3d0)*DEXP((n+15d-1)*mu2)
        Ynp = Ynp* ( (2d0*n+3d0)*(DEXP((2d0*n+3d0)*mu1)+1d0) - w2*DEXP((2d0*n+3d0)*mu1) + w1 )
        Ynp = Ynp/ ( DEXP((2d0*n+3d0)*m12) - 1d0 )
         Zn = DSQRT(8d0)*(2d0*n+1d0)*DEXP((n+05d-1)*mu2) * (DEXP((2d0*n+1d0)*mu1)-1d0) / ( DEXP((2d0*n+1d0)*m12) - 1d0 )
        Znp = DSQRT(8d0)*(2d0*n+3d0)*DEXP((n+15d-1)*mu2) * (DEXP((2d0*n+3d0)*mu1)-1d0) / ( DEXP((2d0*n+3d0)*m12) - 1d0 )
        F2z = 2d0 * DCOS(psi)**2 * Yn / (2d0*n+1d0) * ( Yn - 2d0*DCOSH(mu2)*(n+1d0)/(2d0*n+3d0)*Ynp )
        F2z = F2z + DSIN(psi)**2 * n*(n+1d0)/(2d0*n+1d0)*Zn * ( Zn - 2d0*DCOSH(mu2)*(n+2d0)/(2d0*n+3d0)*Znp )
        F2z = F2zo + F2z * eps0/4d0*(c*E0)**2
        F2x = (n+1d0)/(2d0*n+1d0)/(2d0*n+3d0) * ( (n+2d0)*Znp*Yn - n*Zn*Ynp )
        F2x = F2xo + F2x * eps0/4d0*(c*E0)**2 * DSIN(2d0*psi) * DSINH(mu2)
        rel = DMAX1( DABS(F2z-F2zo)/DABS(F2z), DABS(F2x-F2xo)/DABS(F2x) )
       F2zo = F2z
       F2xo = F2x
          n = n + 1d0
      ENDDO
      F1z = Ez * ( q1 + q2 ) - F2z
      F1x = Ex * ( q1 + q2 ) - F2x
      F1  = F1z
      F2  = F2z

      END SUBROUTINE D64

      FUNCTION S(m,xi,m12,acu)
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER m
       So = 0d0
        S = 0d0
        n = 0d0
      rel = 1d0
      DO WHILE (rel .GT. acu)
         nn1 = 2d0*n + 1d0
           S = So + nn1**m * DEXP(nn1*xi) / ( DEXP(nn1*m12) - 1d0 )
         rel = DABS(S-So)/DABS(S)
          So = S
           n = n + 1d0
      ENDDO
      END FUNCTION

! ======================================================================

! === Khachatourian, Chan, Stace, Bichoutskaia (2014) ==================
!     Khachatourian, A., Chan, H. K., Stace, A. J., & Bichoutskaia, E. (2014). Electrostatic force between a charged sphere and a planar surface: A general solution for dielectric materials. The Journal of chemical physics, 140(7).

      SUBROUTINE KCSB14_T(a1,a2,r,q1,q2,epsr,F12,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

          k0 = 1d0  ! air relative permittivity
          k1 = epsr
          k2 = epsr
          km = 1d0

           h = r - ( a1 + a2 )
         eps =  h / a1
        alam = a2 / a1
      coshal = 1d0 + eps   * ( alam + eps / 2d0 ) / ( 1d0 + alam + eps )
      coshbe = 1d0 + eps/alam*( 1d0 + eps / 2d0 ) / ( 1d0 + alam + eps )
         et1 = DACOSH(coshal)
         et2 = DACOSH(coshbe)
           c = a1 * DSINH(et1)

      F12o = 0d0
!     F21o = 0d0
        iN = 32
1     CALL CPU_TIME(t0)
      ALLOCATE ( T(2*iN+2,8) )
         T = 0d0

      DO i = 1, iN+1
         n = DBLE(i-1)
        fm = DEXP(-(n-05d-1)*(et1+et2))
        fn = DEXP(-(n+05d-1)*(et1+et2))
        fp = DEXP(-(n+15d-1)*(et1+et2))

         IF (i.GT.1) THEN
            T(2*i-1,2) = -5d-1*n*(km+k1)
            T(2*i-1,3) =  5d-1*n*(km-k1)*fm
            T(2*i  ,1) =  5d-1*n*(km-k2)*fm
            T(2*i  ,2) = -5d-1*n*(km+k2)
         ENDIF
            T(2*i-1,4) = (5d-1+n)*DCOSH(et1)*(km+k1) + 5d-1*DSINH(et1) * (km-k1)
            T(2*i-1,5) = (-(5d-1+n)*DCOSH(et1)+5d-1*DSINH(et1))*(km-k1)* fn
            T(2*i  ,3) = (-(5d-1+n)*DCOSH(et2)+5d-1*DSINH(et2))*(km-k2)* fn
            T(2*i  ,4) = (5d-1+n)*DCOSH(et2)*(km+k2) + 5d-1*DSINH(et2) * (km-k2)
         IF (i.LT.iN+1) THEN
            T(2*i-1,6) = -5d-1*(n+1d0)*(km+k1)
            T(2*i-1,7) =  5d-1*(n+1d0)*(km-k1)*fp
            T(2*i  ,5) =  5d-1*(n+1d0)*(km-k2)*fp
            T(2*i  ,6) = -5d-1*(n+1d0)*(km+k2)
         ENDIF
!           K factor not needed in ESUnits:
            T(2*i-1,8) = DSQRT(2d0)*c*DEXP(-(n+5d-1)*et1)*q1/a1**2
            T(2*i  ,8) = DSQRT(2d0)*c*DEXP(-(n+5d-1)*et2)*q2/a2**2
      ENDDO
      CALL CPU_TIME(t1)
      CALL THOMAS(2*iN+2,3,3,T)
      CALL CPU_TIME(t2)
      WRITE(*,*) 1d3*(t1-t0), 1d3*(t2-t1), 1d3*(t2-t0)
!        n = 0:
        fn = DEXP(-5d-1*(et1+et2))
       F12 = fn * 5d-1 * ( -T(1,8) + T(3,8)*DEXP(-et1) ) * T(2,8)
!      F21 = fn * 5d-1 * ( -T(2,8) + T(4,8)*DEXP(-et2) ) * T(1,8)
!        n = N:
         n = DBLE(iN)
        fn = DEXP(-(n+5d-1)*(et1+et2))
       F12 = F12 + fn * ( n/2d0*T(2*iN-1,8)*DEXP( et1)  &
                 - (n+5d-1)*T(2*iN+1,8) ) * T(2*iN+2,8)
!        n = 1, N-1:
      DO i = 1, iN-1
         n = DBLE(i)
        fn = DEXP(-(n+5d-1)*(et1+et2))
       F12 = F12 + fn * ( n/2d0*T(2*i-1,8)*DEXP( et1) - (n+5d-1)*T(2*i+1,8) &
                 +  (n+1d0)/2d0*T(2*i+3,8)*DEXP(-et1) ) * T(2*i+2,8)
!      F21 = F21 + fn * ( n/2d0*T(2*i  ,8)*DEXP( et2) - (n+5d-1)*T(2*i+2,8) &
!                +  (n+1d0)/2d0*T(2*i+4,8)*DEXP(-et2) ) * T(2*i+1,8)
      ENDDO
      DEALLOCATE ( T )

!     Force in SI:
!     F12 = -F12/k
!     F21 = -F21/k
!     Force in dynes (K factor not needed in ESUnits):
      F12 = -F12
!     F21 = -F21
      rel = DABS(F12-F12o)/DABS(F12)
      IF ( iN .LT. 1024 ) THEN
!     IF ( rel .GT. acu ) THEN
         F12o = F12
!        F21o = F21
         iN = INT(2 * FLOAT(iN)) ! 50% increase
         GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === Khachatourian, Chan, Stace, Bichoutskaia (2014) ==================
!     Khachatourian, A., Chan, H. K., Stace, A. J., & Bichoutskaia, E. (2014). Electrostatic force between a charged sphere and a planar surface: A general solution for dielectric materials. The Journal of chemical physics, 140(7).

      SUBROUTINE KCSB14_H(a1,a2,r,q1,q2,epsr,F12,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: a, b, c, d, e, f, g, rhs, x
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

          k0 = 1d0  ! air relative permittivity
          k1 = epsr
          k2 = epsr
          km = 1d0

           h = r - ( a1 + a2 )
         eps =  h / a1
        alam = a2 / a1
      coshal = 1d0 + eps   * ( alam + eps / 2d0 ) / ( 1d0 + alam + eps )
      coshbe = 1d0 + eps/alam*( 1d0 + eps / 2d0 ) / ( 1d0 + alam + eps )
         et1 = DACOSH(coshal)
         et2 = DACOSH(coshbe)
          cc = a1 * DSINH(et1)

      F12o = 0d0
!     F21o = 0d0
        iN = 32
1     CALL CPU_TIME(t0)
      ALLOCATE ( a(2*iN+2), b(2*iN+2), c(2*iN+2), d(2*iN+2), e(2*iN+2), f(2*iN+2), g(2*iN+2), rhs(2*iN+2), x(2*iN+2) )
         a = 0d0
         b = 0d0
         c = 0d0
         d = 0d0
         e = 0d0
         f = 0d0
         g = 0d0
         rhs = 0d0

      DO i = 1, iN+1
         n = DBLE(i-1)
        fm = DEXP(-(n-05d-1)*(et1+et2))
        fn = DEXP(-(n+05d-1)*(et1+et2))
        fp = DEXP(-(n+15d-1)*(et1+et2))

         IF (i.GT.1) THEN
            b(2*i-1) = -5d-1*n*(km+k1)
            c(2*i-1) =  5d-1*n*(km-k1)*fm
            a(2*i  ) =  5d-1*n*(km-k2)*fm
            b(2*i  ) = -5d-1*n*(km+k2)
         ENDIF
            d(2*i-1) = (5d-1+n)*DCOSH(et1)*(km+k1) + 5d-1*DSINH(et1) * (km-k1)
            e(2*i-1) = (-(5d-1+n)*DCOSH(et1)+5d-1*DSINH(et1))*(km-k1)* fn
            c(2*i  ) = (-(5d-1+n)*DCOSH(et2)+5d-1*DSINH(et2))*(km-k2)* fn
            d(2*i  ) = (5d-1+n)*DCOSH(et2)*(km+k2) + 5d-1*DSINH(et2) * (km-k2)
         IF (i.LT.iN+1) THEN
            f(2*i-1) = -5d-1*(n+1d0)*(km+k1)
            g(2*i-1) =  5d-1*(n+1d0)*(km-k1)*fp
            e(2*i  ) =  5d-1*(n+1d0)*(km-k2)*fp
            f(2*i  ) = -5d-1*(n+1d0)*(km+k2)
         ENDIF
!           K factor not needed in ESUnits:
            rhs(2*i-1) = DSQRT(2d0)*cc*DEXP(-(n+5d-1)*et1)*q1/a1**2
            rhs(2*i  ) = DSQRT(2d0)*cc*DEXP(-(n+5d-1)*et2)*q2/a2**2
      ENDDO

      CALL CPU_TIME(t1)
      CALL HDMA(2*iN+2, a, b, c, d, e, f, g, rhs, x)
      CALL CPU_TIME(t2)
      WRITE(*,*) 1d3*(t1-t0), 1d3*(t2-t1), 1d3*(t2-t0)

!        n = 0:
        fn = DEXP(-5d-1*(et1+et2))
       F12 = fn * 5d-1 * ( -x(1) + x(3)*DEXP(-et1) ) * x(2)
!        n = N:
         n = DBLE(iN)
        fn = DEXP(-(n+5d-1)*(et1+et2))
       F12 = F12 + fn * ( n/2d0*x(2*iN-1)*DEXP( et1)  &
                 - (n+5d-1)*x(2*iN+1) ) * x(2*iN+2)
!        n = 1, N-1:
      DO i = 1, iN-1
         n = DBLE(i)
        fn = DEXP(-(n+5d-1)*(et1+et2))
       F12 = F12 + fn * ( n/2d0*x(2*i-1)*DEXP( et1) - (n+5d-1)*x(2*i+1) &
                 +  (n+1d0)/2d0*x(2*i+3)*DEXP(-et1) ) * x(2*i+2)
      ENDDO
      DEALLOCATE ( a, b, c, d, e, f, g, rhs, x )

!     Force in SI:
!     F12 = -F12/k
!     F21 = -F21/k
!     Force in dynes (K factor not needed in ESUnits):
      F12 = -F12
!     F21 = -F21
      rel = DABS(F12-F12o)/DABS(F12)
      IF ( iN .LT. 1024 ) THEN
!     IF ( rel .GT. acu ) THEN
         F12o = F12
!        F21o = F21
         iN = INT(2 * FLOAT(iN)) ! 50% increase
         GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === Khachatourian, Chan, Stace, Bichoutskaia (2014) ==================
!     Khachatourian, A., Chan, H. K., Stace, A. J., & Bichoutskaia, E. (2014). Electrostatic force between a charged sphere and a planar surface: A general solution for dielectric materials. The Journal of chemical physics, 140(7).

      SUBROUTINE KCSB14_G(a1,a2,r,q1,q2,epsr,F12,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T

          k0 = 1d0  ! air relative permittivity
          k1 = epsr
          k2 = epsr
          km = 1d0

           h = r - ( a1 + a2 )
         eps =  h / a1
        alam = a2 / a1
      coshal = 1d0 + eps   * ( alam + eps / 2d0 ) / ( 1d0 + alam + eps )
      coshbe = 1d0 + eps/alam*( 1d0 + eps / 2d0 ) / ( 1d0 + alam + eps )
         et1 = DACOSH(coshal)
         et2 = DACOSH(coshbe)
           c = a1 * DSINH(et1)

      F12o = 0d0
!     F21o = 0d0
        iN = 32
1     CALL CPU_TIME(t0)
      ALLOCATE ( T(2*iN+2,2*iN+3) )
         T = 0d0

      DO i = 1, iN+1
         n = DBLE(i-1)
        fm = DEXP(-(n-05d-1)*(et1+et2))
        fn = DEXP(-(n+05d-1)*(et1+et2))
        fp = DEXP(-(n+15d-1)*(et1+et2))

         IF (i.GT.1) THEN
            T(2*i-1,2*i-3) = -5d-1*n*(km+k1)
            T(2*i-1,2*i-2) =  5d-1*n*(km-k1)*fm
            T(2*i  ,2*i-3) =  5d-1*n*(km-k2)*fm
            T(2*i  ,2*i-2) = -5d-1*n*(km+k2)
         ENDIF
            T(2*i-1,2*i-1) = (5d-1+n)*DCOSH(et1)*(km+k1) + 5d-1*DSINH(et1) * (km-k1)
            T(2*i-1,2*i  ) = (-(5d-1+n)*DCOSH(et1)+5d-1*DSINH(et1))*(km-k1)* fn
            T(2*i  ,2*i-1) = (-(5d-1+n)*DCOSH(et2)+5d-1*DSINH(et2))*(km-k2)* fn
            T(2*i  ,2*i  ) = (5d-1+n)*DCOSH(et2)*(km+k2) + 5d-1*DSINH(et2) * (km-k2)
         IF (i.LT.iN+1) THEN
            T(2*i-1,2*i+1) = -5d-1*(n+1d0)*(km+k1)
            T(2*i-1,2*i+2) =  5d-1*(n+1d0)*(km-k1)*fp
            T(2*i  ,2*i+1) =  5d-1*(n+1d0)*(km-k2)*fp
            T(2*i  ,2*i+2) = -5d-1*(n+1d0)*(km+k2)
         ENDIF
!           K factor not needed in ESUnits:
            T(2*i-1,2*iN+3) = DSQRT(2d0)*c*DEXP(-(n+5d-1)*et1)*q1/a1**2
            T(2*i  ,2*iN+3) = DSQRT(2d0)*c*DEXP(-(n+5d-1)*et2)*q2/a2**2
      ENDDO
      CALL CPU_TIME(t1)
      CALL GAUSSB(2*iN+2,3,3,T)
!     CALL GAUSS(2*iN+2,2*iN+3,T)
      CALL CPU_TIME(t2)
      WRITE(*,*) 1d3*(t1-t0), 1d3*(t2-t1), 1d3*(t2-t0)

!        n = 0:
        fn = DEXP(-5d-1*(et1+et2))
       F12 = fn * 5d-1 * ( -T(1,2*iN+3) + T(3,2*iN+3)*DEXP(-et1) ) * T(2,2*iN+3)
!        n = N:
         n = DBLE(iN)
        fn = DEXP(-(n+5d-1)*(et1+et2))
       F12 = F12 + fn * ( n/2d0*T(2*iN-1,2*iN+3)*DEXP(et1) &
                 -      (n+5d-1)*T(2*iN+1,2*iN+3) ) * T(2*iN+2,2*iN+3)
!        n = 1, N-1:
      DO i = 1, iN-1
         n = DBLE(i)
        fn = DEXP(-(n+5d-1)*(et1+et2))
       F12 = F12 + fn * ( n/2d0*T(2*i-1,2*iN+3)*DEXP( et1) - (n+5d-1)*T(2*i+1,2*iN+3) &
                 +  (n+1d0)/2d0*T(2*i+3,2*iN+3)*DEXP(-et1) ) * T(2*i+2,2*iN+3)
!      F21 = F21 + fn * ( n/2d0*T(2*i  ,8)*DEXP( et2) - (n+5d-1)*T(2*i+2,8) &
!                +  (n+1d0)/2d0*T(2*i+4,8)*DEXP(-et2) ) * T(2*i+1,8)
      ENDDO
      DEALLOCATE ( T )

!     Force in SI:
!     F12 = -F12/k
!     F21 = -F21/k
!     Force in dynes (K factor not needed in ESUnits):
      F12 = -F12
!     F21 = -F21
      rel = DABS(F12-F12o)/DABS(F12)
      IF ( iN .LT. 1024 ) THEN
!     IF ( rel .GT. acu ) THEN
         F12o = F12
!        F21o = F21
         iN = INT(2 * FLOAT(iN)) ! 50% increase
         GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === T H O M A S   A L G O R I T H M   B A N D E D   M A T R I X ======
!     KL = Lower band: No. of sub-diagonals
!     KU = Upper band: No. of super-diagonals
!     If KL = KU = 1 then the solver works
!     similar to TDMA. The system of equations
!     has to be given to the solver in the
!     following compact form:
!     beginning from the left-most column
!     we fill T(:,j) with vectors containing
!     sub-diagonal, diagonal, super-diagonal
!     and finally the RHS (vector b) elements.
!     Example: N = 5, KL = 1, KU = 2
!     2  3  4  0  0 | 5
!     1  2  3  4  0 | 5
!     0  1  2  3  4 | 5
!     0  0  1  2  3 | 5
!     0  0  0  1  2 | 5
!     This system has to be rearranged to:
!     0  2  3  4 | 5
!     1  2  3  4 | 5
!     1  2  3  4 | 5
!     1  2  3  0 | 5
!     1  2  0  0 | 5

      SUBROUTINE THOMAS(N,KL,KU,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N,KL+KU+2)

      DO K = 1, N-1
        NI = K + KL
        IF ( NI .GT. N ) NI = N
        DO I = K+1, NI
           U = T(I, K+KL-I+1) / T(K, KL+1)
           IF ( ABS(T(K, KL+1)) .LT. 1D-15 ) &
           WRITE(*,*) 'Check: diagonal element = 0'
           NJ = K + KL + KU - I + 1
           DO J = K+KL-I+2, NJ
              T(I,J) = T(I,J) - T(K, I+J-K) * U
           ENDDO
           T(I, KL+KU+2) = T(I, KL+KU+2) - T(K, KL+KU+2) * U
        ENDDO
      ENDDO

      DO I = N, 1, -1
         S = 0D0
         DO J = KL+2, KL+KU+1
            K = I + J - ( KL + 1 )
            IF ( K .GT. N ) EXIT
            S = S + T(I,J) * T(K, KL+KU+2)
         ENDDO
         T(I, KL+KU+2) = ( T(I, KL+KU+2) - S ) / T(I, KL+1)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === G A U S S   E L I M I N A T I O N   B A N D E D   M A T R I X ====
!     KL = Lower band: No. of sub-diagonals
!     KU = Upper band: No. of super-diagonals
!     Example: N = 5, KL = 1, KU = 2
!     X  X  X  0  0 | X
!     X  X  X  X  0 | X
!     0  X  X  X  X | X
!     0  0  X  X  X | X
!     0  0  0  X  X | X
      SUBROUTINE GAUSSB(N,KL,KU,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(N,N+1)

      DO K = 1, N-1
        UK = T(K,N+1) / T(K,K)
        IJ = 1
        NI = K + KL
        IF ( NI .GT. N ) NI = N
        DO I = K+1, NI
          UI = T(I,K) / T(K,K)
          NJ = K + KU
          IF ( NJ .GT. N ) NJ = N
          DO J = K+1, NJ
            T(I,J) =  T(I,J)  - T(K,J) * UI
          ENDDO
          T(I,N+1) = T(I,N+1) - T(I,K) * UK
        ENDDO
      ENDDO

      DO I = N, 1, -1
        NJ = I + KU
        IF ( NJ .GT. N ) NJ = N
         S = 0D0
         DO J = I+1, NJ
            S = S + T(I,J) * T(J,N+1)
         ENDDO
         T(I,N+1) = ( T(I,N+1) - S ) / T(I,I)
      ENDDO

      END SUBROUTINE
! ======================================================================

! === G A U S S   E L I M I N A T I O N ================================
!     Written by Alexander Zinchenko (University of Colorado Boulder)

!     For every n, write the system in a matrix form AX=b,
!     where X=(A_n, B_n, C_n, D_n) and b is the RHS vector,
!     and solve this system by Gauss elimination with pivoting.
!     For solving a system of N linear equations for N unknowns with
!     M-N (M minus N) right-side vectors b (so, you can have solutions
!     simutaneously for matrix A and several RHS vectors b). The
!     parameter (NMAX=100, MMAX=200) is just an example to make it
!     suitable for solving with <=100 unknowns with <= 200-100 RHS
!     vectors. In your specific case (4 eqns with one RHS vector),
!     NMAX=4, MMAX=5 would suffice, but you do not have to make that
!     change, just make sure, one way or the other, that dimensions of
!     T are consistent in GAUSS and the calling routine.

!     To use GAUSS for AX=b, put A into the (N,N) part of matrix T,
!     (i.e. T_{iJ}=A_{ij} for i,j <=N) and set
!     T(1,N+1)=b_1, T(2,N+1)=b_2, ... T(N,N+1)=b_N (if you have one
!     RHS vector b). After call GAUSS(N,N+1,T) you will get the
!     solution X_i= T(i,N+1) for i<=N. In your case, N=4.

!     If you want the solutions for several different b-vectors with
!     the same matrix A, Then these vectors must be placed into the
!     columns (N+1), (N+2), ... M of matrix T before calling GAUSS,
!     and the solutions will then be found in those columns after GAUSS.
!     Of course, all calcs must be in DOUBLE PRECISION.

      SUBROUTINE GAUSS(N,M,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     PARAMETER (NMAX=100, MMAX=200)
!     DIMENSION T(NMAX,MMAX)
      DIMENSION T(N,M) ! My modification
      INTEGER S

      DO 1 I=1,N
         S=I
         R=T(I,I)
         DO 2 K=I+1,N
            IF(DABS(T(K,I)).GT.DABS(R)) THEN
              S=K
              R=T(K,I)
            ENDIF
2        CONTINUE
         T(S,I)=T(I,I)
         DO 3 J=I+1,M
            U=T(S,J)/R
            T(S,J)=T(I,J)
            T(I,J)=U
            DO 4 K=I+1,N
4              T(K,J)=T(K,J)-T(K,I)*U
3        CONTINUE
1     CONTINUE
      DO 5 I=N-1,1,-1
         DO 5 J=N+1,M
            DO 5 K=I+1,N
               R=T(I,K)
               T(I,J)=T(I,J)-R*T(K,J)
5     CONTINUE
      RETURN
      END SUBROUTINE
! ======================================================================
!     Hepta-diagonal Matrix Algorithm
      SUBROUTINE HDMA(n, a, b, c, d, e, f, g, r, x)
      DOUBLE PRECISION, DIMENSION(n) :: a, b, c, d, e, f, g, r, x
      INTEGER :: i,n

      ! Forward elimination
      DO i = 2, n-2
         c(i) = c(i) / d(i-1)
         d(i) = d(i) - c(i) * e(i-1)
         e(i) = e(i) - c(i) * f(i-1)
         f(i) = f(i) - c(i) * g(i-1)
         r(i) = r(i) - c(i) * r(i-1)

       b(i+1) = b(i+1) / d(i-1)
       c(i+1) = c(i+1) - b(i+1) * e(i-1)
       d(i+1) = d(i+1) - b(i+1) * f(i-1)
       e(i+1) = e(i+1) - b(i+1) * g(i-1)
       r(i+1) = r(i+1) - b(i+1) * r(i-1)

       a(i+2) = a(i+2) / d(i-1)
       b(i+2) = b(i+2) - a(i+2) * e(i-1)
       c(i+2) = c(i+2) - a(i+2) * f(i-1)
       d(i+2) = d(i+2) - a(i+2) * g(i-1)
       r(i+2) = r(i+2) - a(i+2) * r(i-1)
      ENDDO
       c(n-1) = c(n-1) / d(n-2)
       d(n-1) = d(n-1) - c(n-1) * e(n-2)
       e(n-1) = e(n-1) - c(n-1) * f(n-2)
       r(n-1) = r(n-1) - c(n-1) * r(n-2)

         c(n) = c(n) / d(n-1)
         d(n) = d(n) - c(n) * e(n-1)
         r(n) = r(n) - c(n) * r(n-1)

      ! Backward substitution
         x(n) = r(n) / d(n)
       x(n-1) = (r(n-1) - e(n-1) * x(n)) / d(n-1)
       x(n-2) = (r(n-2) - e(n-2) * x(n-1) - f(n-2) * x(n)) / d(n-2)

      DO i = n-3, 1, -1
         x(i) = (r(i) - e(i) * x(i+1) - f(i) * x(i+2) - g(i) * x(i+3)) / d(i)
      ENDDO

      END SUBROUTINE HDMA
