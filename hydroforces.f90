      PROGRAM HYDROFORCES
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      PARAMETER ( acu = 1d-9 )
      DIMENSION fg(4,4)
      LOGICAL opp, ncl
      INTEGER sample

! ============ I N P U T S ==============

        alam = 1.00d-0 ! Radius ratio
          a1 = 1d+0    ! Larger drop radius [μm]
         mur = 1d+2    ! Viscosity ratio
         ncl = .TRUE.  ! Non-continuum lubrication ON
!        ncl = .FALSE. ! Non-continuum lubrication OFF
         opp = .TRUE.  ! Orientation: opposing
!        opp = .FALSE. ! Orientation: same

! ========= L O G   D I S T R. ========== 
! Logarithmic distribution of normalized
! gap size ξ = s — 2 in JO84 notation:
      xi_min = 1d-3
      xi_max = 1d+2
      sample = 10
      dlt_xi = DLOG ( xi_max / xi_min ) / DBLE(sample-1)
           s = 2d0 + xi_min + 3d-16
! ======================================= 
      DO i = 1, sample
      CALL MAP(s,alam,al,be) ! (s,λ) —> (α,β)

! ============ M E T H O D ==============

!     CALL ONM70R(opp,al,be,fg,F1,F2,T1,T2,acu)
!     CALL ONM70T(opp,al,be,fg,F1,F2,T1,T2,acu)

!     CALL RM74_T(opp,al,be,a1,F1,F2,acu)
      CALL RM74_H(opp,al,be,a1,F1,F2,acu)
!     CALL RM74_G(opp,al,be,a1,F1,F2,acu)

!     CALL RSD22_T(opp,mur,mur,al,be,a1,F1,F2,acu)

! ============ O U T P U T ==============
!     WRITE(*,*) s-2d0, F1, F2, T1, T2
! ======================================= 

      STOP

! ========= L O G   D I S T R. ==========
      s = DLOG ( s - 2d0 ) + dlt_xi
      s = DEXP ( s ) + 2d0
      ENDDO
! ======================================= 

      CLOSE ( 1 )
      END PROGRAM HYDROFORCES

! ======================= S U B R O U T I N E S ========================

! ============ M A P P I N G ============
      SUBROUTINE MAP(s,alam,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
!     Zinchenko explicit transform
!     Here we map the regular JO84 (s,λ)
!     geometrical setting in two sets of
!     spherical polar coordinates (JO84 Sec. 2)
!     system to spherical bipolar, known as
!     bispherical, coordinates (ε,k), 
!     to obtain two spheres interfaces, e.g.
!     (α,β) in Stimson & Jeffery (1926)
!     notation or (ξ₁, ξ₂) in Jeffery (1915)
!     notation, using formulae given by Zinchenko
!     for example in Eq. (2.2) of his paper:
!     doi.org/10.1016/0021-8928(78)90051-5.
!     The methodology to derive these, however,
!     is not straightforward, and he kindly
!     sent it to us in a personal communication.

!     Radii ratio:
!     "k" in Zinchenko (1977) notation
!     "λ = 1/k" in JO84 notation

!     Mapping formulae: 
!     (s,λ) —> (ε,k) —> (α,β)
!     They are modified according to our problem setting.
           k = 1d0 / alam                         ! size ratio
         eps = ( s - 2d0 ) * ( alam + 1d0 ) / 2d0 ! normalized clearance: ε.a₁ = gap
      coshal = 1d0 + eps   * ( alam + eps / 2d0 ) / ( 1d0 + alam + eps )
      coshbe = 1d0 + eps/alam*( 1d0 + eps / 2d0 ) / ( 1d0 + alam + eps )
!     coshal = 1d0 + eps  *   ( 1d0 + k * eps / 2d0 ) / ( 1d0 + k * ( 1d0 + eps ) )
!     coshbe = 1d0 + eps * k**2 * ( 1d0 + eps / 2d0 ) / ( 1d0 + k * ( 1d0 + eps ) )
      al = DACOSH(coshal)
      be =-DACOSH(coshbe)
!     be = asinh(-k*DSINH(al))
!     WRITE(*,*) "ε, k, α, β =",eps,k,al,be

      END SUBROUTINE
! ======================================================================

! === O'Neill & Majumdar (1970) - Rotation =============================
!     O'Neill, M. E., & Majumdar, R. (1970). Asymmetrical slow viscous fluid motions caused by the translation or rotation of two spheres. Part I: The determination of exact solutions for any values of the ratio of radii and separation parameters. Zeitschrift für angewandte Mathematik und Physik ZAMP, 21(2), 164-179.
      SUBROUTINE ONM70R(opp,al,be,fg,F1,F2,T1,T2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: An, Bn
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: T
      DOUBLE PRECISION :: fg(4,4)
      LOGICAL opp

      alam =-DSINH(al)/DSINH(be)
      rlam = 1d0/alam

      F1o= 0d0
      F2o= 0d0
      T1o= 0d0
      T2o= 0d0
      iN = 25
1     ALLOCATE ( T(2*iN,8), An(-1:iN+1), Bn(-1:iN+1) )
       j = 1 ! first droplet
2      T = 0d0
      DO i = 1, iN
         n = DBLE(i)
         fctr1 = (2d0*n-1d0)*alpha(n-1d0,al,be)-(2d0*n-3d0)*alpha(n,al,be)+2d0*kappa(al,be)
         fctr2 = (2d0*n+5d0)*alpha(n,al,be)-(2d0*n+3d0)*alpha(n+1d0,al,be)+2d0*kappa(al,be)
         fctr3 = (2d0*n-1d0)*(gamm(n-1d0,al,be)-dlta(n-1d0,al,be))-(2d0*n-3d0)*(gamm(n,al,be)-dlta(n,al,be))-4d0
         fctr4 = (2d0*n+5d0)*(gamm(n,al,be)-dlta(n,al,be))-(2d0*n+3d0)*(gamm(n+1d0,al,be)-dlta(n+1d0,al,be))+4d0
         IF (i.GT.1) THEN
         T(2*i-1,2) = fctr3 * (n-1d0)/(2d0*n-1d0) ! B_n-1 (39)
         T(2*i-1,3) = fctr1 * (n-1d0)/(2d0*n-1d0) ! A_n-1 (39)
         ENDIF
         T(2*i-1,4) = fctr3 * -n/(2d0*n+1d0) - fctr4 * (n+1d0)/(2d0*n+1d0) ! B_n (39)
         T(2*i-1,5) = fctr1 * -n/(2d0*n+1d0) - fctr2 * (n+1d0)/(2d0*n+1d0) ! A_n (39)
         IF (i.LT.iN) THEN
         T(2*i-1,6) = fctr4 * (n+2d0)/(2d0*n+3d0) ! B_n+1 (39)
         T(2*i-1,7) = fctr2 * (n+2d0)/(2d0*n+3d0) ! A_n+1 (39)
         ENDIF
         T(2*i-1,8) = D_1R(n,al,be)

         fctr1 = (2d0*n-1d0)*(gamm(n-1d0,al,be)+dlta(n-1d0,al,be))-(2d0*n-3d0)*(gamm(n,al,be)+dlta(n,al,be))+4d0
         fctr2 = (2d0*n+5d0)*(gamm(n,al,be)+dlta(n,al,be))-(2d0*n+3d0)*(gamm(n+1d0,al,be)+dlta(n+1d0,al,be))-4d0
         fctr3 = (2d0*n-1d0)*alpha(n-1d0,al,be)-(2d0*n-3d0)*alpha(n,al,be)-2d0*kappa(al,be)
         fctr4 = (2d0*n+5d0)*alpha(n,al,be)-(2d0*n+3d0)*alpha(n+1d0,al,be)-2d0*kappa(al,be)

         IF (i.GT.1) THEN
         T(2*i,1) = fctr3 * (n-1d0)/(2d0*n-1d0) ! B_n-1 (40)
         T(2*i,2) = fctr1 * (n-1d0)/(2d0*n-1d0) ! A_n-1 (40)
         ENDIF
         T(2*i,3) = fctr3 * -n/(2d0*n+1d0) - fctr4 * (n+1d0)/(2d0*n+1d0) ! B_n (40)
         T(2*i,4) = fctr1 * -n/(2d0*n+1d0) - fctr2 * (n+1d0)/(2d0*n+1d0) ! A_n (40)
         IF (i.LT.iN) THEN
         T(2*i,5) = fctr4 * (n+2d0)/(2d0*n+3d0) ! B_n+1 (40)
         T(2*i,6) = fctr2 * (n+2d0)/(2d0*n+3d0) ! A_n+1 (40)
         ENDIF
         T(2*i,8) = D_2R(n,al,be)
      ENDDO

      CALL THOMAS(2*iN,3,3,T)

      An = 0d0
      Bn = 0d0
      DO i = 1, iN
         An(i) = T(2*i  ,8)
         Bn(i) = T(2*i-1,8)
      ENDDO

       f11 = 0d0
       g11 = 0d0
       f12 = 0d0
       g12 = 0d0
      DO i = 0, iN
         n = DBLE(i)
        En =-(kappa(al,be)+alpha(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*An(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1))-(gamm(n,al,be) &
           - dlta(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*Bn(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1))+n*(n-1d0)/(2d0*n-1d0)*Bn(i-1) &
           + (n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1)+(lambda(n,al)/DSINH(al) &
           - DSQRT(2d0)*(2d0*n+1d0)*DEXP(-(n+5d-1)*al))*DSINH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))

        Fn = (gamm(n,al,be)+dlta(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*An(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1))-(kappa(al,be) &
           - alpha(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*Bn(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1))+n*(n-1d0)/(2d0*n-1d0)*An(i-1) &
           + (n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1)-(lambda(n,al)/DSINH(al) &
           - DSQRT(2d0)*(2d0*n+1d0)*DEXP(-(n+5d-1)*al))*DCOSH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))

       f11 = f11 + ( En + Fn )
       f12 = f12 + ( En - Fn )
       g11 = g11 + ( En + Fn ) * ( 2d0*n + 1d0 - 1d0 / DTANH( al) )
       g12 = g12 + ( En - Fn ) * ( 2d0*n + 1d0 - 1d0 / DTANH(-be) )
      ENDDO

      f11 = DSQRT(2d0)/3d0*DSINH( al)**2*f11 ! (5.7)
      f12 = DSQRT(2d0)/3d0*DSINH(-be)**2*f12 ! (5.9)
      g11 = DSQRT(2d0)/4d0*DSINH( al)**3*g11 ! (5.8)
      g12 = DSQRT(2d0)/4d0*DSINH(-be)**3*g12 ! (5.10) typo !

!     IF ( j .EQ. 1 ) WRITE(*,*)"f11,f12,g11,g12=",f11,f12,g11,g12
!     IF ( j .EQ. 2 ) WRITE(*,*)"f12,f11,g12,g11=",f12,f11,g12,g11

      IF ( j .EQ. 1 ) THEN

          YB11 = f11*3d0/2d0
          YB12 = f12*6d0/(1d0+rlam)**2
          YC11 = g11
          YC21 =-g12*8d0/(1d0+rlam)**3
              
       fg(1,1) = f11
       fg(1,2) = f12
       fg(1,3) = g11
       fg(1,4) = g12
           tmp =-al
            al =-be
            be = tmp
             j = 2 ! second droplet
            GOTO 2
      ENDIF
      DEALLOCATE ( T, An, Bn )
      fg(2,1) = f12
      fg(2,2) = f11
      fg(2,3) = g12
      fg(2,4) = g11
           tmp=-al
           al =-be
           be = tmp

         YB22 =-f11*3d0/2d0
         YB21 =-f12*6d0/(1d0+alam)**2
         YC22 = g11
         YC12 =-g12*8d0/(1d0+alam)**3

!     Normalization: —6πμa²Ω  &  —8πμa³Ω
      IF (opp) THEN
         F1 =-fg(1,1) - fg(2,1)
         F2 = fg(1,2) + fg(2,2)
         T1 = fg(1,3) + fg(2,3)
         T2 = fg(1,4) + fg(2,4)
      ELSE
         F1 =-fg(1,1) + fg(2,1)
         F2 =-fg(1,2) + fg(2,2)
         T1 = fg(1,3) - fg(2,3)
         T2 =-fg(1,4) + fg(2,4)
      ENDIF

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2),DABS(T1-T1o)/DABS(T1),DABS(T2-T2o)/DABS(T2))
      IF ( rel .GT. acu ) THEN
         iN = INT(1.5 * FLOAT(iN)) ! 50% increase
        F1o = F1
        F2o = F2
        T1o = T1
        T2o = T2
!       WRITE(*,*) "n_max, F1, F2, T1, T2 = ", iN,F1,F2,T1,T2,rel
        GOTO 1
      ENDIF

!     WRITE(1,*) YB11, YB12, YB21, YB22
!     WRITE(1,*) YC11, YC12, YC21, YC22

      END SUBROUTINE
! ======================================================================

! === O'Neill & Majumdar (1970) - Translation ==========================
!     O'Neill, M. E., & Majumdar, R. (1970). Asymmetrical slow viscous fluid motions caused by the translation or rotation of two spheres. Part I: The determination of exact solutions for any values of the ratio of radii and separation parameters. Zeitschrift für angewandte Mathematik und Physik ZAMP, 21(2), 164-179.
      SUBROUTINE ONM70T(opp,al,be,fg,F1,F2,T1,T2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: An, Bn
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: T
      DOUBLE PRECISION :: fg(4,4)
      LOGICAL opp

      F1o= 0d0
      F2o= 0d0
      T1o= 0d0
      T2o= 0d0
      iN = 25
1     ALLOCATE ( T(2*iN,8), An(-1:iN+1), Bn(-1:iN+1) )
       j = 1 ! first droplet
2      T = 0d0

      DO i = 1, iN
         n = DBLE(i)
         fctr1 = (2d0*n-1d0)*alpha(n-1d0,al,be)-(2d0*n-3d0)*alpha(n,al,be)+2d0*kappa(al,be)
         fctr2 = (2d0*n+5d0)*alpha(n,al,be)-(2d0*n+3d0)*alpha(n+1d0,al,be)+2d0*kappa(al,be)
         fctr3 = (2d0*n-1d0)*(gamm(n-1d0,al,be)-dlta(n-1d0,al,be))-(2d0*n-3d0)*(gamm(n,al,be)-dlta(n,al,be))-4d0
         fctr4 = (2d0*n+5d0)*(gamm(n,al,be)-dlta(n,al,be))-(2d0*n+3d0)*(gamm(n+1d0,al,be)-dlta(n+1d0,al,be))+4d0
         IF (i.GT.1) THEN
         T(2*i-1,2) = fctr3 * (n-1d0)/(2d0*n-1d0) ! B_n-1 (39)
         T(2*i-1,3) = fctr1 * (n-1d0)/(2d0*n-1d0) ! A_n-1 (39)
         ENDIF
         T(2*i-1,4) = fctr3 * -n/(2d0*n+1d0) - fctr4 * (n+1d0)/(2d0*n+1d0) ! B_n (39)
         T(2*i-1,5) = fctr1 * -n/(2d0*n+1d0) - fctr2 * (n+1d0)/(2d0*n+1d0) ! A_n (39)
         IF (i.LT.iN) THEN
         T(2*i-1,6) = fctr4 * (n+2d0)/(2d0*n+3d0) ! B_n+1 (39)
         T(2*i-1,7) = fctr2 * (n+2d0)/(2d0*n+3d0) ! A_n+1 (39)
         ENDIF
         T(2*i-1,8) = D_1T(n,al,be)

         fctr1 = (2d0*n-1d0)*(gamm(n-1d0,al,be)+dlta(n-1d0,al,be))-(2d0*n-3d0)*(gamm(n,al,be)+dlta(n,al,be))+4d0
         fctr2 = (2d0*n+5d0)*(gamm(n,al,be)+dlta(n,al,be))-(2d0*n+3d0)*(gamm(n+1d0,al,be)+dlta(n+1d0,al,be))-4d0
         fctr3 = (2d0*n-1d0)*alpha(n-1d0,al,be)-(2d0*n-3d0)*alpha(n,al,be)-2d0*kappa(al,be)
         fctr4 = (2d0*n+5d0)*alpha(n,al,be)-(2d0*n+3d0)*alpha(n+1d0,al,be)-2d0*kappa(al,be)
         IF (i.GT.1) THEN
         T(2*i,1) = fctr3 * (n-1d0)/(2d0*n-1d0) ! B_n-1 (40)
         T(2*i,2) = fctr1 * (n-1d0)/(2d0*n-1d0) ! A_n-1 (40)
         ENDIF
         T(2*i,3) = fctr3 * -n/(2d0*n+1d0) - fctr4 * (n+1d0)/(2d0*n+1d0) ! B_n (40)
         T(2*i,4) = fctr1 * -n/(2d0*n+1d0) - fctr2 * (n+1d0)/(2d0*n+1d0) ! A_n (40)
         IF (i.LT.iN) THEN
         T(2*i,5) = fctr4 * (n+2d0)/(2d0*n+3d0) ! B_n+1 (40)
         T(2*i,6) = fctr2 * (n+2d0)/(2d0*n+3d0) ! A_n+1 (40)
         ENDIF
         T(2*i,8) = D_2T(n,al,be)
      ENDDO

      CALL THOMAS(2*iN,3,3,T)

      An = 0d0
      Bn = 0d0
      DO i = 1, iN
         An(i) = T(2*i  ,8)
         Bn(i) = T(2*i-1,8)
      ENDDO

       f21 = 0d0
       g21 = 0d0
       f22 = 0d0
       g22 = 0d0
      DO i = 0, iN
         n = DBLE(i)
        En =-(kappa(al,be)+alpha(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*An(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1))-(gamm(n,al,be) &
           - dlta(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*Bn(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1))+n*(n-1d0)/(2d0*n-1d0)*Bn(i-1) &
           + (n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1)-DSQRT(8d0)*DEXP(-(n+5d-1)*al) &
           * DSINH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))

        Fn = (gamm(n,al,be)+dlta(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*An(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1))-(kappa(al,be) &
           - alpha(n,al,be))/2d0*(n*(n-1d0)/(2d0*n-1d0)*Bn(i-1)-(n+1d0)*(n+2d0)/(2d0*n+3d0)*Bn(i+1))+n*(n-1d0)/(2d0*n-1d0)*An(i-1) &
           + (n+1d0)*(n+2d0)/(2d0*n+3d0)*An(i+1)+DSQRT(8d0)*DEXP(-(n+5d-1)*al) &
           * DCOSH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))

       f21 = f21 + ( En + Fn )
       f22 = f22 + ( En - Fn )
       g21 = g21 + ( En + Fn ) * ( 2d0*n + 1d0 - 1d0 / DTANH( al) )
       g22 = g22 + ( En - Fn ) * ( 2d0*n + 1d0 - 1d0 / DTANH(-be) )
      ENDDO

      f21 = DSQRT(2d0)/3d0*DSINH( al)  * f21 ! (5.11)
      f22 = DSQRT(2d0)/3d0*DSINH(-be)  * f22 ! (5.11)
      g21 = DSQRT(2d0)/4d0*DSINH( al)**2*g21 ! (5.12)
      g22 = DSQRT(2d0)/4d0*DSINH(-be)**2*g22 ! (5.13)

!     IF ( j .EQ. 1 ) WRITE(*,*)"f21,f22,g21,g22=",f21,f22,g21,g22
!     IF ( j .EQ. 2 ) WRITE(*,*)"f22,f21,g22,g21=",f22,f21,g22,g21

      IF ( j .EQ. 1 ) THEN ! first droplet
         fg(3,1) = f21
         fg(3,2) = f22
         fg(3,3) = g21
         fg(3,4) = g22
         IF (opp) THEN
         F1 = f21
         F2 =-f22
         T1 =-g21
         T2 =-g22
         ELSE
         F1 = f21
         F2 = f22
         T1 =-g21
         T2 = g22
         ENDIF
        tmp =-al
         al =-be
         be = tmp
          j = 2 ! second droplet
         GOTO 2
      ENDIF
      DEALLOCATE ( T, An, Bn )
      tmp=-al
      al =-be
      be = tmp
      fg(4,1) = f22
      fg(4,2) = f21
      fg(4,3) = g22
      fg(4,4) = g21

!     Normalization: —6πμaV  &  —8πμa²V
      IF (opp) THEN
         F1 = F1 - f22
         F2 = F2 + f21
         T1 = T1 + g22
         T2 = T2 + g21
      ELSE
         F1 = F1 + f22
         F2 = F2 + f21
         T1 = T1 - g22
         T2 = T2 + g21
      ENDIF

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2),DABS(T1-T1o)/DABS(T1),DABS(T2-T2o)/DABS(T2))
      IF ( rel .GT. acu ) THEN
         iN = INT(1.3 * FLOAT(iN)) ! 30% increase
        F1o = F1
        F2o = F2
        T1o = T1
        T2o = T2
!       WRITE(*,*) "n_max, F1, F2, T1, T2 = ", iN,F1,F2,T1,T2,rel
        GOTO 1
      ENDIF

      END SUBROUTINE
! ======================================================================

! === F U N C T I O N S :  O M (1970) ==================================

      FUNCTION alpha(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      alpha = DSINH(al-be)/DSINH(al)/DSINH(be)/DSINH((n+5d-1)*(al-be))
      alpha = alpha * DSINH((n+5d-1)*(al+be))
      END

      FUNCTION gamm(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      gamm = DSINH(al-be)/DSINH(al)/DSINH(be)/DSINH((n+5d-1)*(al-be))
      gamm = gamm * DCOSH((n+5d-1)*(al+be))
      END

      FUNCTION dlta(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      dlta = DSINH(al-be)/DSINH(al)/DSINH(be)/DSINH((n+5d-1)*(al-be))
      dlta = dlta * DCOSH((n+5d-1)*(al-be))
      END

      FUNCTION kappa(al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      kappa = DSINH(al+be)/DSINH(al)/DSINH(be)
      END

      FUNCTION lambda(n,al)
      IMPLICIT DOUBLE PRECISION (A-Z)
      lambda = DEXP(-(n+5d-1)*al)/DSQRT(2d0)
      lambda = lambda * ( DEXP(-al)/(2d0*n+3d0) - DEXP(al)/(2d0*n-1d0) )
      END

      FUNCTION D_1R(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D_1R =-(2d0*n+1d0)**2*DSINH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))
      D_1R = D_1R * ( DEXP(al)/(2d0*n-1d0) + DEXP(-al)/(2d0*n+3d0) )
      D_1R = D_1R + (2d0*n-1d0)*DSINH((n-05d-1)*be)/DSINH((n-05d-1)*(al-be))
      D_1R = D_1R + (2d0*n+3d0)*DSINH((n+15d-1)*be)/DSINH((n+15d-1)*(al-be))
      D_1R = D_1R * DSQRT(8d0) *DEXP(-(n+05d-1)*al)/(2d0*n+1d0)/DSINH(al)
      END

      FUNCTION D_2R(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D_2R =-(2d0*n+1d0)**2*DCOSH((n+5d-1)*be)/DSINH((n+5d-1)*(al-be))
      D_2R = D_2R * ( DEXP(al)/(2d0*n-1d0) + DEXP(-al)/(2d0*n+3d0) )
      D_2R = D_2R + (2d0*n-1d0)*DCOSH((n-05d-1)*be)/DSINH((n-05d-1)*(al-be))
      D_2R = D_2R + (2d0*n+3d0)*DCOSH((n+15d-1)*be)/DSINH((n+15d-1)*(al-be))
      D_2R = D_2R * DSQRT(8d0) *DEXP(-(n+05d-1)*al)/(2d0*n+1d0)/DSINH(al)
      END

      FUNCTION D_1T(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D_1T =  2d0 * DEXP(-(n+05d-1)*al)*DSINH((n+05d-1)*be)/DSINH((n+05d-1)*(al-be))
      D_1T = D_1T - DEXP(-(n-05d-1)*al)*DSINH((n-05d-1)*be)/DSINH((n-05d-1)*(al-be))
      D_1T = D_1T - DEXP(-(n+15d-1)*al)*DSINH((n+15d-1)*be)/DSINH((n+15d-1)*(al-be))
      D_1T = D_1T * DSQRT(8d0)
      END

      FUNCTION D_2T(n,al,be)
      IMPLICIT DOUBLE PRECISION (A-Z)
      D_2T =  2d0 * DEXP(-(n+05d-1)*al)*DCOSH((n+05d-1)*be)/DSINH((n+05d-1)*(al-be))
      D_2T = D_2T - DEXP(-(n-05d-1)*al)*DCOSH((n-05d-1)*be)/DSINH((n-05d-1)*(al-be))
      D_2T = D_2T - DEXP(-(n+15d-1)*al)*DCOSH((n+15d-1)*be)/DSINH((n+15d-1)*(al-be))
      D_2T = D_2T * DSQRT(8d0)
      END

! ======================================================================

! === Reed & Morrison (1974) ===========================================
!     Reed, L. D., & Morrison Jr, F. A. (1974). Particle interactions in viscous flow at small values of Knudsen number. Journal of Aerosol Science, 5(2), 175-189.
      SUBROUTINE RM74_T(opp,xi1,xi2,a1,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T
      LOGICAL opp

      Cm = 1.0d0 ! Momentum accommodation (assumed)
      mfp= 0.1d0 ! Air mean free path
       c = a1 * DSINH(xi1)
      a2 = a1 * DSINH(xi1) /-DSINH(xi2)
      Clc= Cm * mfp / c
      C1 = Cm * mfp / a1
      C2 = Cm * mfp / a2

      F1o= 0d0
      F2o= 0d0
      iN = 32
1     CALL CPU_TIME(t0)
      ALLOCATE ( T(4*iN,14) )
       T = 0d0

      DO i = 1, iN
         n = DBLE(i)
         k = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)

         ! (19): ξ₁
         T(4*i-3,08) = DCOSH((n-05d-1)*xi1)
         T(4*i-3,09) = DSINH((n-05d-1)*xi1)
         T(4*i-3,10) = DCOSH((n+15d-1)*xi1)
         T(4*i-3,11) = DSINH((n+15d-1)*xi1)
         T(4*i-3,14)=-k*( (2d0*n+3d0)*DEXP(-(n-05d-1)*xi1) - (2d0*n-1d0)*DEXP(-(n+15d-1)*xi1) )
         ! the factors "Uc^2/sq(2)" are ignored due to (22) & (23)

         ! (19): ξ₂
         T(4*i-2,07) = DCOSH((n-05d-1)*xi2)
         T(4*i-2,08) = DSINH((n-05d-1)*xi2)
         T(4*i-2,09) = DCOSH((n+15d-1)*xi2)
         T(4*i-2,10) = DSINH((n+15d-1)*xi2)
         T(4*i-2,14) =-k*( (2d0*n+3d0)*DEXP((n-05d-1)*xi2) - (2d0*n-1d0)*DEXP((n+15d-1)*xi2) )
         IF (opp) T(4*i-2,14) = -T(4*i-2,14)
         
         ! (20): ξ₁
         IF (i.GT.1) THEN
         T(4*i-1,02) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DCOSH((n-15d-1)*xi1)
         T(4*i-1,03) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DSINH((n-15d-1)*xi1)
         T(4*i-1,04) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi1)
         T(4*i-1,05) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi1)
         ENDIF
         T(4*i-1,06) = 2d0*(2d0*n-1d0)*DSINH((n-05d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n-1d0)**2*DCOSH((n-05d-1)*xi1)
         T(4*i-1,07) = 2d0*(2d0*n-1d0)*DCOSH((n-05d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n-1d0)**2*DSINH((n-05d-1)*xi1)
         T(4*i-1,08) = 2d0*(2d0*n+3d0)*DSINH((n+15d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n+3d0)**2*DCOSH((n+15d-1)*xi1)
         T(4*i-1,09) = 2d0*(2d0*n+3d0)*DCOSH((n+15d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n+3d0)**2*DSINH((n+15d-1)*xi1)
         IF (i.LT.iN) THEN
         T(4*i-1,10) =-Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi1)
         T(4*i-1,11) =-Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi1)
         T(4*i-1,12) =-Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DCOSH((n+25d-1)*xi1)
         T(4*i-1,13) =-Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DSINH((n+25d-1)*xi1)
         ENDIF
         T(4*i-1,14) =-Clc*n*(n+1d0)*( DCOSH(xi1)* &
         ((2d0*n-1d0)*DEXP(-(n-05d-1)*xi1)-(2d0*n+3d0)*DEXP(-(n+15d-1)*xi1)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP(-(n-15d-1)*xi1)-(2d0*n+1d0)*DEXP(-(n+05d-1)*xi1)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP(-(n+05d-1)*xi1)-(2d0*n+5d0)*DEXP(-(n+25d-1)*xi1)))&
         +2d0*n*(n+1d0)      *            (DEXP(-(n-05d-1)*xi1)      -      DEXP(-(n+15d-1)*xi1) )

         ! (20): ξ₂
         IF (i.GT.1) THEN
         T(4*i,01) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DCOSH((n-15d-1)*xi2)
         T(4*i,02) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DSINH((n-15d-1)*xi2)
         T(4*i,03) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi2)
         T(4*i,04) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi2)
         ENDIF
         T(4*i,05) = 2d0*(2d0*n-1d0)*DSINH((n-05d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n-1d0)**2*DCOSH((n-05d-1)*xi2)
         T(4*i,06) = 2d0*(2d0*n-1d0)*DCOSH((n-05d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n-1d0)**2*DSINH((n-05d-1)*xi2)
         T(4*i,07) = 2d0*(2d0*n+3d0)*DSINH((n+15d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n+3d0)**2*DCOSH((n+15d-1)*xi2)
         T(4*i,08) = 2d0*(2d0*n+3d0)*DCOSH((n+15d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n+3d0)**2*DSINH((n+15d-1)*xi2)
         IF (i.LT.iN) THEN
         T(4*i,09) =+Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi2)
         T(4*i,10) =+Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi2)
         T(4*i,11) =+Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DCOSH((n+25d-1)*xi2)
         T(4*i,12) =+Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DSINH((n+25d-1)*xi2)
         ENDIF
         IF (opp) THEN
         T(4*i,14) =-Clc*n*(n+1d0)*( DCOSH(xi2)* &
         ((2d0*n-1d0)*DEXP((n-05d-1)*xi2)-(2d0*n+3d0)*DEXP((n+15d-1)*xi2)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP((n-15d-1)*xi2)-(2d0*n+1d0)*DEXP((n+05d-1)*xi2)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP((n+05d-1)*xi2)-(2d0*n+5d0)*DEXP((n+25d-1)*xi2)))&
         +2d0*n*(n+1d0)      *            (DEXP((n-05d-1)*xi2)      -      DEXP((n+15d-1)*xi2) )
         ELSE
         T(4*i,14) =+Clc*n*(n+1d0)*( DCOSH(xi2)* &
         ((2d0*n-1d0)*DEXP((n-05d-1)*xi2)-(2d0*n+3d0)*DEXP((n+15d-1)*xi2)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP((n-15d-1)*xi2)-(2d0*n+1d0)*DEXP((n+05d-1)*xi2)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP((n+05d-1)*xi2)-(2d0*n+5d0)*DEXP((n+25d-1)*xi2)))&
         -2d0*n*(n+1d0)      *            (DEXP((n-05d-1)*xi2)      -      DEXP((n+15d-1)*xi2) )
         ENDIF
      ENDDO
      CALL CPU_TIME(t1)
      CALL THOMAS(4*iN,7,5,T)
      CALL CPU_TIME(t2)
      WRITE(*,*) 1d3*(t1-t0), 1d3*(t2-t1), 1d3*(t2-t0)

      F1 = SUM(T(:,14))
      F2 = 0d0
      DO i = 1, iN
         an = T(4*i-3,14)
         bn = T(4*i-2,14)
         cn = T(4*i-1,14)
         dn = T(4*i  ,14)
         F2 = F2 + an - bn + cn - dn
      ENDDO
      DEALLOCATE ( T )

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
      IF ( iN .LT. 16384 ) THEN
         F1o = F1
         F2o = F2
          iN = INT(2 * FLOAT(iN)) ! 50% increase
         GOTO 1
      ENDIF

!     Normalized by (1): —6πμaᵢVᵢ(1+2Cl/aᵢ)/(1+3Cl/aᵢ)
      F1 = DSINH(xi1)/3d0 * (1d0+3d0*C1)/(1d0+2d0*C1) * DABS(F1)  ! (22)
      F2 =-DSINH(xi2)/3d0 * (1d0+3d0*C2)/(1d0+2d0*C2) * DABS(F2)  ! (23)

      END SUBROUTINE
! ======================================================================                 

! === Reed & Morrison (1974) ===========================================
!     Reed, L. D., & Morrison Jr, F. A. (1974). Particle interactions in viscous flow at small values of Knudsen number. Journal of Aerosol Science, 5(2), 175-189.
      SUBROUTINE RM74_H(opp,xi1,xi2,a1,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: a, b, c, d, e, f, g, h, o, p, q, r, s, rhs, x
      LOGICAL opp

      Cm = 1.0d0 ! Momentum accommodation (assumed)
      mfp= 0.1d0 ! Air mean free path
      cc = a1 * DSINH(xi1)
      a2 = a1 * DSINH(xi1) /-DSINH(xi2)
      Clc= Cm * mfp / cc
      C1 = Cm * mfp / a1
      C2 = Cm * mfp / a2

      F1o= 0d0
      F2o= 0d0
      iN = 32
1     CALL CPU_TIME(t0)
      ALLOCATE ( a(4*iN), b(4*iN), c(4*iN), d(4*iN), e(4*iN), f(4*iN), g(4*iN), h(4*iN) )
      ALLOCATE ( o(4*iN), p(4*iN), q(4*iN), r(4*iN), s(4*iN), rhs(4*iN), x(4*iN) )
         a = 0d0
         b = 0d0
         c = 0d0
         d = 0d0
         e = 0d0
         f = 0d0
         g = 0d0
         h = 0d0
         o = 0d0
         p = 0d0
         q = 0d0
         r = 0d0
         s = 0d0
         rhs = 0d0

      DO i = 1, iN
         n = DBLE(i)
         k = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)

         ! (19): ξ₁
         h(4*i-3) = DCOSH((n-05d-1)*xi1)
         o(4*i-3) = DSINH((n-05d-1)*xi1)
         p(4*i-3) = DCOSH((n+15d-1)*xi1)
         q(4*i-3) = DSINH((n+15d-1)*xi1)
         rhs(4*i-3)=-k*( (2d0*n+3d0)*DEXP(-(n-05d-1)*xi1) - (2d0*n-1d0)*DEXP(-(n+15d-1)*xi1) )
         ! the factors "Uc^2/sq(2)" are ignored due to (22) & (23)

         ! (19): ξ₂
         g(4*i-2) = DCOSH((n-05d-1)*xi2)
         h(4*i-2) = DSINH((n-05d-1)*xi2)
         o(4*i-2) = DCOSH((n+15d-1)*xi2)
         p(4*i-2) = DSINH((n+15d-1)*xi2)
         rhs(4*i-2) =-k*( (2d0*n+3d0)*DEXP((n-05d-1)*xi2) - (2d0*n-1d0)*DEXP((n+15d-1)*xi2) )
         IF (opp) rhs(4*i-2) = -rhs(4*i-2)

         ! (20): ξ₁
         IF (i.GT.1) THEN
         b(4*i-1) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DCOSH((n-15d-1)*xi1)
         c(4*i-1) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DSINH((n-15d-1)*xi1)
         d(4*i-1) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi1)
         e(4*i-1) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi1)
         ENDIF
         f(4*i-1) = 2d0*(2d0*n-1d0)*DSINH((n-05d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n-1d0)**2*DCOSH((n-05d-1)*xi1)
         g(4*i-1) = 2d0*(2d0*n-1d0)*DCOSH((n-05d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n-1d0)**2*DSINH((n-05d-1)*xi1)
         h(4*i-1) = 2d0*(2d0*n+3d0)*DSINH((n+15d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n+3d0)**2*DCOSH((n+15d-1)*xi1)
         o(4*i-1) = 2d0*(2d0*n+3d0)*DCOSH((n+15d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n+3d0)**2*DSINH((n+15d-1)*xi1)
         IF (i.LT.iN) THEN
         p(4*i-1) =-Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi1)
         q(4*i-1) =-Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi1)
         r(4*i-1) =-Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DCOSH((n+25d-1)*xi1)
         s(4*i-1) =-Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DSINH((n+25d-1)*xi1)
         ENDIF
         rhs(4*i-1) =-Clc*n*(n+1d0)*( DCOSH(xi1)* &
         ((2d0*n-1d0)*DEXP(-(n-05d-1)*xi1)-(2d0*n+3d0)*DEXP(-(n+15d-1)*xi1)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP(-(n-15d-1)*xi1)-(2d0*n+1d0)*DEXP(-(n+05d-1)*xi1)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP(-(n+05d-1)*xi1)-(2d0*n+5d0)*DEXP(-(n+25d-1)*xi1)))&
         +2d0*n*(n+1d0)      *            (DEXP(-(n-05d-1)*xi1)      -      DEXP(-(n+15d-1)*xi1) )

         ! (20): ξ₂
         IF (i.GT.1) THEN
         a(4*i) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DCOSH((n-15d-1)*xi2)
         b(4*i) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DSINH((n-15d-1)*xi2)
         c(4*i) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi2)
         d(4*i) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi2)
         ENDIF
         e(4*i) = 2d0*(2d0*n-1d0)*DSINH((n-05d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n-1d0)**2*DCOSH((n-05d-1)*xi2)
         f(4*i) = 2d0*(2d0*n-1d0)*DCOSH((n-05d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n-1d0)**2*DSINH((n-05d-1)*xi2)
         g(4*i) = 2d0*(2d0*n+3d0)*DSINH((n+15d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n+3d0)**2*DCOSH((n+15d-1)*xi2)
         h(4*i) = 2d0*(2d0*n+3d0)*DCOSH((n+15d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n+3d0)**2*DSINH((n+15d-1)*xi2)
         IF (i.LT.iN) THEN
         o(4*i) =+Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi2)
         p(4*i) =+Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi2)
         q(4*i) =+Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DCOSH((n+25d-1)*xi2)
         r(4*i) =+Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DSINH((n+25d-1)*xi2)
         ENDIF
         IF (opp) THEN
         rhs(4*i) =-Clc*n*(n+1d0)*( DCOSH(xi2)* &
         ((2d0*n-1d0)*DEXP((n-05d-1)*xi2)-(2d0*n+3d0)*DEXP((n+15d-1)*xi2)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP((n-15d-1)*xi2)-(2d0*n+1d0)*DEXP((n+05d-1)*xi2)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP((n+05d-1)*xi2)-(2d0*n+5d0)*DEXP((n+25d-1)*xi2)))&
         +2d0*n*(n+1d0)      *            (DEXP((n-05d-1)*xi2)      -      DEXP((n+15d-1)*xi2) )
         ELSE
         rhs(4*i) =+Clc*n*(n+1d0)*( DCOSH(xi2)* &
         ((2d0*n-1d0)*DEXP((n-05d-1)*xi2)-(2d0*n+3d0)*DEXP((n+15d-1)*xi2)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP((n-15d-1)*xi2)-(2d0*n+1d0)*DEXP((n+05d-1)*xi2)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP((n+05d-1)*xi2)-(2d0*n+5d0)*DEXP((n+25d-1)*xi2)))&
         -2d0*n*(n+1d0)      *            (DEXP((n-05d-1)*xi2)      -      DEXP((n+15d-1)*xi2) )
         ENDIF
      ENDDO
      CALL CPU_TIME(t1)
      CALL TA_SEVEN_FIVE(4*iN, a, b, c, d, e, f, g, h, o, p, q, r, s, rhs, x)
      CALL CPU_TIME(t2)
      WRITE(*,*) 1d3*(t1-t0), 1d3*(t2-t1), 1d3*(t2-t0)

      F1 = SUM(x(:))
      F2 = 0d0
      DO i = 1, iN
         an = x(4*i-3)
         bn = x(4*i-2)
         cn = x(4*i-1)
         dn = x(4*i  )
         F2 = F2 + an - bn + cn - dn
      ENDDO
      DEALLOCATE ( a, b, c, d, e, f, g, h, o, p, q, r, s, rhs, x )

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
      IF ( iN .LT. 1024 ) THEN
         F1o = F1
         F2o = F2
          iN = INT(2 * FLOAT(iN)) ! 50% increase
         GOTO 1
      ENDIF

!     Normalized by (1): —6πμaᵢVᵢ(1+2Cl/aᵢ)/(1+3Cl/aᵢ)
      F1 = DSINH(xi1)/3d0 * (1d0+3d0*C1)/(1d0+2d0*C1) * DABS(F1)  ! (22)
      F2 =-DSINH(xi2)/3d0 * (1d0+3d0*C2)/(1d0+2d0*C2) * DABS(F2)  ! (23)

      END SUBROUTINE
! ======================================================================                 

! === Reed & Morrison (1974) ===========================================
!     Reed, L. D., & Morrison Jr, F. A. (1974). Particle interactions in viscous flow at small values of Knudsen number. Journal of Aerosol Science, 5(2), 175-189.
      SUBROUTINE RM74_G(opp,xi1,xi2,a1,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T
      LOGICAL opp

      Cm = 1.0d0 ! Momentum accommodation (assumed)
      mfp= 0.1d0 ! Air mean free path
       c = a1 * DSINH(xi1)
      a2 = a1 * DSINH(xi1) /-DSINH(xi2)
      Clc= Cm * mfp / c
      C1 = Cm * mfp / a1
      C2 = Cm * mfp / a2

      F1o= 0d0
      F2o= 0d0
      iN = 32
1     CALL CPU_TIME(t0)
      ALLOCATE ( T(4*iN, 4*iN+2) )
       T = 0d0

      DO i = 1, iN
         n = DBLE(i)
         k = n*(n+1d0)/(2d0*n-1d0)/(2d0*n+3d0)

         ! (19): ξ₁
         T(4*i-3,4*i-3) = DCOSH((n-05d-1)*xi1)
         T(4*i-3,4*i-2) = DSINH((n-05d-1)*xi1)
         T(4*i-3,4*i-1) = DCOSH((n+15d-1)*xi1)
         T(4*i-3,4*i  ) = DSINH((n+15d-1)*xi1)
         T(4*i-3,4*iN+1)=-k*( (2d0*n+3d0)*DEXP(-(n-05d-1)*xi1) - (2d0*n-1d0)*DEXP(-(n+15d-1)*xi1) )
         ! the factors "Uc^2/sq(2)" are ignored due to (22) & (23)

         ! (19): ξ₂
         T(4*i-2,4*i-3) = DCOSH((n-05d-1)*xi2)
         T(4*i-2,4*i-2) = DSINH((n-05d-1)*xi2)
         T(4*i-2,4*i-1) = DCOSH((n+15d-1)*xi2)
         T(4*i-2,4*i  ) = DSINH((n+15d-1)*xi2)
         T(4*i-2,4*iN+1)=-k*( (2d0*n+3d0)*DEXP((n-05d-1)*xi2) - (2d0*n-1d0)*DEXP((n+15d-1)*xi2) )
         IF (opp) T(4*i-2,4*iN+1) = -T(4*i-2,4*iN+1)

         ! (20): ξ₁
         IF (i.GT.1) THEN
         T(4*i-1,4*i-7) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DCOSH((n-15d-1)*xi1)
         T(4*i-1,4*i-6) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DSINH((n-15d-1)*xi1)
         T(4*i-1,4*i-5) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi1)
         T(4*i-1,4*i-4) =-Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi1)
         ENDIF
         T(4*i-1,4*i-3) = 2d0*(2d0*n-1d0)*DSINH((n-05d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n-1d0)**2*DCOSH((n-05d-1)*xi1)
         T(4*i-1,4*i-2) = 2d0*(2d0*n-1d0)*DCOSH((n-05d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n-1d0)**2*DSINH((n-05d-1)*xi1)
         T(4*i-1,4*i-1) = 2d0*(2d0*n+3d0)*DSINH((n+15d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n+3d0)**2*DCOSH((n+15d-1)*xi1)
         T(4*i-1,4*i  ) = 2d0*(2d0*n+3d0)*DCOSH((n+15d-1)*xi1) + Clc*DCOSH(xi1)*(2d0*n+3d0)**2*DSINH((n+15d-1)*xi1)
         IF (i.LT.iN) THEN
         T(4*i-1,4*i+1) =-Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi1)
         T(4*i-1,4*i+2) =-Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi1)
         T(4*i-1,4*i+3) =-Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DCOSH((n+25d-1)*xi1)
         T(4*i-1,4*i+4) =-Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DSINH((n+25d-1)*xi1)
         ENDIF
         T(4*i-1,4*iN+1)=-Clc*n*(n+1d0)*( DCOSH(xi1)* &
         ((2d0*n-1d0)*DEXP(-(n-05d-1)*xi1)-(2d0*n+3d0)*DEXP(-(n+15d-1)*xi1)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP(-(n-15d-1)*xi1)-(2d0*n+1d0)*DEXP(-(n+05d-1)*xi1)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP(-(n+05d-1)*xi1)-(2d0*n+5d0)*DEXP(-(n+25d-1)*xi1)))&
         +2d0*n*(n+1d0)      *            (DEXP(-(n-05d-1)*xi1)      -      DEXP(-(n+15d-1)*xi1) )

         ! (20): ξ₂
         IF (i.GT.1) THEN
         T(4*i,4*i-7) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DCOSH((n-15d-1)*xi2)
         T(4*i,4*i-6) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n-3d0)**2*DSINH((n-15d-1)*xi2)
         T(4*i,4*i-5) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi2)
         T(4*i,4*i-4) =+Clc*(n+1d0)/(2d0*n-1d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi2)
         ENDIF
         T(4*i,4*i-3) = 2d0*(2d0*n-1d0)*DSINH((n-05d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n-1d0)**2*DCOSH((n-05d-1)*xi2)
         T(4*i,4*i-2) = 2d0*(2d0*n-1d0)*DCOSH((n-05d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n-1d0)**2*DSINH((n-05d-1)*xi2)
         T(4*i,4*i-1) = 2d0*(2d0*n+3d0)*DSINH((n+15d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n+3d0)**2*DCOSH((n+15d-1)*xi2)
         T(4*i,4*i  ) = 2d0*(2d0*n+3d0)*DCOSH((n+15d-1)*xi2) - Clc*DCOSH(xi2)*(2d0*n+3d0)**2*DSINH((n+15d-1)*xi2)
         IF (i.LT.iN) THEN
         T(4*i,4*i+1) =+Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DCOSH((n+05d-1)*xi2)
         T(4*i,4*i+2) =+Clc*n/(2d0*n+3d0)*(2d0*n+1d0)**2*DSINH((n+05d-1)*xi2)
         T(4*i,4*i+3) =+Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DCOSH((n+25d-1)*xi2)
         T(4*i,4*i+4) =+Clc*n/(2d0*n+3d0)*(2d0*n+5d0)**2*DSINH((n+25d-1)*xi2)
         ENDIF
         IF (opp) THEN
         T(4*i,4*iN+1)=-Clc*n*(n+1d0)*( DCOSH(xi2)* &
         ((2d0*n-1d0)*DEXP((n-05d-1)*xi2)-(2d0*n+3d0)*DEXP((n+15d-1)*xi2)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP((n-15d-1)*xi2)-(2d0*n+1d0)*DEXP((n+05d-1)*xi2)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP((n+05d-1)*xi2)-(2d0*n+5d0)*DEXP((n+25d-1)*xi2)))&
         +2d0*n*(n+1d0)      *            (DEXP((n-05d-1)*xi2)      -      DEXP((n+15d-1)*xi2) )
         ELSE
         T(4*i,4*iN+1)=+Clc*n*(n+1d0)*( DCOSH(xi2)* &
         ((2d0*n-1d0)*DEXP((n-05d-1)*xi2)-(2d0*n+3d0)*DEXP((n+15d-1)*xi2)) &
         -(n-1d0)/(2d0*n-1d0)*((2d0*n-3d0)*DEXP((n-15d-1)*xi2)-(2d0*n+1d0)*DEXP((n+05d-1)*xi2)) &
         -(n+2d0)/(2d0*n+3d0)*((2d0*n+1d0)*DEXP((n+05d-1)*xi2)-(2d0*n+5d0)*DEXP((n+25d-1)*xi2)))&
         -2d0*n*(n+1d0)      *            (DEXP((n-05d-1)*xi2)      -      DEXP((n+15d-1)*xi2) )
         ENDIF
      ENDDO
      CALL CPU_TIME(t1)
      CALL GAUSS(4*iN,4*iN+1,T)
!     CALL GAUSSB(4*iN,7,5,T)
!     CALL GAUSS_SEIDEL(4*iN,T)
      CALL CPU_TIME(t2)
      WRITE(*,*) 1d3*(t1-t0), 1d3*(t2-t1), 1d3*(t2-t0)

      F1 = SUM(T(:,4*iN+1))
      F2 = 0d0
      DO i = 1, iN
         an = T(4*i-3,4*iN+1)
         bn = T(4*i-2,4*iN+1)
         cn = T(4*i-1,4*iN+1)
         dn = T(4*i  ,4*iN+1)
         F2 = F2 + an - bn + cn - dn
      ENDDO
      DEALLOCATE ( T )

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
      IF ( iN .LT. 1024 ) THEN
         F1o = F1
         F2o = F2
          iN = INT( 2 * FLOAT(iN)) ! 1000% increase
         GOTO 1
      ENDIF

!     Normalized by (1): —6πμaᵢVᵢ(1+2Cl/aᵢ)/(1+3Cl/aᵢ)
      F1 = DSINH(xi1)/3d0 * (1d0+3d0*C1)/(1d0+2d0*C1) * DABS(F1)  ! (22)
      F2 =-DSINH(xi2)/3d0 * (1d0+3d0*C2)/(1d0+2d0*C2) * DABS(F2)  ! (23)

      END SUBROUTINE
! ======================================================================                 

! === Rother, Stark, Davis (2022) ======================================
!     Rother, M. A., Stark, J. K., & Davis, R. H. (2022). Gravitational collision efficiencies of small viscous drops at finite Stokes numbers and low Reynolds numbers. International Journal of Multiphase Flow, 146, 103876.
      SUBROUTINE RSD22_T(opp,mu1,mu2,eta1,eta2,a1,F1,F2,acu)
      IMPLICIT DOUBLE PRECISION (A-H,L-Z)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T
      LOGICAL opp

      IF ( mu1 .EQ. 0d0 ) mu1 = 1d-9
      IF ( mu2 .EQ. 0d0 ) mu2 = 1d-9
      CS = 1.0d0 ! Momentum accommodation (assumed)
      mfp= 0.1d0 ! Air mean free path [μm]
       c = a1 * DSINH(eta1)
      a2 = a1 * DSINH(eta1) /-DSINH(eta2)
      Clc= CS * mfp / c
      C1 = CS * mfp / a1
      C2 = CS * mfp / a2
      C1 = (1d0+2d0/3d0/mu1+2d0*C1)/(1d0+1d0/mu1+3d0*C1)
      C2 = (1d0+2d0/3d0/mu2+2d0*C2)/(1d0+1d0/mu2+3d0*C2)

      F1o= 0d0
      F2o= 0d0
      iN = 50
1     ALLOCATE ( T(4*iN,14) )
       T = 0d0

      DO i = 1, iN
         n = DBLE(i)

         ! (A.11)
         T(4*i-3,08) = 1d0
         T(4*i-3,09) = DEXP((n-05d-1)*(eta2-eta1))
         T(4*i-3,10) = 1d0
         T(4*i-3,11) = DEXP((n+15d-1)*(eta2-eta1))
         T(4*i-3,14) = DEXP((n-05d-1)*-eta1)/(2d0*n-1d0)-DEXP((n+15d-1)*-eta1)/(2d0*n+3d0) ! typo

         ! (A.12)
         T(4*i-2,07) = DEXP((n-05d-1)*(eta2-eta1))
         T(4*i-2,08) = 1d0
         T(4*i-2,09) = DEXP((n+15d-1)*(eta2-eta1))
         T(4*i-2,10) = 1d0
         T(4*i-2,14) = DEXP((n-05d-1)*eta2)/(2d0*n-1d0)-DEXP((n+15d-1)*eta2)/(2d0*n+3d0) ! typo
         IF (opp) T(4*i-2,14) = -T(4*i-2,14)

         ! (A.13)
         fct1 =-Clc*(n-1d0)/(2d0*n-1d0)
         IF (i.GT.1) THEN
         T(4*i-1,02) = fct1 * (n-15d-1)**2
         T(4*i-1,03) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*(eta2-eta1))
         T(4*i-1,04) = fct1 * (n+05d-1)**2
         T(4*i-1,05) = fct1 * (n+05d-1)**2 * DEXP((n+05d-1)*(eta2-eta1))
         ENDIF

         fct2 = Clc*DCOSH(eta1) + (mu1*(2d0*n+1d0))**-1
         T(4*i-1,06) = fct2 * (n-05d-1)**2 + (n-05d-1)
         T(4*i-1,07) = fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*(eta2-eta1)) - (n-05d-1) * DEXP((n-05d-1)*(eta2-eta1))
         T(4*i-1,08) = fct2 * (n+15d-1)**2 + (n+15d-1)
         T(4*i-1,09) = fct2 * (n+15d-1)**2 * DEXP((n+15d-1)*(eta2-eta1)) - (n+15d-1) * DEXP((n+15d-1)*(eta2-eta1))

         fct3 =-Clc*(n+2d0)/(2d0*n+3d0)
         IF (i.LT.iN) THEN
         T(4*i-1,10) = fct3 * (n+05d-1)**2
         T(4*i-1,11) = fct3 * (n+05d-1)**2 * DEXP((n+05d-1)*(eta2-eta1))
         T(4*i-1,12) = fct3 * (n+25d-1)**2
         T(4*i-1,13) = fct3 * (n+25d-1)**2 * DEXP((n+25d-1)*(eta2-eta1))
         ENDIF

         T(4*i-1,14) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*-eta1)/(2d0*n-3d0) &
                     - fct1 * (n+05d-1)**2 * DEXP((n+05d-1)*-eta1)/(2d0*n+1d0) &
                     + fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*-eta1)/(2d0*n-1d0) &
                     - fct2 * (n+15d-1)**2 * DEXP((n+15d-1)*-eta1)/(2d0*n+3d0) &
                     + fct3 * (n+05d-1)**2 * DEXP((n+05d-1)*-eta1)/(2d0*n+1d0) &
                     - fct3 * (n+25d-1)**2 * DEXP((n+25d-1)*-eta1)/(2d0*n+5d0) &
                     +        (n+15d-1)    * DEXP((n+15d-1)*-eta1)/(2d0*n+3d0) &
                     -        (n-05d-1)    * DEXP((n-05d-1)*-eta1)/(2d0*n-1d0)

         ! (A.14)
         fct1 =-Clc*(n-1d0)/(2d0*n-1d0)
         IF (i.GT.1) THEN
!        T(4*i,01) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*eta2)
         T(4*i,01) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*(eta2-eta1)) ! typo
         T(4*i,02) = fct1 * (n-15d-1)**2
         T(4*i,03) = fct1 * (n+05d-1)**2 * DEXP((n+05d-1)*(eta2-eta1))
         T(4*i,04) = fct1 * (n+05d-1)**2
         ENDIF

         fct2 = Clc*DCOSH(eta2) + (mu2*(2d0*n+1d0))**-1
!        T(4*i,05) = fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*eta2) - (n-05d-1) * DEXP((n-05d-1)*(eta2-eta1))
         T(4*i,05) = fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*(eta2-eta1)) - (n-05d-1) * DEXP((n-05d-1)*(eta2-eta1)) ! typo
         T(4*i,06) = fct2 * (n-05d-1)**2 + (n-05d-1)
         T(4*i,07) = fct2 * (n+15d-1)**2 * DEXP((n+15d-1)*(eta2-eta1)) - (n+15d-1) * DEXP((n+15d-1)*(eta2-eta1))
         T(4*i,08) = fct2 * (n+15d-1)**2 + (n+15d-1)

         fct3 =-Clc*(n+2d0)/(2d0*n+3d0)
         IF (i.LT.iN) THEN
         T(4*i,09) = fct3 * (n+05d-1)**2 * DEXP((n+05d-1)*(eta2-eta1))
         T(4*i,10) = fct3 * (n+05d-1)**2
         T(4*i,11) = fct3 * (n+25d-1)**2 * DEXP((n+25d-1)*(eta2-eta1))
         T(4*i,12) = fct3 * (n+25d-1)**2
         ENDIF

         T(4*i,14) = fct1 * (n-15d-1)**2 * DEXP((n-15d-1)*eta2)/(2d0*n-3d0) &
                   - fct1 * (n+05d-1)**2 * DEXP((n+05d-1)*eta2)/(2d0*n+1d0) &
                   + fct2 * (n-05d-1)**2 * DEXP((n-05d-1)*eta2)/(2d0*n-1d0) &
                   - fct2 * (n+15d-1)**2 * DEXP((n+15d-1)*eta2)/(2d0*n+3d0) &
                   + fct3 * (n+05d-1)**2 * DEXP((n+05d-1)*eta2)/(2d0*n+1d0) & ! typo
                   - fct3 * (n+25d-1)**2 * DEXP((n+25d-1)*eta2)/(2d0*n+5d0) &
                   +        (n+15d-1)    * DEXP((n+15d-1)*eta2)/(2d0*n+3d0) &
                   -        (n-05d-1)    * DEXP((n-05d-1)*eta2)/(2d0*n-1d0)

         IF (opp) T(4*i,14) = -T(4*i,14)
      ENDDO

      CALL THOMAS(4*iN,7,5,T)

      F1 = 0d0
      F2 = 0d0
      DO i = 1, iN
         n = DBLE(i)
        En = T(4*i-3,14)
        Fn = T(4*i-2,14)
        Gn = T(4*i-1,14)
        Hn = T(4*i  ,14)
        F1 = F1 + n*(n+1d0) * ( En * DEXP(-(n-05d-1)*eta1) + Gn * DEXP(-(n+15d-1)*eta1) )
        F2 = F2 + n*(n+1d0) * ( Fn * DEXP( (n-05d-1)*eta2) + Hn * DEXP( (n+15d-1)*eta2) )
      ENDDO
      DEALLOCATE ( T )

      rel = MAX(DABS(F1-F1o)/DABS(F1),DABS(F2-F2o)/DABS(F2))
      IF ( rel .GT. acu ) THEN
          iN = INT(1.5 * FLOAT(iN)) ! 50% increase
         F1o = F1
         F2o = F2
!        WRITE(*,*) 'n_max, F1, F2, rel = ', iN,F1,F2,rel
         GOTO 1
      ENDIF

      F1 = 2d0/3d0 * DSINH(eta1) / C1 * DABS(F1)
      F2 =-2d0/3d0 * DSINH(eta2) / C2 * DABS(F2)

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
      SUBROUTINE TA_SEVEN_FIVE(n, a, b, c, d, e, f, g, h, o, p, q, r, s, rhs, x)
      DOUBLE PRECISION, DIMENSION(n) :: a, b, c, d, e, f, g, h, o, p, q, r, s, rhs, x
      INTEGER :: i,n

      ! Forward elimination
      DO i = 2, n-6
         g(i) = g(i) / h(i-1)
         h(i) = h(i) - g(i) * o(i-1)
         o(i) = o(i) - g(i) * p(i-1)
         p(i) = p(i) - g(i) * q(i-1)
         q(i) = q(i) - g(i) * r(i-1)
         r(i) = r(i) - g(i) * s(i-1)
       rhs(i) = rhs(i) - g(i) * rhs(i-1)

       f(i+1) = f(i+1) / h(i-1)
       g(i+1) = g(i+1) - f(i+1) * o(i-1)
       h(i+1) = h(i+1) - f(i+1) * p(i-1)
       o(i+1) = o(i+1) - f(i+1) * q(i-1)
       p(i+1) = p(i+1) - f(i+1) * r(i-1)
       q(i+1) = q(i+1) - f(i+1) * s(i-1)
     rhs(i+1) = rhs(i+1) - f(i+1) * rhs(i-1)

       e(i+2) = e(i+2) / h(i-1)
       f(i+2) = f(i+2) - e(i+2) * o(i-1)
       g(i+2) = g(i+2) - e(i+2) * p(i-1)
       h(i+2) = h(i+2) - e(i+2) * q(i-1)
       o(i+2) = o(i+2) - e(i+2) * r(i-1)
       p(i+2) = p(i+2) - e(i+2) * s(i-1)
     rhs(i+2) = rhs(i+2) - e(i+2) * rhs(i-1)

       d(i+3) = d(i+3) / h(i-1)
       e(i+3) = e(i+3) - d(i+3) * o(i-1)
       f(i+3) = f(i+3) - d(i+3) * p(i-1)
       g(i+3) = g(i+3) - d(i+3) * q(i-1)
       h(i+3) = h(i+3) - d(i+3) * r(i-1)
       o(i+3) = o(i+3) - d(i+3) * s(i-1)
     rhs(i+3) = rhs(i+3) - d(i+3) * rhs(i-1)

       c(i+4) = c(i+4) / h(i-1)
       d(i+4) = d(i+4) - c(i+4) * o(i-1)
       e(i+4) = e(i+4) - c(i+4) * p(i-1)
       f(i+4) = f(i+4) - c(i+4) * q(i-1)
       g(i+4) = g(i+4) - c(i+4) * r(i-1)
       h(i+4) = h(i+4) - c(i+4) * s(i-1)
     rhs(i+4) = rhs(i+4) - c(i+4) * rhs(i-1)

       b(i+5) = b(i+5) / h(i-1)
       c(i+5) = c(i+5) - b(i+5) * o(i-1)
       d(i+5) = d(i+5) - b(i+5) * p(i-1)
       e(i+5) = e(i+5) - b(i+5) * q(i-1)
       f(i+5) = f(i+5) - b(i+5) * r(i-1)
       g(i+5) = g(i+5) - b(i+5) * s(i-1)
     rhs(i+5) = rhs(i+5) - b(i+5) * rhs(i-1)

       a(i+6) = a(i+6) / h(i-1)
       b(i+6) = b(i+6) - a(i+6) * o(i-1)
       c(i+6) = c(i+6) - a(i+6) * p(i-1)
       d(i+6) = d(i+6) - a(i+6) * q(i-1)
       e(i+6) = e(i+6) - a(i+6) * r(i-1)
       f(i+6) = f(i+6) - a(i+6) * s(i-1)
     rhs(i+6) = rhs(i+6) - a(i+6) * rhs(i-1)
      ENDDO
       g(n-5) = g(n-5) / h(n-6)
       h(n-5) = h(n-5) - g(n-5) * o(n-6)
       o(n-5) = o(n-5) - g(n-5) * p(n-6)
       p(n-5) = p(n-5) - g(n-5) * q(n-6)
       q(n-5) = q(n-5) - g(n-5) * r(n-6)
       r(n-5) = r(n-5) - g(n-5) * s(n-6)
       rhs(n-5) = rhs(n-5) - g(n-5) * rhs(n-6)
       f(n-4) = f(n-4) / h(n-6)
       g(n-4) = g(n-4) - f(n-4) * o(n-6)
       h(n-4) = h(n-4) - f(n-4) * p(n-6)
       o(n-4) = o(n-4) - f(n-4) * q(n-6)
       p(n-4) = p(n-4) - f(n-4) * r(n-6)
       q(n-4) = q(n-4) - f(n-4) * s(n-6)
       rhs(n-4) = rhs(n-4) - f(n-4) * rhs(n-6)
       e(n-3) = e(n-3) / h(n-6)
       f(n-3) = f(n-3) - e(n-3) * o(n-6)
       g(n-3) = g(n-3) - e(n-3) * p(n-6)
       h(n-3) = h(n-3) - e(n-3) * q(n-6)
       o(n-3) = o(n-3) - e(n-3) * r(n-6)
       p(n-3) = p(n-3) - e(n-3) * s(n-6)
       rhs(n-3) = rhs(n-3) - e(n-3) * rhs(n-6)
       d(n-2) = d(n-2) / h(n-6)
       e(n-2) = e(n-2) - d(n-2) * o(n-6)
       f(n-2) = f(n-2) - d(n-2) * p(n-6)
       g(n-2) = g(n-2) - d(n-2) * q(n-6)
       h(n-2) = h(n-2) - d(n-2) * r(n-6)
       o(n-2) = o(n-2) - d(n-2) * s(n-6)
       rhs(n-2) = rhs(n-2) - d(n-2) * rhs(n-6)
       c(n-1) = c(n-1) / h(n-6)
       d(n-1) = d(n-1) - c(n-1) * o(n-6)
       e(n-1) = e(n-1) - c(n-1) * p(n-6)
       f(n-1) = f(n-1) - c(n-1) * q(n-6)
       g(n-1) = g(n-1) - c(n-1) * r(n-6)
       h(n-1) = h(n-1) - c(n-1) * s(n-6)
       rhs(n-1) = rhs(n-1) - c(n-1) * rhs(n-6)
       b( n ) = b( n ) / h(n-6)
       c( n ) = c( n ) - b( n ) * o(n-6)
       d( n ) = d( n ) - b( n ) * p(n-6)
       e( n ) = e( n ) - b( n ) * q(n-6)
       f( n ) = f( n ) - b( n ) * r(n-6)
       g( n ) = g( n ) - b( n ) * s(n-6)
       rhs( n ) = rhs( n ) - b( n ) * rhs(n-6)

       g(n-4) = g(n-4) / h(n-5)
       h(n-4) = h(n-4) - g(n-4) * o(n-5)
       o(n-4) = o(n-4) - g(n-4) * p(n-5)
       p(n-4) = p(n-4) - g(n-4) * q(n-5)
       q(n-4) = q(n-4) - g(n-4) * r(n-5)
       r(n-4) = r(n-4) - g(n-4) * s(n-5)
       rhs(n-4) = rhs(n-4) - g(n-4) * rhs(n-5)
       f(n-3) = f(n-3) / h(n-5)
       g(n-3) = g(n-3) - f(n-3) * o(n-5)
       h(n-3) = h(n-3) - f(n-3) * p(n-5)
       o(n-3) = o(n-3) - f(n-3) * q(n-5)
       p(n-3) = p(n-3) - f(n-3) * r(n-5)
       q(n-3) = q(n-3) - f(n-3) * s(n-5)
       rhs(n-3) = rhs(n-3) - f(n-3) * rhs(n-5)
       e(n-2) = e(n-2) / h(n-5)
       f(n-2) = f(n-2) - e(n-2) * o(n-5)
       g(n-2) = g(n-2) - e(n-2) * p(n-5)
       h(n-2) = h(n-2) - e(n-2) * q(n-5)
       o(n-2) = o(n-2) - e(n-2) * r(n-5)
       p(n-2) = p(n-2) - e(n-2) * s(n-5)
       rhs(n-2) = rhs(n-2) - e(n-2) * rhs(n-5)
       d(n-1) = d(n-1) / h(n-5)
       e(n-1) = e(n-1) - d(n-1) * o(n-5)
       f(n-1) = f(n-1) - d(n-1) * p(n-5)
       g(n-1) = g(n-1) - d(n-1) * q(n-5)
       h(n-1) = h(n-1) - d(n-1) * r(n-5)
       o(n-1) = o(n-1) - d(n-1) * s(n-5)
       rhs(n-1) = rhs(n-1) - d(n-1) * rhs(n-5)
       c( n ) = c( n ) / h(n-5)
       d( n ) = d( n ) - c( n ) * o(n-5)
       e( n ) = e( n ) - c( n ) * p(n-5)
       f( n ) = f( n ) - c( n ) * q(n-5)
       g( n ) = g( n ) - c( n ) * r(n-5)
       h( n ) = h( n ) - c( n ) * s(n-5)
       rhs( n ) = rhs( n ) - c( n ) * rhs(n-5)

       g(n-3) = g(n-3) / h(n-4)
       h(n-3) = h(n-3) - g(n-3) * o(n-4)
       o(n-3) = o(n-3) - g(n-3) * p(n-4)
       p(n-3) = p(n-3) - g(n-3) * q(n-4)
       q(n-3) = q(n-3) - g(n-3) * r(n-4)
       rhs(n-3) = rhs(n-3) - g(n-3) * rhs(n-4)
       f(n-2) = f(n-2) / h(n-4)
       g(n-2) = g(n-2) - f(n-2) * o(n-4)
       h(n-2) = h(n-2) - f(n-2) * p(n-4)
       o(n-2) = o(n-2) - f(n-2) * q(n-4)
       p(n-2) = p(n-2) - f(n-2) * r(n-4)
       rhs(n-2) = rhs(n-2) - f(n-2) * rhs(n-4)
       e(n-1) = e(n-1) / h(n-4)
       f(n-1) = f(n-1) - e(n-1) * o(n-4)
       g(n-1) = g(n-1) - e(n-1) * p(n-4)
       h(n-1) = h(n-1) - e(n-1) * q(n-4)
       o(n-1) = o(n-1) - e(n-1) * r(n-4)
       rhs(n-1) = rhs(n-1) - e(n-1) * rhs(n-4)
       d( n ) = d( n ) / h(n-4)
       e( n ) = e( n ) - d( n ) * o(n-4)
       f( n ) = f( n ) - d( n ) * p(n-4)
       g( n ) = g( n ) - d( n ) * q(n-4)
       h( n ) = h( n ) - d( n ) * r(n-4)
       rhs( n ) = rhs( n ) - d( n ) * rhs(n-4)
       
       g(n-2) = g(n-2) / h(n-3)
       h(n-2) = h(n-2) - g(n-2) * o(n-3)
       o(n-2) = o(n-2) - g(n-2) * p(n-3)
       p(n-2) = p(n-2) - g(n-2) * q(n-3)
       rhs(n-2) = rhs(n-2) - g(n-2) * rhs(n-3)
       f(n-1) = f(n-1) / h(n-3)
       g(n-1) = g(n-1) - f(n-1) * o(n-3)
       h(n-1) = h(n-1) - f(n-1) * p(n-3)
       o(n-1) = o(n-1) - f(n-1) * q(n-3)
       rhs(n-1) = rhs(n-1) - f(n-1) * rhs(n-3)
       e( n ) = e( n ) / h(n-3)
       f( n ) = f( n ) - e( n ) * o(n-3)
       g( n ) = g( n ) - e( n ) * p(n-3)
       h( n ) = h( n ) - e( n ) * q(n-3)
       rhs( n ) = rhs( n ) - e( n ) * rhs(n-3)

       g(n-1) = g(n-1) / h(n-2)
       h(n-1) = h(n-1) - g(n-1) * o(n-2)
       o(n-1) = o(n-1) - g(n-1) * p(n-2)
       rhs(n-1) = rhs(n-1) - g(n-1) * rhs(n-2)
       f( n ) = f( n ) / h(n-2)
       g( n ) = g( n ) - f( n ) * o(n-2)
       h( n ) = h( n ) - f( n ) * p(n-2)
       rhs( n ) = rhs( n ) - f( n ) * rhs(n-2)

       g(n) = g(n) / h(n-1)
       h(n) = h(n) - g(n) * o(n-1)
       rhs(n) = rhs(n) - g(n) * rhs(n-1)

      ! Backward substitution
         x(n) = rhs(n) / h(n)
       x(n-1) = (rhs(n-1) - o(n-1) * x(n)) / h(n-1)
       x(n-2) = (rhs(n-2) - o(n-2) * x(n-1) - p(n-2) * x(n)) / h(n-2)
       x(n-3) = (rhs(n-3) - o(n-3) * x(n-2) - p(n-3) * x(n-1) - q(n-3) * x(n)) / h(n-3)
       x(n-4) = (rhs(n-4) - o(n-4) * x(n-3) - p(n-4) * x(n-2) - q(n-4) * x(n-1) - r(n-4) * x(n)) / h(n-4)
       x(n-5) = (rhs(n-5) - o(n-5) * x(n-4) - p(n-5) * x(n-3) - q(n-5) * x(n-2) - r(n-5) * x(n-1) - s(n-5) * x(n)) / h(n-5)
       x(n-6) = (rhs(n-6) - o(n-6) * x(n-5) - p(n-6) * x(n-4) - q(n-6) * x(n-3) - r(n-6) * x(n-2) - s(n-6) * x(n-1)) / h(n-6)

      DO i = n-7, 1, -1
         x(i) = (rhs(i) - o(i) * x(i+1) - p(i) * x(i+2) - q(i) * x(i+3) - r(i) * x(i+4) - s(i) * x(i+5)) / h(i)
      ENDDO

      END SUBROUTINE TA_SEVEN_FIVE
