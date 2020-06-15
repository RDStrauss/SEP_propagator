 program SEP_propagator

! -----------------------------------------------------------------------------------
! Solves the focused transport equation for one spatial dimension, i.e. the Roelof equation.
! Developed and described by:
! *Strauss & Fichnter (2015): https://ui.adsabs.harvard.edu/abs/2015ApJ...801...29S/abstract
! *Strauss et al. (2018): https://ui.adsabs.harvard.edu/abs/2017SoPh..292...51S/abstract
! *MSc thesis of PKN Heita "Numerical investigation of solar energetic particle transport between 
! the Sun, Earth, and Mars". North-West University, South Africa.

! Email: dutoit.strauss@nwu.ac.za / dutoit.strauss@gmail.com
! https://github.com/RDStrauss/SEP_propagator

! To add in future:

! -----------------------------------------------------------------------------------   
 IMPLICIT NONE
  
 ! Initiate parameter
 ! N = number of grid cells in Z-direction, default, N = 200
 ! M = number of grd cells in pitch-angle space, default, M = 99 (has to be an odd number)
 ! Z_max = maximum extent of the Z-axis
 
  INTEGER, PARAMETER :: N = 100, M = 49
  INTEGER :: i,j, zlimiter, mulimiter, z_index, species, writer(5), injection_swtich
  REAL :: Delta_t, Delta_z, Delta_mu, CFL_coeff, V_sw, D_mumu_max, D_mumu(N,M), fin, fout
  REAL :: Z(N), MU(M), time, totaltime, L_max, anisotropy, D_mumu_dmu(N,M), Z_max = 3
  REAL :: f(N,M), f0(N,M), f00(N,M), L(N), f_omni_direc, time_printer, acceleration_time
  REAL :: A(M), B(N,M), speed, left_lim, right_lim, limiter, B_max, r_printer, escape_time
  REAL :: energy, lambda, times(5), pitch_angle_distribution(5,M), r_position

! lambda = value of the effective radial mean-free-path, lambda_rr
! energy = energy of the mono-energetic SEP distribution in MeV
  lambda = 0.12 !in AU
  energy = 0.08 !MeV
! species = a way to switch between electrons (species = 1) and protons (species = 2)  
  species = 1
! r_position is the radial position at which the solution is needed
  r_position = 1. !AU
! The total time in [h] that the code should compute
  totaltime = 10. !in hours
! V_sw is the solar wind speed in km/s
  V_sw = 400. !in km/s
! The different times (in hours) where the pitch-angle distribution should be printed out at 1AU
  times(1) = 2.5
  times(2) = 3.5
  times(3) = 4.5
  times(4) = 5.5
  times(5) = 6.5
! For a delta-function in time, injection_swtich = 1. For time-dependent injection, injection_swtich = 2
! Time-dependent injection, a Reid-Axford profile is used with the following acceleration and escape times
  injection_swtich = 2
  acceleration_time = 0.1 !in hrs
  escape_time = 1. !in hrs
! setting up the write flag
  DO i = 1, 5
    writer(i) = 0
  END DO

! Output files
! grid_z.txt ontains the values of the z grid
! grid_mu.txt contains the values of the mu-grid
! output_at_1AU.txt contains Z(z_index), time, f_omni_direc. anisotropy, where
! z_index = value of Z-grid correspnding to Earth, R = 1 AU, time = simulation time in hours, and
! f_omni_direc = omni-directional intensity at z = z_index, anisotropy = calulated anisotropy at 1AU

OPEN(500,file='grid_z.txt',status='unknown')
OPEN(600,file='grid_mu.txt',status='unknown')
OPEN(700,file='output.txt',status='unknown')
OPEN(400,file='pitch_angle_distribution.txt',status='unknown')
OPEN(666,file='model_setup.txt',status='unknown')

 ! initiate the grids  
 Delta_z = Z_max/(REAL(N) - 1)
 Delta_mu = 2./REAL(M)

 Z(1) = 0.05
 
   WRITE(500,"(1(ES18.8))") Z(1)
   
 DO i = 2, N, 1
 
  Z(i) = Z(i - 1) + Delta_z
  
     WRITE(500,"(1(ES18.8))") Z(i)
 
 END DO

  MU(1) = -1. + Delta_mu/2.
  
     WRITE(600,"(1(ES18.8))") MU(1)
     
 DO i = 2, M, 1
 
  MU(i) = MU(i - 1) + Delta_mu
  
    WRITE(600,"(1(ES18.8))") MU(i)
    
 END DO
 
 !--------------------------------------------------------
 ! Apply the initial conditions for time-independent injection

   DO i = 1, N
 
    DO j = 1, M
    
    f(i,j) = 0.0
    f00(i,j) = 0.0
     
    END DO
    
   END DO
 
 IF (injection_swtich.EQ.1) THEN
 
 DO i = 1, N
 
    DO j = 1, M

    ! This small Gaussian approximates a delta-like injecion in space
    f0(i,j) = 100.*exp(-(Z(i) - 0.05)*(Z(i) - 0.05)/0.0005)
    
    END DO
    
  END DO
 
 END IF
!--------------------------------------------------------
! Define the coefficients and time step.
! The time step is calculated from the CFL condition to keep the numerical scheme
! stable and therefore changes when any transport parameters are changed.
 
! Coefficient for focussing is B(i,j) and is of course z and mu dependent
! D_mumu is the pitch-angle diffusion coefficient and is z and mu dependent
! The units of D_mumu is [1/h]
 
 CALL DEF_COEFFICIENTS(speed,N,L,Z,M,MU,B,D_mumu,D_mumu_dmu,A,z_index,energy,lambda,species,V_sw,r_position,r_printer)
 
 L_max = MAXVAL(L)
 
 D_mumu_max = MAXVAL(D_mumu)!maximum value of D_mumu and used for CFL calculation
 
 B_max = speed/2./L_max!maximum value of the focussing coefficient, used in CFL calculation
 
! The CFL cofficients is a parameter that must be between (0,1)
 CFL_coeff = 1.

Delta_t = MIN(abs(CFL_coeff*Delta_z/speed),abs(CFL_coeff*Delta_mu/B_max),abs(CFL_coeff*Delta_mu*Delta_mu/D_mumu_max/2.))

 WRITE(*,*) 'Delta time (CFL):', Delta_t
 
!-------------------------------------------------------- 
! Only play around here when you truly understand the numerics!
!Choose a flux limiter to use: 
!0 = no limiter
!1 = minmod
!2 = van leer
!3 = superbee
 zlimiter = 2
 mulimiter = 2

 !-------------------------------------------------------
 ! Write code set-up to file
 WRITE(666,*) lambda 
 WRITE(666,*) energy
 WRITE(666,*) species
 WRITE(666,*) injection_swtich
 WRITE(666,*) V_sw
 
  IF (injection_swtich.EQ.2) THEN
 
   WRITE(666,*) acceleration_time
   WRITE(666,*) escape_time
  
  ELSE
  
   WRITE(666,*) 0.
   WRITE(666,*) 0.
  
  ENDIF
  
 WRITE(666,*) N 
 WRITE(666,*) M
 WRITE(666,*) Z_max
 WRITE(666,*) mulimiter
 WRITE(666,*) zlimiter
 !--------------------------------------------------------
 ! Do the iteration

 WRITE(*,*) 'Start the iteration...'
 
 time = 0.
 time_printer = 0.

 DO WHILE (time.LT.totaltime)
 
 time = time + Delta_t
 !======================================================
 !The option of a time-dependent injection profiles
 IF (injection_swtich.EQ.2) THEN
 
 
   DO i = 1, N
 
    DO j = 1, M
    
    ! This small Gaussian approximates a delta-like injecion in space
    f0(i,j) = f0(i,j) + 100.*exp(-(Z(i) - 0.05)*(Z(i) - 0.05)/0.0005)/time*exp(-acceleration_time/time - time/escape_time)
      
    END DO
    
   END DO
   
  END IF
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 ! First step: Diffusion - A simple time forward Euler type numerical scheme
 
  DO i = 2, N - 1
  
    DO j = 2, M - 1
    
      f(i,j) = f0(i,j) + D_mumu_dmu(i,j)*Delta_t/2./Delta_mu*(f0(i,j+1) - f0(i,j-1)) + &
      D_mumu(i,j)*Delta_t/Delta_mu/Delta_mu*(f0(i,j+1) - 2.*f0(i,j) + f0(i,j-1))
                
    END DO
    
  END DO

!Update the older f's and boundary conditions
   DO i = 1, N

! Conservation of diffusive flux as boundary condition
  f(i,1) = f0(i,1) - Delta_t/Delta_mu*(-(D_mumu(i,1)+D_mumu(i,2))/2.*(f0(i,2) - f0(i,1))/Delta_mu)
  f(i,M) = f0(i,M) + Delta_t/Delta_mu*(-(D_mumu(i,M)+D_mumu(i,M-1))/2.*(f0(i,M) - f0(i,M-1))/Delta_mu)
  
  END DO
  
  DO i = 1, N
  
    DO j = 1, M
  
        f0(i,j) = f(i,j)
  
    END DO
    
  END DO

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 ! Second step: Convection in z - upwind scheme + flux delimiter
 ! NOTE that the advection speed is multiplied by mu - it therefore changes sign!
 
  DO j = 1, M
  
  IF (A(j).GE.0.) THEN ! If the convection speed is positive
  
    DO i = 2, N - 1

      !Do a half-time step using the upwind scheme
      f(i,j) = f0(i,j) - Delta_t/Delta_z*(A(j)*f0(i,j) - A(j)*f0(i-1,j))/2.
	      
    END DO
     
    DO i = 2, N - 1 !Apply the flux limiter on the half-step fluxes
        
	  left_lim = (A(j)*f(i,j) - A(j)*f(i-1,j))/2.
	  right_lim = (A(j)*f(i+1,j) - A(j)*f(i,j))/2.

        IF (zlimiter.EQ.1) THEN
        limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(left_lim),abs(right_lim))
        END IF
        
        IF (zlimiter.EQ.2) THEN
        limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
        ENDIF
        
        IF (zlimiter.EQ.3) THEN
        IF (abs(left_lim).GE.abs(right_lim)) THEN
        
        	  limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(left_lim),abs(2.*right_lim))
        	  
	ELSE
	
		  limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(2.*left_lim),abs(right_lim))
	
	END IF
	ENDIF
        	  
	IF (zlimiter.EQ.0) THEN
	limiter = 0. 
	ENDIF
	  
	IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN and then no limiter is applied
	  
	IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs (i.e. gradients have different signs) are different
	  
	f00(i,j) = A(j)*f(i,j) + limiter !apply the limiter, i.e. correct the flux !NOTE: here f00 is a flux, while f and f0 are intensities
	        
     END DO
    
     DO i = 2, N - 1
            
	f(i,j) = f0(i,j) - Delta_t/Delta_z*(f00(i,j) - f00(i-1,j)) !Do a full time step usin the upwind scheme and the flux corrected values
          
     END DO
  
  ELSE !If convection speed is negative, basically the same of positive except that the spatial indices change
  
      DO i = 2, N - 1

      !Do a half-time step using the upwind scheme
      f(i,j) = f0(i,j) - Delta_t/Delta_z*(A(j)*f0(i+1,j) - A(j)*f0(i,j))/2.
	      
    END DO
     
    DO i = 2, N - 1 !Apply the flux limiter on the hafl-step fluxes
        
	  left_lim = (A(j)*f(i-1,j) - A(j)*f(i,j))/2.
	  right_lim = (A(j)*f(i,j) - A(j)*f(i+1,j))/2.
        
        IF (zlimiter.EQ.1) THEN
        limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(left_lim),abs(right_lim))
        END IF
        
        IF (zlimiter.EQ.2) THEN
        limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
        ENDIF
        
        IF (zlimiter.EQ.3) THEN
        IF (abs(left_lim).GE.abs(right_lim)) THEN
        
        	  limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(left_lim),abs(2.*right_lim))
        	  
	ELSE
	
		  limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(2.*left_lim),abs(right_lim))
	
	END IF
	ENDIF
        	  
	IF (zlimiter.EQ.0) THEN
	limiter = 0. 
	ENDIF
	  
	IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN
	  
	IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs are different
	  
	f00(i,j) = A(j)*f(i,j) + limiter !apply the limiter
	        
     END DO
    
     DO i = 2, N - 1
            
	f(i,j) = f0(i,j) - Delta_t/Delta_z*(f00(i+1,j) - f00(i,j)) !Do the full time step
          
     END DO
  
  END IF !If for sign of convection speed test
    
  END DO !the loop over mu

!Update the older f's and boundary conditions
   DO i = 1, N

    DO j = 1, M
    
    !Boundaries for z
  f(1,j) = f(2,j)
  f(N,j) = f(N - 1,j)
    
    f0(i,j) = f(i,j) !Update f^t-1
            
    END DO
    
  END DO 
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 ! Third step: Focussing - upwind scheme + flux delimiter
 ! Same numerical scheme as used to advect in z-direction
  
    DO i = 2, N - 1
  
      DO j = 2, M - 1
      
	  !Do a half-time step using the upwind scheme
	  f(i,j) = f0(i,j) - Delta_t/Delta_mu*(B(i,j)*f0(i,j) - B(i,j-1)*f0(i,j-1))/2.
      
    END DO
    
    f(i,1) = f0(i,1) - Delta_t/Delta_mu/2.*f0(i,1)*B(i,1)
    f(i,M) = f0(i,M) + Delta_t/Delta_mu/2.*f0(i,M-1)*B(i,M-1)

    DO j = 2, M - 1 !Apply the flux limiter on the hafl-step fluxes
      
	  left_lim = (B(i,j)*f(i,j) - B(i,j-1)*f(i,j-1))/2.
	  right_lim = (B(i,j+1)*f(i,j+1) - B(i,j)*f(i,j))/2.
        
        IF (mulimiter.EQ.1) THEN
        limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(left_lim),abs(right_lim))
        END IF
        
        IF (mulimiter.EQ.2) THEN
        limiter = 2.*left_lim*right_lim/(left_lim + right_lim)
        ENDIF
        
        IF (mulimiter.EQ.3) THEN
        IF (abs(left_lim).GE.abs(right_lim)) THEN
        
        	  limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(left_lim),abs(2.*right_lim))
        	  
	ELSE
	
		  limiter = 0.5*(left_lim/abs(left_lim) + right_lim/abs(right_lim))*min(abs(2.*left_lim),abs(right_lim))
	
	END IF
	ENDIF
        	  
	IF (mulimiter.EQ.0) THEN
	limiter = 0. 
	ENDIF
	  
	IF (limiter.NE.limiter) limiter = 0. !sometimes the limiter becomes a NaN
	  
	IF(left_lim*right_lim.LT.0.) limiter = 0. !limiter not applied near extrema where signs are different
	  
	f00(i,j) = B(i,j)*f(i,j) + limiter !apply the limiter
	        
     END DO
     
    f00(i,1) = B(i,1)*f(i,1) + 0.
    f00(i,M) = B(i,M)*f(i,M) + 0.
    
     DO j = 2, M - 1
      
	f(i,j) = f0(i,j) - Delta_t/Delta_mu*(f00(i,j) - f00(i,j-1)) !Do the full time step

    END DO
    
    f(i,1) = f0(i,1) - Delta_t/Delta_mu*f00(i,1)
    f(i,M) = f0(i,M) + Delta_t/Delta_mu*f00(i,M-1)
  
END DO !the loop over mu

!Update the older f's and boundary conditions
   DO j = 1, M

    DO i = 1, N
    
    f0(i,j) = f(i,j) !Update f^t-1
        
    END DO
    
  END DO 

  time_printer = time_printer + Delta_t

! Print some code output every ~1 min = 0.01667 hours..
!-----------------------------------------------------------------------
  IF(time_printer.GT.0.01667) THEN
  !Write the time-dependent omni-directional output at 1 AU, i.e. i = 28
  
  WRITE(*,*) 'Time: ', time, ' of ', totaltime, ' hrs'
  
  fin = 0.
  fout = 0.
  
  
    DO j = 1, M/2.

! integrate forward moving distribution
    fin = fin + f(z_index,j)+ 0.0001

   END DO
   
   
    DO j = M/2, M

! integrate forward moving distribution
    fout = fout + f(z_index,j)+ 0.0001

   END DO
  
  
  
 time_printer = 0.
 f_omni_direc = 0.
    
    DO j = 1, M

! integrate the distribution function to find the omni-directional intensity
    f_omni_direc = f_omni_direc + f(z_index,j) + 0.0001

   END DO
 
   anisotropy = 0.

     DO j = 1, M
     
! integrate the distribution function * mu to find the omni-directional intensity    
    anisotropy = anisotropy + 3.*(f(z_index,j) + 0.0001)*MU(j)

   END DO
   
   anisotropy = anisotropy/f_omni_direc
   
   WRITE(700,"(7(ES18.8))") Z(z_index),r_position, time, f_omni_direc*Delta_mu, anisotropy, fin, fout
      
    DO i = 1, 5
      
      IF ((time.GT.times(i)).AND.(writer(i).EQ.0)) THEN
    
	DO j = 1, M
	
	  pitch_angle_distribution(i,j) = f(z_index,j)
	  writer(i) = 1
      
      END DO
    
      ENDIF
    
    END DO
      
  ENDIF
  
  !--------------------------------------------------------
  ! the time loop... 
  END DO
  
  DO j = 1, M
  
    WRITE(400,"(11(ES18.8))") MU(j), times(1), pitch_angle_distribution(1,j), times(2), pitch_angle_distribution(2,j), &
    times(3), pitch_angle_distribution(3,j), times(4), pitch_angle_distribution(4,j), times(5), pitch_angle_distribution(5,j)
  
  END DO

 WRITE(*,*) 'Total time [h]:', time
 
! Close all open files! 
 CLOSE(500)
 CLOSE(600)
 CLOSE(700)
 CLOSE(400)
 CLOSE(666) 
 END
!--------------------------------------------------------
SUBROUTINE DEF_COEFFICIENTS(speed,N,L,Z,M,MU,B,D_mumu,D_mumu_dmu,A,z_index,energy,lambda,species, V_sw,r_position,r_printer)
 
 IMPLICIT NONE
 
 INTEGER :: i,j, z_index, species, N, M
 REAL :: speed, L(N), Z(N), D_0, MU(M), B(N,M), D_mumu(N,M), D_mumu_dmu(N,M), A(M)
 REAL :: V_sw, Omega, R(N), SunBeta, r0, delta_r, test, integral, psi, beta
 REAL :: energy, lambda, rigidity, r_position, r_printer, rest_energy
 REAL, PARAMETER :: speed_of_light = 7.2 !in units of AU/hr
 
 V_sw = V_sw/1.5d8*60.*60 !solar wind speed in AU/hr
 Omega = 2.*3.14/25.4/24. !solar rotation rate in /AU
 
 SunBeta = Omega/V_sw
  
 IF (species.EQ.1) THEN
 
  rest_energy = 0.51 !Electron rest mass in MeV
 
 ELSE
 
  rest_energy = 938.27 !Proton rest mass in MV
 
 ENDIF
 
 rigidity = SQRT(energy*(energy + 2.*rest_energy)) ! in MV
 beta = rigidity/(energy + rest_energy)
 speed = beta*speed_of_light
 
 delta_r = 0.0001
 
! Calculate the r-value correspnding to every z-value
  DO i = 1, N
  
      r0 = 0.
      test = 0.
  
    DO WHILE (test.LE.Z(i))
    
      test = test + SQRT(1. + SunBeta*SunBeta*r0*r0)*delta_r
      
      r0 = r0 + delta_r
    
    END DO
  
  IF ((r0.GT.r_position).AND.((r0 - 0.01).LT.r_position)) THEN
! Find the closest z values to r = r_position
  z_index = i
  
  ENDIF
  
  R(i) = r0
  
  END DO
  
  r_printer = R(z_index)
  
! L is the focussing length in units of [AU]
    DO i = 1, N
    
      L(i) = R(i)*(1. + SunBeta*SunBeta*R(i)*R(i))**(3./2.)/(2. + SunBeta*SunBeta*R(i)*R(i))
 
 END DO
 
!D_0 is just the constant in front of D_mumu, in units of [1/h], and is
! determined from the normalization to lambda
 D_0 = 1.
 
  DO j = 1, M
 
   DO i = 1, N
   
      B(i,j) = speed/2./L(i)*(1. - mu(j)*mu(j))
      
! The form fo Dmumu used. q = 1.67 is related to the turbulence spectra and H = 0.005
! gives finite scattering through mu = 0
      D_mumu(i,j) = D_0*(1. - mu(j)*mu(j))*((ABS(MU(j)))**(1.67 - 1.) + 0.05)
       
   END DO
 
 END DO
 
!Do the integral of D_uu to calculate normalization constant
  integral = 0.
 
    DO j = 1, M

      integral = integral + (1. - MU(j)*MU(j))*(1. - MU(j)*MU(j))/D_mumu(1,j)
    
    END DO
    
   integral = integral*ABS(MU(1) - MU(2))
    
    DO j = 1, M
        
    DO i = 1, N
    
	psi = ATAN(SunBeta*R(i))
	
	D_mumu(i,j) = D_mumu(i,j)*COS(psi)*COS(psi)*integral*3.*speed/8./lambda
    
    END DO
    
    END DO
    
 ! DO the derivative of D_mumu numerically
  DO i = 1, N
 
   DO j = 2, M - 1
   
        D_mumu_dmu(i,j) = (D_mumu(i,j+1) - D_mumu(i,j-1))/2./ABS(MU(2) - MU(1))
 
   END DO
    
     D_mumu_dmu(i,1) = (D_mumu(i,1) - D_mumu(i,2))/ABS(MU(2) - MU(1))
     D_mumu_dmu(i,M-1) = (D_mumu(i,M) - D_mumu(i,j-1))/ABS(MU(2) - MU(1))
 
  
 END DO

! Coefficient for z-convection is A(j), which is only mu dependent (i.e. speed is constant)
 DO j = 1, M
 
   A(j) = speed*MU(j)
 
 END DO
 
 RETURN
 
 END
!----------------------------------------------------
