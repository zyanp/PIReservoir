program makeLorenzData
  
  implicit none
  
  integer     (8), parameter :: NN = 100000, every = 100, NNN = NN * every
  real        (8), parameter :: sigma = 10.0d0, rho = 28.0d0, beta = 8.0d0/3.0d0, dt = 0.0001d0
  real        (8) :: xx ( 3, 1:NN )
  real        (8) :: u ( 3 ), UU ( 3 )
  integer     (8) :: t, case
  character (200) :: filepath
  
  !!! calculate numerical soltion of Lorenz system !!!
  ! set initial conditions
  u = 1.0d0
  
  ! cut the transient
  do t = 1, 1000000
    call RK4 ( u, UU, sigma, rho, beta, dt )
    u = UU
  end do
  
  ! collect timeseries
  do t = 1, NNN
    call RK4 ( u, UU, sigma, rho, beta, dt )
    u = UU
    if ( mod ( t, every ) == 0 ) then
      xx ( 1:3, t/every ) = u ( 1:3 )
    end if
  end do
  
  !!! output !!!
  write ( filepath, '( A )' ) '../data/Target/data1_xyz_lorenz1e-4_every1e2.bin'
  open  ( 81, file = filepath, form = 'unformatted', access = 'direct', recl = 8*3*NN )
  write ( 81, rec = 1 ) xx
  close ( 81 )
  
  write ( filepath, '( A )' ) '../data/Target/data1_xyz_lorenz1e-4_every1e2.dat'
  open  ( 21, file = filepath )
  do t = 1, NN
    write ( 21, '( F9.4, 3F15.8 )' ) dble(t*every)*dt, xx(1:3,t)
  end do
  close ( 21 )
  
  !!! for test data !!!
  ! cut the transient
  do t = 1, 1000000
    call RK4 ( u, UU, sigma, rho, beta, dt )
    u = UU
  end do
  
  ! collect timeseries
  do t = 1, NNN
    call RK4 ( u, UU, sigma, rho, beta, dt )
    u = UU
    if ( mod ( t, every ) == 0 ) then
      xx ( 1:3, t/every ) = u ( 1:3 )
    end if
  end do
  
  !!! output !!!
  write ( filepath, '( A )' ) '../data/Target/data2_xyz_lorenz1e-4_every1e2.bin'
  open  ( 81, file = filepath, form = 'unformatted', access = 'direct', recl = 8*3*NN )
  write ( 81, rec = 1 ) xx
  close ( 81 )

  write ( filepath, '( A )' ) '../data/Target/data2_xyz_lorenz1e-4_every1e2.dat'
  open  ( 21, file = filepath )
  do t = 1, NN
    write ( 21, '( F9.4, 3F15.8 )' ) dble(t*every)*dt, xx(1:3,t)
  end do
  close ( 21 )
  
  
  print '( A )', 'all done'
end program makeLorenzData

!!!!!!!!!!!!!!!!!

subroutine rhs( x, XX, sigma, rho, beta )
  implicit none
  real (8), intent (in)  :: x ( 3 ), sigma, rho, beta
  real (8), intent (out) :: XX ( 3 )
  
  XX(1) = sigma * (x(2) - x(1))
  XX(2) = x(1) * (rho - x(3) ) - x(2)
  XX(3) = x(1) * x(2) - beta * x(3)

end subroutine rhs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RK4 ( x, XX, sigma, rho, beta, dt )
  
  implicit none
  real    (8), intent (in)  :: x  ( 3 ), sigma, rho, beta, dt
  real    (8), intent (out) :: XX ( 3 )
  real    (8) :: temp ( 3 ), k1 ( 3 ), k2 ( 3 ), k3 ( 3 ), k4 ( 3 )
  
  ! k1
  call rhs ( x, k1, sigma, rho, beta )
  temp = x + dt * k1 / 2.0d0
  
  ! k2
  call rhs ( temp, k2, sigma, rho, beta )
  temp = x + dt * k2 / 2.0d0
  
  ! k3
  call rhs ( temp, k3, sigma, rho, beta )
  temp = x + dt * k3
  
  ! k4
  call rhs ( temp, k4, sigma, rho, beta )
  
  ! RK4
  XX = x + dt * ( k1 + 2.0d0*k2 + 2.0d0*k3 + k4 ) / 6.0d0
  
end subroutine RK4




