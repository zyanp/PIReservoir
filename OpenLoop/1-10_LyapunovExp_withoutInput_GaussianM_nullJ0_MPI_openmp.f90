program runreservoir
  
  implicit none
  include 'mpif.h'
  integer     (8) :: PETOT, myrank, ierr
  integer     (8), parameter :: NN = 10000
  integer     (8), parameter :: N_res = 100 ! N_res = 100, 500, 1000
  integer     (8), parameter :: N_ens = 100
  integer     (8), parameter :: Nb = 50
  real        (8), parameter :: pi = dacos(-1.0d0), alpha = 0.2d0
  real        (8), parameter :: db = 0.01d0
  real        (8) :: M ( N_res, N_res ), II ( N_res, N_res ), JMbase ( N_res, N_res ), JM ( N_res, N_res )
  real        (8) :: g1 ( N_res*N_res/2 ), g2 ( N_res*N_res/2 )
  real        (8) :: g3 ( N_res*N_res/2 ), g4 ( N_res*N_res/2 )
  real        (8) :: r ( N_res ), z ( N_res ), u ( N_res )
  real        (8) :: mu, std, J0, J, unorm, LEs ( N_ens )
  integer     (8) :: ti, i_ens, iab, ia, ib, i, k, ith
  character (200) :: filepath, linebuf
  integer (4), allocatable :: seed ( : )
  integer (4) :: clock ( 8 ), seedsize
  integer (4) :: time1, time2, t_rate, t_max, t_diff
  
  
  !!!! MPI !!!!
  call MPI_INIT (ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, myrank, ierr)
  
  
  ! set random seed
  call random_seed ( size = seedsize )
  allocate ( seed ( seedsize ) )
  call random_seed ( get = seed )
  call date_and_time ( values = clock )
  seed ( seedsize ) = clock (8)
  seed ( 1 ) = clock ( 8 ) * clock ( 7 ) * clock ( 6 )
  call random_seed ( put = seed )
  deallocate ( seed )
  
  
  ! identity matrix
  II = 0.0d0
  do i = 1, N_res
    II ( i, i ) = 1.0d0
  end do
  
  !!! Echo-State Network  !!!
  do iab = 1, Nb
    ib = iab + myrank*50
    print '(A, I3, A, I3, A)', 'i_b=', ib, ' (', myrank, ')'
    
    J = 1.0d0 / ( dble(ib) * db )
    J0 = 0.0d0
    
    mu = J0 / dble(N_res)
    std = J / dsqrt( dble(N_res) )
   
    LEs = 0.0d0
    
    !$omp parallel
    !$omp do private(ith,i_ens,g1,i,k,g2,g3,g4,M,r,z,u,ti,JMbase,JM,unorm)
    do i_ens = 1, N_ens
      ith = omp_get_thread_num()
      
      ! box-muller
      call random_number( g1 )
      do i = 1, size(g1)
        if ( g1(i)==0.0d0 ) then
          call random_number( g1(i) )
        end if
      end do
      call random_number( g2 )
      g3 = dsqrt( -2.0d0*dlog(g1)) * dcos(2.0d0*pi*g2)
      g4 = dsqrt( -2.0d0*dlog(g1)) * dsin(2.0d0*pi*g2)
      
      ! N(mu, std^2)
      M(1:N_res/2,       :) = reshape( g3, (/ N_res/2, N_res /) )
      M(N_res/2+1:N_res, :) = reshape( g4, (/ N_res/2, N_res /) )
      M = mu + M * std
      
      JMbase = alpha * M + (1.0d0-alpha) * II
      
      ! initialise r and u
      call random_number( r )
      call random_number( u )
      
      !!! run reservoir !!!
      do ti = 1, NN
        z = dtanh( r )
        
        ! Jacobi matrix
        do i = 1, N_res
          do k = 1, N_res
            JM ( i, k ) = alpha * M ( i, k ) * ( 1.0d0 - z(k)*z(k) ) + (1.0d0-alpha) * II ( i, k )
          end do
        end do
        
        ! norm of u
        unorm = dsqrt( sum(u*u) )
        
        ! calculation for Lyapunov exponent
        LEs( i_ens ) = LEs ( i_ens ) + dlog( unorm )
        
        ! update u
        u = matmul( JM, u/unorm )
        
        ! update r
        r = alpha * (matmul(M, z)) + (1.0d0-alpha) * r
        
      end do
      
    end do
    !$omp end do
    !$omp end parallel
    
    LEs = LEs / dble(NN)
    
    ! output Lyapunov exponent
    write ( filepath, '( A, I4.4, A, I3.3, A )' ) &
      '../data/LyapunovExp_period1e4_withoutInput_GaussM_Nres', N_res, '_withoutBias_nullJ0Linear/ib', ib, '.dat'
    open  ( 20, file = filepath )
    do i_ens = 1, N_ens
      write ( 20, '( E21.8E3 )' ) LEs(i_ens)
    end do
    close ( 20 )
    
  end do
  
  print ' ( A )', 'all done'
  
  call MPI_FINALIZE (ierr)
  
end program runreservoir

