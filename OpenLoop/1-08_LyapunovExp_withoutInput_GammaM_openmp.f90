program runreservoir
  use mkl_vsl_type
  use mkl_vsl
  integer     (8), parameter :: NN = 10000
  integer     (8), parameter :: N_tran = 1000
  integer     (8), parameter :: N_res = 500, N_resSQ = N_res*N_res, N_in = 1, N_out = 3
  integer     (8), parameter :: N_ens = 100
  integer     (8), parameter :: Na = 40, Nb = 40
  real        (8), parameter :: pi = dacos(-1.0d0), alpha = 0.2d0
  real        (8), parameter :: da = 0.05d0, db = 0.05d0
  real        (8) :: II ( N_res, N_res ), JMbase ( N_res, N_res ), JM ( N_res, N_res )
  real        (8), allocatable :: M ( :, : ), Mtemp ( : )
  real        (8) :: r ( N_res ), z ( N_res ), u ( N_res )
  real        (8) :: shapek, theta, aGamma = 0.0d0, xx, J0, J, unorm, LEs ( N_ens )
  integer     (8) :: ti, i_ens, iab, ia, ib, i, k, ith
  character (200) :: filepath, linebuf
  integer (4), allocatable :: seed ( : )
  integer (4) :: clock ( 8 ), seedsize
  integer     (4) :: N, LDA, LDB, NRHS, LWORK, INFO
  integer     (4), allocatable :: IPIV ( : )
  real        (8), allocatable :: WORK ( : )
  integer (4) :: time1, time2, t_rate, t_max, t_diff
  integer (4) :: method = VSL_RNG_METHOD_GAMMA_GNORM
  integer (4) :: intseed
  type (vsl_stream_state) :: newstr
  
  
  ! set random seed
  call random_seed ( size = seedsize )
  allocate ( seed ( seedsize ) )
  call random_seed ( get = seed )
  call date_and_time ( values = clock )
  seed ( seedsize ) = clock (8)
  seed ( 1 ) = clock ( 8 ) * clock ( 7 ) * clock ( 6 )
  call random_seed ( put = seed )
  deallocate ( seed )
  
  ! vslGamma
  intseed = 3698
  
  ! dsysv preparation
  N = N_res
  LDA = N
  LDB = N
  NRHS = N_out
  LWORK = (N+1)*(N+1)
  
  ! identity matrix
  II = 0.0d0
  do i = 1, N_res
    II ( i, i ) = 1.0d0
  end do
  
  !!! Echo-State Network  !!!
  do iab = 1, Na*Nb
    ia = mod(iab-1, Na) + 1
    ib = (iab-1) / Na + 1
    print '(A, I3, A, I3, A)', 'i_a=', ia, ': i_b=', ib
    
    J = 1.0d0 / ( dble(ib) * db )
    J0 = J * dble(ia) * da
    
    shapek = J0*J0/J/J / dble(N_res)
    theta  = J*J / J0

    call random_number ( xx )
    intseed = int( xx * 10000 ) + 1
    state = vslnewstream ( newstr, VSL_BRNG_MCG31, intseed )
   
    LEs = 0.0d0
     
    !$omp parallel
    !$omp do private(ith,i_ens,i,k,M,Mtemp,r,z,u,ti,JMbase,JM,unorm)
    do i_ens = 1, N_ens
      ith = omp_get_thread_num()
      allocate ( M ( N_res, N_res ), Mtemp ( N_resSQ ) )
      state = vdrnggamma ( method, newstr, N_resSQ, Mtemp, shapek, aGamma, theta )
      M = reshape( Mtemp, (/ N_res, N_res /) )
      
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
      deallocate ( M, Mtemp )
      
    end do
    !$omp end do
    !$omp end parallel

    LEs = LEs / dble(NN)
    
    ! output Lyapunov exponent
    write ( filepath, '( A, I4.4, A, I3.3, A, I3.3, A )' ) &
      '../data/LyapunovExp_period1e4_withoutInput_GammaM_Nres', N_res, '_withoutBias_J0JLinear/ia', ia, '_ib', ib, '.dat'
    open  ( 20, file = filepath )
    do i_ens = 1, N_ens
      write ( 20, '( E21.8E3 )' ) LEs(i_ens)
    end do
    close ( 20 )
    
  end do
  
  print ' ( A )', 'all done'
  
end program runreservoir

