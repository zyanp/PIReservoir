program runreservoir
  !$use omp_lib
  use mkl_vsl_type
  use mkl_vsl
  ! implicit none
  integer     (8), parameter :: NN = 100000
  integer     (8), parameter :: N_tran = 1000
  integer     (8), parameter :: N_res = 500 ! N_res = 100, 500, 1000
  integer     (8), parameter :: N_resSQ = N_res * N_res
  integer     (8), parameter :: N_ens = 100, N_in = 1, N_out = 3
  integer     (8), parameter :: lambexp = -2
  integer     (8), parameter :: Na = 40, Nb = 10, myrank = 0 ! myrank = 0, 1, 2, 3
  real        (8), parameter :: pi = dacos(-1.0d0), alpha = 0.2d0, lamb = dble(10**lambexp)
  real        (8), parameter :: da = 0.05d0, db = 0.05d0
  real        (8) :: X1 ( 3, NN ), X2 ( 3, NN ), W_in ( N_res )
  real        (8) :: MSEtrain ( 3, N_ens ), MSEtest ( 3, N_ens )
  real        (8), allocatable :: r1s ( :, : ), r2s ( :, : ), M ( :, : ), Mtemp ( : ), RR ( :, : ), RY ( :, : )
  real        (8), allocatable :: X1p ( :, : ), X2p ( :, : ), sgnM ( :, : )
  real        (8) :: r1 ( N_res ), r2 ( N_res ), Z1 ( N_res ), Z2 ( N_res )
  real        (8) :: shapek, theta, aGamma=0.0d0, J0, J, mean1 ( 3 ), var1 ( 3 ), mean2 ( 3 ), var2 ( 3 )
  integer     (8) :: ti, i_ens, iab, ia, ib, i, k, ith
  character (200) :: filepath, linebuf
  integer (4), allocatable :: seed ( : )
  integer (4) :: clock ( 8 ), seedsize
  integer     (4) :: N, LDA, LDB, NRHS, LWORK, INFO
  integer     (4), allocatable :: IPIV ( : )
  real        (8), allocatable :: WORK ( : )
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
  
  
  ! read input data X1 and X2
  write ( filepath, '( A )' ) '../data/Target/data1_xyz_lorenz1e-4_every1e2.bin'
  open  ( 81, file = filepath, form = 'unformatted', access = 'direct', recl = 8*3*NN )
  read  ( 81, rec = 1 ) X1
  close ( 81 )

  write ( filepath, '( A )' ) '../data/Target/data2_xyz_lorenz1e-4_every1e2.bin'
  open  ( 81, file = filepath, form = 'unformatted', access = 'direct', recl = 8*3*NN )
  read  ( 81, rec = 1 ) X2
  close ( 81 )
  
  ! normalisation
  mean1 = 0.0d0
  mean2 = 0.0d0
  
  do i = 1, 3
    mean1(i) = sum(X1(i, :)) / dble(NN)
    mean2(i) = sum(X2(i, :)) / dble(NN)
  end do
  
  var1 = 0.0d0
  var2 = 0.0d0
  do i = 1, 3
    var1(i) = sum( (X1(i, :)-mean1(i))*(X1(i, :)-mean1(i)) ) / dble(NN)
    var2(i) = sum( (X2(i, :)-mean2(i))*(X2(i, :)-mean2(i)) ) / dble(NN)
  end do
  
  do i = 1, 3
    X1(i, :) = ( X1(i, :) - mean1(i) ) / dsqrt(var1(i))
    X2(i, :) = ( X2(i, :) - mean2(i) ) / dsqrt(var2(i))
  end do
  
  ! dsysv preparation
  N = N_res
  LDA = N
  LDB = N
  NRHS = N_out
  LWORK = (N+1)*(N+1)
  
  !!! Echo-State Network  !!!
  do iab = 1, Na*Nb
    ia = mod(iab-1, Na) + 1
    ib = (iab-1) / Na + 1 + myrank*10
    print '(A, I3, A, I3, A)', 'i_a=', ia, ': i_b=', ib
    
    J = 1.0d0 / ( dble(ib) * db )
    J0 = J * dble(ia) * da
    
    shapek = J0*J0/J/J / dble(N_res)
    theta  = J*J / J0
    
    ! ! preparation for Gamma
    call random_number ( xx )
    intseed = int( xx * 10000 ) + 1
    state = vslnewstream ( newstr, VSL_BRNG_MCG31, intseed )
    
    !$omp parallel
    !$omp do private(ith,W_in,i_ens,i,k,M,Mtemp,sgnM,r1,r2,r1s,r2s,Z1,Z2,RR,RY,IPIV,WORK,INFO,ti,X1p,X2p)
    do i_ens = 1, N_ens
      ith = omp_get_thread_num()
      
      !!! reservoir setup !!!
      ! set M
      allocate ( M ( N_res, N_res ), Mtemp ( N_resSQ ) )
      state = vdrnggamma ( method, newstr, N_resSQ, Mtemp, shapek, aGamma, theta )
      M = reshape( Mtemp, (/ N_res, N_res /) )
      deallocate ( Mtemp )
      
      allocate ( sgnM ( N_res, N_res ) )
      call random_number ( sgnM )
      sgnM = 2.0d0 * dble( floor( 2.0d0*sgnM ) ) - 1.0d0
      M = M * sgnM
      deallocate ( sgnM )
      
      ! set W_in
      call random_number( W_in )
      W_in = 2.0d0 * (W_in-0.5d0) * w0
      
      ! initialise r
      call random_number( r1 )
      call random_number( r2 )
      allocate( r1s ( N_res, NN ), r2s ( N_res, NN ) )
      r1s = 0.0d0
      r2s = 0.0d0
      
      !!! run reservoir !!!
      do ti = 1, NN
        Z1 = matmul( M, dtanh(r1) ) + W_in*X1(1, ti)
        Z2 = matmul( M, dtanh(r2) ) + W_in*X2(1, ti)
        r1 = alpha * Z1 + (1.0d0-alpha) * r1
        r2 = alpha * Z2 + (1.0d0-alpha) * r2
        r1s(:, ti) = r1
        r2s(:, ti) = r2
      end do
      deallocate ( M )
      r1s = dtanh(r1s)
      r2s = dtanh(r2s)
       
      ! estimate W_out for 1-step prediction of X-variable of Lorenz
      allocate ( RR ( N_res, N_res ), RY ( N_res, 3 ) )
      RR = matmul( r1s(:, N_tran:NN), transpose(r1s(:, N_tran:NN)) ) / dble(NN-N_tran+1)
      RY = matmul( r1s(:, N_tran:NN), transpose(X1 (:, N_tran:NN)) ) / dble(NN-N_tran+1)
      
      ! regularisation if needed
      ! do i = 1, N_res
      !   RR(i, i) = RR(i, i) + lamb
      ! end do
      
      allocate ( IPIV ( N ), WORK ( LWORK ) )
      call dsysv ( 'U', N, NRHS, RR, LDA, IPIV, RY, LDB, WORK, LWORK, INFO )
      deallocate ( IPIV, WORK, RR )
      
      allocate ( X1p ( 3, NN ), X2p ( 3, NN ) )
      X1p = matmul( transpose(RY), r1s )
      X2p = matmul( transpose(RY), r2s )
      deallocate ( r1s, r2s, RY )
      
      ! MSE
      do i = 1, 3
        MSEtrain(i, i_ens) = sum( ( X1p(i, N_tran:NN)-X1(i, N_tran:NN) )*( X1p(i, N_tran:NN)-X1(i, N_tran:NN) ))
        MSEtest (i, i_ens) = sum( ( X2p(i, N_tran:NN)-X2(i, N_tran:NN) )*( X2p(i, N_tran:NN)-X2(i, N_tran:NN) ))
      end do
      
      deallocate ( X1p, X2p )
      
    end do
    !$omp end do
    !$omp end parallel
    
    MSEtrain = MSEtrain / (NN-N_tran+1)
    MSEtest  = MSEtest  / (NN-N_tran+1)
    
    ! output MSE
    write ( filepath, '( A, I4.4, A, I3.3, A, I3.3, A )' ) &
      '../data/MSE_OLP_Nres', N_res, '_sGammaM_tanhr_J0J_Linear/ia', ia, '_ib', ib, '.dat'
    open  ( 20, file = filepath )
    do i_ens = 1, N_ens
      write ( 20, '( 6E15.6E3 )' ) MSEtrain(1:3, i_ens), MSEtest(1:3, i_ens)
    end do
    close ( 20 )
    
  end do
  
  print ' ( A )', 'all done'
  
  
end program runreservoir

