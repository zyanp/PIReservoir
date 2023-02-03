program runreservoir
  
  implicit none
  include 'mpif.h'
  integer     (8) :: PETOT, myrank, ierr
  integer     (8), parameter :: NN = 100000
  integer     (8), parameter :: N_tran = 1000, N_pred = 1110
  integer     (8), parameter :: N_res = 500 ! N_res = 100, 500, 1000
  integer     (8), parameter :: N_ens = 100, N_in = 1, N_out = 1
  integer     (8), parameter :: lambexp = -1
  integer     (8), parameter :: Na = 40, Nb = 10
  real        (8), parameter :: pi = dacos(-1.0d0), alpha = 0.2d0, lamb = dble(10**lambexp)
  real        (8), parameter :: da = 0.05d0, db = 0.05d0
  real        (8) :: X1 ( 3, NN ), X2 ( 3, NN ), X1p ( N_pred ), X2p ( N_pred )
  real        (8) :: W_in ( N_res ), M ( N_res, N_res )
  real        (8) :: g1 ( N_res*N_res/2 ), g2 ( N_res*N_res/2 )
  real        (8) :: g3 ( N_res*N_res/2 ), g4 ( N_res*N_res/2 )
  real        (8) :: r1 ( N_res ), r2 ( N_res ), Z1 ( N_res ), Z2 ( N_res )
  real        (8), allocatable :: r1s ( :, : )
  real        (8) :: u ( N_in, 1 ), RR ( N_res, N_res ), RY ( N_res )
  real        (8) :: mu, std, J0, J, mean1 ( 3 ), var1 ( 3 ), mean2 ( 3 ), var2 ( 3 )
  real        (8) :: meanX, meanXp, varX, varXp
  real        (8) :: MSEtrain ( N_ens ), MSEtest ( N_ens ), CRCtrain ( N_ens ), CRCtest ( N_ens )
  integer     (8) :: ti, i_ens, iab, ia, ib, i, ith
  character (200) :: filepath, linebuf
  integer (4), allocatable :: seed ( : )
  integer (4) :: clock ( 8 ), seedsize
  integer     (4) :: N, LDA, LDB, NRHS, LWORK, INFO
  integer     (4), allocatable :: IPIV ( : )
  real        (8), allocatable :: WORK ( : )
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
    ib = (iab-1)/Na + 1 + myrank*10
    print '(A, I3, A, I3, A, I3, A )', 'i_a=', ia, ': i_b=', ib, ' (', myrank, ')'
    
    J = 1.0d0 / ( dble(ib)*db )
    J0 = J * dble(ia)*da
    
    mu = J0 / dble(N_res)
    std = J / dsqrt( dble(N_res) )
    
    MSEtrain = 0.0d0
    MSEtest  = 0.0d0
    CRCtrain = 0.0d0
    CRCtest  = 0.0d0
    
    !$omp parallel
    !$omp do private(ith,W_in,g1,i,g2,g3,g4,M,r1,r2,r1s,ti,Z1,Z2,RR,RY,IPIV,WORK,INFO,X1p,X2p,meanX,varX,meanXp,varXp,filepath)
    do i_ens = 1, N_ens
      ith = omp_get_thread_num()
      
      !!! reservoir setup !!!
      call random_number( W_in )
      ! call random_number( bias )
      W_in = 2.0d0 * (W_in-0.5d0) * w0
      ! bias = 2.0d0 * (bias-0.5d0) * b0
      
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
      
      ! initialise r
      call random_number( r1 )
      allocate( r1s ( N_res, NN ) )
      r1s = 0.0d0
      r1s(:, 1) = r1
      
      !!! run reservoir !!!
      do ti = 2, NN
        Z1 = matmul( M, dtanh(r1) ) + W_in*X1(1, ti-1)
        r1 = alpha * Z1 + (1.0d0-alpha) * r1
        r1s(:, ti) = r1
      end do
      
      ! estimate W_out for 1-step prediction of X-variable of Lorenz
      RR = matmul( r1s(:, N_tran:NN), transpose(r1s(:, N_tran:NN)) ) / dble(NN-N_tran+1)
      RY = matmul( r1s(:, N_tran:NN), X1 (1, N_tran:NN) )  / dble(NN-N_tran+1)
     
      ! regularisation if needed
      do i = 1, N_res
        RR(i, i) = RR(i, i) + lamb
      end do
       
      allocate ( IPIV ( N ), WORK ( LWORK ) )
      call dsysv ( 'U', N, NRHS, RR, LDA, IPIV, RY, LDB, WORK, LWORK, INFO )
      
      
      !!! closed loop prediction !!!
      ! initialise r
      call random_number( r1 )
      call random_number( r2 )
      
      X1p ( 1 ) = X1 ( 1, 1 )
      X2p ( 1 ) = X2 ( 1, 1 )
      ! run reservoir
      do ti = 2, N_pred
        Z1 = matmul( M, dtanh(r1) ) + W_in*X1p(ti-1)
        Z2 = matmul( M, dtanh(r2) ) + W_in*X2p(ti-1)
        r1 = alpha * Z1 + (1.0d0-alpha) * r1
        r2 = alpha * Z2 + (1.0d0-alpha) * r2
        
        !!! set input !!!
        ! feed correct time-series until reservoir state is fit to the time series
        if ( ti < N_tran ) then
          ! listening phase
          X1p(ti) = X1(1, ti)
          X2p(ti) = X2(1, ti)
        else
          ! predicting phase
          X1p(ti) = sum(RY*r1)
          X2p(ti) = sum(RY*r2) 
        end if
      end do
      
      
      ! MSE
      MSEtrain(i_ens) = sum( ( X1p(N_tran:N_pred)-X1(1, N_tran:N_pred) )*( X1p(N_tran:N_pred)-X1(1, N_tran:N_pred) ))
      MSEtest (i_ens) = sum( ( X2p(N_tran:N_pred)-X2(1, N_tran:N_pred) )*( X2p(N_tran:N_pred)-X2(1, N_tran:N_pred) ))
      deallocate ( IPIV, WORK, r1s ) !, r2s )
      
      ! Correlation Coeffcient
      meanX = sum( X1(1, N_tran:N_pred) ) / dble(N_pred-N_tran+1)
      varX  = sum( (X1(1, N_tran:N_pred)-meanX)*(X1(1, N_tran:N_pred)-meanX) ) / dble(N_pred-N_tran+1)
      meanXp = sum( X1p( N_tran:N_pred ) ) / dble(N_pred-N_tran+1)
      varXp  = sum( (X1p( N_tran:N_pred )-meanXp)*(X1p( N_tran:N_pred )-meanXp) ) / dble(N_pred-N_tran+1)
      CRCtrain(i_ens) = sum( (X1(1, N_tran:N_pred)-meanX)*(X1p(N_tran:N_pred)-meanXp) ) / dble(N_pred-N_tran+1) / dsqrt(varX) / dsqrt(varXp)
      
      meanX = sum( X2(1, N_tran:N_pred) ) / dble(N_pred-N_tran+1)
      varX  = sum( (X2(1, N_tran:N_pred)-meanX)*(X2(1, N_tran:N_pred)-meanX) ) / dble(N_pred-N_tran+1)
      meanXp = sum( X2p( N_tran:N_pred ) ) / dble(N_pred-N_tran+1)
      varXp  = sum( (X2p( N_tran:N_pred )-meanXp)*(X2p( N_tran:N_pred )-meanXp) ) / dble(N_pred-N_tran+1)
      CRCtest (i_ens) = sum( (X2(1, N_tran:N_pred)-meanX)*(X2p(N_tran:N_pred)-meanXp) ) / dble(N_pred-N_tran+1) / dsqrt(varX) / dsqrt(varXp)
      
    end do
    !$omp end do
    !$omp end parallel
    
    MSEtrain = MSEtrain / (N_pred-N_tran)
    MSEtest  = MSEtest  / (N_pred-N_tran)
    
    ! output MSE and CRC
    write ( filepath, '( A, I4.4, A, I3.3, A, I3.3, A )' ) &
      & '../data/Nres', N_res, '_regp1e-1_CLP_SoA_GaussianM_J0J_Linear/MSE_CRC_NormalisedLorenz_ia', ia, '_ib', ib, '_traintest.dat'
    print *, filepath
    open  ( 20, file = filepath )
    do i_ens = 1, N_ens
      write ( 20, '( 4E15.6E3 )' ) MSEtrain(i_ens), MSEtest(i_ens), CRCtrain(i_ens), CRCtest(i_ens)
    end do
    close ( 20 )
    
    print *, maxval(MSEtrain), maxval(MSEtest), maxval(CRCtrain), maxval(CRCtest)
  end do
  
  print ' ( A )', 'all done'
  
  call MPI_FINALIZE (ierr)
  
end program runreservoir

