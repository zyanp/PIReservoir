program runreservoir
  !$use omp_lib
  ! use mkl_vsl_type
  ! use mkl_vsl
  ! implicit none
  integer     (8), parameter :: NN = 100000
  integer     (8), parameter :: N_tran = 1000, N_lyap = 10000
  integer     (8), parameter :: N_res = 500, N_in = 1, N_out = 1
  integer     (8), parameter :: N_ens = 100
  integer     (8), parameter :: lambexp = -2
  integer     (8), parameter :: Na = 40, Nb = 10, myrank = 0
  real        (8), parameter :: pi = dacos(-1.0d0), alpha = 0.2d0, lamb = dble(10**lambexp)
  real        (8), parameter :: da = 0.05d0, db = 0.05d0
  real        (8) :: X1 ( NN ), W_in ( N_res ), tanhr ( N_res ), X1p
  real        (8), allocatable :: rs ( :, : ), M ( :, : ), Mtemp ( : ), RR ( :, : ), RY ( : )
  real        (8), allocatable :: g1 ( : ), g2 ( : ), g3 ( : ), g4 ( : ), JM ( :, : ), X1temp ( :, : )
  real        (8) :: II ( N_res, N_res )
  real        (8) :: r ( N_res ), z ( N_res ), u ( N_res )
  real        (8) :: mu, std, J0, J, unorm, LEs ( N_ens ), mean1, var1 !, mean2 ( 3 ), var2 ( 3 )
  integer     (8) :: ti, i_ens, iab, ia, ib, i, k, ith
  character (200) :: filepath, linebuf
  integer     (4), allocatable :: seed ( : )
  integer     (4) :: clock ( 8 ), seedsize
  integer     (4) :: N, LDA, LDB, NRHS, LWORK, INFO
  integer     (4), allocatable :: IPIV ( : )
  real        (8), allocatable :: WORK ( : )
  
  
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
  allocate ( X1temp ( 3, NN ) )
  write ( filepath, '( A )' ) '../../NormalisedLorenz/data/Target/data1_xyz_lorenz1e-4_every1e2.bin'
  open  ( 81, file = filepath, form = 'unformatted', access = 'direct', recl = 8*3*NN )
  read  ( 81, rec = 1 ) X1temp
  close ( 81 )
  X1 = X1temp ( 1, : )
  deallocate ( X1temp )
  
  ! normalisation
  mean1 = 0.0d0
  mean1 = sum ( X1 ) / dble( NN )
  var1 = 0.0d0
  var1 = sum ( (X1-mean1)*(X1-mean1) ) / dble(NN)
  X1 = ( X1 - mean1 ) / dsqrt( var1 )
  
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
    ib = (iab-1) / Na + 1 + myrank*10
    print '(A, I3, A, I3, A)', 'i_a=', ia, ': i_b=', ib
    
    J = 1.0d0 / ( dble(ib) * db )
    J0 = J * dble(ia) * da
    
    mu = J0 / dble(N_res)
    std = J / dsqrt( dble(N_res) )
    
    LEs = 0.0d0
    
    !$omp parallel
    !$omp do private(ith,W_in,i_ens,i,k,M,g1,g2,g3,g4,r,rs,Z,RR,RY,IPIV,WORK,INFO,u,ti,JM,unorm,tanhr,X1p)
    do i_ens = 1, N_ens
      ith = omp_get_thread_num()
      
      !!! reservoir setup !!!
      call random_number( W_in )
      W_in = 2.0d0 * (W_in-0.5d0)
      
      ! set M
      allocate ( M ( N_res, N_res ) , g1 ( N_res*N_res/2 ), g2 ( N_res*N_res/2 ), &
                                    & g3 ( N_res*N_res/2 ), g4 ( N_res*N_res/2 ) )
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
      deallocate ( g1, g2, g3, g4 )
      
      ! initialise r
      call random_number( r )
      allocate( rs ( N_res, NN ) )
      rs = 0.0d0
      rs(:, 1) = r
      
      !!! run reservoir !!!
      do ti = 2, NN
        Z = matmul( M, dtanh(r) ) + W_in*X1(ti-1)
        r = alpha * Z + (1.0d0-alpha) * r
        rs(:, ti) = r
      end do
      rs = dtanh( rs )
      
      ! estimate W_out for 1-step prediction of X-variable of Lorenz
      allocate ( RR ( N_res, N_res ), RY ( N_res ) )
      RR = matmul( rs(:, N_tran:NN), transpose(rs(:, N_tran:NN)) ) / dble(NN-N_tran+1)
      RY = matmul( rs(:, N_tran:NN), X1 (N_tran:NN) )  / dble(NN-N_tran+1)
      deallocate ( rs )
      
      ! regularisation if needed
      do i = 1, N_res
        RR(i, i) = RR(i, i) + lamb
      end do
      
      allocate ( IPIV ( N ), WORK ( LWORK ) )
      call dsysv ( 'U', N, NRHS, RR, LDA, IPIV, RY, LDB, WORK, LWORK, INFO )
      deallocate ( IPIV, WORK, RR )
      
      !!! closed loop !!!
      
      !!! run reservoir !!!
      ! initialise r
      call random_number( r )
      
      allocate ( JM ( N_res, N_res ) ) ! ,X1p ( N_tran+N_lyap) )
      X1p = X1 ( 1 )
      ! run reservoir
      do ti = 2, N_tran + N_lyap
        tanhr = dtanh ( r )
        Z = matmul( M, tanhr ) + W_in*X1p
        r = alpha * Z + (1.0d0-alpha) * r
        
        !!! set input !!!
        ! feed correct time-series until reservoir state is fit to the time series
        if ( ti < N_tran ) then
          ! listening phase
          X1p = X1(ti)
        else
          if ( ti == N_tran ) then
            call random_number( u )
          end if
          
          ! predicting phase
          X1p = sum(RY*dtanh(r))
          
          ! Jacobi matrix
          do i = 1, N_res
            do k = 1, N_res
              JM ( i, k ) = alpha * ( M(i,k) + W_in(i)*RY(k) ) * (1.0d0-tanhr(k)*tanhr(k)) + (1.0d0-alpha) * II(i,k)
            end do
          end do
          
          ! norm of u
          unorm = dsqrt( sum(u*u) )
          
          ! calculation for Lyapunov exponent
          LEs( i_ens ) = LEs ( i_ens ) + dlog( unorm )
          
          ! update u
          u = matmul( JM, u/unorm )
        end if
      end do
      
      deallocate ( M, RY, JM ) ! X1p
      
    end do
    !$omp end do
    !$omp end parallel
    
    LEs = LEs / dble(N_lyap)
    
    ! output Lyapunov exponent
    write ( filepath, '( A, I4.4, A, I3.3, A, I3.3, A )' ) &
      '../data/LyapunovExp_period1e4_CLP_Nres', N_res, '_regpexp-2_GaussianM_tanhr_J0J_Linear/ia', ia, '_ib', ib, '.dat'
    open  ( 20, file = filepath )
    do i_ens = 1, N_ens
      write ( 20, '( E21.8E3 )' ) LEs(i_ens)
    end do
    close ( 20 )
    
  end do
  
  print ' ( A )', 'all done'
  
  
end program runreservoir

