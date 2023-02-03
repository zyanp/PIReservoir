program runreservoir
  
  implicit none
  include 'mpif.h'
  integer     (8) :: PETOT, myrank, ierr
  integer     (8), parameter :: N_res = 500
  integer     (8), parameter :: N_ens = 100
  integer     (8), parameter :: Na = 40, Nb = 10
  real        (8), parameter :: pi = dacos(-1.0d0), alpha = 0.2d0
  real        (8), parameter :: da = 0.05d0, db = 0.05d0
  real        (8) :: M ( N_res, N_res )
  real        (8) :: g1 ( N_res*N_res/2 ), g2 ( N_res*N_res/2 )
  real        (8) :: g3 ( N_res*N_res/2 ), g4 ( N_res*N_res/2 )
  real        (8) :: mu, std, J0, J, SRs ( N_ens )
  integer     (8) :: i_ens, iab, ia, ib, i, ith
  character (200) :: filepath, linebuf
  integer (4), allocatable :: seed ( : )
  integer (4) :: clock ( 8 ), seedsize
  integer     (4) :: N, LDA, LDVL, LDVR, LWORK, INFO
  real        (8), allocatable :: WR ( : ), WI ( : ), VL ( :, : ), VR ( :, : ), WORK ( : )
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
  
  ! dsysv preparation
  N = N_res
  LDA = N
  LDVR = N
  LDVL = N
  LWORK = 4*N
  allocate ( WR(N), WI(N), VL(LDVL, N), VR(LDVR, N), WORK(LWORK) )
  
  !!! Echo-State Network  !!!
  !$omp parallel
  !$omp do private(iab,ith,ia,ib,J,J0,mu,std,i_ens,g1,i,g2,g3,g4,M,WR,WI,VL,VR,WORK,INFO,filepath,SRs)
  do iab = 1, Na*Nb
    ith = omp_get_thread_num()
    ia = mod(iab-1, Na) + 1
    ib = (iab-1) / Na + 1 + myrank*10
    print '(A, I3, A, I3, A, I3, A, I3, A)', 'i_a=', ia, ': i_b=', ib, ' (', myrank, ', ', ith, ')'
    
    J = 1.0d0 / ( dble(ib) * db )
    J0 = J * dble(ia) * da
    
    mu = J0 / dble(N_res)
    std = J / dsqrt( dble(N_res) )
   
    SRs = 0.0d0
     
    do i_ens = 1, N_ens
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
      
      ! calc eigenvalues
      call dgeev ( 'N', 'V', N, M, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      
      SRs ( i_ens ) = dsqrt ( maxval ( WR*WR + WI*WI ) )
    end do
    
    ! output Lyapunov exponent
    write ( filepath, '( A, I4.4, A, I3.3, A, I3.3, A )' ) &
      '../data/SpectrumRadius_GaussM_Nres', N_res, '_J0JLinear/ia', ia, '_ib', ib, '.dat'
    open  ( 20+ith, file = filepath )
    do i_ens = 1, N_ens
      write ( 20+ith, '( E21.8E3 )' ) SRs(i_ens)
    end do
    close ( 20+ith )
    
  end do
  !$omp end do
  !$omp end parallel
  
  deallocate ( WR, WI, VL, VR, WORK )
  print ' ( A )', 'all done'
  
  call MPI_FINALIZE (ierr)
  
end program runreservoir

