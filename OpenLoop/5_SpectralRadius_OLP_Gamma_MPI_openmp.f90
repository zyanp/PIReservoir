program runreservoir
  use mkl_vsl_type
  use mkl_vsl
  !$ use omp_lib
  !implicit none
  include 'mpif.h'
  integer     (8) :: PETOT, myrank, ierr
  integer     (8), parameter :: N_res = 500, N_resSQ = N_res*N_res
  integer     (8), parameter :: N_ens = 100
  integer     (8), parameter :: Na = 40, Nb = 10
  real        (8), parameter :: pi = dacos(-1.0d0), alpha = 0.2d0
  real        (8), parameter :: da = 0.05d0, db = 0.05d0
  real        (8), allocatable :: M ( :, : ), Mtemp ( : )
  real        (8) :: g1 ( N_res*N_res/2 ), g2 ( N_res*N_res/2 )
  real        (8) :: g3 ( N_res*N_res/2 ), g4 ( N_res*N_res/2 )
  real        (8) :: shapek, theta, aGamma = 0.0d0, J0, J, SRs ( N_ens ), xx
  integer     (8) :: i_ens, iab, ia, ib, i, ith
  character (200) :: filepath, linebuf
  integer (4), allocatable :: seed ( : )
  integer (4) :: clock ( 8 ), seedsize
  integer     (4) :: N, LDA, LDVL, LDVR, LWORK, INFO
  real        (8), allocatable :: WR ( : ), WI ( : ), VL ( :, : ), VR ( :, : ), WORK ( : )
  integer (4) :: time1, time2, t_rate, t_max, t_diff
  integer (4) :: method = VSL_RNG_METHOD_GAMMA_GNORM
  integer (4) :: intseed
  type (vsl_stream_state) :: newstr
  
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
  
  ! vslGamma
  intseed = 3698
  
  ! dsysv preparation
  N = N_res
  LDA = N
  LDVR = N
  LDVL = N
  LWORK = 4*N
  allocate ( WR(N), WI(N), VL(LDVL, N), VR(LDVR, N), WORK(LWORK) )
  
  !!! Echo-State Network  !!!
  !$omp parallel
  !$omp do private(iab,ith,ia,ib,J,J0,newstr,method,shapek,theta,i_ens,Mtemp,M,WR,WI,VL,VR,WORK,INFO,filepath,SRs)
  do iab = 1, Na*Nb
    ith = omp_get_thread_num()
    
    ia = mod(iab-1, Na) + 1
    ib = (iab-1) / Na + 1 + myrank*10
    print '(A, I3, A, I3, A, I3, A, I3, A)', 'i_a=', ia, ': i_b=', ib, ' (', myrank, ', ', ith, ')'
    
    J = 1.0d0 / ( dble(ib) * db )
    J0 = J * dble(ia) * da
    
    shapek = J0*J0/J/J / dble(N_res)
    theta  = J*J / J0
    
    call random_number ( xx )
    intseed = int( xx * 10000 ) + 1
    state = vslnewstream ( newstr, VSL_BRNG_MCG31, intseed )
    SRs = 0.0d0
    
    allocate( M ( N_res, N_res ), Mtemp ( N_res*N_res ) )
     
    do i_ens = 1, N_ens
      state = vdrnggamma ( method, newstr, N_resSQ, Mtemp, shapek, aGamma, theta )
      M = reshape( Mtemp, (/ N_res, N_res /) )
      
      ! calc eigenvalues
      call dgeev ( 'N', 'V', N, M, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      
      SRs ( i_ens ) = dsqrt ( maxval ( WR*WR + WI*WI ) )
    end do
    
    ! output Lyapunov exponent
    write ( filepath, '( A, I4.4, A, I3.3, A, I3.3, A )' ) &
      '../data/SpectrumRadius_GammaM_Nres', N_res, '_J0JLinear/ia', ia, '_ib', ib, '.dat'
    print *, filepath
    open  ( 20+ith, file = filepath )
    do i_ens = 1, N_ens
      write ( 20+ith, '( E21.8E3 )' ) SRs(i_ens)
    end do
    close ( 20+ith )
    
    deallocate ( M, Mtemp )
  end do
  !$omp end do
  !$omp end parallel
  
  deallocate ( WR, WI, VL, VR, WORK )
  print ' ( A )', 'all done'
  
  call MPI_FINALIZE (ierr)
  
end program runreservoir

