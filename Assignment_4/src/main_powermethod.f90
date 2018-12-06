! Subroutine to generate PRN using Intel MKL 
  include 'Random.f90'
  

  ! Main program
  program matrix_eigenvalue
    use globals
    use auxiliary
    implicit none
    logical :: ok
    integer :: ierr, STATUS(MPI_STATUS_SIZE), seed_slaves
    double precision :: wtime, eps, error, shift
    double precision, allocatable, dimension(:,:) :: A,C
    double precision, allocatable, dimension(:,:) :: DUMMY
    double precision, allocatable, dimension(:) :: G, v, v_next, v_error

    ! Initialize parameters (n = Matrix Size (nxn); Ni = Desired number of matrices/max eigenvalues to be collected by the Master)
    n=1000
    Ni=10000

    ! File 'Frequency_distribution' contains a list of maximum eigenvalues
    open(14,file='Frequency_distribution')

    call start_MPI(ok)
    
    ! Flags to follow account for the jobs done and left.
    flag=1 ! Slaves work when flag==1
    job_count=0 ! Count of the eigenvalues received by master
    seed_slaves=0 ! Slave Identifier        

    if(ok) then
       wtime=mpi_wtime()
			
       !Slaves
       if(proc_num.ne.0) then
          allocate(A(n,n))
          allocate(v(n))
          allocate(v_next(n))
          allocate(v_error(n))
          
          ! Initialize power iteration parameters.	
          shift=10d0               ! Shifting parameter used in the power iteration
          eps=1d-5                 ! Tolerance
          v=1d0                    ! Eigenvector for first iteration
          
          do while (flag.eq.1)
             
             call open_random_stream

             ! Make random matrix (NOTE: Here all the entries will be N(0,1))
             ierr=vdrnggaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,n*n,A,0d0,1d0) ! Draw (n+1) numbers from standard normal distribution.

             ! Generate a symmetric matrix using the PRN generator
             do i =1,n
                do j=i,n
                   A(j,i) = A(i,j)
                   if (i.eq.j) then
                      A(i,j) = 2*A(i,j) + shift
                   end if
                end do
             end do
             
             ! Implementation of shifted power iterations using LAPACK subroutine dgemv for matrix-vector multiplication
             error = eps+1d0 ! Initial error
             
             do while (error.gt.eps)
                call dgemv('n',n,n, 1d0,A,n,v,1,0d0,v_next,1)
                v_next = v_next/dnorm2(v_next)
                v_error = v_next - v
                error=dnorm2(v_error)
                v = v_next
             end do
             
             call dgemv('n',n,n, 1d0,A,n,v,1,0d0,v_next,1) ! Calculate shift Eigenvalue from Eigenvector
            
             EV = (ddot_product(v_next,v)/dnorm2(v)) - shift ! Maximum Eigenvalue obtained by shifting back
             
             call mpi_send(EV,1,MPI_DOUBLE_PRECISION,proc_num.eq.0,DEFAULT_TAG,MPI_COMM_WORLD,ierr) ! Send Eigenvalue to Master
             
             seed_slaves=seed_slaves+1
             call mpi_recv(flag,1,MPI_DOUBLE_PRECISION,0,DEFAULT_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr) ! Receive instruction from Master
          end do
          deallocate (A, v, v_next, v_error)
       end if

       ! Master
       if(proc_num.eq.0) then
          allocate(G(Ni))
          ith = 1 ! Pointer
          
          ! Save received eigenvalues, and send flag if desired number of eigenvalues obtained
          do while (ith.le.Ni) 
             call mpi_recv(G(ith),1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,DEFAULT_TAG,MPI_COMM_WORLD,STATUS,ierr) ! Receive eigenvalues from slaves
             
             job_count = job_count+1 ! Update number of receives
             
             if (job_count.le.Ni-num_procs+1) then
                flag=1
                call mpi_send(flag,1,MPI_DOUBLE_PRECISION,STATUS(3),DEFAULT_TAG,MPI_COMM_WORLD,ierr) ! Send continue flag to slave
             else 
                flag=0
                call mpi_send(flag,1,MPI_DOUBLE_PRECISION,STATUS(3),DEFAULT_TAG,MPI_COMM_WORLD,ierr)
             end if
             ith = ith+1
          end do
          
          ! Write maximum eigenvalues to a file
          do i=1,Ni
             write (14,*) G(i)
          end do
          
       end if
    end if

    if(proc_num.eq.0) then
       wtime=mpi_wtime()-wtime
       print *, 'proc_num=', proc_num, 'wtime=',wtime
    end if
    
    call stop_MPI
    

  contains  

    !functions used at PM iteration
    function DDOT_PRODUCT(v1,v2) result(dot)
      double precision :: v1(:),v2(:)
      double precision :: dot
      integer :: n, i
      n=size(v1)
      dot =0d0
      do i = 1,n
         dot = dot + v1(i)*v2(i)
      end do
    end function DDOT_PRODUCT

    function DNORM2(v) result(r)
      double precision :: v(:)
      double precision :: r
      integer :: n,i
      n = size(v)
      r =0d0
      do i=1,n
         r = r+v(i)**2
      end do
      r = sqrt(r)
    end function DNORM2

  end program matrix_eigenvalue

  
  subroutine start_MPI(ok)
    use globals
    implicit none
    ! MPI initialization
    logical :: my_ok=.true.,ok
    integer :: ierr

    call mpi_init(ierr)

    if(ierr.ne.0) then
       print *,'MPI_init failed!'
       my_ok=.false.
    end if

    if(my_ok) then
       call mpi_comm_size(MPI_COMM_WORLD,num_procs,ierr)
       if(ierr.ne.0) then
          print *,'MPI_comm_size failed!'
          my_ok=.false.
       end if
       
       if(my_ok) then
          call mpi_comm_rank(MPI_COMM_WORLD,proc_num,ierr)
          if(ierr.ne.0) then
             print *,'MPI_comm_rank failed!'
             my_ok=.false.
          end if
       end if
    end if
       
    ! Check if everyone is ok.
    call mpi_allreduce(my_ok,ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
  end subroutine start_MPI


  subroutine stop_MPI
    use globals
    implicit none
    logical :: init
    integer :: ierr
    ! Wait until everybody is done. One process "finalize"ing while others are still working can cause ugly crashes.
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    ! Check if MPI has been initialized.
    call mpi_initialized(init,ierr)
    ! If it is, call finalize.
    if(init) call mpi_finalize(ierr)
  end subroutine stop_MPI
  
