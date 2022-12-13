	program linear_solver
	use MKL_SPBLAS
	implicit none
	integer,parameter :: short =4
	integer,parameter :: double=8
	integer :: n,msize
	real(kind=double),allocatable,dimension(:) :: a_val,x,b,y
	integer(kind=short),allocatable,dimension(:) :: a_col,a_row,a_ptB,a_ptE !csr-format
	type(SPARSE_MATRIX_T) :: a_matrix
	type(MATRIX_DESCR) :: descrA
	integer :: create_matrix
	integer(kind=double) :: pt(64)
	integer :: maxfct,mnum,mtype,phase,nrhs,iparm(64),msglvl,mv_product,operation,error
	integer,allocatable,dimension(:) :: perm
	integer :: i

	msize=4
	allocate(x(msize))
	allocate(b(msize))
	allocate(y(msize))
	n=16
	allocate(a_val(n))
	allocate(a_col(n))
	allocate(a_row(msize+1))
	allocate(a_ptB(msize))
	allocate(a_ptE(msize))
	allocate(perm(msize))

	!load matrix
	a_val(1 )=1
	a_val(2 )=2
	a_val(3 )=3
	a_val(4 )=0
	a_val(5 )=5
	a_val(6 )=6
	a_val(7 )=0
	a_val(8 )=1
	a_val(9 )=9
	a_val(10)=10
	a_val(11)=1
	a_val(12)=0
	a_val(13)=13
	a_val(14)=14
	a_val(15)=1
	a_val(16)=1
	do i=1,4
		b(i)=i+16
	enddo

	do i=1,13,4
		a_col(i+0)=1+0
		a_col(i+1)=1+1
		a_col(i+2)=1+2
		a_col(i+3)=1+3
	enddo

	a_row(1)=1
	a_row(2)=5
	a_row(3)=9
	a_row(4)=13
	a_row(5)=17

	a_ptB(1)=1
	a_ptB(2)=5
	a_ptB(3)=9
	a_ptB(4)=13

	a_ptE(1)=5
	a_ptE(2)=9
	a_ptE(3)=13
	a_ptE(4)=17

	! create matrix
	create_matrix=mkl_sparse_d_create_csr(a_matrix, sparse_index_base_one, msize, msize, a_ptB, a_ptE, a_col, a_val)
	if(create_matrix==sparse_status_success)then
		! pardiso solver
		pt=0
		maxfct=1
		mnum=1
		mtype=11
		phase=13
		perm=0
		nrhs=1
		iparm(1)=0
		msglvl=0
		call pardiso(pt, maxfct, mnum, mtype, phase, msize, a_val, a_row, a_col, perm, nrhs, iparm, msglvl, b, x, error)
		if(error==0)then
			write(*,*) "successfully solved..."
		endif
	endif
	! matrix-vector product
	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
	descrA.mode = SPARSE_FILL_MODE_UPPER;
	descrA.diag = SPARSE_DIAG_NON_UNIT;
	mv_product=mkl_sparse_d_mv(sparse_operation_non_transpose,1.0d0,a_matrix,descrA,x,0.0d0,y)

	! check solution
	write(*,*) y
	write(*,*) b
	write(*,*) create_matrix
	pause
	end program linear_solver



	!program mkl
	!    implicit none
	!    integer, parameter        :: n = 4
	!    real(kind=8)              ::  x(n), y(n),A(n)
	!    integer(kind=4)           ::i(5),j(10)
	!    real(kind=8)           ::csr(10)
	!    character(len=1)          :: uplo = 'u'
	!    integer(kind=8)                 :: ipt(64)
	!    integer(kind=8)                 :: idum(n)
	!    integer(kind=4)                 :: maxfct, mnum, mtype, phase, nrhs, error, msglvl
	!    integer(kind=4)                 :: iparm(64)
	!    x=(/1,1,1,1/)
	!    i=(/1,5,8,10,11/)
	!    j=(/1,2,3,4,2,3,4,3,4,4/)
	!    csr=(/6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25/)
	!    call mkl_dcsrsymv( uplo, n, csr, i, j, x, y )
	!    ipt    = 0
	!    maxfct = 1
	!    mnum   = 1
	!    mtype  = -2
	!    phase  = 13
	!    nrhs   = 1
	!    msglvl = 0
	!    error  = 0
	!    iparm  = 0
	!    idum=0
	!    call pardiso(ipt, maxfct, mnum, mtype, phase, n, csr, i, j, idum, nrhs, iparm, msglvl,y, A, error)
	!end program mkl
