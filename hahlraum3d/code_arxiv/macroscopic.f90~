!--------------------------------------------------
!> macroscopic
!--------------------------------------------------
module macroscopic
    use global_data
    use MKL_SPBLAS
    implicit none
contains

    !--------------------------------------------------
    !> update cell averaged values
    !--------------------------------------------------
    subroutine source_iteration()
        real(kind=double),allocatable,dimension(:) :: rhophi
        real(kind=double),allocatable,dimension(:) :: RHS
        real(kind=double),allocatable,dimension(:) :: solution
        real(kind=double),allocatable,dimension(:) :: Tk,Tk1
        real(kind=double),allocatable,dimension(:) :: phik,phik1
        real(kind=double),allocatable,dimension(:) :: rhok,rhok1
        real(kind=double),allocatable,dimension(:) :: flux_limiter
        real(kind=double) :: loop_res,loop_eps
        integer(kind=short) :: loop,loop_max
        integer(kind=short) :: matrix_type ! 1=general 2=sparse
        integer(kind=short) :: i,j
        ! general matrix
        real(kind=double),allocatable,dimension(:,:) :: matrix_explicit
        real(kind=double),allocatable,dimension(:,:) :: matrix_implicit
        real(kind=double),allocatable,dimension(:,:) :: matrix_out
        integer(kind=short),allocatable,dimension(:) :: IPIV
        integer(kind=short) :: matrix_dim,leading_dimA,leading_dimB
        integer(kind=short) :: INFO
        ! sparse matrix variables
        type(sparse_matrix_t) :: ex_matrix
        type(matrix_descr)    :: ex_descr
        real(kind=double),allocatable,dimension(:)   :: ex_val
        real(kind=double),allocatable,dimension(:)   :: im_val
        integer(kind=short),allocatable,dimension(:) :: ex_col,ex_row,ex_ptB,ex_ptE
        integer(kind=short),allocatable,dimension(:) :: im_col,im_row,im_ptB,im_ptE
        integer(kind=short),allocatable,dimension(:) :: perm
        integer(kind=short)  :: ex_nval,ex_nrow,ex_ncol
        integer(kind=short)  :: im_nval,im_nrow,im_ncol
        integer(kind=short)  :: flag_create,flag_product
        integer(kind=double) :: pt(64)
        integer  :: maxfct,mnum,mtype,phase,nrhs,msglvl,operation,error,iparm(64)
        !debug
        real(kind=double),allocatable,dimension(:) :: solution2

        matrix_type=2
        allocate(rhophi(2*cell_number))
        allocate(RHS(2*cell_number_inner))
        allocate(solution(2*cell_number_inner))
        allocate(Tk(cell_number_inner))
        allocate(Tk1(cell_number_inner))
        allocate(phik(cell_number_inner))
        allocate(phik1(cell_number_inner))
        allocate(rhok(cell_number_inner))
        allocate(rhok1(cell_number_inner))
        allocate(flux_limiter(cell_number))
        !if(matrix_type==1)then
            ! general matrix
            allocate(matrix_explicit(2*cell_number_inner,2*cell_number))
            allocate(matrix_implicit(2*cell_number_inner,2*cell_number_inner))
            allocate(matrix_out(2*cell_number_inner,2*cell_number_inner))
            allocate(IPIV(2*cell_number_inner))
        !elseif(matrix_type==2)then
            ! allocate sparse matrix
            im_nrow=2*cell_number_inner
            im_ncol=2*cell_number_inner
            ex_nrow=2*cell_number_inner
            ex_ncol=2*cell_number
            allocate(perm(im_nrow))
            allocate(im_ptB(ex_nrow))
            allocate(im_ptE(ex_nrow))
            allocate(im_row(ex_nrow+1))
            allocate(ex_ptB(ex_nrow))
            allocate(ex_ptE(ex_nrow))
            allocate(ex_row(ex_nrow+1))
            call pre_loadmatrix(im_nrow,im_ncol,im_nval,im_ptB,im_ptE,im_row,&
                ex_nrow,ex_ncol,ex_nval,ex_ptB,ex_ptE,ex_row)
            allocate(im_val(im_nval))
            allocate(im_col(im_nval))
            allocate(ex_val(ex_nval))
            allocate(ex_col(ex_nval))
        !endif
        allocate(solution2(2*cell_number_inner))

        !set init value
        loop=0
        loop_max=200
        loop_eps=1.0d-8
        Tk=ctr(1:cell_number_inner)%T
        rhok=ctr(1:cell_number_inner)%rho
        do i=cell_number_inner+1,cell_number
            if(ctr(i)%group==cell01)then
                !sigma
                ctr(i)%absorbtion=2.0d-3
                !capacity
                ctr(i)%capacity=1.0d-1
                !kappa=c/3*sigma
                ctr(i)%kappa=lightspeed/(3.0d0*ctr(i)%absorbtion)
                !absorbtion_scaled=c*sigma/epsilon^2
                ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                !capacity_scaled=4acT^3/cv
                ctr(i)%capacity_scaled=max(4.0d0*radconst*ctr(i)%T**3.0d0/ctr(i)%capacity,1.0d-8)
                !!effective diffusive coefficient
                !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                !    2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))
                !!!UGKWP
                !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                !                2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))*&
                !                (1.0d0-exp(-dt*ctr(i)%absorbtion_scaled))
                !ctr(i)%kappa_effective=ctr(i)%kappa*flux_limiter(i)
                ctr(i)%kappa_effective=ctr(i)%kappa
            elseif(ctr(i)%group==cell02)then
                !sigma
                ctr(i)%absorbtion=2.0d-3
                !capacity
                ctr(i)%capacity=1.0d-1
                !kappa=c/3*sigma
                ctr(i)%kappa=lightspeed/(3.0d0*ctr(i)%absorbtion)
                !absorbtion_scaled=c*sigma/epsilon^2
                ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                !capacity_scaled=4acT^3/cv
                ctr(i)%capacity_scaled=max(4.0d0*radconst*ctr(i)%T**3.0d0/ctr(i)%capacity,1.0d-8)
                !!effective diffusive coefficient
                !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                !    2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))
                !!!UGKWP
                !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                !                2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))*&
                !                (1.0d0-exp(-dt*ctr(i)%absorbtion_scaled))
                !ctr(i)%kappa_effective=ctr(i)%kappa*flux_limiter(i)
                ctr(i)%kappa_effective=ctr(i)%kappa
            else
                write(*,*) "error: ghost cell group ..."
                pause
            endif
        enddo

        !source iteration
        do while(.true.)
            loop=loop+1
            !-------------------------------------
            !> diffusive coefficient  : kappa
            !> absorbtion coefficient : absorbtion_scaled
            !> capacity coefficient   : capacity_scaled
            !-------------------------------------
            do i=1,cell_number_inner
                if(ctr(i)%group==cell11 .or. ctr(i)%group==cell13)then
                    !sigma
                    ctr(i)%absorbtion=2.0d-3
                    !capacity
                    ctr(i)%capacity=1.0d-1
                    !kappa=c/3*sigma
                    ctr(i)%kappa=lightspeed/(3.0d0*ctr(i)%absorbtion)
                    !absorbtion_scaled=c*sigma/epsilon^2
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                    !capacity_scaled=4acT^3/cv
                    ctr(i)%capacity_scaled=max(4.0d0*radconst*Tk(i)**3.0d0/ctr(i)%capacity,1.0d-8)
                    !!effective diffusive coefficient
                    !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                    !    2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))
                    !!!UGKWP
                    !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                    !                2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))*&
                    !                (1.0d0-exp(-dt*ctr(i)%absorbtion_scaled))
                    !ctr(i)%kappa_effective=ctr(i)%kappa*flux_limiter(i)
                    ctr(i)%kappa_effective=ctr(i)%kappa
                elseif(ctr(i)%group==cell12 .or. ctr(i)%group==cell14)then
                    !sigma
                    ctr(i)%absorbtion=2.0d4
                    !capacity
                    ctr(i)%capacity=1.0d-1
                    !kappa=c/3*sigma
                    ctr(i)%kappa=lightspeed/(3.0d0*ctr(i)%absorbtion)
                    !absorbtion_scaled=c*sigma/epsilon^2
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                    !capacity_scaled=4acT^3/cv
                    ctr(i)%capacity_scaled=max(4.0d0*radconst*Tk(i)**3.0d0/ctr(i)%capacity,1.0d-8)
                    !!effective diffusive coefficient
                    !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                    !    2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))
                    !!!UGKWP
                    !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                    !                2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))*&
                    !                (1.0d0-exp(-dt*ctr(i)%absorbtion_scaled))
                    !ctr(i)%kappa_effective=ctr(i)%kappa*flux_limiter(i)
                    ctr(i)%kappa_effective=ctr(i)%kappa
                else
                    write(*,*) "error: source iteration cross section error..."
                    pause
                endif
            enddo
            do i=1,face_number
                if(face(i)%group==face1)then
                    !face(i)%kappa=2.0d0*(ctr(face(i)%cell_id(1))%kappa_effective*ctr(face(i)%cell_id(2))%kappa_effective)/&
                    !    sum(ctr(face(i)%cell_id(1:2))%kappa_effective)
                    face(i)%kappa=max(ctr(face(i)%cell_id(1))%kappa_effective,ctr(face(i)%cell_id(2))%kappa_effective)
                endif
            enddo

            !------------------------------------------
            !>solve linear equation
            !------------------------------------------
            !if(matrix_type==1)then !>general case
                !load matrix
                call loadmatrix_general(matrix_implicit,matrix_explicit,cell_number_inner,cell_number)
                !right hand side
                rhophi(1:cell_number)=ctr(1:cell_number)%rho
                rhophi(cell_number+1:2*cell_number)=ctr(1:cell_number)%phi
                do i=1,cell_number_inner*2
                    RHS(i)=sum(matrix_explicit(i,:)*rhophi(:))
                enddo
                !solve linear system
                NRHS=1
                matrix_dim=cell_number_inner*2
                leading_dimA=cell_number_inner*2
                leading_dimB=cell_number_inner*2
                matrix_out=matrix_implicit
                solution=RHS
                ! gauss elimination
                call dgesv(matrix_dim,NRHS,matrix_out,leading_dimA,IPIV,solution,leading_dimB,INFO)
                if(INFO/=0)then
                    write(*,*) "linear solver error: dgesv..."
                    pause
                endif
            !elseif(matrix_type==2)then !>sparse case
                !load matrix
                call loadmatrix_sparse(im_nval,im_nrow,im_ncol,im_val,im_col,ex_nval,ex_nrow,ex_ncol,ex_val,ex_col)
                !right hand side explicit*rhophi
                flag_create=mkl_sparse_d_create_csr(ex_matrix,sparse_index_base_one,ex_nrow,ex_ncol,ex_ptB,ex_ptE,ex_col,ex_val)
                if(flag_create/=sparse_status_success)then
                    write(*,*) "linear solver error: mkl_sparse_d_create_csr..."
                    pause
                endif
                rhophi(1:cell_number)=ctr(1:cell_number)%rho
                rhophi(cell_number+1:2*cell_number)=ctr(1:cell_number)%phi
                ex_descr.type = SPARSE_MATRIX_TYPE_GENERAL;
                ex_descr.mode = SPARSE_FILL_MODE_LOWER;
                ex_descr.diag = SPARSE_DIAG_NON_UNIT;
                flag_product=mkl_sparse_d_mv(sparse_operation_non_transpose,1.0d0,ex_matrix,ex_descr,rhophi,0.0d0,RHS)
                if(flag_product/=sparse_status_success)then
                    write(*,*) "linear solver error: mkl_sparse_d_mv..."
                    pause
                endif
                ! solve implicit*solution=RHS
                pt=0
                maxfct=1
                mnum=1
                mtype=11
                phase=13
                perm=0
                nrhs=1
                iparm(1)=0
                msglvl=0
                call pardiso(pt,maxfct,mnum,mtype,phase,im_nrow,im_val,im_row,im_col,perm,nrhs,iparm,msglvl,RHS,solution2,error)
                if(error/=0)then
                    write(*,*) "linear solver error: pardiso..."
                    pause
                endif
            !endif
            if(sum((solution-solution2)**2.0d0) > 1.0d-10)then
                write(*,*) "error: linear solver..."
                pause
            endif
            pause

            !-------------------------------------
            !>check convergence
            !-------------------------------------
            rhok1=solution(1:cell_number_inner)
            phik1=solution(cell_number_inner+1:2*cell_number_inner)
            Tk1=(abs(phik1/radconst/lightspeed))**0.25d0
            do i=1,cell_number_inner
                if(isnan(Tk1(i)))then
                    write(*,*) "Tk1 is nan",solution(cell_number_inner+i)
                endif
            enddo
            loop_res=sqrt(sum((rhok1-rhok)**2.0d0)+sum((Tk1-Tk)**2.0d0))/sqrt(sum((rhok)**2.0d0)+sum((Tk)**2.0d0))
            if(loop_res<loop_eps)then
                ctr(1:cell_number_inner)%T=Tk1
                ctr(1:cell_number_inner)%rho=rhok1
                ctr(1:cell_number_inner)%phi=radconst*lightspeed*ctr(1:cell_number_inner)%T**4.0d0
                exit
            elseif(loop>loop_max)then
                ctr(1:cell_number_inner)%T=Tk1
                ctr(1:cell_number_inner)%rho=rhok1
                ctr(1:cell_number_inner)%phi=radconst*lightspeed*ctr(1:cell_number_inner)%T**4.0d0
                write(*,*) "Source iteration not converge, res=", loop_res
                exit
            else
                Tk=Tk1
                rhok=rhok1
                !write(*,*) "source_iteration_res=",loop_res
            endif
        enddo
    end subroutine source_iteration

    !--------------------------------------------------------
    !> nonzeros of sparse matrix
    !> implicit and explicit: nval, pointB, pointE, rowIndex
    !--------------------------------------------------------
    subroutine pre_loadmatrix(im_nrow,im_ncol,im_nval,im_ptB,im_ptE,im_row,&
                              ex_nrow,ex_ncol,ex_nval,ex_ptB,ex_ptE,ex_row)
        integer(kind=short),intent(in)  :: im_nrow,im_ncol
        integer(kind=short),intent(in)  :: ex_nrow,ex_ncol
        integer(kind=short),intent(out) :: im_nval
        integer(kind=short),intent(out) :: ex_nval
        integer(kind=short),intent(out) :: im_ptB(im_nrow)
        integer(kind=short),intent(out) :: im_ptE(im_nrow)
        integer(kind=short),intent(out) :: im_row(im_nrow+1)
        integer(kind=short),intent(out) :: ex_ptB(ex_nrow)
        integer(kind=short),intent(out) :: ex_ptE(ex_nrow)
        integer(kind=short),intent(out) :: ex_row(ex_nrow+1)
        integer(kind=short) :: face_id,node_id,cell_id
        integer(kind=short) :: flag_im(im_ncol)
        integer(kind=short) :: flag_ex(ex_ncol)
        integer(kind=short) :: i,j,k

        im_nval=0
        ex_nval=0
        im_ptB =1
        im_ptE =1
        im_row =1
        ex_ptB =1
        ex_ptE =1
        ex_row =1

        !> rho: loop row
        do i=1,cell_number_inner
            flag_im=0
            flag_ex=0
            !> ------------------------------------------------
            !> time derivative
            !> ------------------------------------------------
            !> matrix_implicit(i,i)=matrix_implicit(i,i)+1.0d0
            if(flag_im(i)==0)then
                im_nval=im_nval+1
                flag_im(i)=im_nval
                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
            endif
            !> matrix_explicit(i,i)=matrix_explicit(i,i)+1.0d0
            if(flag_ex(i)==0)then
                ex_nval=ex_nval+1
                flag_ex(i)=ex_nval
                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
            endif
            !> ------------------------------------------------
            !> source term
            !> ------------------------------------------------
            !> matrix_implicit(i,i)=matrix_implicit(i,i)+&
            !>     dt*ctr(i)%absorbtion_scaled
            if(flag_im(i)==0)then
                im_nval=im_nval+1
                flag_im(i)=im_nval
                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
            endif
            !> matrix_implicit(i,i+cell_number_inner)=matrix_implicit(i,i+cell_number_inner)-&
            !>     dt*ctr(i)%absorbtion_scaled
            cell_id=i+cell_number_inner
            if(flag_im(cell_id)==0)then
                im_nval=im_nval+1
                flag_im(cell_id)=im_nval
                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
            endif
            !> ------------------------------------------------
            !> diffusive flux
            !> ------------------------------------------------
            do j=1,ctr(i)%face_number
                face_id=ctr(i)%face_id(j)
                if(face(face_id)%cell_id(1)==i)then ! inner norm
                    !nomal derivative
                    cell_id=face(face_id)%cell_id(2)
                    if(cell_id/=0)then
                        !matrix_implicit(i,i)=matrix_implicit(i,i)+&
                        !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        if(flag_im(i)==0)then
                            im_nval=im_nval+1
                            flag_im(i)=im_nval
                            im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                        endif
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    endif
                    !tangent derivative
                    !> node1-node2
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    enddo
                    node_id=face(face_id)%node_id(2)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    enddo
                    !tangent derivative
                    !> node1-node3
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    enddo
                    node_id=face(face_id)%node_id(3)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    enddo
                else !> outer norm
                    !nomal derivative
                    cell_id=face(face_id)%cell_id(1)
                    if(cell_id/=0)then
                        !matrix_implicit(i,i)=matrix_implicit(i,i)+&
                        !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        if(flag_im(i)==0)then
                            im_nval=im_nval+1
                            flag_im(i)=im_nval
                            im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                        endif
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    endif
                    !tangent derivative
                    !> node1-node2
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    enddo
                    node_id=face(face_id)%node_id(2)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    enddo
                    !tangent derivative
                    !> node1-node3
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    enddo
                    node_id=face(face_id)%node_id(3)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval=im_nval+1
                                flag_im(cell_id)=im_nval
                                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
                            endif
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval=ex_nval+1
                                flag_ex(cell_id)=ex_nval
                                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
                            endif
                        endif
                    enddo
                endif
            enddo
        enddo
        !> phi
        do i=cell_number_inner+1,cell_number_inner*2
            flag_im=0
            flag_ex=0
            !> ------------------------------------------------
            !> time derivative
            !> ------------------------------------------------
            !matrix_implicit(i,i)=matrix_implicit(i,i)+1.0d0
            if(flag_im(i)==0)then
                im_nval=im_nval+1
                flag_im(i)=im_nval
                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
            endif
            !matrix_explicit(i,i+(cell_number-cell_number_inner))=matrix_explicit(i,i+(cell_number-cell_number_inner))+1.0d0
            cell_id=i+(cell_number-cell_number_inner)
            if(flag_ex(cell_id)==0)then
                ex_nval=ex_nval+1
                flag_ex(cell_id)=ex_nval
                ex_ptE(i:im_nrow)=ex_ptE(i:im_nrow)+1
            endif
            !> ------------------------------------------------
            !> source term
            !> ------------------------------------------------
            !matrix_implicit(i,i)=matrix_implicit(i,i)+&
            !    dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
            if(flag_im(i)==0)then
                im_nval=im_nval+1
                flag_im(i)=im_nval
                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
            endif
            !matrix_implicit(i,i-cell_number_inner)=matrix_implicit(i,i-cell_number_inner)-&
            !    dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
            cell_id=i-cell_number_inner
            if(flag_im(cell_id)==0)then
                im_nval=im_nval+1
                flag_im(cell_id)=im_nval
                im_ptE(i:im_nrow)=im_ptE(i:im_nrow)+1
            endif
        enddo
        !> calculate pointB and rowIndex
        im_ptB(1)=1
        im_ptB(2:im_nrow)=im_ptE(1:im_nrow-1)
        im_row(1)=1
        im_row(2:im_nrow+1)=im_ptE(1:im_nrow)
        ex_ptB(1)=1
        ex_ptB(2:ex_nrow)=ex_ptE(1:ex_nrow-1)
        ex_row(1)=1
        ex_row(2:ex_nrow+1)=ex_ptE(1:ex_nrow)
    end subroutine pre_loadmatrix

    !--------------------------------------------------
    !> load matrix
    !> implicit and explicit entries
    !--------------------------------------------------
    subroutine loadmatrix_sparse(im_nval,im_nrow,im_ncol,im_val,im_col,& 
            ex_nval,ex_nrow,ex_ncol,ex_val,ex_col)
        integer(kind=short),intent(in)  :: im_nval,im_nrow,im_ncol
        integer(kind=short),intent(in)  :: ex_nval,ex_nrow,ex_ncol
        real(kind=double),intent(out)   :: ex_val(ex_nval)
        integer(kind=short),intent(out) :: ex_col(ex_nval)
        real(kind=double),intent(out)   :: im_val(im_nval)
        integer(kind=short),intent(out) :: im_col(im_nval)
        integer(kind=short) :: flag_im(im_ncol)
        integer(kind=short) :: flag_ex(ex_ncol)
        integer(kind=short) :: face_id,node_id,cell_id
        integer(kind=short) :: im_nval_temp,ex_nval_temp
        integer(kind=short) :: i,j,k

        !> initialization
        ex_val=0
        ex_col=0
        im_val=0
        im_col=0
        im_nval_temp=0
        ex_nval_temp=0

        !> rho
        do i=1,cell_number_inner
            flag_im=0
            flag_ex=0
            !> ------------------------------------------------
            !> time derivative
            !> ------------------------------------------------
            !> matrix_implicit(i,i)=matrix_implicit(i,i)+1.0d0
            if(flag_im(i)==0)then
                im_nval_temp=im_nval_temp+1
                flag_im(i)=im_nval_temp
            endif
            im_val(flag_im(i))=im_val(flag_im(i))+1.0d0
            im_col(flag_im(i))=i
            !> matrix_explicit(i,i)=matrix_explicit(i,i)+1.0d0
            if(flag_ex(i)==0)then
                ex_nval_temp=ex_nval_temp+1
                flag_ex(i)=ex_nval_temp
            endif
            ex_val(flag_ex(i))=ex_val(flag_ex(i))+1.0d0
            ex_col(flag_ex(i))=i
            !> ------------------------------------------------
            !> source term
            !> ------------------------------------------------
            !> matrix_implicit(i,i)=matrix_implicit(i,i)+&
            !>     dt*ctr(i)%absorbtion_scaled
            if(flag_im(i)==0)then
                im_nval_temp=im_nval_temp+1
                flag_im(i)=im_nval_temp
            endif
            im_val(flag_im(i))=im_val(flag_im(i))+&
                dt*ctr(i)%absorbtion_scaled
            im_col(flag_im(i))=i
            !> matrix_implicit(i,i+cell_number_inner)=matrix_implicit(i,i+cell_number_inner)-&
            !>     dt*ctr(i)%absorbtion_scaled
            cell_id=i+cell_number_inner
            if(flag_im(cell_id)==0)then
                im_nval_temp=im_nval_temp+1
                flag_im(cell_id)=im_nval_temp
            endif
            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))-&
                dt*ctr(i)%absorbtion_scaled
            im_col(flag_im(cell_id))=cell_id
            !> ------------------------------------------------
            !> diffusive flux
            !> ------------------------------------------------
            do j=1,ctr(i)%face_number
                face_id=ctr(i)%face_id(j)
                if(face(face_id)%cell_id(1)==i)then ! inner norm
                    !nomal derivative
                    cell_id=face(face_id)%cell_id(2)
                    if(cell_id/=0)then
                        !matrix_implicit(i,i)=matrix_implicit(i,i)+&
                        !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        if(flag_im(i)==0)then
                            im_nval_temp=im_nval_temp+1
                            flag_im(i)=im_nval_temp
                        endif
                        im_val(flag_im(i))=im_val(flag_im(i))+&
                            dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        im_col(flag_im(i))=i
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    endif
                    !tangent derivative
                    !> node1-node2
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    enddo
                    node_id=face(face_id)%node_id(2)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    enddo
                    !tangent derivative
                    !> node1-node3
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    enddo
                    node_id=face(face_id)%node_id(3)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    enddo
                else !> outer norm
                    !nomal derivative
                    cell_id=face(face_id)%cell_id(1)
                    if(cell_id/=0)then
                        !matrix_implicit(i,i)=matrix_implicit(i,i)+&
                        !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        if(flag_im(i)==0)then
                            im_nval_temp=im_nval_temp+1
                            flag_im(i)=im_nval_temp
                        endif
                        im_val(flag_im(i))=im_val(flag_im(i))+&
                            dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        im_col(flag_im(i))=i
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    endif
                    !tangent derivative
                    !> node1-node2
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    enddo
                    node_id=face(face_id)%node_id(2)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    enddo
                    !tangent derivative
                    !> node1-node3
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    enddo
                    node_id=face(face_id)%node_id(3)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                            ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            !matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_im(cell_id)==0)then
                                im_nval_temp=im_nval_temp+1
                                flag_im(cell_id)=im_nval_temp
                            endif
                            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            im_col(flag_im(cell_id))=cell_id
                        else
                            !matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                            !    dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            if(flag_ex(cell_id)==0)then
                                ex_nval_temp=ex_nval_temp+1
                                flag_ex(cell_id)=ex_nval_temp
                            endif
                            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                            ex_col(flag_ex(cell_id))=cell_id
                        endif
                    enddo
                endif
            enddo
        enddo
        !> phi
        do i=cell_number_inner+1,cell_number_inner*2
            flag_im=0
            flag_ex=0
            !> ------------------------------------------------
            !> time derivative
            !> ------------------------------------------------
            !matrix_implicit(i,i)=matrix_implicit(i,i)+1.0d0
            if(flag_im(i)==0)then
                im_nval_temp=im_nval_temp+1
                flag_im(i)=im_nval_temp
            endif
            im_val(flag_im(i))=im_val(flag_im(i))+1.0d0
            im_col(flag_im(i))=i
            !matrix_explicit(i,i+(cell_number-cell_number_inner))=matrix_explicit(i,i+(cell_number-cell_number_inner))+1.0d0
            cell_id=i+(cell_number-cell_number_inner)
            if(flag_ex(cell_id)==0)then
                ex_nval_temp=ex_nval_temp+1
                flag_ex(cell_id)=ex_nval_temp
            endif
            ex_val(flag_ex(cell_id))=ex_val(flag_ex(cell_id))+1.0d0
            ex_col(flag_ex(cell_id))=cell_id
            !> ------------------------------------------------
            !> source term
            !> ------------------------------------------------
            !matrix_implicit(i,i)=matrix_implicit(i,i)+&
            !    dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
            if(flag_im(i)==0)then
                im_nval_temp=im_nval_temp+1
                flag_im(i)=im_nval_temp
            endif
            im_val(flag_im(i))=im_val(flag_im(i))+&
                dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
            im_col(flag_im(i))=i
            !matrix_implicit(i,i-cell_number_inner)=matrix_implicit(i,i-cell_number_inner)-&
            !    dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
            cell_id=i-cell_number_inner
            if(flag_im(cell_id)==0)then
                im_nval_temp=im_nval_temp+1
                flag_im(cell_id)=im_nval_temp
            endif
            im_val(flag_im(cell_id))=im_val(flag_im(cell_id))-&
                dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
            im_col(flag_im(cell_id))=cell_id
        enddo
    end subroutine loadmatrix_sparse

    !--------------------------------------------------
    !> load matrix
    !> implicit and explicit entries
    !--------------------------------------------------
    subroutine loadmatrix_general(matrix_implicit,matrix_explicit,cell_number_inner,cell_number)
        integer(kind=short) :: cell_number_inner,cell_number
        real(kind=double),intent(out) :: matrix_implicit(2*cell_number_inner,2*cell_number_inner)
        real(kind=double),intent(out) :: matrix_explicit(2*cell_number_inner,2*cell_number)
        integer(kind=short) :: face_id,node_id,cell_id
        integer(kind=short) :: i,j,k

        matrix_implicit=0
        matrix_explicit=0

        !> time derivative
        do i=1,cell_number_inner
            matrix_implicit(i,i)=matrix_implicit(i,i)+1.0d0
            matrix_explicit(i,i)=matrix_explicit(i,i)+1.0d0
        enddo
        do i=cell_number_inner+1,2*cell_number_inner
            matrix_implicit(i,i)=matrix_implicit(i,i)+1.0d0
            matrix_explicit(i,i+(cell_number-cell_number_inner))=matrix_explicit(i,i+(cell_number-cell_number_inner))+1.0d0
        enddo

        !> source term
        do i=1,cell_number_inner
            matrix_implicit(i,i)=matrix_implicit(i,i)+&
                dt*ctr(i)%absorbtion_scaled
            matrix_implicit(i,i+cell_number_inner)=matrix_implicit(i,i+cell_number_inner)-&
                dt*ctr(i)%absorbtion_scaled
        enddo
        do i=cell_number_inner+1,cell_number_inner*2
            matrix_implicit(i,i)=matrix_implicit(i,i)+&
                dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
            matrix_implicit(i,i-cell_number_inner)=matrix_implicit(i,i-cell_number_inner)-&
                dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
        enddo

        !> diffusive flux
        do i=1,cell_number_inner
            do j=1,ctr(i)%face_number
                face_id=ctr(i)%face_id(j)
                if(face(face_id)%cell_id(1)==i)then ! inner norm
                    !nomal derivative
                    cell_id=face(face_id)%cell_id(2)
                    if(cell_id/=0)then
                        matrix_implicit(i,i)=matrix_implicit(i,i)+&
                            dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        endif
                    endif
                    !tangent derivative
                    !> node1-node2
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                    node_id=face(face_id)%node_id(2)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                    !tangent derivative
                    !> node1-node3
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                    node_id=face(face_id)%node_id(3)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                else !> outer norm
                    !nomal derivative
                    cell_id=face(face_id)%cell_id(1)
                    if(cell_id/=0)then
                        matrix_implicit(i,i)=matrix_implicit(i,i)+&
                            dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(1)*face(face_id)%kappa
                        endif
                    endif
                    !tangent derivative
                    !> node1-node2
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                    node_id=face(face_id)%node_id(2)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                    !tangent derivative
                    !> node1-node3
                    node_id=face(face_id)%node_id(1)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                    node_id=face(face_id)%node_id(3)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                                dt/ctr(i)%volume*face(face_id)%area*face(face_id)%weight(3)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                endif
            enddo
        enddo
    end subroutine loadmatrix_general

    !-------------------------------------
    !> Gauss-Seidel iteration
    !-------------------------------------
    subroutine GaussSeidel(A,b,x,Adim)
        real(kind=double),parameter :: gs_eps=1.0d-30
        integer(kind=short),parameter :: max_iter=1000
        integer(kind=short),intent(in) :: Adim
        real(kind=double),intent(in) :: A(Adim,Adim),b(Adim)
        real(kind=double),intent(out) :: x(Adim)
        real(kind=double) :: x_new(Adim),error(Adim)
        integer :: i,j,iter

        !initial variables
        x=0
        !iteration
        do iter=1,max_iter
            !calculate x_new
            x_new(1)=1.0d0/A(1,1)*(b(1)-sum(A(1,2:Adim)*x(2:Adim)))
            do i=2,Adim
                x_new(i)=1.0d0/A(i,i)*(b(i)-sum(A(i,1:i-1)*x_new(1:i-1))-sum(A(i,i+1:Adim)*x(i+1:Adim)))
            enddo
            !check convergence
            do i=1,Adim
                error(i)=b(i)-sum(A(i,:)*x)
            enddo
            if(sqrt(sum(error**2.0d0))<gs_eps)then
                x=x_new
                exit
            else
                x=x_new
            endif
        enddo
    end subroutine GaussSeidel
end module macroscopic
