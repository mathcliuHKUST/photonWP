!--------------------------------------------------
!> macroscopic
!--------------------------------------------------
module macroscopic
    use global_data
    use lapack95
    use f95_precision
    implicit none
contains
    
    !--------------------------------------------------
    !>update cell averaged values
    !--------------------------------------------------
    subroutine source_iteration()
        real(kind=double),allocatable,dimension(:,:) :: matrix_explicit
        real(kind=double),allocatable,dimension(:,:) :: matrix_implicit
        real(kind=double),allocatable,dimension(:,:) :: matrix_out
        real(kind=double),allocatable,dimension(:) :: RHS
        real(kind=double),allocatable,dimension(:) :: solution
        real(kind=double),allocatable,dimension(:) :: Tk,Tk1
        real(kind=double),allocatable,dimension(:) :: phik,phik1
        real(kind=double),allocatable,dimension(:) :: rhok,rhok1
        real(kind=double),allocatable,dimension(:) :: flux_limiter
        real(kind=double) :: source_iteration_res,source_iteration_eps
        integer(kind=short),allocatable,dimension(:) :: IPIV
        integer(kind=short) :: loop,source_iteration_max
        integer(kind=short) :: matrix_dim,NRHS,leading_dimA,leading_dimB
        integer(kind=short) :: INFO
        integer(kind=short) :: i,j

        allocate(matrix_explicit(2*cell_number_inner,2*cell_number))
        allocate(matrix_implicit(2*cell_number_inner,2*cell_number_inner))
        allocate(matrix_out(2*cell_number_inner,2*cell_number_inner))
        allocate(solution(2*cell_number_inner))
        allocate(IPIV(2*cell_number_inner))
        allocate(RHS(2*cell_number))
        allocate(Tk(cell_number_inner))
        allocate(Tk1(cell_number_inner))
        allocate(phik(cell_number_inner))
        allocate(phik1(cell_number_inner))
        allocate(rhok(cell_number_inner))
        allocate(rhok1(cell_number_inner))
        allocate(flux_limiter(cell_number))

        !set init value
        loop=0
        source_iteration_max=200
        source_iteration_eps=1.0d-8
        Tk=ctr(1:cell_number_inner)%T
        rhok=ctr(1:cell_number_inner)%rho
        do i=cell_number_inner+1,cell_number
            if(ctr(i)%coords(2)>0.5)then
                !sigma
                ctr(i)%absorbtion=2.0d4
                !capacity
                ctr(i)%capacity=1.0d-1
                !kappa=ac/3*sigma
                ctr(i)%kappa=radconst*lightspeed/(3.0d0*ctr(i)%absorbtion)
                !absorbtion_scaled=c*sigma/epsilon^2
                ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                !capacity_scaled=4acT^3/cv
                ctr(i)%capacity_scaled=max(4.0d0*radconst*ctr(i)%T**3.0d0/ctr(i)%capacity,1.0d-8)
                !effective diffusive coefficient
                flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                    2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))
                !!!UGKWP
                !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                !                2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))*&
                !                (1.0d0-exp(-dt*ctr(i)%absorbtion_scaled))
                !ctr(i)%kappa_effective=ctr(i)%kappa*flux_limiter(i)
                ctr(i)%kappa_effective=ctr(i)%kappa
            else
                !sigma
                ctr(i)%absorbtion=2.0d-3
                !capacity
                ctr(i)%capacity=1.0d-1
                !kappa=ac/3*sigma
                ctr(i)%kappa=radconst*lightspeed/(3.0d0*ctr(i)%absorbtion)
                !absorbtion_scaled=c*sigma/epsilon^2
                ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                !capacity_scaled=4acT^3/cv
                ctr(i)%capacity_scaled=max(4.0d0*radconst*ctr(i)%T**3.0d0/ctr(i)%capacity,1.0d-8)
                !effective diffusive coefficient
                flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                    2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))
                !!!UGKWP
                !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                !                2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))*&
                !                (1.0d0-exp(-dt*ctr(i)%absorbtion_scaled))
                !ctr(i)%kappa_effective=ctr(i)%kappa*flux_limiter(i)
                ctr(i)%kappa_effective=ctr(i)%kappa
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
                    !kappa=ac/3*sigma
                    ctr(i)%kappa=radconst*lightspeed/(3.0d0*ctr(i)%absorbtion)
                    !absorbtion_scaled=c*sigma/epsilon^2
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                    !capacity_scaled=4acT^3/cv
                    ctr(i)%capacity_scaled=1 !max(4.0d0*radconst*Tk(i)**3.0d0/ctr(i)%capacity,1.0d-8)
                    !effective diffusive coefficient
                    flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                        2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))
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
                    !kappa=ac/3*sigma
                    ctr(i)%kappa=radconst*lightspeed/(3.0d0*ctr(i)%absorbtion)
                    !absorbtion_scaled=c*sigma/epsilon^2
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                    !capacity_scaled=4acT^3/cv
                    ctr(i)%capacity_scaled=1 !max(4.0d0*radconst*Tk(i)**3.0d0/ctr(i)%capacity,1.0d-8)
                    !effective diffusive coefficient
                    flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                        2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))
                    !!!UGKWP
                    !flux_limiter(i)=(1.0d0-2.0d0/(dt*ctr(i)%absorbtion_scaled)+exp(-dt*ctr(i)%absorbtion_scaled)+&
                    !                2.0d0*exp(-dt*ctr(i)%absorbtion_scaled)/(dt*ctr(i)%absorbtion_scaled))*&
                    !                (1.0d0-exp(-dt*ctr(i)%absorbtion_scaled))
                    !ctr(i)%kappa_effective=ctr(i)%kappa*flux_limiter(i)
                    ctr(i)%kappa_effective=ctr(i)%kappa
                else
                    write(*,*) "source iteration cross section error..."
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

            !-------------------------------------
            !>solve linear equation
            !-------------------------------------
            !load matrix
            !write(*,*) "loading matrix..."
            call loadmatrix(matrix_implicit,matrix_explicit,cell_number_inner,cell_number)
            !right hand side
            RHS(1:cell_number)=ctr(1:cell_number)%rho
            RHS(cell_number+1:2*cell_number)=ctr(1:cell_number)%phi
            do i=1,cell_number_inner*2
                solution(i)=sum(matrix_explicit(i,:)*RHS(:))
                if(isnan(solution(i)))then
                    write(*,*) "solution is nan"
                endif
                do j=1,cell_number_inner*2
                    if(isnan(matrix_implicit(i,j)))then
                        write(*,*) i,j,matrix_implicit(i,j),dt
                        write(*,*) "matrix is nan"
                        pause
                    endif
                enddo
            enddo
            !solve linear system
            NRHS=1
            matrix_dim=cell_number_inner*2
            leading_dimA=cell_number_inner*2
            leading_dimB=cell_number_inner*2
            matrix_out=matrix_implicit
            !write(*,*) "solving matrix..."
            ! gauss elimination
            call dgesv(matrix_dim,NRHS,matrix_out,leading_dimA,IPIV,solution,leading_dimB,INFO)
            if(INFO/=0)then
                write(*,*) "linear solver dgesv error..."
                pause
            endif
            ! gauss-seidel iteration
            !call GaussSeidel(matrix_implicit,RHS,solution,matrix_dim)
            do i=1,cell_number_inner*2
                if(isnan(solution(i)))then
                    write(*,*) "solution is nan"
                endif
            enddo

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
            source_iteration_res=sqrt(sum((rhok1-rhok)**2.0d0)+sum((Tk1-Tk)**2.0d0))/sqrt(sum((rhok)**2.0d0)+sum((Tk)**2.0d0))
            if(source_iteration_res<source_iteration_eps)then
                ctr(1:cell_number_inner)%T=Tk1
                ctr(1:cell_number_inner)%rho=rhok1
                ctr(1:cell_number_inner)%phi=radconst*lightspeed*ctr(1:cell_number_inner)%T**4.0d0
                exit
            elseif(loop>source_iteration_max)then
                ctr(1:cell_number_inner)%T=Tk1
                ctr(1:cell_number_inner)%rho=rhok1
                ctr(1:cell_number_inner)%phi=radconst*lightspeed*ctr(1:cell_number_inner)%T**4.0d0
                write(*,*) "Source iteration not converge, res=", source_iteration_res
                exit
            else
                Tk=Tk1
                rhok=rhok1
                !write(*,*) "source_iteration_res=",source_iteration_res
            endif
        enddo
    end subroutine source_iteration

    !--------------------------------------------------
    !> load matrix
    !> implicit and explicit entries
    !--------------------------------------------------
    subroutine loadmatrix(matrix_implicit,matrix_explicit,cell_number_inner,cell_number)
        integer(kind=short) :: cell_number_inner,cell_number
        real(kind=double),intent(out) :: matrix_implicit(2*cell_number_inner,2*cell_number_inner)
        real(kind=double),intent(out) :: matrix_explicit(2*cell_number_inner,2*cell_number)
        integer(kind=short) :: face_id,node_id,cell_id
        integer(kind=short) :: i,j,k

        matrix_implicit=0
        matrix_explicit=0

        !time derivative
        do i=1,cell_number_inner
            matrix_implicit(i,i)=matrix_implicit(i,i)+1.0d0
            matrix_explicit(i,i)=matrix_explicit(i,i)+1.0d0
        enddo
        do i=cell_number_inner+1,2*cell_number_inner
            matrix_implicit(i,i)=matrix_implicit(i,i)+1.0d0
            matrix_explicit(i,i+(cell_number-cell_number_inner))=matrix_explicit(i,i+(cell_number-cell_number_inner))+1.0d0
        enddo

        !source term
        do i=1,cell_number_inner
            matrix_implicit(i,i)=matrix_implicit(i,i)+dt*ctr(i)%absorbtion_scaled
            matrix_implicit(i,i+cell_number_inner)=matrix_implicit(i,i+cell_number_inner)-&
                dt*ctr(i)%absorbtion_scaled
        enddo
        do i=cell_number_inner+1,cell_number_inner*2
            matrix_implicit(i,i)=matrix_implicit(i,i)+dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
            matrix_implicit(i,i-cell_number_inner)=matrix_implicit(i,i-cell_number_inner)-&
                dt*ctr(i-cell_number_inner)%absorbtion_scaled*ctr(i-cell_number_inner)%capacity_scaled
        enddo

        do i=1,cell_number_inner
            !diffusive flux
            do j=1,ctr(i)%face_number
                face_id=ctr(i)%face_id(j)
                if(face(face_id)%cell_id(1)==i)then !inner normal direction
                    !nomal derivative
                    cell_id=face(face_id)%cell_id(2)
                    if(cell_id/=0)then
                        matrix_implicit(i,i)=matrix_implicit(i,i)+&
                            dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(1)*face(face_id)%kappa
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(1)*face(face_id)%kappa
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(1)*face(face_id)%kappa
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
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                    node_id=face(face_id)%node_id(2)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                else !> outer normal
                    !nomal derivative
                    cell_id=face(face_id)%cell_id(1)
                    if(cell_id/=0)then
                        matrix_implicit(i,i)=matrix_implicit(i,i)+&
                            dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(1)*face(face_id)%kappa
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)-&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(1)*face(face_id)%kappa
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(1)*face(face_id)%kappa
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
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)+&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                    node_id=face(face_id)%node_id(2)
                    do k=1,node(node_id)%cell_number
                        cell_id=node(node_id)%cell_id(k)
                        if(ctr(cell_id)%group==cell11 .or. ctr(cell_id)%group==cell12 .or.&
                           ctr(cell_id)%group==cell13 .or. ctr(cell_id)%group==cell14)then
                            matrix_implicit(i,cell_id)=matrix_implicit(i,cell_id)+&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        else
                            matrix_explicit(i,cell_id)=matrix_explicit(i,cell_id)-&
                                dt/ctr(i)%area*face(face_id)%length*face(face_id)%weight(2)*face(face_id)%kappa/node(node_id)%cell_number
                        endif
                    enddo
                endif
            enddo
        enddo
    end subroutine loadmatrix

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
