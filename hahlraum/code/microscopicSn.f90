!--------------------------------------------------
!>macroscopic_sn
!--------------------------------------------------
module microscopic_sn
    use global_data
    use tools
    implicit none
contains

    !--------------------------------------------------
    !>calculate the slope of distribution function
    !--------------------------------------------------
    subroutine interpolation()
        integer :: i,j

        !no interpolation if first order (already set to zero slope when initializing)
        if (method_interp==FIRST_ORDER) return

        !$omp parallel
        !$omp do
        do i=1,cell_number
            call interp_inner(ctr(i))
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine interpolation

    !--------------------------------------------------
    !>calculate the flux across the interfaces
    !--------------------------------------------------
    subroutine evolution()
        integer :: i,j

        !$omp parallel
        !$omp do
        do i=1,face_number
            if(face(i)%group==face1)then
                call calc_flux(face(i))
            elseif(face(i)%group==face0)then
                call calc_flux_boundary(face(i))
            endif
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine evolution

    !--------------------------------------------------
    !>update cell averaged values
    !--------------------------------------------------
    subroutine update_flux_sn()
        integer :: i,j

        !----------------------------------------------------------
        !> calculate cross section
        !----------------------------------------------------------
        do i=1,cell_number
            if(ctr(i)%group==cell01 .or. ctr(i)%group==cell02)then
                if(ctr(i)%coords(2)>0.5)then
                    ctr(i)%absorbtion=2.0d3
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                else
                    ctr(i)%absorbtion=2.0d-1
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                endif
            else
                if(ctr(i)%group==cell11 .or. ctr(i)%group==cell13)then
                    ctr(i)%absorbtion=2.0d-1
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                elseif(ctr(i)%group==cell12 .or. ctr(i)%group==cell14)then
                    ctr(i)%absorbtion=2.0d3
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                endif
            endif
        enddo
        
        do i=1,face_number
            if(face(i)%group==face1)then
                !face(i)%kappa=2.0d0*(ctr(face(i)%cell_id(1))%kappa_effective*ctr(face(i)%cell_id(2))%kappa_effective)/&
                !    sum(ctr(face(i)%cell_id(1:2))%kappa_effective)
                face(i)%absorbtion_scaled=min(ctr(face(i)%cell_id(1))%absorbtion_scaled,ctr(face(i)%cell_id(2))%absorbtion_scaled)
            else
                face(i)%absorbtion_scaled=ctr(face(i)%cell_id(1))%absorbtion_scaled
            endif
            face(i)%flux_h=face(i)%flux_h*exp(-dt*face(i)%absorbtion_scaled)
        enddo

        !$omp parallel
        !$omp do
        do i=1,cell_number
            if(ctr(i)%group==cell11 .or. ctr(i)%group==cell12 .or.&
               ctr(i)%group==cell13 .or. ctr(i)%group==cell14)then
                call update_flux(ctr(i),i)
            endif
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine update_flux_sn

    !--------------------------------------------------
    !>update cell averaged values
    !--------------------------------------------------
    subroutine update_source_sn()
        integer :: i,j

        !$omp parallel
        !$omp do
        do i=1,cell_number
            if(ctr(i)%group==cell11 .or. ctr(i)%group==cell12 .or.&
               ctr(i)%group==cell13 .or. ctr(i)%group==cell14)then
                call update_source(ctr(i))
            endif
        end do
        !$omp end do nowait
        !$omp end parallel
    end subroutine update_source_sn

    !--------------------------------------------------
    !>interpolation of the inner cells
    !>@param[in]    cell_L :the left cell
    !>@param[inout] cell_N :the target cell
    !>@param[in]    cell_R :the right cell
    !>@param[in]    idx    :the index indicating i or j direction
    !--------------------------------------------------
    subroutine interp_inner(cell)
        type(cell_type),intent(inout) :: cell

    end subroutine interp_inner

    !--------------------------------------------------
    !>update flux
    !--------------------------------------------------
    subroutine update_flux(cell,cell_id)
        type(cell_type),intent(inout) :: cell
        type(face_type) :: face_local(3)
        integer,intent(in) :: cell_id
        real(kind=double) :: temp_rho
        integer :: direction(3)
        integer :: i,j

        do i=1,3
            face_local(i)=face(cell%face_id(i))
            if(cell_id==face_local(i)%cell_id(1))then
                direction(i)=1
            elseif(cell_id==face_local(i)%cell_id(2))then
                direction(i)=-1
            else
                write(*,*) "check update_flux..."
            endif
        enddo
        !update distribution function
        cell%h=cell%h+(direction(1)*face_local(1)%flux_h+direction(2)*face_local(2)%flux_h+direction(3)*face_local(3)%flux_h)/cell%area
        cell%rho=sum(weight*cell%h)
        do i=1,unum
            do j=1,vnum
                cell%h(i,j)=max(cell%h(i,j),0.0d0)
            enddo
        enddo
        temp_rho=sum(weight*cell%h)
        cell%h=cell%h*cell%rho/temp_rho
    end subroutine update_flux

    !--------------------------------------------------
    !>update flux
    !--------------------------------------------------
    subroutine update_source(cell)
        type(cell_type),intent(inout) :: cell
        real(kind=double),dimension(0:4):: coefficient
        complex(kind=double),dimension(0:3):: solution
        real(kind=double) :: realpart,absimaginepart
        real(kind=double) :: temp_alpha,temp_beta,temp_rho
        integer :: i

        !temp coefficients
        temp_alpha=dt*lightspeed*cell%absorbtion/kn**2.0d0
        !temp_beta=dt*cell%absorbtion/cell%capacity/kn**2.0d0

        !!update T
        !coefficient(4)=temp_beta/(1+temp_alpha)*radconst*lightspeed
        !coefficient(3)=0.0d0
        !coefficient(2)=0.0d0
        !coefficient(1)=1.0d0
        !coefficient(0)=-cell%T-temp_beta/(1+temp_alpha)*cell%rho
        !CALL QuarticRoots(coefficient,solution)
        !realpart=0.0d0
        !absimaginepart=1.0d8
        !do i=0,3
        !    if(dble(solution(i))>realpart .and. abs(aimag(solution(i)))<absimaginepart)then
        !        realpart=dble(solution(i))
        !        absimaginepart=abs(aimag(solution(i)))
        !    endif
        !enddo
        !cell%T=realpart

        !update h
        !cell%h=1.0d0/(1.0d0+temp_alpha)*(cell%h+temp_alpha*radconst*lightspeed/(4.0d0*pi)*cell%T**4.0d0)
        cell%h=exp(-temp_alpha)*abs(cell%h)+(1-exp(-temp_alpha))*cell%rho/(4.0d0*pi)
        temp_rho=sum(weight*cell%h)
        cell%h=cell%h*cell%rho/temp_rho
    end subroutine update_source

    !--------------------------------------------------
    !>calculate flux of inner interface
    !>@param[in]    cell_L :cell left to the target interface
    !>@param[inout] face   :the target interface
    !>@param[in]    cell_R :cell right to the target interface
    !>@param[in]    idx    :index indicating i or j direction
    !--------------------------------------------------
    subroutine calc_flux(face)
        type(face_type),intent(inout) :: face
        real(kind=double),allocatable,dimension(:,:) :: vn,vt !normal and tangential micro velocity
        real(kind=double),allocatable,dimension(:,:) :: h !distribution function at the interface
        real(kind=double),allocatable,dimension(:,:) :: sh !slope of distribution function at the interface
        integer,allocatable,dimension(:,:) :: delta !Heaviside step function
        integer :: i,j

        !allocate array
        allocate(vn(unum,vnum))
        allocate(vt(unum,vnum))
        allocate(delta(unum,vnum))
        allocate(h(unum,vnum))
        allocate(sh(unum,vnum))

        !convert the micro velocity to local frame,right hand
        vn = uspace*face%norm(1)+vspace*face%norm(2)
        vt =-uspace*face%norm(2)+vspace*face%norm(1)

        !Heaviside step function
        delta = (sign(UP,vn)+1)/2

        !--------------------------------------------------
        !>reconstruct initial distribution
        !>cell2->cell1
        !--------------------------------------------------
        face%h = ctr(face%cell_id(2))%h*delta+ctr(face%cell_id(1))%h*(1-delta)

        !--------------------------------------------------
        !calculate flux of distribution function
        !--------------------------------------------------
        face%flux_h = dt*face%length*vn*face%h
    end subroutine calc_flux

    !--------------------------------------------------
    !>calculate flux of inner interface
    !>@param[in]    cell_L :cell left to the target interface
    !>@param[inout] face   :the target interface
    !>@param[in]    cell_R :cell right to the target interface
    !>@param[in]    idx    :index indicating i or j direction
    !--------------------------------------------------
    subroutine calc_flux_boundary(face)
        type(face_type),intent(inout) :: face
        real(kind=double),allocatable,dimension(:,:) :: vn,vt !normal and tangential micro velocity
        real(kind=double),allocatable,dimension(:,:) :: h !distribution function at the interface
        real(kind=double),allocatable,dimension(:,:) :: sh !slope of distribution function at the interface
        integer,allocatable,dimension(:,:) :: delta !Heaviside step function
        integer :: i,j

        !allocate array
        allocate(vn(unum,vnum))
        allocate(vt(unum,vnum))
        allocate(delta(unum,vnum))
        allocate(h(unum,vnum))
        allocate(sh(unum,vnum))

        !convert the micro velocity to local frame,right hand
        vn = uspace*face%norm(1)+vspace*face%norm(2)
        vt =-uspace*face%norm(2)+vspace*face%norm(1)

        !Heaviside step function
        delta = (sign(UP,vn)+1)/2

        !--------------------------------------------------
        !>reconstruct initial distribution
        !>cell2->cell1
        !--------------------------------------------------
        do i=1,unum
            do j=1,vnum
                face%h(i,j) = ctr(face%cell_id(1))%h(i,vnum-j+1)*delta(i,j)+ctr(face%cell_id(1))%h(i,j)*(1-delta(i,j))
            enddo
        enddo

        !--------------------------------------------------
        !calculate flux of distribution function
        !--------------------------------------------------
        face%flux_h = dt*face%length*vn*face%h

    end subroutine calc_flux_boundary
end module microscopic_sn
