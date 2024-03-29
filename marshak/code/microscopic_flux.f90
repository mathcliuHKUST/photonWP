!--------------------------------------------------
!>microscopic
!--------------------------------------------------
module microscopic
    use global_data
    implicit none
contains
    !--------------------------------------------------
    !>update particle
    !--------------------------------------------------
    subroutine update_particle()
        real(kind=double),allocatable,dimension(:) :: rand_num
        integer :: i,j,k

        !----------------------------------------------------------
        !> calculate cross section
        !----------------------------------------------------------
        do i=1,cell_number
            if(ctr(i)%group==cell01 .or. ctr(i)%group==cell02)then
                if(ctr(i)%coords(2)>0.5)then
                    ctr(i)%absorbtion=min(30.0d0/ctr(i)%T**3.0d0,1.0d8)
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                else
                    ctr(i)%absorbtion=min(30.0d0/ctr(i)%T**3.0d0,1.0d8)
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                endif
            else
                if(ctr(i)%group==cell11 .or. ctr(i)%group==cell13)then
                    ctr(i)%absorbtion=min(30.0d0/ctr(i)%T**3.0d0,1.0d8)
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                elseif(ctr(i)%group==cell12 .or. ctr(i)%group==cell14)then
                    ctr(i)%absorbtion=min(30.0d0/ctr(i)%T**3.0d0,1.0d8)
                    ctr(i)%absorbtion_scaled=lightspeed*ctr(i)%absorbtion/Kn**2.0d0
                endif
            endif
        enddo

        !----------------------------------------------------------
        !> stream particles
        !----------------------------------------------------------
        allocate(rand_num(particle_number))
        do i=1,face_number
            face(i)%flux_particle=0
        enddo
        do while(.true.)
            ! random number
            do i=1,particle_number
                if(particle(i)%flag_track==1)then
                    call random_number(rand_num(i))
                endif
            enddo
            !track particle
            !$omp parallel
            !$omp do
            do i=1,particle_number
                if(particle(i)%flag_track==1)then
                    call stream(i,rand_num(i))
                endif
            enddo
            !$omp end do nowait
            !$omp end parallel
            if(sum(particle(1:particle_number)%flag_track)==0) exit
        enddo

        !----------------------------------------------------------
        !> microscopic flux
        !----------------------------------------------------------
        do i=1,face_number
            if(face(i)%group==face1)then
                !face(i)%kappa=2.0d0*(ctr(face(i)%cell_id(1))%kappa_effective*ctr(face(i)%cell_id(2))%kappa_effective)/&
                !    sum(ctr(face(i)%cell_id(1:2))%kappa_effective)
                face(i)%absorbtion_scaled=max(ctr(face(i)%cell_id(1))%absorbtion_scaled,ctr(face(i)%cell_id(2))%absorbtion_scaled)
            else
                face(i)%absorbtion_scaled=ctr(face(i)%cell_id(1))%absorbtion_scaled
            endif
        enddo
        do i=1,face_number
            ctr(face(i)%cell_id(1))%rho=ctr(face(i)%cell_id(1))%rho-face(i)%flux_particle*exp(-dt*face(i)%absorbtion_scaled)/ctr(face(i)%cell_id(1))%area
            if(face(i)%cell_id(2) /=0)then
                ctr(face(i)%cell_id(2))%rho=ctr(face(i)%cell_id(2))%rho+face(i)%flux_particle*exp(-dt*face(i)%absorbtion_scaled)/ctr(face(i)%cell_id(2))%area
            endif
            !ctr(face(i)%cell_id(1))%rho=ctr(face(i)%cell_id(1))%rho-face(i)%flux_particle/ctr(face(i)%cell_id(1))%area
            !if(face(i)%cell_id(2) /=0)then
            !    ctr(face(i)%cell_id(2))%rho=ctr(face(i)%cell_id(2))%rho+face(i)%flux_particle/ctr(face(i)%cell_id(2))%area
            !endif
        enddo
    end subroutine update_particle

    !--------------------------------------------------
    !>stream particle
    !--------------------------------------------------
    subroutine stream(particle_id,rand_num)
        real(kind=double),intent(in) :: rand_num
        real(kind=double) :: velocity(2)
        real(kind=double) :: vector(2)
        real(kind=double) :: point(2)
        real(kind=double) :: coefficient(3,2)
        real(kind=double) :: crosspoint(2,3)
        real(kind=double) :: vnorm
        real(kind=double) :: tb(3)
        integer,intent(in) :: particle_id
        integer :: cell_id,face_id,node_id(2)
        integer :: i,j,k

        !--------------------------------------------------
        ! sample collision time
        !--------------------------------------------------
        cell_id=particle(particle_id)%cell
        if(ctr(cell_id)%group==cell01 .or. ctr(cell_id)%group==cell02)then
            particle(particle_id)%ta=maxnumber
        else
            particle(particle_id)%ta=maxnumber
            !particle(particle_id)%ta=-log(rand_num)/ctr(cell_id)%absorbtion_scaled
        endif

        !--------------------------------------------------
        !> calculate cross points
        !--------------------------------------------------
        ! particle trajectory
        velocity(1)=particle(particle_id)%v(1)
        velocity(2)=particle(particle_id)%v(2)
        point(1)=particle(particle_id)%x(1)
        point(2)=particle(particle_id)%x(2)
        coefficient(1,2)=velocity(2)
        coefficient(2,2)=-velocity(1)
        coefficient(3,2)=velocity(1)*point(2)-velocity(2)*point(1)

        ! boundary face
        do i=1,ctr(cell_id)%face_number
            face_id=ctr(cell_id)%face_id(i)
            if(face_id==particle(particle_id)%face_in)then
                tb(i)=maxnumber
            else
                node_id(1:2)=face(face_id)%node_id(1:2)
                vector(1)=node(node_id(2))%coords(1)-node(node_id(1))%coords(1)
                vector(2)=node(node_id(2))%coords(2)-node(node_id(1))%coords(2)
                point(1)=node(node_id(1))%coords(1)
                point(2)=node(node_id(1))%coords(2)
                coefficient(1,1)=vector(2)
                coefficient(2,1)=-vector(1)
                coefficient(3,1)=vector(1)*point(2)-vector(2)*point(1)
                ! cross point i
                crosspoint(1,i)=(coefficient(2,1)*coefficient(3,2)-coefficient(2,2)*coefficient(3,1))/&
                    (coefficient(1,1)*coefficient(2,2)-coefficient(1,2)*coefficient(2,1)+SMV)
                crosspoint(2,i)=(coefficient(3,1)*coefficient(1,2)-coefficient(3,2)*coefficient(1,1))/&
                    (coefficient(1,1)*coefficient(2,2)-coefficient(1,2)*coefficient(2,1)+SMV)
                tb(i)=((crosspoint(1,i)-particle(particle_id)%x(1))*velocity(1)+(crosspoint(2,i)-particle(particle_id)%x(2))*velocity(2))/sum(velocity(1:2)**2.0d0)
            endif
        enddo

        ! cross point and time
        particle(particle_id)%tb=maxnumber
        ! pick 1.cross point 2.facein 3.flag
        do i=1,ctr(cell_id)%face_number
            if(tb(i)>0 .and. particle(particle_id)%tb>tb(i)) then
                particle(particle_id)%tb=tb(i)
                particle(particle_id)%x_new=crosspoint(:,i)
                if(cell_id==face(ctr(cell_id)%face_id(i))%cell_id(1))then
                    particle(particle_id)%cell_new=face(ctr(cell_id)%face_id(i))%cell_id(2)
                    if(particle(particle_id)%cell_new/=0)then
                        do k=1,ctr(particle(particle_id)%cell_new)%face_number
                            if(face(ctr(particle(particle_id)%cell_new)%face_id(k))%cell_id(1)==cell_id)then
                                particle(particle_id)%face_in_new=ctr(particle(particle_id)%cell_new)%face_id(k)
                                if(ctr(cell_id)%face_id(i) /= particle(particle_id)%face_in_new)then
                                    write(*,*) "face id error..."
                                endif
                            endif
                        enddo
                        particle(particle_id)%face_in_new=ctr(cell_id)%face_id(i)
                    elseif(ctr(particle(particle_id)%cell)%group==cell13 .or. ctr(particle(particle_id)%cell)%group==cell14)then
                        particle(particle_id)%face_in_new=ctr(cell_id)%face_id(i)
                    else
                        particle(particle_id)%flag_delete=1
                        particle(particle_id)%flag_track=0
                        particle(particle_id)%face_in_new=ctr(cell_id)%face_id(i)
                    endif
                else
                    particle(particle_id)%cell_new=face(ctr(cell_id)%face_id(i))%cell_id(1)
                    if(particle(particle_id)%cell_new/=0)then
                        do k=1,ctr(particle(particle_id)%cell_new)%face_number
                            if(face(ctr(particle(particle_id)%cell_new)%face_id(k))%cell_id(2)==cell_id)then
                                particle(particle_id)%face_in_new=ctr(particle(particle_id)%cell_new)%face_id(k)
                                if(ctr(cell_id)%face_id(i) /= particle(particle_id)%face_in_new)then
                                    write(*,*) "face id error..."
                                endif
                            endif
                        enddo
                        particle(particle_id)%face_in_new=ctr(cell_id)%face_id(i)
                    elseif(ctr(particle(particle_id)%cell)%group==cell13 .or. ctr(particle(particle_id)%cell)%group==cell14)then
                        particle(particle_id)%face_in_new=ctr(cell_id)%face_id(i)
                    else
                        particle(particle_id)%flag_delete=1
                        particle(particle_id)%flag_track=0
                        particle(particle_id)%face_in_new=ctr(cell_id)%face_id(i)
                    endif
                endif
            endif
        enddo

        ! update particle position & free stream flux
        if(particle(particle_id)%tb<particle(particle_id)%ta)then
            ! free stream
            if(particle(particle_id)%t+particle(particle_id)%tb<dt)then
                particle(particle_id)%t=particle(particle_id)%t+particle(particle_id)%tb
                !ctr(particle(particle_id)%cell)%rho=ctr(particle(particle_id)%cell)%rho-particle(particle_id)%weight/ctr(particle(particle_id)%cell)%area
                if(particle(particle_id)%cell_new/=0)then
                    particle(particle_id)%x=particle(particle_id)%x_new
                    particle(particle_id)%cell=particle(particle_id)%cell_new
                    particle(particle_id)%face_in=particle(particle_id)%face_in_new
                    if(particle(particle_id)%cell==face(particle(particle_id)%face_in)%cell_id(1))then
                        face(particle(particle_id)%face_in)%flux_particle=face(particle(particle_id)%face_in)%flux_particle-particle(particle_id)%weight
                    else
                        face(particle(particle_id)%face_in)%flux_particle=face(particle(particle_id)%face_in)%flux_particle+particle(particle_id)%weight
                    endif
                    !ctr(particle(particle_id)%cell)%rho=ctr(particle(particle_id)%cell)%rho+particle(particle_id)%weight/ctr(particle(particle_id)%cell)%area
                elseif(ctr(particle(particle_id)%cell)%group==cell13 .or. ctr(particle(particle_id)%cell)%group==cell14)then
                    particle(particle_id)%v(2)=-particle(particle_id)%v(2)
                    particle(particle_id)%x=particle(particle_id)%x_new
                    particle(particle_id)%cell=particle(particle_id)%cell
                    particle(particle_id)%face_in=particle(particle_id)%face_in_new
                    !ctr(particle(particle_id)%cell)%rho=ctr(particle(particle_id)%cell)%rho+particle(particle_id)%weight/ctr(particle(particle_id)%cell)%area
                else
                    particle(particle_id)%flag_track=0
                    particle(particle_id)%flag_delete=1
                    particle(particle_id)%face_in=particle(particle_id)%face_in_new
                    face(particle(particle_id)%face_in)%flux_particle=face(particle(particle_id)%face_in)%flux_particle+particle(particle_id)%weight
                endif
            else
                particle(particle_id)%x=particle(particle_id)%x+particle(particle_id)%v*(dt-particle(particle_id)%t)
                particle(particle_id)%t=dt
                particle(particle_id)%flag_track=0
            endif
        else
            ! absorbed
            particle(particle_id)%flag_track=0
            particle(particle_id)%flag_delete=1
        endif
    endsubroutine stream

    !--------------------------------------------------
    !>resample particle
    !--------------------------------------------------
    subroutine reinit_particle()
        real(kind=double),allocatable,dimension(:,:) :: rand_num
        real(kind=double),allocatable,dimension(:) :: rand_ta
        integer :: particle_number_survive
        integer :: cell_id
        integer :: i,j,k

        !----------------------------------------------------------------------------
        !> absorb photon
        !----------------------------------------------------------------------------
        allocate(rand_ta(particle_number))
        do i=1,particle_number
            ! sample collision time
            cell_id=particle(i)%cell
            if(ctr(cell_id)%group==cell01 .or. ctr(cell_id)%group==cell02)then
                particle(i)%flag_delete=1
            else
                call random_number(rand_ta(i))
                particle(i)%ta=-log(rand_ta(i))/ctr(cell_id)%absorbtion_scaled
                if(particle(i)%ta<dt)then
                    particle(i)%flag_delete=1
                endif
            endif
        enddo

        ! particle number
        do i=1,cell_number
            ctr(i)%particle_number=0
            ctr(i)%rho_particle=0
        enddo
        do i=1,particle_number
            if(particle(i)%flag_delete==0 .and.&
              (ctr(particle(i)%cell)%group==cell11 .or.&
               ctr(particle(i)%cell)%group==cell12 .or.&
               ctr(particle(i)%cell)%group==cell13 .or.&
               ctr(particle(i)%cell)%group==cell14))then
                cell_id=particle(i)%cell
                ctr(cell_id)%particle_number=ctr(cell_id)%particle_number+1
                ctr(cell_id)%rho_particle=ctr(cell_id)%rho_particle+particle(i)%weight/ctr(cell_id)%area
            endif
        enddo

        !----------------------------------------------------------------------------
        !> calculate rho_wave
        !----------------------------------------------------------------------------
        do i=1,cell_number
            if(ctr(i)%group==cell01 .and. ctr(i)%coords(2)>0.5)then
                ctr(i)%T          = 1.0d0
                ctr(i)%phi        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho_wave   = ctr(i)%rho
                ctr(i)%WaveParticleRatio=0.0d0
            elseif(ctr(i)%group==cell01 .and. ctr(i)%coords(2)<0.5)then
                ctr(i)%T          = 1.0d0
                ctr(i)%phi        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho_wave   = ctr(i)%rho
                ctr(i)%WaveParticleRatio=0.0d0
            elseif(ctr(i)%group==cell02)then
                ctr(i)%T          = 1.0d-6
                ctr(i)%phi        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho_wave   = ctr(i)%rho
                ctr(i)%WaveParticleRatio=0.0d0
            else
                !!UGKP
                !method 1
                ctr(i)%rho_wave=max(ctr(i)%rho-ctr(i)%rho_particle,0.0d0)
                ctr(i)%WaveParticleRatio=ctr(i)%rho_particle/(ctr(i)%rho+1.0d-10)
                !!method 2
                !ctr(i)%rho_wave=ctr(i)%rho*(1-exp(-dt*ctr(i)%absorbtion_scaled))

                !!!UGKWP
                !!method 1
                !!!ctr(i)%rho_wave=(ctr(i)%rho*(1-exp(-dt*ctr(i)%absorbtion_scaled)))*exp(-dt*ctr(i)%absorbtion_scaled)
                !!method 1
                !!!ctr(i)%rho_wave=ctr(i)%rho*(1-exp(-dt*ctr(i)%absorbtion_scaled))*exp(-dt*ctr(i)%absorbtion_scaled)
            endif
        enddo

        do i=1,cell_number
            ctr(i)%particle_number_wave=max(int(ctr(i)%rho_wave*ctr(i)%area/ref_weight),0)
            ctr(i)%particle_number=ctr(i)%particle_number+ctr(i)%particle_number_wave
        enddo
        particle_number_temp=sum(ctr(1:cell_number)%particle_number)
        allocate(particle_temp(particle_number_temp))

        ! particle->particle_temp->particle
        particle_number_survive=0
        do i=1,particle_number
            if(particle(i)%flag_delete==0 .and.&
              (ctr(particle(i)%cell)%group==cell11 .or.&
               ctr(particle(i)%cell)%group==cell12 .or.&
               ctr(particle(i)%cell)%group==cell13 .or.&
               ctr(particle(i)%cell)%group==cell14))then
                particle_number_survive=particle_number_survive+1
                particle_temp(particle_number_survive)=particle(i)
                particle_temp(particle_number_survive)%t=0
                particle_temp(particle_number_survive)%flag_track=1
            endif
        enddo

        deallocate(particle)
        allocate(particle(particle_number_temp))
        particle(1:particle_number_survive)=particle_temp(1:particle_number_survive)
        particle_number=particle_number_temp
        particle_number_wave=particle_number_temp-particle_number_survive
        deallocate(particle_temp)

        ! resample particle
        !random number
        allocate(rand_num(particle_number_wave,4))
        do i=1,particle_number_wave
            do k=1,4
                call random_number(rand_num(i,k))
            enddo
        enddo

        ! particle weight, velocity, position
        !$omp parallel
        !$omp do
        do i=1,cell_number
            call resample_particle_cell(i,particle_number_wave,rand_num,particle_number_survive)
        enddo
        !$omp end do nowait
        !$omp end parallel
    end subroutine reinit_particle

    !--------------------------------------------------
    !>sample particle in cell
    !--------------------------------------------------
    subroutine resample_particle_cell(cell_number,particle_number_wave,rand_num,particle_number_survive)
        integer,intent(in) :: cell_number
        integer,intent(in) :: particle_number_wave
        integer,intent(in) :: particle_number_survive
        real(kind=double),intent(in) :: rand_num(particle_number_wave,4)
        real(kind=double) :: theta,phi
        real(kind=double) :: x1_temp,x2_temp
        real(kind=double) :: transMatrix(2,2)
        integer :: particle_id_start,particle_id_end
        integer :: k

        ! transport matrix (1->2,1->3)
        transMatrix(1,1)=node(ctr(cell_number)%node_id(2))%coords(1)-node(ctr(cell_number)%node_id(1))%coords(1)
        transMatrix(2,1)=node(ctr(cell_number)%node_id(2))%coords(2)-node(ctr(cell_number)%node_id(1))%coords(2)
        transMatrix(1,2)=node(ctr(cell_number)%node_id(3))%coords(1)-node(ctr(cell_number)%node_id(1))%coords(1)
        transMatrix(2,2)=node(ctr(cell_number)%node_id(3))%coords(2)-node(ctr(cell_number)%node_id(1))%coords(2)

        ! sample particle
        particle_id_start=particle_number_survive+sum(ctr(1:cell_number-1)%particle_number_wave)+1
        particle_id_end=particle_id_start-1+ctr(cell_number)%particle_number_wave
        do k=particle_id_start,particle_id_end
            !cell number
            particle(k)%cell=cell_number
            particle(k)%cell_new=0
            particle(k)%face_in=0
            particle(k)%face_in_new=0
            !flag
            particle(k)%flag_track=1
            particle(k)%flag_delete=0
            !time
            particle(k)%t=0
            particle(k)%ta=0
            particle(k)%tb=0
            !weight
            particle(k)%weight=ctr(cell_number)%rho_wave*ctr(cell_number)%area/ctr(cell_number)%particle_number_wave
            !velocity
            theta=acos(2.0d0*rand_num(k-particle_number_survive,1)-1)
            phi=2.0d0*pi*rand_num(k-particle_number_survive,2)
            particle(k)%v(1)=lightspeed*sin(theta)*cos(phi)
            particle(k)%v(2)=lightspeed*sin(theta)*sin(phi)
            !position
            x1_temp=rand_num(k-particle_number_survive,3)
            x2_temp=rand_num(k-particle_number_survive,4)
            if(x2_temp>(1.0d0-x1_temp))then
                x1_temp=1.0d0-rand_num(k-particle_number_survive,3)
                x2_temp=1.0d0-rand_num(k-particle_number_survive,4)
            endif
            particle(k)%x(1)=transMatrix(1,1)*x1_temp+transMatrix(1,2)*x2_temp+node(ctr(cell_number)%node_id(1))%coords(1)
            particle(k)%x(2)=transMatrix(2,1)*x1_temp+transMatrix(2,2)*x2_temp+node(ctr(cell_number)%node_id(1))%coords(2)
        enddo
    end subroutine resample_particle_cell
end module microscopic
