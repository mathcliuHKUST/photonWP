!--------------------------------------------------
!>input and output
!--------------------------------------------------
module io
    use global_data
    implicit none
contains
    !--------------------------------------------------
    !>main initialization subroutine
    !--------------------------------------------------
    subroutine init()

        !control
        iteration      = 1
        cfl            = 5.0d-1
        dt             = 1.0d-1/30.0d0
        dt_ex          = 1.0d99
        implicit_factor= 1
        sim_time       = 0.0d0
        end_time       = 1.0d3
        kn             = 1.0d0
        lightspeed     = 29.98d0
        radconst       = 0.01372d0
        ref_number     = 5000
        ref_weight     = 1.0d0
        inputfilename  = 'tophat.msh'
        outputfilename = 'tophat.plt'

        !probes
        probe = 0

        probe_x(1) = 0.25d0
        probe_x(2) = 2.75d0
        probe_x(3) = 3.5d0
        probe_x(4) = 4.25d0
        probe_x(5) = 6.75d0

        probe_y(1) = 0.0d0
        probe_y(2) = 0.0d0
        probe_y(3) = 1.25d0
        probe_y(4) = 0.0d0
        probe_y(5) = 0.0d0

        !method
        method_scheme = SN
        method_scheme = UGKWP
        method_interp = SECOND_ORDER !second order interpolation

        call init_geometry() !initialize the geometry
        call init_velocity() !initialize discrete velocity space
        call init_flow_field() !set the initial condition
        call init_particle() !initialize discrete velocity space
    end subroutine init

    !--------------------------------------------------
    !>initialize the geometry
    !>@param[in] xnum,ynum       :number of cells in x and y direction
    !>@param[in] xlength,ylength :domain length in x and y direction
    !--------------------------------------------------
    subroutine init_geometry()
        type(cell_type),allocatable,dimension(:) :: ctr_temp
        real(kind=double) :: tangent(2),normal(2),center(2)
        real(kind=double) :: tempcosine,innerproduct
        real(kind=double) :: cross_ab,cross_ap,cross_pb
        integer :: boundary_face_number
        integer :: element_number,element_id,element_type,group
        integer :: node_id(3),node_id_temp
        integer :: swap_id
        integer :: temp,tempid(3)
        integer :: sort
        integer :: fileunit
        integer :: cell_number2,boundary_face_number2
        integer :: i,j,k,sort1,sort2
        character(30) :: cline

        !node information
        fileunit=20
        open(unit = fileunit, file = inputfilename, status = 'old')
        do while(.true.)
            read(fileunit,*) cline
            cline = trim(cline)
            if (cline == '$Nodes') exit
        enddo
        read(fileunit,*) node_number
        allocate(node(1:node_number))
        do i = 1, node_number
            node(i)%group=node1
            read(fileunit,*) temp,node(i)%coords(1),node(i)%coords(2)
            node(i)%cell_number=0
        end do
        close(unit = 20)

        !cell number, boundary number
        open(unit = fileunit, file = inputfilename, status = 'old')
        do while(.true.)
            read(fileunit,*) cline
            cline = trim(cline)
            if (cline == '$Elements') exit
        enddo
        read(fileunit,*) element_number
        cell_number=0
        cell_number_inner=0
        boundary_face_number=0
        do i = 1, element_number
            read(fileunit,*) element_id,element_type,temp,group
            if(element_type==2)then
                cell_number=cell_number+1
                if(group==311 .or. group==312)then ! 311-thin, 312-thick
                    cell_number_inner=cell_number_inner+1
                endif
            elseif(element_type==1)then
                boundary_face_number=boundary_face_number+1
            endif
        end do
        close(unit = 20)

        !---------------------------------------------------
        !>read cell information
        !---------------------------------------------------
        !cell information
        allocate(ctr_temp(cell_number))
        do i=1,cell_number
            ctr_temp(i)%group=0
            ctr_temp(i)%node_id=0
            ctr_temp(i)%face_id=0
            ctr_temp(i)%face_number=0
            !geometry
            ctr_temp(i)%coords=0
            ctr_temp(i)%area=0
            ctr_temp(i)%length=0
            ctr_temp(i)%rho=0
            ctr_temp(i)%T=0
            ctr_temp(i)%capacity=0
            ctr_temp(i)%absorbtion=0
        enddo
        open(unit = 20, file = inputfilename, status = 'old')
        do while(.true.)
            read(fileunit,*) cline
            cline = trim(cline)
            if (cline == '$Elements') exit
        enddo
        read(fileunit,*) element_number
        j=0
        do i = 1, element_number
            if(i>(element_number-cell_number))then
                j=j+1
                read(fileunit,*) element_id,element_type,temp,group,temp,node_id(1),node_id(2),node_id(3)
                if(group==311)then     ! inner thin region
                    ctr_temp(j)%group=cell11
                elseif(group==312)then ! inner thick region
                    ctr_temp(j)%group=cell12
                elseif(group==301)then ! left ghost cell
                    ctr_temp(j)%group=cell01
                elseif(group==302)then ! right ghost cell
                    ctr_temp(j)%group=cell02
                else
                    write(*,*) "reading cell information error..."
                    pause
                endif

                !sort node 1<2<3
                do sort1=1,2
                    do sort2=1,3-sort1
                        if(node_id(sort2)<node_id(sort2+1))then
                            swap_id=node_id(sort2)
                            node_id(sort2)=node_id(sort2+1)
                            node_id(sort2+1)=swap_id
                        endif
                    enddo
                enddo
                ctr_temp(j)%node_id=node_id

                if(element_type/=2)then
                    write(*,*) "error: mesh file type error..."
                    pause
                endif
            else
                read(fileunit,*)
            endif
        end do
        close(unit = 20)

        !cell volume
        do i=1,cell_number
            ctr_temp(i)%node_number=3
            tempid(1)=ctr_temp(i)%node_id(1)
            tempid(2)=ctr_temp(i)%node_id(2)
            tempid(3)=ctr_temp(i)%node_id(3)
            ctr_temp(i)%area=abs(0.5d0*((node(tempid(2))%coords(1)-node(tempid(1))%coords(1))*&
                (node(tempid(3))%coords(2)-node(tempid(1))%coords(2))-&
                (node(tempid(2))%coords(2)-node(tempid(1))%coords(2))*&
                (node(tempid(3))%coords(1)-node(tempid(1))%coords(1))))
            ctr_temp(i)%coords=(node(tempid(1))%coords+node(tempid(2))%coords+node(tempid(3))%coords)/3.0d0
        enddo

        !time step
        do i=1,cell_number
            dt_ex=min(dt_ex,cfl*sqrt(ctr_temp(i)%area)/lightspeed)
        enddo
        implicit_factor=max(floor(dt/dt_ex),1)
        dt=implicit_factor*dt_ex
        ref_weight=1.0d0*ctr_temp(1)%area/ref_number

        !-------------------------------------------------
        !> read face information
        !-------------------------------------------------
        face_number=(cell_number*3+boundary_face_number)/2
        allocate(face(face_number))
        !cell interface
        do i=1,face_number
            face(i)%group=face1
            face(i)%node_id=0
            face(i)%cell_id=0
            face(i)%coords=0
            face(i)%length=0
            face(i)%flux=0
            face(i)%weight=0
        end do
        do i=1,face_number
            face(i)%cell_id=0
            face(i)%node_id=0
        enddo
        write(*,*) "geometry information ..."
        j=0
        ! cell id 1
        i=1
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr_temp(i)%node_id(1)
        face(j)%node_id(2)=ctr_temp(i)%node_id(2)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr_temp(i)%node_id(2)
        face(j)%node_id(2)=ctr_temp(i)%node_id(3)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr_temp(i)%node_id(1)
        face(j)%node_id(2)=ctr_temp(i)%node_id(3)
        ! cell id 2,3,...
        do i = 2, cell_number
            !boudary 1 : 1->2
            temp=0
            do k=1,j
                if((face(k)%node_id(1)==ctr_temp(i)%node_id(1) .and. face(k)%node_id(2)==ctr_temp(i)%node_id(2)))then
                    face(k)%cell_id(2)=i
                    temp=1
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr_temp(i)%node_id(1)
                face(j)%node_id(2)=ctr_temp(i)%node_id(2)
            endif
            !boudary 2 : 2->3
            temp=0
            do k=1,j
                if((face(k)%node_id(1)==ctr_temp(i)%node_id(2) .and. face(k)%node_id(2)==ctr_temp(i)%node_id(3)))then
                    face(k)%cell_id(2)=i
                    temp=1
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr_temp(i)%node_id(2)
                face(j)%node_id(2)=ctr_temp(i)%node_id(3)
            endif
            !boudary 3 : 1->3
            temp=0
            do k=1,j
                if((face(k)%node_id(1)==ctr_temp(i)%node_id(1) .and. face(k)%node_id(2)==ctr_temp(i)%node_id(3)))then
                    face(k)%cell_id(2)=i
                    temp=1
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr_temp(i)%node_id(1)
                face(j)%node_id(2)=ctr_temp(i)%node_id(3)
            endif
        enddo
        write(*,*) "geometry information finished, boundary number = ", j

        !> pick boundary
        do i=1,face_number
            if(face(i)%cell_id(2)==0)then
                face(i)%group=face0
                node(face(i)%node_id(1:2))%group=node0
                if(ctr_temp(face(i)%cell_id(1))%group==cell11)then
                    ctr_temp(face(i)%cell_id(1))%group=cell13
                elseif(ctr_temp(face(i)%cell_id(1))%group==cell12)then
                    ctr_temp(face(i)%cell_id(1))%group=cell14
                endif
            endif
        enddo

        !---------------------------------------------------
        !> resort cell
        !---------------------------------------------------
        !cell information
        sort=0
        allocate(ctr(cell_number))
        do i=1,cell_number
            if(ctr_temp(i)%group==cell11)then
                sort=sort+1
                ctr(sort)=ctr_temp(i)
            endif
        enddo
        do i=1,cell_number
            if(ctr_temp(i)%group==cell12)then
                sort=sort+1
                ctr(sort)=ctr_temp(i)
            endif
        enddo
        do i=1,cell_number
            if(ctr_temp(i)%group==cell13)then
                sort=sort+1
                ctr(sort)=ctr_temp(i)
            endif
        enddo
        do i=1,cell_number
            if(ctr_temp(i)%group==cell14)then
                sort=sort+1
                ctr(sort)=ctr_temp(i)
            endif
        enddo
        if(sort/=cell_number_inner)then
            write(*,*) "error: sort/=cell_number_inner..."
            write(*,*) sort,cell_number_inner
            pause
        endif
        do i=1,cell_number
            if(ctr_temp(i)%group==cell01)then
                sort=sort+1
                ctr(sort)=ctr_temp(i)
            endif
        enddo
        do i=1,cell_number
            if(ctr_temp(i)%group==cell02)then
                sort=sort+1
                ctr(sort)=ctr_temp(i)
            endif
        enddo
        if(sort/=cell_number)then
            write(*,*) "error: sort/=cell_number..."
            pause
        endif

        ! node's cell
        do i=1,cell_number
            do j=1,ctr(i)%node_number
                node_id_temp=ctr(i)%node_id(j)
                node(node_id_temp)%cell_number=node(node_id_temp)%cell_number+1
                node(node_id_temp)%cell_id(node(node_id_temp)%cell_number)=i
                !> reflection boundary condition
                !if(node(node_id)%group==node0)then
                !    node(node_id)%cell_number=node(node_id)%cell_number+1
                !    node(node_id)%cell_id(node(node_id)%cell_number)=i
                !endif
            enddo
        enddo

        !face information
        !cell interface
        do i=1,face_number
            face(i)%group=face1
            face(i)%node_id=0
            face(i)%cell_id=0
            face(i)%coords=0
            face(i)%length=0
            face(i)%flux=0
            face(i)%weight=0
        end do
        do i=1,face_number
            face(i)%cell_id=0
            face(i)%node_id=0
        enddo
        write(*,*) "geometry information ..."
        j=0
        ! cell id 1
        i=1
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr(i)%node_id(1)
        face(j)%node_id(2)=ctr(i)%node_id(2)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr(i)%node_id(2)
        face(j)%node_id(2)=ctr(i)%node_id(3)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr(i)%node_id(1)
        face(j)%node_id(2)=ctr(i)%node_id(3)
        ! cell id 2,3,...
        do i = 2, cell_number
            !boundary 1
            temp=0
            do k=1,j
                if((face(k)%node_id(1)==ctr(i)%node_id(1) .and. face(k)%node_id(2)==ctr(i)%node_id(2)))then
                    face(k)%cell_id(2)=i
                    temp=1
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr(i)%node_id(1)
                face(j)%node_id(2)=ctr(i)%node_id(2)
            endif
            !boundary 2
            temp=0
            do k=1,j
                if((face(k)%node_id(1)==ctr(i)%node_id(2) .and. face(k)%node_id(2)==ctr(i)%node_id(3)))then
                    face(k)%cell_id(2)=i
                    temp=1
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr(i)%node_id(2)
                face(j)%node_id(2)=ctr(i)%node_id(3)
            endif
            !boundary 3
            temp=0
            do k=1,j
                if((face(k)%node_id(1)==ctr(i)%node_id(1) .and. face(k)%node_id(2)==ctr(i)%node_id(3)))then
                    face(k)%cell_id(2)=i
                    temp=1
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr(i)%node_id(1)
                face(j)%node_id(2)=ctr(i)%node_id(3)
            endif
        enddo
        write(*,*) "geometry information finished..."

        !> check boundary
        do i=1,face_number
            if(face(i)%cell_id(2)==0)then
                face(i)%group=face0
                if(node(face(i)%node_id(1))%group/=node0 .or.&
                    node(face(i)%node_id(2))%group/=node0)then
                    write(*,*) "error: face boundary 2 ..."
                endif
                if(ctr(face(i)%cell_id(1))%group/= cell13 .and. &
                    ctr(face(i)%cell_id(1))%group/= cell14 .and. &
                    ctr(face(i)%cell_id(1))%group/= cell01 .and. &
                    ctr(face(i)%cell_id(1))%group/= cell02)then
                    write(*,*) "error: face boundary 3 ..."
                endif
            endif
        enddo

        do i=1,face_number
            !> cell's faceid
            ctr(face(i)%cell_id(1))%face_number=ctr(face(i)%cell_id(1))%face_number+1
            ctr(face(i)%cell_id(1))%face_id(ctr(face(i)%cell_id(1))%face_number)=i
            if(face(i)%cell_id(2)/=0)then
                ctr(face(i)%cell_id(2))%face_number=ctr(face(i)%cell_id(2))%face_number+1
                ctr(face(i)%cell_id(2))%face_id(ctr(face(i)%cell_id(2))%face_number)=i
                face(i)%coords=(node(face(i)%node_id(1))%coords+node(face(i)%node_id(2))%coords)/2.0d0
                !> tangent direction node2-->node1
                tangent(1)=node(face(i)%node_id(1))%coords(1)-node(face(i)%node_id(2))%coords(1)
                tangent(2)=node(face(i)%node_id(1))%coords(2)-node(face(i)%node_id(2))%coords(2)
                face(i)%length=sqrt(sum(tangent(:)**2.0d0))
                tangent=tangent/face(i)%length
                !> normal direction center2-->center1
                normal(1)=tangent(2)
                normal(2)=-tangent(1)
                innerproduct=sum(normal(1:2)*(ctr(face(i)%cell_id(1))%coords(1:2)-face(i)%coords(1:2)))
                if(innerproduct<0)then
                    face(i)%norm(1:2) = -normal(1:2)
                else
                    face(i)%norm(1:2) = normal(1:2)
                endif
                !calculate cosine
                ! center vector cell2-->cell1
                center(1:2)=(ctr(face(i)%cell_id(1))%coords(1:2)-ctr(face(i)%cell_id(2))%coords(1:2))/&
                    sqrt((ctr(face(i)%cell_id(1))%coords(1)-ctr(face(i)%cell_id(2))%coords(1))**2.0d0+&
                    (ctr(face(i)%cell_id(1))%coords(2)-ctr(face(i)%cell_id(2))%coords(2))**2.0d0)
                tempcosine=abs(sum(face(i)%norm(1:2)*center(1:2)))
                ! sec(theta_ij)/|cj-ci|
                face(i)%weight(1)=1.0d0/tempcosine/sqrt(sum((ctr(face(i)%cell_id(1))%coords(1:2)-ctr(face(i)%cell_id(2))%coords(1:2))**2.0d0))
                !tan(theta_ij)/|vj-vi|
                face(i)%weight(2)=sum((face(i)%norm-center/tempcosine)*tangent)/face(i)%length
                !check bug
                if(((face(i)%norm(1)-center(1)/tempcosine)*tangent(2)-(face(i)%norm(2)-center(2)/tempcosine)*tangent(1))>1.0d-10)then
                    write(*,*) "error in weight 2: not parallel", &
                        ((face(i)%norm(1)-center(1)/tempcosine)*tangent(2)-(face(i)%norm(2)-center(2)/tempcosine)*tangent(1)),&
                        (face(i)%norm-center/tempcosine),&
                        tangent
                    pause
                endif
                !check bug
                if(abs(tempcosine)>(1.0d0+1.0d-10))then
                    write(*,*) "check face weights..."
                    pause
                endif
            else
                face(i)%coords=(node(face(i)%node_id(1))%coords+node(face(i)%node_id(2))%coords)/2.0d0
                !> tangent direction node2-->node1
                tangent(1)=node(face(i)%node_id(1))%coords(1)-node(face(i)%node_id(2))%coords(1)
                tangent(2)=node(face(i)%node_id(1))%coords(2)-node(face(i)%node_id(2))%coords(2)
                face(i)%length=sqrt(sum(tangent(:)**2.0d0))
                tangent=tangent/face(i)%length
                !> normal direction center2-->center1
                normal(1)=tangent(2)
                normal(2)=-tangent(1)
                innerproduct=sum(normal(1:2)*(ctr(face(i)%cell_id(1))%coords(1:2)-face(i)%coords(1:2)))
                if(innerproduct<0)then
                    face(i)%norm(1:2) = -normal(1:2)
                else
                    face(i)%norm(1:2) = normal(1:2)
                endif
            endif
        enddo

        !----------------------------------------------------------
        !>check bug
        !----------------------------------------------------------
        do i=1,cell_number
            if(ctr(i)%face_number/=3)then
                write(*,*) "check geometry1..."
                pause
            endif
        enddo
        boundary_face_number2=0
        do i=1,face_number
            if(face(i)%group==face0)then
                boundary_face_number2=boundary_face_number2+1
            endif
        enddo
        if(boundary_face_number/=boundary_face_number2)then
            write(*,*) "check geometry2.."
            pause
        endif
        cell_number2=0
        do i=1,node_number
            if(node(i)%group==node1)then
                cell_number2=cell_number2+node(i)%cell_number
            elseif(node(i)%group==node0)then
                cell_number2=cell_number2+node(i)%cell_number
            else
                write(*,*) "node type unfinished..."
                pause
            endif
        enddo
        if(cell_number/=cell_number2/3.0d0)then
            write(*,*) cell_number,cell_number2
            write(*,*) "check geometry3.."
            pause
        endif

        !----------------------------------------------------------
        !>find probes' location
        !----------------------------------------------------------
        do i = 1, 5
            do j = 1, cell_number_inner
                !check if the probe is the cell's node
                do k = 1, 3
                    if ( abs(node(ctr(j)%node_id(k))%coords(1) - probe_x(i)) < 1.0d-10 .and. &
                        abs(node(ctr(j)%node_id(k))%coords(2) - probe_y(i)) < 1.0d-10 ) then
                        probe(i) = j
                        exit
                    end if
                end do

                ! if find the cell, exit
                if (probe(i) /= 0) exit

                !check if the probe is in the cell's boundary
                do k = 1, 3
                    !if y = const
                    if (abs(node(face(ctr(j)%face_id(k))%node_id(1))%coords(2) - node(face(ctr(j)%face_id(k))%node_id(2))%coords(2)) < 1.0d-10) then
                        !if y = probe_y & x_1 < probe_x < x_2
                        if (abs(node(face(ctr(j)%face_id(k))%node_id(1))%coords(2) - probe_y(i)) < 1.0d-10) then
                            if ((node(face(ctr(j)%face_id(k))%node_id(1))%coords(1) < probe_x(i) .and. probe_x(i) < node(face(ctr(j)%face_id(k))%node_id(2))%coords(1)) .or. &
                                (node(face(ctr(j)%face_id(k))%node_id(2))%coords(1) < probe_x(i) .and. probe_x(i) < node(face(ctr(j)%face_id(k))%node_id(1))%coords(1)) ) then
                                probe(i) = j
                                exit  
                            end if
                        end if	      
                    end if

                    !if x = const
                    if (abs(node(face(ctr(j)%face_id(k))%node_id(1))%coords(1) - node(face(ctr(j)%face_id(k))%node_id(2))%coords(1)) < 1.0d-10) then
                        !if x = probe_x & y_1 < probe_y < y_2
                        if (abs(node(face(ctr(j)%face_id(k))%node_id(1))%coords(1) - probe_x(i)) < 1.0d-10) then
                            if ((node(face(ctr(j)%face_id(k))%node_id(1))%coords(2) < probe_y(i) .and. probe_y(i) < node(face(ctr(j)%face_id(k))%node_id(2))%coords(2)) .or. &
                                (node(face(ctr(j)%face_id(k))%node_id(2))%coords(2) < probe_y(i) .and. probe_y(i) < node(face(ctr(j)%face_id(k))%node_id(1))%coords(2)) ) then
                                probe(i) = j
                                exit  
                            end if
                        end if	      
                    end if
                end do

                ! if find the cell, exit
                if (probe(i) /= 0) exit

                !check if the probe is inside the cell
                cross_ab = 0.0d0
                do k = 1, 3
                    !calculate the area of the triangle consist of probe and the nodes of face k	      
                    cross_ab = cross_ab + 0.5d0* &
                        abs( node(face(ctr(j)%face_id(k))%node_id(1))%coords(1)*(node(face(ctr(j)%face_id(k))%node_id(2))%coords(2) - probe_y(i)) + &
                        node(face(ctr(j)%face_id(k))%node_id(2))%coords(1)*(probe_y(i) - node(face(ctr(j)%face_id(k))%node_id(1))%coords(2)) + &
                        probe_x(i)*(node(face(ctr(j)%face_id(k))%node_id(1))%coords(2) - node(face(ctr(j)%face_id(k))%node_id(2))%coords(2)) )	    
                end do
                ! if the sum of the 3 triangle's area is equal to the cell's area, found the location
                if (abs(cross_ab - ctr(j)%area) < 1.0d-10) probe(i) = j

            end do

            ! if still not find the cell, error
            if (probe(i) == 0) then
                write(*,'(A,I4)') "error: can not find probe ", i
                pause
            end if	  
            write(*,'(A,I4,A,I9)') "probe ",i," is in cell ",probe(i)
            write(*,'(A,D15.6,D15.6,D15.6)') "node 1: ",node(ctr(probe(i))%node_id(1))%coords(1),node(ctr(probe(i))%node_id(1))%coords(2)
            write(*,'(A,D15.6,D15.6,D15.6)') "node 2: ",node(ctr(probe(i))%node_id(2))%coords(1),node(ctr(probe(i))%node_id(2))%coords(2)
            write(*,'(A,D15.6,D15.6,D15.6)') "node 3: ",node(ctr(probe(i))%node_id(3))%coords(1),node(ctr(probe(i))%node_id(3))%coords(2)
        end do
    end subroutine init_geometry

    !--------------------------------------------------
    !>set discrete velocity space using Newton¨CCotes formulas
    !>@param[inout] num_u,num_v :number of velocity points
    !>@param[in]    min_u,min_v :smallest discrete velocity
    !>@param[in]    max_u,max_v :largest discrete velocity
    !--------------------------------------------------
    subroutine init_velocity()
        real(kind=double)  :: min_u,max_u,min_v,max_v
        real(kind=double)  :: du,dv !spacing in u and v velocity space
        integer :: i,j

        !velocity number
        unum=32 !polar
        vnum=32 !rotation

        !allocate array
        allocate(uspace(unum,vnum))
        allocate(vspace(unum,vnum))
        allocate(weight(unum,vnum))

        !spacing in u and v velocity space
        min_u =-1.0d0
        max_u = 1.0d0
        min_v = 0.0d0
        max_v = 2.0d0*pi
        du = (max_u-min_u)/unum
        dv = (max_v-min_v)/vnum

        !velocity space
        do i=1,unum
            do j=1,vnum
                uspace(i,j) = lightspeed*sin(acos(min_u+(i-0.5d0)*du))*cos(min_v+(j-0.5d0)*dv)
                vspace(i,j) = lightspeed*sin(acos(min_u+(i-0.5d0)*du))*sin(min_v+(j-0.5d0)*dv)
                weight(i,j) = du*dv
            enddo
        enddo
    end subroutine init_velocity


    !--------------------------------------------------
    !>set the initial condition
    !>@param[in] init_gas :initial condition
    !--------------------------------------------------
    subroutine init_flow_field()
        integer :: i,j

        do i=1,cell_number
            allocate(ctr(i)%h(unum,vnum))
            allocate(ctr(i)%sh(unum,vnum,2))
            ctr(i)%h=0
            ctr(i)%sh=0
        enddo

        do i=1,face_number
            allocate(face(i)%h(unum,vnum))
            allocate(face(i)%sh(unum,vnum,2))
            allocate(face(i)%flux_h(unum,vnum))
            face(i)%h=0
            face(i)%sh=0
            face(i)%flux_h=0
        end do

        !initial condition
        do i=1,cell_number
            if(ctr(i)%group==cell11 .or. ctr(i)%group==cell12 .or.&
               ctr(i)%group==cell13 .or. ctr(i)%group==cell14 .or.&
               ctr(i)%group==cell02)then
                ctr(i)%T          = 5.0d-2
                ctr(i)%phi        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho_wave   = ctr(i)%rho
                ctr(i)%h          = ctr(i)%rho/(4.0d0*pi)
                ctr(i)%sh         = 0.0d0
                ctr(i)%WaveParticleRatio=0.0d0
            elseif(ctr(i)%group==cell01 .and. ctr(i)%coords(2)>0.5)then
                ctr(i)%T          = 5.0d-2
                ctr(i)%phi        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho_wave   = ctr(i)%rho
                ctr(i)%h          = ctr(i)%rho/(4.0d0*pi)
                ctr(i)%sh         = 0.0d0
                ctr(i)%WaveParticleRatio=0.0d0
            elseif(ctr(i)%group==cell01 .and. ctr(i)%coords(2)<0.5)then
                ctr(i)%T          = 5.0d-1
                ctr(i)%phi        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho_wave   = ctr(i)%rho
                ctr(i)%h          = ctr(i)%rho/(4.0d0*pi)
                ctr(i)%sh         = 0.0d0
                ctr(i)%WaveParticleRatio=0.0d0
            endif
        end do

        !cell interface
        do i=1,face_number
            face(i)%h      = 0.0d0
            face(i)%sh     = 0.0d0
            face(i)%flux_h = 0.0d0
        end do
    end subroutine init_flow_field

    !--------------------------------------------------
    !>initialize particle
    !--------------------------------------------------
    subroutine init_particle()
        real(kind=double),allocatable,dimension(:,:) :: rand_num
        integer :: i,j,k

        ! particle number
        particle_number=0
        do i=1,cell_number
            ctr(i)%particle_number=max(int(ctr(i)%rho*ctr(i)%area/ref_weight),0)
            particle_number=particle_number+ctr(i)%particle_number
        enddo
        allocate(particle(particle_number))
        allocate(rand_num(particle_number,4))

        !random number
        do i=1,particle_number
            do k=1,4
                call random_number(rand_num(i,k))
            enddo
        enddo

        ! particle weight, velocity, position
        !$omp parallel
        !$omp do
        do i=1,cell_number
            call sample_particle_cell(i,particle_number,rand_num)
        enddo
        !$omp end do nowait
        !$omp end parallel

        !check bug
        do i=1,particle_number
            !ctr(particle(i)%cell)%rho_wave=ctr(particle(i)%cell)%rho_wave+particle(i)%weight/ctr(particle(i)%cell)%area
            if(particle(i)%x(1)>max(node(ctr(particle(i)%cell)%node_id(1))%coords(1),node(ctr(particle(i)%cell)%node_id(2))%coords(1),node(ctr(particle(i)%cell)%node_id(3))%coords(1)) .or. &
                particle(i)%x(1)<min(node(ctr(particle(i)%cell)%node_id(1))%coords(1),node(ctr(particle(i)%cell)%node_id(2))%coords(1),node(ctr(particle(i)%cell)%node_id(3))%coords(1)) .or. &
                particle(i)%x(2)>max(node(ctr(particle(i)%cell)%node_id(1))%coords(2),node(ctr(particle(i)%cell)%node_id(2))%coords(2),node(ctr(particle(i)%cell)%node_id(3))%coords(2)) .or. &
                particle(i)%x(2)<min(node(ctr(particle(i)%cell)%node_id(1))%coords(2),node(ctr(particle(i)%cell)%node_id(2))%coords(2),node(ctr(particle(i)%cell)%node_id(3))%coords(2)))then
                write(*,*) "particle not in cell..."
                write(*,*) particle(i)%x(1)>max(node(ctr(particle(i)%cell)%node_id(1))%coords(1),node(ctr(particle(i)%cell)%node_id(2))%coords(1),node(ctr(particle(i)%cell)%node_id(3))%coords(1))
                write(*,*) particle(i)%x(1)<min(node(ctr(particle(i)%cell)%node_id(1))%coords(1),node(ctr(particle(i)%cell)%node_id(2))%coords(1),node(ctr(particle(i)%cell)%node_id(3))%coords(1))
                write(*,*) particle(i)%x(2)>max(node(ctr(particle(i)%cell)%node_id(1))%coords(2),node(ctr(particle(i)%cell)%node_id(2))%coords(2),node(ctr(particle(i)%cell)%node_id(3))%coords(2))
                write(*,*) particle(i)%x(2)<min(node(ctr(particle(i)%cell)%node_id(1))%coords(2),node(ctr(particle(i)%cell)%node_id(2))%coords(2),node(ctr(particle(i)%cell)%node_id(3))%coords(2))
                write(*,*) node(ctr(particle(i)%cell)%node_id(1))%coords
                write(*,*) node(ctr(particle(i)%cell)%node_id(2))%coords
                write(*,*) node(ctr(particle(i)%cell)%node_id(3))%coords
                pause
            endif
        enddo
    end subroutine init_particle

    !--------------------------------------------------
    !>sample particle in cell
    !--------------------------------------------------
    subroutine sample_particle_cell(cell_number,particle_number,rand_num)
        integer,intent(in) :: cell_number
        integer,intent(in) :: particle_number
        real(kind=double),intent(in) :: rand_num(particle_number,4)
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
        particle_id_start=sum(ctr(1:cell_number-1)%particle_number)+1
        particle_id_end=particle_id_start-1+ctr(cell_number)%particle_number
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
            particle(k)%weight=ctr(cell_number)%rho*ctr(cell_number)%area/ctr(cell_number)%particle_number
            !velocity
            theta=acos(2.0d0*rand_num(k,1)-1)
            phi=2.0d0*pi*rand_num(k,2)
            particle(k)%v(1)=lightspeed*sin(theta)*cos(phi)
            particle(k)%v(2)=lightspeed*sin(theta)*sin(phi)
            !position
            x1_temp=rand_num(k,3)
            x2_temp=rand_num(k,4)
            if(x2_temp>(1.0d0-x1_temp))then
                x1_temp=1.0d0-rand_num(k,3)
                x2_temp=1.0d0-rand_num(k,4)
            endif
            particle(k)%x(1)=transMatrix(1,1)*x1_temp+transMatrix(1,2)*x2_temp+node(ctr(cell_number)%node_id(1))%coords(1)
            particle(k)%x(2)=transMatrix(2,1)*x1_temp+transMatrix(2,2)*x2_temp+node(ctr(cell_number)%node_id(1))%coords(2)
        enddo
    end subroutine sample_particle_cell

    !--------------------------------------------------
    !>write result
    !--------------------------------------------------
    subroutine output()
        character(len=10) :: outfilename0
        character(len=100) :: outfilename1
        integer :: i,j

        write(outfilename0,'(I10.10)') iteration
        outfilename1=trim(outfilename0)//trim(outputfilename)
        open(unit = 21, file = outfilename1, status = 'replace')
        write(21,'(a)') 'VARIABLES=X,Y,rho,T,kappa,ParticlePercentage'
        write(21,*) 'ZONE N=',node_number,',E=',cell_number_inner,',datapacking=block'
        write(21,'(a)') 'varlocation=([3-6]=cellcentered)'
        write(21,'(a)') 'ZONETYPE=FETRIANGLE'

        !node data
        do i = 1, node_number
            write(21,*) node(i)%coords(1)
        end do

        do i = 1, node_number
            write(21,*) node(i)%coords(2)
        end do

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%rho
        enddo

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%T
        enddo

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%kappa_effective
        enddo

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%WaveParticleRatio
        enddo

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%node_id(1),ctr(i)%node_id(2),ctr(i)%node_id(3)
        enddo
        close(unit = 21)
    end subroutine output

    !--------------------------------------------------
    !>write result
    !--------------------------------------------------
    subroutine output2()
        character(len=10) :: outfilename0
        character(len=100) :: outfilename1
        integer :: i,j

        write(outfilename0,'(I10.10)') iteration
        outfilename1=trim(outfilename0)//trim(outputfilename)
        open(unit = 21, file = outfilename1, status = 'replace')
        write(21,'(a)') 'VARIABLES=X,Y,rho,Tr,Te,kappa,ParticlePercentage'
        write(21,*) 'ZONE N=',node_number*2,',E=',cell_number*2,',datapacking=block'
        write(21,'(a)') 'varlocation=([3-7]=cellcentered)'
        write(21,'(a)') 'ZONETYPE=FETRIANGLE'

        !node data
        do i = 1, node_number
            write(21,*) node(i)%coords(1)
        end do
        do i = 1, node_number
            write(21,*) node(i)%coords(1)
        end do

        do i = 1, node_number
            write(21,*) node(i)%coords(2)
        end do
        do i = 1, node_number
            write(21,*) -node(i)%coords(2)
        end do

        do i = 1, cell_number
            write(21,*) ctr(i)%rho
        enddo
        do i = 1, cell_number
            write(21,*) ctr(i)%rho
        enddo

        do i = 1, cell_number
            write(21,*) (max(ctr(i)%rho/radconst/lightspeed,1.0d-10))**0.25d0
        enddo
        do i = 1, cell_number
            write(21,*) (max(ctr(i)%rho/radconst/lightspeed,1.0d-10))**0.25d0
        enddo

        do i = 1, cell_number
            write(21,*) ctr(i)%T
        enddo
        do i = 1, cell_number
            write(21,*) ctr(i)%T
        enddo

        do i = 1, cell_number
            write(21,*) ctr(i)%kappa_effective
        enddo
        do i = 1, cell_number
            write(21,*) ctr(i)%kappa_effective
        enddo

        do i = 1, cell_number
            write(21,*) ctr(i)%WaveParticleRatio
        enddo
        do i = 1, cell_number
            write(21,*) ctr(i)%WaveParticleRatio
        enddo

        do i = 1, cell_number
            write(21,*) ctr(i)%node_id(1),ctr(i)%node_id(2),ctr(i)%node_id(3)
        enddo
        do i = 1, cell_number
            write(21,*) ctr(i)%node_id(1)+node_number,ctr(i)%node_id(2)+node_number,ctr(i)%node_id(3)+node_number
        enddo
        close(unit = 21)
    end subroutine output2
end module io
