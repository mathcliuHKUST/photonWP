!--------------------------------------------------
!>input and output
!--------------------------------------------------
module io
    use global_data
    use lapack95
    use f95_precision
    implicit none
contains
    !--------------------------------------------------
    !>main initialization subroutine
    !--------------------------------------------------
    subroutine init()

        !control
        iteration     = 1
        dt            = 1.0d99
        sim_time      = 0.0d0
        end_time      = 1.0d0
        cfl           = 1.0d-2
        kn            = 1.0d0
        lightspeed    = 29.98d0
        radconst      = 1.372d-2
        ref_number    = 50000
        ref_weight    = 1.0d0
        inputfilename = 'hahlraum3d_f.msh'
        outputfilename= 'hahlraum3d'

        !method
        method_scheme = SN
        method_scheme = UGKWP
        method_interp = SECOND_ORDER !second order interpolation

        call init_geometry() !initialize the geometry
        !call init_velocity() !initialize discrete velocity space
        call init_flow_field() !set the initial condition
        !call init_particle() !initialize discrete velocity space
    end subroutine init

    !--------------------------------------------------
    !>initialize the geometry
    !--------------------------------------------------
    subroutine init_geometry()
        type(cell_type),allocatable,dimension(:) :: ctr_temp
        real(kind=double)   :: local(3,3),local_norm(3)
        real(kind=double)   :: norm(3),norm_norm
        real(kind=double)   :: innerproduct
        real(kind=double)   :: v1(3),v2(3),v3(3)
        real(kind=double)   :: n1(3),n2(3),n3(3),n4(3)
        real(kind=double)   :: matrix2(2,2),matrix3(3,3)
        real(kind=double)   :: RHS2(2),RHS3(3)
        real(kind=double)   :: IPIV2(2),IPIV3(3)
        real(kind=double)   :: volume
        integer(kind=short) :: boundary_face_number
        integer(kind=short) :: element_number,element_id,element_type,group
        integer(kind=short) :: node_id(4),swap_id
        integer(kind=short) :: temp_id
        integer(kind=short) :: temp
        integer(kind=short) :: tempid(4)
        integer(kind=short) :: i,j,k,sort1,sort2
        integer(kind=short) :: sort
        integer(kind=short) :: fileunit
        integer(kind=short) :: cell_number2,boundary_face_number2
        integer(kind=short) :: INFO
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
            read(fileunit,*) temp,node(i)%coords(1),node(i)%coords(2),node(i)%coords(3)
            node(i)%cell_number=0
            if(temp /=i) then
                write(*,*) "error: node information ..."
                pause
            endif
        end do
        close(unit = fileunit)

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
            read(fileunit,*) element_id,element_type,temp,group,temp
            if(element_type==4)then
                cell_number=cell_number+1
                if(group==311 .or. group==312 .or. group==313)then
                    cell_number_inner=cell_number_inner+1
                endif
            elseif(element_type==2)then
                boundary_face_number=boundary_face_number+1
            endif
        end do
        close(unit = fileunit)

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
            ctr_temp(i)%volume=0
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
        do i = 1,element_number-cell_number
            read(fileunit,*) element_id,element_type
            if(element_type==4)then
                write(*,*) "error: cell nodes reading..."
            endif
        enddo
        do i = element_number-cell_number+1, element_number
            read(fileunit,*) element_id,element_type,temp,group,temp,node_id(1),node_id(2),node_id(3),node_id(4)
            if(element_type==4)then
                j=j+1
                if(group==311)then
                    ctr_temp(j)%group=cell11
                elseif(group==312)then
                    ctr_temp(j)%group=cell12
                elseif(group==313)then
                    ctr_temp(j)%group=cell15
                elseif(group==301)then
                    ctr_temp(j)%group=cell01
                elseif(group==302)then
                    ctr_temp(j)%group=cell02
                else
                    write(*,*) "error: reading cell information..."
                    pause
                endif
                !sort node 1<2<3<4
                do sort1=1,3
                    do sort2=1,4-sort1
                        if(node_id(sort2)<node_id(sort2+1))then
                            swap_id=node_id(sort2)
                            node_id(sort2)=node_id(sort2+1)
                            node_id(sort2+1)=swap_id
                        endif
                    enddo
                enddo
                ctr_temp(j)%node_number=4
                ctr_temp(j)%node_id=node_id
            else
                !read(fileunit,*)
                write(*,*) "error: cell nodes reading..."
            endif
        end do
        close(unit = fileunit)

        !cell volume
        do i=1,cell_number
            !node
            tempid=ctr_temp(i)%node_id
            n1(1:3)=node(tempid(1))%coords(1:3)
            n2(1:3)=node(tempid(2))%coords(1:3)
            n3(1:3)=node(tempid(3))%coords(1:3)
            n4(1:3)=node(tempid(4))%coords(1:3)
            !vector
            v1(1:3)=n2(1:3)-n1(1:3)
            v2(1:3)=n3(1:3)-n1(1:3)
            v3(1:3)=n4(1:3)-n1(1:3)
            ctr_temp(i)%volume=abs(-n1(1)*n2(2)*n3(3)+n1(1)*n2(2)*n4(3)+n1(1)*n3(2)*n2(3)-&
                n1(1)*n3(2)*n4(3)-n1(1)*n4(2)*n2(3)+n1(1)*n4(2)*n3(3)+&
                n2(1)*n1(2)*n3(3)-n2(1)*n1(2)*n4(3)-n2(1)*n3(2)*n1(3)+&
                n2(1)*n3(2)*n4(3)+n2(1)*n4(2)*n1(3)-n2(1)*n4(2)*n3(3)-&
                n3(1)*n1(2)*n2(3)+n3(1)*n1(2)*n4(3)+n3(1)*n2(2)*n1(3)-&
                n3(1)*n2(2)*n4(3)-n3(1)*n4(2)*n1(3)+n3(1)*n4(2)*n2(3)+&
                n4(1)*n1(2)*n2(3)-n4(1)*n1(2)*n3(3)-n4(1)*n2(2)*n1(3)+&
                n4(1)*n2(2)*n3(3)+n4(1)*n3(2)*n1(3)-n4(1)*n3(2)*n2(3))/6.0d0
            volume=abs(v3(3)*(v1(1)*v2(2)-v1(2)*v2(1))+&
                v3(2)*(v1(3)*v2(1)-v1(1)*v2(3))+&
                v3(1)*(v1(2)*v2(3)-v1(3)*v2(2)))/6.0d0
            if(abs(ctr_temp(i)%volume-volume)/volume>1.0d-4)then
                write(*,*) "error: volume calculation... "
                write(*,*) ctr_temp(i)%volume ,volume
                pause
            endif
            ctr_temp(i)%coords=(node(tempid(1))%coords+node(tempid(2))%coords+&
                node(tempid(3))%coords+node(tempid(4))%coords)/4.0d0
        enddo

        !time step
        do i=1,cell_number
            dt=min(dt,cfl*(ctr_temp(i)%volume)**0.33/lightspeed)
        enddo
        ref_weight=1.0d0*ctr_temp(1)%volume/ref_number

        !-------------------------------------------------
        !> read face information
        !-------------------------------------------------
        write(*,*) "geometry information ..."
        face_number=(cell_number*4+boundary_face_number)/2
        allocate(face(face_number))
        !cell interface
        do i=1,face_number
            face(i)%group=face1
            face(i)%node_id=0
            face(i)%cell_id=0
            face(i)%coords=0
            face(i)%area=0
            face(i)%flux=0
            face(i)%weight=0
        enddo
        j=0
        ! cell id 1
        i=1
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr_temp(i)%node_id(1)
        face(j)%node_id(2)=ctr_temp(i)%node_id(2)
        face(j)%node_id(3)=ctr_temp(i)%node_id(3)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr_temp(i)%node_id(1)
        face(j)%node_id(2)=ctr_temp(i)%node_id(2)
        face(j)%node_id(3)=ctr_temp(i)%node_id(4)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr_temp(i)%node_id(1)
        face(j)%node_id(2)=ctr_temp(i)%node_id(3)
        face(j)%node_id(3)=ctr_temp(i)%node_id(4)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr_temp(i)%node_id(2)
        face(j)%node_id(2)=ctr_temp(i)%node_id(3)
        face(j)%node_id(3)=ctr_temp(i)%node_id(4)
        ! cell id 2,3,...
        do i = 2, cell_number
            !boudary 1 = 1,2,3
            temp=0
            do k=1,j
                if(face(k)%cell_id(2)==0)then
                    if((face(k)%node_id(1)==ctr_temp(i)%node_id(1) .and. &
                        face(k)%node_id(2)==ctr_temp(i)%node_id(2) .and. &
                        face(k)%node_id(3)==ctr_temp(i)%node_id(3)))then
                        face(k)%cell_id(2)=i
                        temp=1
                    endif
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr_temp(i)%node_id(1)
                face(j)%node_id(2)=ctr_temp(i)%node_id(2)
                face(j)%node_id(3)=ctr_temp(i)%node_id(3)
            endif
            !boudary 2 = 1,2,4
            temp=0
            do k=1,j
                if(face(k)%cell_id(2)==0)then
                    if((face(k)%node_id(1)==ctr_temp(i)%node_id(1) .and. &
                        face(k)%node_id(2)==ctr_temp(i)%node_id(2) .and. &
                        face(k)%node_id(3)==ctr_temp(i)%node_id(4)))then
                        face(k)%cell_id(2)=i
                        temp=1
                    endif
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr_temp(i)%node_id(1)
                face(j)%node_id(2)=ctr_temp(i)%node_id(2)
                face(j)%node_id(3)=ctr_temp(i)%node_id(4)
            endif
            !boudary 3 = 1,3,4
            temp=0
            do k=1,j
                if(face(k)%cell_id(2)==0)then
                    if((face(k)%node_id(1)==ctr_temp(i)%node_id(1) .and. &
                        face(k)%node_id(2)==ctr_temp(i)%node_id(3) .and. &
                        face(k)%node_id(3)==ctr_temp(i)%node_id(4)))then
                        face(k)%cell_id(2)=i
                        temp=1
                    endif
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr_temp(i)%node_id(1)
                face(j)%node_id(2)=ctr_temp(i)%node_id(3)
                face(j)%node_id(3)=ctr_temp(i)%node_id(4)
            endif
            !boudary 4 = 2,3,4
            temp=0
            do k=1,j
                if(face(k)%cell_id(2)==0)then
                    if((face(k)%node_id(1)==ctr_temp(i)%node_id(2) .and. &
                        face(k)%node_id(2)==ctr_temp(i)%node_id(3) .and. &
                        face(k)%node_id(3)==ctr_temp(i)%node_id(4)))then
                        face(k)%cell_id(2)=i
                        temp=1
                    endif
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr_temp(i)%node_id(2)
                face(j)%node_id(2)=ctr_temp(i)%node_id(3)
                face(j)%node_id(3)=ctr_temp(i)%node_id(4)
            endif
        enddo !end cell loop
        write(*,*) "geometry information finished: boundary number = ", j

        !> pick boundary face, boundary node, boundary cell
        do i=1,face_number
            if(face(i)%cell_id(2)==0)then
                face(i)%group=face0
                node(face(i)%node_id(1:3))%group=node0
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
        do i=1,cell_number
            if(ctr_temp(i)%group==cell15)then
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
                temp_id=ctr(i)%node_id(j)
                node(temp_id)%cell_number=node(temp_id)%cell_number+1
                node(temp_id)%cell_id(node(temp_id)%cell_number)=i
                !if(node(temp_id)%group==node0)then
                !    node(temp_id)%cell_number=node(temp_id)%cell_number+1
                !    node(temp_id)%cell_id(node(temp_id)%cell_number)=i
                !endif
            enddo
        enddo

        !face information
        write(*,*) "geometry information ..."
        do i=1,face_number
            face(i)%group=face1
            face(i)%node_id=0
            face(i)%cell_id=0
            face(i)%coords=0
            face(i)%area=0
            face(i)%flux=0
            face(i)%weight=0
        enddo
        j=0
        ! cell id 1
        i=1
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr(i)%node_id(1)
        face(j)%node_id(2)=ctr(i)%node_id(2)
        face(j)%node_id(3)=ctr(i)%node_id(3)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr(i)%node_id(1)
        face(j)%node_id(2)=ctr(i)%node_id(2)
        face(j)%node_id(3)=ctr(i)%node_id(4)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr(i)%node_id(1)
        face(j)%node_id(2)=ctr(i)%node_id(3)
        face(j)%node_id(3)=ctr(i)%node_id(4)
        j=j+1
        face(j)%cell_id(1)=i
        face(j)%node_id(1)=ctr(i)%node_id(2)
        face(j)%node_id(2)=ctr(i)%node_id(3)
        face(j)%node_id(3)=ctr(i)%node_id(4)
        ! cell id 2,3,...
        do i = 2, cell_number
            !boudary 1 = 1,2,3
            temp=0
            do k=1,j
                if(face(k)%cell_id(2)==0)then
                    if((face(k)%node_id(1)==ctr(i)%node_id(1) .and. &
                        face(k)%node_id(2)==ctr(i)%node_id(2) .and. &
                        face(k)%node_id(3)==ctr(i)%node_id(3)))then
                        face(k)%cell_id(2)=i
                        temp=1
                    endif
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr(i)%node_id(1)
                face(j)%node_id(2)=ctr(i)%node_id(2)
                face(j)%node_id(3)=ctr(i)%node_id(3)
            endif
            !boudary 2 = 1,2,4
            temp=0
            do k=1,j
                if(face(k)%cell_id(2)==0)then
                    if((face(k)%node_id(1)==ctr(i)%node_id(1) .and. &
                        face(k)%node_id(2)==ctr(i)%node_id(2) .and. &
                        face(k)%node_id(3)==ctr(i)%node_id(4)))then
                        face(k)%cell_id(2)=i
                        temp=1
                    endif
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr(i)%node_id(1)
                face(j)%node_id(2)=ctr(i)%node_id(2)
                face(j)%node_id(3)=ctr(i)%node_id(4)
            endif
            !boudary 3 = 1,3,4
            temp=0
            do k=1,j
                if(face(k)%cell_id(2)==0)then
                    if((face(k)%node_id(1)==ctr(i)%node_id(1) .and. &
                        face(k)%node_id(2)==ctr(i)%node_id(3) .and. &
                        face(k)%node_id(3)==ctr(i)%node_id(4)))then
                        face(k)%cell_id(2)=i
                        temp=1
                    endif
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr(i)%node_id(1)
                face(j)%node_id(2)=ctr(i)%node_id(3)
                face(j)%node_id(3)=ctr(i)%node_id(4)
            endif
            !boudary 4 = 2,3,4
            temp=0
            do k=1,j
                if(face(k)%cell_id(2)==0)then
                    if((face(k)%node_id(1)==ctr(i)%node_id(2) .and. &
                        face(k)%node_id(2)==ctr(i)%node_id(3) .and. &
                        face(k)%node_id(3)==ctr(i)%node_id(4)))then
                        face(k)%cell_id(2)=i
                        temp=1
                    endif
                endif
            enddo
            if(temp==0)then
                j=j+1
                face(j)%cell_id(1)=i
                face(j)%node_id(1)=ctr(i)%node_id(2)
                face(j)%node_id(2)=ctr(i)%node_id(3)
                face(j)%node_id(3)=ctr(i)%node_id(4)
            endif
        enddo !end cell loop
        write(*,*) "geometry information finished: boundary number = ", j

        !> check boundary
        do i=1,face_number
            if(face(i)%cell_id(2)==0)then
                face(i)%group=face0
                if(node(face(i)%node_id(1))%group/=node0 .or.&
                    node(face(i)%node_id(2))%group/=node0 .or.&
                    node(face(i)%node_id(3))%group/=node0)then
                    write(*,*) "error: face boundary 2 ..."
                    pause
                endif
                if(ctr(face(i)%cell_id(1))%group/= cell13 .and. &
                    ctr(face(i)%cell_id(1))%group/= cell14 .and. &
                    ctr(face(i)%cell_id(1))%group/= cell15 .and. &
                    ctr(face(i)%cell_id(1))%group/= cell01 .and. &
                    ctr(face(i)%cell_id(1))%group/= cell02)then
                    write(*,*) "error: face boundary 3 ..."
                    pause
                endif
            endif
        enddo

        !> 1.cell's faceid
        !> 2.face area
        !> 3.face norm
        !> 4.face weights
        do i=1,face_number
            !> cell's faceid
            ctr(face(i)%cell_id(1))%face_number=ctr(face(i)%cell_id(1))%face_number+1
            ctr(face(i)%cell_id(1))%face_id(ctr(face(i)%cell_id(1))%face_number)=i
            if(face(i)%cell_id(2)/=0)then
                ctr(face(i)%cell_id(2))%face_number=ctr(face(i)%cell_id(2))%face_number+1
                ctr(face(i)%cell_id(2))%face_id(ctr(face(i)%cell_id(2))%face_number)=i
                face(i)%coords=(node(face(i)%node_id(1))%coords+node(face(i)%node_id(2))%coords+node(face(i)%node_id(3))%coords)/3.0d0
                !> local direction: c2->c1, node2-->node1, node3-->node1
                local(:,1)=ctr(face(i)%cell_id(1))%coords-ctr(face(i)%cell_id(2))%coords
                local(:,2)=node(face(i)%node_id(1))%coords-node(face(i)%node_id(2))%coords
                local(:,3)=node(face(i)%node_id(1))%coords-node(face(i)%node_id(3))%coords
                face(i)%area=sqrt((local(1,2)*local(2,3)-local(2,2)*local(1,3))**2.0d0+&
                    (local(3,2)*local(1,3)-local(1,2)*local(3,3))**2.0d0+&
                    (local(2,2)*local(3,3)-local(3,2)*local(2,3))**2.0d0)/2.0d0
                if(abs(local(1,2)*local(2,3)-local(2,2)*local(1,3))>1.0d-8)then
                    matrix2(1,1)= local(1,2)
                    matrix2(1,2)= local(2,2)
                    RHS2(1)     =-local(3,2)
                    matrix2(2,1)= local(1,3)
                    matrix2(2,2)= local(2,3)
                    RHS2(2)     =-local(3,3)
                    call dgesv(2,1,matrix2,2,IPIV2,RHS2,2,INFO)
                    if(INFO/=0)then
                        write(*,*) "linear solver dgesv error..."
                        pause
                    endif
                    norm(1)=RHS2(1)
                    norm(2)=RHS2(2)
                    norm(3)=1
                elseif(abs(local(1,2)*local(3,3)-local(3,2)*local(1,3))>1.0d-8)then
                    matrix2(1,1)= local(1,2)
                    RHS2(1)     =-local(2,2)
                    matrix2(1,2)= local(3,2)
                    matrix2(2,1)= local(1,3)
                    RHS2(2)     =-local(2,3)
                    matrix2(2,2)= local(3,3)
                    call dgesv(2,1,matrix2,2,IPIV2,RHS2,2,INFO)
                    if(INFO/=0)then
                        write(*,*) "linear solver dgesv error..."
                        pause
                    endif
                    norm(1)=RHS2(1)
                    norm(2)=1
                    norm(3)=RHS2(2)
                elseif(abs(local(2,2)*local(3,3)-local(3,2)*local(2,3))>1.0d-8)then
                    RHS2(1)     =-local(1,2)
                    matrix2(1,1)= local(2,2)
                    matrix2(1,2)= local(3,2)
                    RHS2(2)     =-local(1,3)
                    matrix2(2,1)= local(2,3)
                    matrix2(2,2)= local(3,3)
                    call dgesv(2,1,matrix2,2,IPIV2,RHS2,2,INFO)
                    if(INFO/=0)then
                        write(*,*) "linear solver dgesv error..."
                        pause
                    endif
                    norm(1)=1
                    norm(2)=RHS2(1)
                    norm(3)=RHS2(2)
                else
                    write(*,*) "error: calculate normal vector..."
                    pause
                endif
                !nomalize local and normal
                do j=1,3
                    local_norm(j)=sqrt(sum(local(:,j)**2.0d0))
                    local(:,j)=local(:,j)/local_norm(j)
                enddo
                norm_norm=sqrt(sum(norm**2.0d0))
                norm=norm/norm_norm
                !> normal direction parallel to local_1 c2->c1
                innerproduct=sum(norm*local(:,1))
                if(innerproduct<0)then
                    norm = -norm
                endif
                face(i)%norm = norm
                if(abs(sum(norm*local(:,2)))>1.0d-10 .or. abs(sum(norm*local(:,3)))>1.0d-10)then
                    write(*,*) "error: normal direction..."
                    write(*,*) abs(sum(norm*local(:,2))),abs(sum(norm*local(:,3)))
                    pause
                endif
                !> calculate weights
                !> center vector cell2-->cell1
                matrix3(1:3,1)=local(:,1) ! center2->center1
                matrix3(1:3,2)=local(:,2) ! node2->node1
                matrix3(1:3,3)=local(:,3) ! node3->node1
                RHS3(1:3)=norm
                call dgesv(3,1,matrix3,3,IPIV3,RHS3,3,INFO)
                if(INFO/=0)then
                    write(*,*) "linear solver dgesv error..."
                    pause
                endif
                if(abs((RHS3(1)*local(1,1)+RHS3(2)*local(1,2)+RHS3(3)*local(1,3))-norm(1))>1.0d-10 .or.&
                    abs((RHS3(1)*local(2,1)+RHS3(2)*local(2,2)+RHS3(3)*local(2,3))-norm(2))>1.0d-10 .or.&
                    abs((RHS3(1)*local(3,1)+RHS3(2)*local(3,2)+RHS3(3)*local(3,3))-norm(3))>1.0d-10)then
                    write(*,*) "error: face weights..."
                    pause
                else
                    face(i)%weight(1)=RHS3(1)/local_norm(1)
                    face(i)%weight(2)=RHS3(2)/local_norm(2)
                    face(i)%weight(3)=RHS3(3)/local_norm(3)
                endif
            else
                face(i)%coords=(node(face(i)%node_id(1))%coords+node(face(i)%node_id(2))%coords+node(face(i)%node_id(3))%coords)/3.0d0
                !> local direction: c2->c1, node2-->node1, node3-->node1
                local(:,1)=ctr(face(i)%cell_id(1))%coords-face(i)%coords
                local(:,2)=node(face(i)%node_id(1))%coords-node(face(i)%node_id(2))%coords
                local(:,3)=node(face(i)%node_id(1))%coords-node(face(i)%node_id(3))%coords
                face(i)%area=sqrt((local(1,2)*local(2,3)-local(2,2)*local(1,3))**2.0d0+&
                    (local(3,2)*local(1,3)-local(1,2)*local(3,3))**2.0d0+&
                    (local(2,2)*local(3,3)-local(3,2)*local(2,3))**2.0d0)/2.0d0
                if(abs(local(1,2)*local(2,3)-local(2,2)*local(1,3))>1.0d-8)then
                    matrix2(1,1)= local(1,2)
                    matrix2(1,2)= local(2,2)
                    RHS2(1)     =-local(3,2)
                    matrix2(2,1)= local(1,3)
                    matrix2(2,2)= local(2,3)
                    RHS2(2)     =-local(3,3)
                    call dgesv(2,1,matrix2,2,IPIV2,RHS2,2,INFO)
                    if(INFO/=0)then
                        write(*,*) "linear solver dgesv error..."
                        pause
                    endif
                    norm(1)=RHS2(1)
                    norm(2)=RHS2(2)
                    norm(3)=1
                elseif(abs(local(1,2)*local(3,3)-local(3,2)*local(1,3))>1.0d-8)then
                    matrix2(1,1)= local(1,2)
                    RHS2(1)     =-local(2,2)
                    matrix2(1,2)= local(3,2)
                    matrix2(2,1)= local(1,3)
                    RHS2(2)     =-local(2,3)
                    matrix2(2,2)= local(3,3)
                    call dgesv(2,1,matrix2,2,IPIV2,RHS2,2,INFO)
                    if(INFO/=0)then
                        write(*,*) "linear solver dgesv error..."
                        pause
                    endif
                    norm(1)=RHS2(1)
                    norm(2)=1
                    norm(3)=RHS2(2)
                elseif(abs(local(2,2)*local(3,3)-local(3,2)*local(2,3))>1.0d-8)then
                    RHS2(1)     =-local(1,2)
                    matrix2(1,1)= local(2,2)
                    matrix2(1,2)= local(3,2)
                    RHS2(2)     =-local(1,3)
                    matrix2(2,1)= local(2,3)
                    matrix2(2,2)= local(3,3)
                    call dgesv(2,1,matrix2,2,IPIV2,RHS2,2,INFO)
                    if(INFO/=0)then
                        write(*,*) "linear solver dgesv error..."
                        pause
                    endif
                    norm(1)=1
                    norm(2)=RHS2(1)
                    norm(3)=RHS2(2)
                else
                    write(*,*) "error: calculate normal vector..."
                    pause
                endif
                !nomalize local and normal
                do j=1,3
                    local_norm(j)=sqrt(sum(local(:,j)**2.0d0))
                    local(:,j)=local(:,j)/local_norm(j)
                enddo
                norm_norm=sqrt(sum(norm**2.0d0))
                norm=norm/norm_norm
                !> normal direction parallel to local_1 c2->c1
                innerproduct=sum(norm*local(:,1))
                if(innerproduct<0)then
                    norm = -norm
                endif
                face(i)%norm = norm
                if(abs(sum(norm*local(:,2)))>1.0d-10 .or. abs(sum(norm*local(:,3)))>1.0d-10)then
                    write(*,*) "error: normal direction..."
                    pause
                endif
                !> calculate weights
                face(i)%weight(1)=1
                face(i)%weight(2:3)=0
            endif
        enddo

        !----------------------------------------------------------
        !>check bugs
        !----------------------------------------------------------
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
                !cell_number2=cell_number2+node(i)%cell_number/2.0d0
                cell_number2=cell_number2+node(i)%cell_number
            else
                write(*,*) "node type unfinished..."
                pause
            endif
        enddo
        if(cell_number/=cell_number2/4.0d0)then
            write(*,*) cell_number,cell_number2
            write(*,*) "check geometry3.."
            pause
        endif
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
        unum=16 !polar
        vnum=16 !rotation

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
                ctr(i)%T          = 1.0d-2
                ctr(i)%phi        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho        = radconst*lightspeed*ctr(i)%T**4.0d0
                ctr(i)%rho_wave   = ctr(i)%rho
                ctr(i)%h          = ctr(i)%rho/(4.0d0*pi)
                ctr(i)%sh         = 0.0d0
                ctr(i)%WaveParticleRatio=0.0d0
            elseif(ctr(i)%group==cell01)then
                ctr(i)%T          = 1.0d0
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
            ctr(i)%particle_number=max(int(ctr(i)%rho*ctr(i)%volume/ref_weight),0)
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
            !ctr(particle(i)%cell)%rho_wave=ctr(particle(i)%cell)%rho_wave+particle(i)%weight/ctr(particle(i)%cell)%volume
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
            particle(k)%weight=ctr(cell_number)%rho*ctr(cell_number)%volume/ctr(cell_number)%particle_number
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
        write(21,'(a)') "# vtk DataFile Version 2.0"
        write(21,'(a)') "Unstructured Grid Example"
        write(21,'(a)') "ASCII"
        write(21,'(a)') "DATASET UNSTRUCTURED_GRID"
        write(21,*) "POINTS",node_number,"float"
        do i = 1, node_number
            write(21,*) node(i)%coords(1),node(i)%coords(2),node(i)%coords(3)
        end do
        write(21,*) "CELLS",cell_number,cell_number*5
        do i = 1, cell_number
            write(21,*) "4 ",ctr(i)%node_id(1)-1,ctr(i)%node_id(2)-1,ctr(i)%node_id(3)-1,ctr(i)%node_id(4)-1
        end do
        write(21,*) "CELL_TYPES", cell_number
        do i = 1, cell_number
            write(21,*) "10"
        end do

        write(21,*) "CELL_DATA", cell_number
        write(21,*) "SCALARS rho float", 1
        write(21,*) "LOOKUP_TABLE default"
        do i = 1, cell_number
            write(21,*) ctr(i)%rho
        end do

        write(21,*) "SCALARS T float", 1
        write(21,*) "LOOKUP_TABLE default"
        do i = 1, cell_number
            write(21,*) ctr(i)%T
        end do
        close(unit = 21)
    end subroutine output

    !--------------------------------------------------
    !>write result
    !--------------------------------------------------
    subroutine output_2pi()
        character(len=10)  :: outfilename0
        character(len=100) :: outfilename1
        type(node_type),allocatable,dimension(:) :: node_out,node_outz
        type(cell_type),allocatable,dimension(:) :: ctr_out,ctr_outz
        integer :: cell_number_out,cell_number_outz
        integer :: node_number_out,node_number_outz
        integer :: N2pi,Nfold
        integer :: temp_2,temp_1,temp0,temp1,temp_check
        integer :: temp_node,temp_cell,temp_cell0
        integer :: i,j,k

        !>-----------------------------------------------------------------
        !> theta direction
        !>-----------------------------------------------------------------
        Nfold=5
        N2pi=2**5
        !> project to r-theta
        do i=1,node_number
            if(abs(node(i)%coords(1)-0.0d0)<1.0d-15 .and. abs(node(i)%coords(2)-0.0d0)<1.0d-15)then
                node(i)%coords_c(2)=0
            else
                node(i)%coords_c(1)=sqrt(sum(node(i)%coords(1:2)**2.0d0)) !r
                node(i)%coords_c(2)=acos(node(i)%coords(1)/node(i)%coords_c(1)) !theta 0->2pi
                node(i)%coords_c(3)=node(i)%coords(3) !z
            endif
        enddo


        !> count ref_b number
        temp_2=0
        temp_1=0
        temp0 =0
        temp1 =0
        do i=1,node_number
            node(i)%ref_c=0
            node(i)%ref_id=0
            if(abs(node(i)%coords(1)-0.0d0)<1.0d-15 .and. abs(node(i)%coords(2)-0.0d0)<1.0d-15)then
                node(i)%ref_b=-2
                temp_2=temp_2+1
            elseif(abs(node(i)%coords_c(2)-0.0d0)<1.0d-15)then
                node(i)%ref_b=-1
                temp_1=temp_1+1
            elseif(abs(node(i)%coords_c(2)-2.0d0*pi/N2pi)<1.0d-15)then
                node(i)%ref_b=1
                temp1=temp1+1
            else
                node(i)%ref_b=0
                temp0=temp0+1
            endif
        enddo
        if((temp_2+temp_1+temp0+temp1)/=node_number)then
            write(*,*) "error: node out number..."
            pause
        endif

        !> allocate memory
        node_number_out=temp_2+2**Nfold*(temp_1+temp0+temp1)-2**(Nfold-1)*temp1
        cell_number_out=2**Nfold*cell_number
        do i=Nfold-2,0,-1
            node_number_out=node_number_out-2**i*temp_1
        enddo
        node_number_out=node_number_out-temp_1
        allocate(node_out(node_number_out))
        allocate(ctr_out (cell_number_out))

        temp_node=0
        temp_cell=0
        i=0
        temp_check=0
        do j=1,cell_number
            temp_cell=temp_cell+1
            ctr_out(temp_cell)%rho=ctr(j)%rho
            ctr_out(temp_cell)%T=ctr(j)%T
            ctr_out(temp_cell)%node_number=ctr(j)%node_number
            ! node
            do k=1,ctr(j)%node_number
                if(node(ctr(j)%node_id(k))%ref_c==0)then
                    temp_node=temp_node+1
                    node(ctr(j)%node_id(k))%ref_id=temp_node
                    node_out(temp_node)%coords_c(1)=node(ctr(j)%node_id(k))%coords_c(1)
                    node_out(temp_node)%coords_c(2)=node(ctr(j)%node_id(k))%coords_c(2)
                    node_out(temp_node)%coords_c(3)=node(ctr(j)%node_id(k))%coords_c(3)
                    node_out(temp_node)%coords(1)=node(ctr(j)%node_id(k))%coords(1)
                    node_out(temp_node)%coords(2)=node(ctr(j)%node_id(k))%coords(2)
                    node_out(temp_node)%coords(3)=node(ctr(j)%node_id(k))%coords(3)
                    ctr_out(temp_cell)%node_id(k)=temp_node
                    node(ctr(j)%node_id(k))%ref_c=1
                    temp_check=temp_check+1
                    if(isnan(node_out(temp_node)%coords(1)))then
                        write(*,*) "error: node coords is NaN"
                        pause
                    endif
                else
                    ctr_out(temp_cell)%node_id(k)=node(ctr(j)%node_id(k))%ref_id
                endif
            enddo
        enddo
        if(temp_check/=(temp_2+temp_1+temp0+temp1))then
            write(*,*) "error: node reflection..."
            write(*,*) temp_check,(temp_2+temp_1+temp0+temp1)
            pause
        endif

        do i=1,Nfold-1
            !> reset
            temp_check=0
            temp_2=0
            temp_1=0
            temp0=0
            temp1=0
            do j=1,temp_node
                node_out(j)%ref_c=0
                node_out(j)%ref_id=0
                if(abs(node_out(j)%coords(1)-0.0d0)<1.0d-15 .and. abs(node_out(j)%coords(2)-0.0d0)<1.0d-15)then
                    node_out(j)%ref_b=1
                    temp_2=temp_2+1
                elseif(abs(node_out(j)%coords_c(2)-0.0d0)<1.0d-15)then
                    node_out(j)%ref_b=0
                    temp_1=temp_1+1
                elseif(abs(node_out(j)%coords_c(2)-((2.0d0*pi/N2pi)*2**(i-1)))<1.0d-15)then
                    node_out(j)%ref_b=1
                    temp1=temp1+1
                else
                    node_out(j)%ref_b=0
                    temp0=temp0+1
                endif
            enddo
            temp_cell0=temp_cell
            do j=1,temp_cell0
                temp_cell=temp_cell+1
                ctr_out(temp_cell)%rho=ctr_out(j)%rho
                ctr_out(temp_cell)%T=ctr_out(j)%T
                ctr_out(temp_cell)%node_number=ctr_out(j)%node_number
                ! node
                do k=1,ctr_out(j)%node_number
                    if(node_out(ctr_out(j)%node_id(k))%ref_b==1)then
                        ctr_out(temp_cell)%node_id(k)=ctr_out(j)%node_id(k)
                        node_out(ctr_out(j)%node_id(k))%ref_c=1
                        !temp_check=temp_check+1
                    elseif(node_out(ctr_out(j)%node_id(k))%ref_c==0)then
                        temp_node=temp_node+1
                        node_out(ctr_out(j)%node_id(k))%ref_id=temp_node
                        node_out(temp_node)%coords_c(1)=node_out(ctr_out(j)%node_id(k))%coords_c(1)
                        node_out(temp_node)%coords_c(2)=2*((2.0d0*pi/N2pi)*2**(i-1))-node_out(ctr_out(j)%node_id(k))%coords_c(2)
                        node_out(temp_node)%coords_c(3)=node_out(ctr_out(j)%node_id(k))%coords_c(3)
                        node_out(temp_node)%coords(1)=node_out(temp_node)%coords_c(1)*cos(node_out(temp_node)%coords_c(2))
                        node_out(temp_node)%coords(2)=node_out(temp_node)%coords_c(1)*sin(node_out(temp_node)%coords_c(2))
                        node_out(temp_node)%coords(3)=node_out(temp_node)%coords_c(3)
                        ctr_out(temp_cell)%node_id(k)=temp_node
                        node_out(ctr_out(j)%node_id(k))%ref_c=1
                        temp_check=temp_check+1
                    else
                        ctr_out(temp_cell)%node_id(k)=node_out(ctr_out(j)%node_id(k))%ref_id
                    endif
                enddo
            enddo
            if(temp_check/=(temp0+temp_1))then
                write(*,*) "error: node reflection..."
                pause
            endif
        enddo

        !> reset
        i=Nfold
        temp_check=0
        temp_2=0
        temp_1=0
        temp0=0
        temp1=0
        do j=1,temp_node
            node_out(j)%ref_c=0
            node_out(j)%ref_id=0
            if(abs(node_out(j)%coords(1)-0.0d0)<1.0d-15 .and. abs(node_out(j)%coords(2)-0.0d0)<1.0d-15)then
                node_out(j)%ref_b=1
                temp_2=temp_2+1
            elseif(abs(node_out(j)%coords_c(2)-0.0d0)<1.0d-15)then
                node_out(j)%ref_b=1
                temp_1=temp_1+1
            elseif(abs(node_out(j)%coords_c(2)-((2.0d0*pi/N2pi)*2**(i-1)))<1.0d-15)then
                node_out(j)%ref_b=1
                temp1=temp1+1
            else
                node_out(j)%ref_b=0
                temp0=temp0+1
            endif
        enddo
        temp_cell0=temp_cell
        do j=1,temp_cell0
            temp_cell=temp_cell+1
            ctr_out(temp_cell)%rho=ctr_out(j)%rho
            ctr_out(temp_cell)%T=ctr_out(j)%T
            ctr_out(temp_cell)%node_number=ctr_out(j)%node_number
            ! node
            do k=1,ctr_out(j)%node_number
                if(node_out(ctr_out(j)%node_id(k))%ref_b==1)then
                    ctr_out(temp_cell)%node_id(k)=ctr_out(j)%node_id(k)
                    node_out(ctr_out(j)%node_id(k))%ref_c=1
                    !temp_check=temp_check+1
                elseif(node_out(ctr_out(j)%node_id(k))%ref_c==0)then
                    temp_node=temp_node+1
                    node_out(ctr_out(j)%node_id(k))%ref_id=temp_node
                    node_out(temp_node)%coords_c(1)=node_out(ctr_out(j)%node_id(k))%coords_c(1)
                    node_out(temp_node)%coords_c(2)=2*((2.0d0*pi/N2pi)*2**(i-1))-node_out(ctr_out(j)%node_id(k))%coords_c(2)
                    node_out(temp_node)%coords_c(3)=node_out(ctr_out(j)%node_id(k))%coords_c(3)
                    node_out(temp_node)%coords(1)=node_out(temp_node)%coords_c(1)*cos(node_out(temp_node)%coords_c(2))
                    node_out(temp_node)%coords(2)=node_out(temp_node)%coords_c(1)*sin(node_out(temp_node)%coords_c(2))
                    node_out(temp_node)%coords(3)=node_out(temp_node)%coords_c(3)
                    ctr_out(temp_cell)%node_id(k)=temp_node
                    node_out(ctr_out(j)%node_id(k))%ref_c=1
                    temp_check=temp_check+1
                else
                    ctr_out(temp_cell)%node_id(k)=node_out(ctr_out(j)%node_id(k))%ref_id
                endif
            enddo
        enddo
        if(temp_check/=temp0)then
            write(*,*) "error: node reflection..."
            pause
        endif

        !>-----------------------------------------------------------------
        !> z direction
        !>-----------------------------------------------------------------
        !> count ref_b number
        temp_1=0
        temp0 =0
        do i=1,node_number_out
            node_out(i)%ref_c=0
            node_out(i)%ref_id=0
            if(abs(node_out(i)%coords(3)-0.0d0)<1.0d-15)then
                node_out(i)%ref_b=1
                temp_1=temp_1+1
            else
                node_out(i)%ref_b=0
                temp0=temp0+1
            endif
        enddo
        if((temp_1+temp0)/=node_number_out)then
            write(*,*) "error: node out number..."
            pause
        endif

        !> allocate memory
        node_number_outz=temp_1+2*temp0
        cell_number_outz=2*cell_number_out
        allocate(node_outz(node_number_outz))
        allocate(ctr_outz (cell_number_outz))

        temp_node=0
        temp_cell=0
        temp_check=0
        do j=1,cell_number_out
            temp_cell=temp_cell+1
            ctr_outz(temp_cell)%rho=ctr_out(j)%rho
            ctr_outz(temp_cell)%T=ctr_out(j)%T
            ctr_outz(temp_cell)%node_number=ctr_out(j)%node_number
            ! node
            do k=1,ctr_out(j)%node_number
                if(node_out(ctr_out(j)%node_id(k))%ref_c==0)then
                    temp_node=temp_node+1
                    node_out(ctr_out(j)%node_id(k))%ref_id=temp_node
                    node_outz(temp_node)%coords(1)=node_out(ctr_out(j)%node_id(k))%coords(1)
                    node_outz(temp_node)%coords(2)=node_out(ctr_out(j)%node_id(k))%coords(2)
                    node_outz(temp_node)%coords(3)=node_out(ctr_out(j)%node_id(k))%coords(3)
                    ctr_outz(temp_cell)%node_id(k)=temp_node
                    node_out(ctr_out(j)%node_id(k))%ref_c=1
                    temp_check=temp_check+1
                else
                    ctr_outz(temp_cell)%node_id(k)=node_out(ctr_out(j)%node_id(k))%ref_id
                endif
            enddo
        enddo
        if(temp_check/=(temp_1+temp0))then
            write(*,*) "error: node reflection..."
            pause
        endif

        !> reset
        temp_check=0
        do j=1,temp_node
            node_outz(j)%ref_c=0
            node_outz(j)%ref_id=0
            if(node_outz(j)%coords(3)<1.0d-20)then
                node_outz(j)%ref_b=1
            else
                node_outz(j)%ref_b=0
            endif
        enddo
        temp_cell0=temp_cell
        do j=1,temp_cell0
            temp_cell=temp_cell+1
            ctr_outz(temp_cell)%rho=ctr_outz(j)%rho
            ctr_outz(temp_cell)%T=ctr_outz(j)%T
            ctr_outz(temp_cell)%node_number=ctr_outz(j)%node_number
            ! node
            do k=1,ctr_outz(j)%node_number
                if(node_outz(ctr_outz(j)%node_id(k))%ref_b==1)then
                    ctr_outz(temp_cell)%node_id(k)=ctr_outz(j)%node_id(k)
                    node_outz(ctr_outz(j)%node_id(k))%ref_c=1
                elseif(node_outz(ctr_outz(j)%node_id(k))%ref_c==0)then
                    temp_node=temp_node+1
                    node_outz(ctr_outz(j)%node_id(k))%ref_id=temp_node
                    node_outz(temp_node)%coords(1)= node_outz(ctr_outz(j)%node_id(k))%coords(1)
                    node_outz(temp_node)%coords(2)= node_outz(ctr_outz(j)%node_id(k))%coords(2)
                    node_outz(temp_node)%coords(3)=-node_outz(ctr_outz(j)%node_id(k))%coords(3)
                    ctr_outz(temp_cell)%node_id(k)=temp_node
                    node_outz(ctr_out(j)%node_id(k))%ref_c=1
                    temp_check=temp_check+1
                else
                    ctr_outz(temp_cell)%node_id(k)=node_outz(ctr_outz(j)%node_id(k))%ref_id
                endif
            enddo
        enddo
        if(temp_check/=temp0)then
            write(*,*) "error: node reflection..."
            pause
        endif


        !> output
        write(outfilename0,'(I10.10)') iteration
        outfilename1=trim(outputfilename)//trim(outfilename0)//trim('.vtk')
        open(unit = 21, file = outfilename1, status = 'replace')
        write(21,'(a)') "# vtk DataFile Version 2.0"
        write(21,'(a)') "Unstructured Grid Example"
        write(21,'(a)') "ASCII"
        write(21,'(a)') "DATASET UNSTRUCTURED_GRID"
        write(21,*) "POINTS",node_number_outz,"float"
        do i = 1, node_number_outz
            write(21,*) node_outz(i)%coords(1),node_outz(i)%coords(2),node_outz(i)%coords(3)
        end do
        write(21,*) "CELLS",cell_number_outz,cell_number_outz*5
        do i = 1, cell_number_outz
            write(21,*) "4 ",ctr_outz(i)%node_id(1)-1,ctr_outz(i)%node_id(2)-1,ctr_outz(i)%node_id(3)-1,ctr_outz(i)%node_id(4)-1
        end do
        write(21,*) "CELL_TYPES", cell_number_outz
        do i = 1, cell_number_outz
            write(21,*) "10"
        end do

        write(21,*) "CELL_DATA", cell_number_outz
        write(21,*) "SCALARS rho float", 1
        write(21,*) "LOOKUP_TABLE default"
        do i = 1, cell_number_outz
            write(21,*) ctr_outz(i)%rho
        end do

        write(21,*) "SCALARS T float", 1
        write(21,*) "LOOKUP_TABLE default"
        do i = 1, cell_number_outz
            write(21,*) ctr_outz(i)%T
        end do
        close(unit = 21)
    end subroutine output_2pi

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
        write(21,'(a)') 'VARIABLES=X,Y,rho,T,kappa,ParticlePercentage'
        write(21,*) 'ZONE N=',node_number*2,',E=',cell_number_inner*2,',datapacking=block'
        write(21,'(a)') 'varlocation=([3-6]=cellcentered)'
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

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%rho
        enddo
        do i = 1, cell_number_inner
            write(21,*) ctr(i)%rho
        enddo

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%T
        enddo
        do i = 1, cell_number_inner
            write(21,*) ctr(i)%T
        enddo

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%kappa_effective
        enddo
        do i = 1, cell_number_inner
            write(21,*) ctr(i)%kappa_effective
        enddo

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%WaveParticleRatio
        enddo
        do i = 1, cell_number_inner
            write(21,*) ctr(i)%WaveParticleRatio
        enddo

        do i = 1, cell_number_inner
            write(21,*) ctr(i)%node_id(1),ctr(i)%node_id(2),ctr(i)%node_id(3)
        enddo
        do i = 1, cell_number_inner
            write(21,*) ctr(i)%node_id(1)+node_number,ctr(i)%node_id(2)+node_number,ctr(i)%node_id(3)+node_number
        enddo
        close(unit = 21)
    end subroutine output2
end module io
