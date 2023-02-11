!--------------------------------------------------
!>main program
!--------------------------------------------------
program main
    use global_data
    use macroscopic
    use microscopic
    use microscopic_sn
    use io
    implicit none
    integer :: i

    !initialization
    call init()
    call output2()

    !open file and write header
    open(unit=HSTFILE, file = HSTFILENAME, status = "replace", action="write") 
    write(HSTFILE,*) "VARIABLES = sim_time, T1, E1, T2, E2, T3, E3, T4, E4, T5, E5" !write header

    !iteration
    do while(.true.)
        if(method_scheme==UGKWP)then
            !----------------------------------------------------------------------
            !> UGKWP method
            !----------------------------------------------------------------------
            dt=dt/implicit_factor
            do i=1,implicit_factor
                !> stream particle: free stream flux
                call update_particle()

                !> resample particle: microscopic closure
                call reinit_particle()
            enddo

            dt=dt*implicit_factor
            !> source iteration: diffusive flux
            call source_iteration()
        elseif(method_scheme==SN)then
            !----------------------------------------------------------------------
            !> SN method
            !----------------------------------------------------------------------
            !> calculate flux
            call evolution()

            !> update
            call update_flux_sn()

            !> source iteration
            call source_iteration()

            !> source term
            call update_source_sn()
        endif

        !----------------------------------------------------------------------
        !> output
        !----------------------------------------------------------------------
        if (mod(iteration,10)==0) then
            write(*,"(A18,I15,2E15.7,I15,I15)") "iter,sim_time,dt:",iteration,sim_time,dt,particle_number,implicit_factor
            !write(HSTFILE,"(I15,10E15.7)") sim_time,ctr(probe(1))%matT,ctr(probe(1))%radT,ctr(probe(2))%matT,ctr(probe(2))%radT,&
            !ctr(probe(3))%matT,ctr(probe(3))%radT,ctr(probe(4))%matT,ctr(probe(4))%radT,ctr(probe(5))%matT,ctr(probe(5))%radT
        end if
        !if (mod(iteration,100)==0) then
        !    call output()
        !end if
        if (sim_time>=0.2 .and. sim_time<=0.2+dt) then
            call output()
        end if
        if (sim_time>=0.4 .and. sim_time<=0.4+dt) then
            call output()
        end if
        if (sim_time>=0.6 .and. sim_time<=0.6+dt) then
            call output()
        end if
        if (sim_time>=0.8 .and. sim_time<=0.8+dt) then
            call output()
        end if
        if (sim_time>=1.0 .and. sim_time<=1.0+dt) then
            call output()
            exit
        end if

        !----------------------------------------------------------------------
        !> simulation time
        !----------------------------------------------------------------------
        iteration = iteration+1
        sim_time  = sim_time+dt
        if (sim_time >= end_time) exit
    end do

    !close history file
    close(HSTFILE)

    !output solution
    call output()
end program main
