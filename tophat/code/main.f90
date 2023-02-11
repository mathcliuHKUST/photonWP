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
        implicit_factor=1 !min(floor(sim_time/200*50)+1,50)
        dt=implicit_factor*dt_ex
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
            dt=dt/implicit_factor
            do i=1,implicit_factor
                !> calculate flux
                call evolution()

                !> update
                call update_flux_sn()

                !> source term
                call update_source_sn()
            enddo


            dt=dt*implicit_factor
            !> source iteration: diffusive flux
            call source_iteration()
        endif

        !----------------------------------------------------------------------
        !> output
        !----------------------------------------------------------------------
        if (mod(iteration,10)==0) then
            write(*,"(A18,I15,2E15.7,I15,I15)") "iter,sim_time,dt:",iteration,sim_time,dt,particle_number,implicit_factor
            write(HSTFILE,"(11E15.7)") sim_time,&
            ctr(probe(1))%T,(max(ctr(probe(1))%rho/radconst/lightspeed,1.0d-10))**0.25d0,&
            ctr(probe(2))%T,(max(ctr(probe(2))%rho/radconst/lightspeed,1.0d-10))**0.25d0,&
            ctr(probe(3))%T,(max(ctr(probe(3))%rho/radconst/lightspeed,1.0d-10))**0.25d0,&
            ctr(probe(4))%T,(max(ctr(probe(4))%rho/radconst/lightspeed,1.0d-10))**0.25d0,&
            ctr(probe(5))%T,(max(ctr(probe(5))%rho/radconst/lightspeed,1.0d-10))**0.25d0
        end if
        if(sim_time<1.5d0)then
            if (mod(iteration,200)==0) then
                call output2()
            end if
        else
            if (mod(iteration,200)==0) then
                call output2()
            end if
        endif

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
    call output2()
end program main
