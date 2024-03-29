!--------------------------------------------------
!>store the global variables
!--------------------------------------------------
module global_data
    !--------------------------------------------------
    !kind selection
    !--------------------------------------------------
    integer,parameter :: short  = 4
    integer,parameter :: double = 8

    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(kind=double),parameter :: PI = 4.0d0*atan(1.0d0) !Pi
    real(kind=double),parameter :: SMV = tiny(real(1.0,8)) !small value to avoid 0/0
    real(kind=double),parameter :: maxnumber=1.0d20
    real(kind=double),parameter :: UP = 1.0d0 !used in sign() function
    real(kind=double) :: cfl !global CFL number
    real(kind=double) :: dt,dt_ex !global time step
    real(kind=double) :: sim_time !current simulation time
    real(kind=double) :: end_time !current simulation time
    integer :: iteration !iteration
    integer :: method_scheme ! 1 wp, 2 sn
    integer :: method_interp !interpolation method
    integer :: method_output !output as cell centered or point value
    integer :: implicit_factor
    character(100) :: inputfilename
    character(100) :: outputfilename

    !--------------------------------------------------
    !probes
    !--------------------------------------------------
    integer(kind=double) :: probe(5)
    real(kind=double) :: probe_x(5),probe_y(5)

    !--------------------------------------------------
    !gas properties
    !--------------------------------------------------
    real(kind=double) :: Kn !Knudsen number
    real(kind=double) :: lightspeed ! speed of light
    real(kind=double) :: radconst ! radiation constant
    real(kind=double) :: ref_weight
    integer :: ref_number

    !--------------------------------------------------
    !macros for a readable code
    !--------------------------------------------------
    !file name
    character(len=20),parameter :: HSTFILENAME   = "photonET.plt" !history file name
    character(len=20),parameter :: RSTFILENAME   = "photonwp.plt" !result file name
    character(len=20),parameter :: PROBEFILENAME = "probe.plt" !result file name
    !direction
    integer,parameter :: IDIRC = 1 !i direction
    integer,parameter :: JDIRC = 2 !j direction
    !rotation
    integer,parameter :: RN = 1 !no frame rotation
    integer,parameter :: RY =-1 !with frame rotation
    !I/O
    integer,parameter :: HSTFILE   = 20 !history file ID
    integer,parameter :: RSTFILE   = 21 !result file ID
    integer,parameter :: PROBEFILE = 22 !history file ID
    !method
    integer,parameter :: UGKWP = 1
    integer,parameter :: SN    = 2
    !interpolation
    integer,parameter :: FIRST_ORDER  = 0 !first order interpolation
    integer,parameter :: SECOND_ORDER = 1 !second order interpolation
    !group
    integer,parameter :: cell01 = -1
    integer,parameter :: cell02 = -2
    integer,parameter :: cell11 = -3
    integer,parameter :: cell12 = -4
    integer,parameter :: cell13 = -5
    integer,parameter :: cell14 = -6
    integer,parameter :: face0  = -7
    integer,parameter :: face1  = -8
    integer,parameter :: node0  = -9
    integer,parameter :: node1  = -10

    !--------------------------------------------------
    !basic derived type
    !--------------------------------------------------
    !grid geometry(node coordinates)
    type :: node_type
        integer :: group
        integer :: cell_number
        integer :: cell_id(1:10)
        real(kind=double) :: coords(1:2)
    end type node_type

    type :: face_type
        integer :: group
        integer :: node_id(1:2)
        integer :: cell_id(1:2)
        !geometry
        real(kind=double) :: coords(1:2) !face center coordinates
        real(kind=double) :: length !length of cell interface
        !normal direction
        real(kind=double) :: norm(1:2)
        real(kind=double) :: weight(1:2) !projection of norm direction
        !field
        real(kind=double) :: rho !density
        real(kind=double) :: T !temperature
        real(kind=double) :: kappa !diffusion coefficient
        real(kind=double) :: absorbtion_scaled
        !flow flux
        real(kind=double) :: flux !mass flux
        real(kind=double) :: flux_particle
        real(kind=double),allocatable,dimension(:,:) :: flux_h !flux of distribution function
        real(kind=double),allocatable,dimension(:,:) :: h !distribution function
        real(kind=double),allocatable,dimension(:,:,:) :: sh !slope of distribution function in i and j direction
    end type face_type

    !cell center
    type :: cell_type
        integer :: group
        integer :: node_id(1:3)
        integer :: face_id(1:3)
        integer :: node_number
        integer :: face_number
        !geometry
        real(kind=double) :: area !cell area
        real(kind=double) :: coords(2) !cell center coordinates
        real(kind=double) :: length(2) !length in i and j direction
        !flow field
        real(kind=double) :: rho               ! density \rho
        real(kind=double) :: rho_wave          ! wave part
        real(kind=double) :: rho_particle      ! particle part
        real(kind=double) :: phi               ! emission energy flux \phi=acT^4
        real(kind=double) :: T                 ! material temperature T
        real(kind=double) :: kappa             ! diffusion coefficient
        real(kind=double) :: kappa_effective   ! effective diffusion coefficient
        real(kind=double) :: capacity          ! heat capacity
        real(kind=double) :: capacity_scaled   ! beta
        real(kind=double) :: absorbtion        ! absorbtion coefficient
        real(kind=double) :: absorbtion_scaled ! nu
        real(kind=double) :: WaveParticleRatio !index
        !particle
        integer :: particle_number
        integer :: particle_number_wave
        !velocity
        real(kind=double),allocatable,dimension(:,:) :: h ! distribution function
        real(kind=double),allocatable,dimension(:,:,:) :: sh ! slope of distribution function in i and j direction
    end type cell_type

    !particle
    type :: particle_type
        integer :: cell
        integer :: cell_new
        integer :: face_in
        integer :: face_in_new
        integer :: flag_track ! 1 need to track
        integer :: flag_delete ! 1 need to delete
        real(kind=double) :: t  !time
        real(kind=double) :: ta !absorbtion time
        real(kind=double) :: tb !boundary time
        real(kind=double) :: weight ! energy flux as h
        real(kind=double) :: x(1:2) !position
        real(kind=double) :: x_new(1:2) !position
        real(kind=double) :: v(1:2) !velocity
    end type particle_type

    !--------------------------------------------------
    !flow field
    !--------------------------------------------------
    integer :: node_number
    integer :: face_number
    integer :: cell_number,cell_number_inner
    integer :: particle_number,particle_number_temp,particle_number_wave
    type(node_type),allocatable,dimension(:) :: node !geometry (node coordinates)
    type(cell_type),allocatable,dimension(:) :: ctr !cell centers
    type(face_type),allocatable,dimension(:) :: face !vertical and horizontal interfaces
    type(particle_type),allocatable,dimension(:) :: particle,particle_temp !particle

    !--------------------------------------------------
    !discrete velocity space
    !--------------------------------------------------
    integer :: unum,vnum !number of velocity points for u and v
    real(kind=double),allocatable,dimension(:,:) :: uspace,vspace !u and v discrete velocity space
    real(kind=double),allocatable,dimension(:,:) :: weight !weight at velocity u_k and v_l
end module global_data
