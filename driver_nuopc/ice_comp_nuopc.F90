module ice_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for CICE
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                  , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                  , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC                  , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC                  , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model            , only : model_routine_SS           => SetServices
  use NUOPC_Model            , only : model_label_Advance        => label_Advance
  use NUOPC_Model            , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model            , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model            , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model            , only : NUOPC_ModelGet, SetVM
  use ice_import_export     , only : ice_advertise_fields, ice_realize_fields
  use ice_import_export     , only : ice_import, ice_export, ice_cpl_dt !, tlast_coupled
  use shr_kind_mod           , only : cl=>shr_kind_cl, cs=>shr_kind_cs
  !r8 slips in through another module
  use shr_file_mod           
  use shr_orb_mod            , only : shr_orb_decl, shr_orb_params, SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT
  use shr_cal_mod            , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use perf_mod               , only : t_startf, t_stopf, t_barrierf
  use nuopc_shr_methods      , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use nuopc_shr_methods      , only : set_component_logging, get_component_instance, log_clock_advance

  use mpas_derived_types
  use mpas_timekeeping
  use mpas_stream_manager
  use mpas_abort
  use mpas_pool_routines
  use mpas_framework
  use mpas_timer
  
  use seaice_column, only : seaice_column_reinitialize_fluxes
  use seaice_forcing, only : post_atmospheric_coupling, post_oceanic_coupling,   &
                             seaice_forcing_get, seaice_forcing_write_restart_times
  use seaice_analysis_driver, only : seaice_analysis_precompute, seaice_analysis_compute, &
                             seaice_analysis_restart, seaice_analysis_write
  use seaice_initialize, only : seaice_init_post_clock_advance
  use seaice_core_interface, only : seaice_setup_core, seaice_setup_domain
  use seaice_time_integration, only : seaice_timestep
  !use mpas_seaice_constants, only : coupleAlarmID

!$ use OMP_LIB               , only : omp_set_num_threads
  implicit none

  public  :: SetServices
  public  :: SetVM

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize

  character(len=CL) :: flds_scalar_name = ''
  integer           :: flds_scalar_num = 0
  integer           :: flds_scalar_index_nx = 0
  integer           :: flds_scalar_index_ny = 0
  integer           :: flds_scalar_index_nextsw_cday = 0

  integer     , parameter :: dbug = 0
  character(*), parameter :: modName =  "(ice_comp_nuopc)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

  character(len=CL)      :: orb_mode        ! attribute - orbital mode
  integer                :: orb_iyear       ! attribute - orbital year
  integer                :: orb_iyear_align ! attribute - associated with model year
  real(R8)               :: orb_obliq       ! attribute - obliquity in degrees
  real(R8)               :: orb_mvelp       ! attribute - moving vernal equinox longitude
  real(R8)               :: orb_eccen       ! attribute and update-  orbital eccentricity

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

  integer                 :: nthrds   ! Number of threads to use in this component

  type (core_type), pointer :: corelist => null()
  type (dm_info), pointer :: dminfo
  type (domain_type), pointer :: domain_ptr
  type (iosystem_desc_t), pointer :: io_system 

  integer, private ::   &
      my_task, lmpicom, iam
  integer, parameter  :: stdout = 6
   integer :: iceLogUnit ! unit number for ice log
   character(len=StrKIND) :: runtype, coupleTimeStamp
  character (len=*), parameter :: coupleAlarmID = 'coupling'
   integer :: itimestep

!=======================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !--------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local variables
    character(len=CL)  :: cvalue
    character(len=CL)  :: logmsg
    logical            :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write(iceLogUnit,*) subname//'Need to set attribute ScalarFieldName'
       call mpas_dmpar_abort(domain_ptr % dminfo)
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write(iceLogUnit,*) subname//'Need to set attribute ScalarFieldCount'
       call mpas_dmpar_abort(domain_ptr % dminfo)
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write(iceLogUnit,*) subname//'Need to set attribute ScalarFieldIdxGridNX'
       call mpas_dmpar_abort(domain_ptr % dminfo)
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write(iceLogUnit,*) subname//'Need to set attribute ScalarFieldIdxGridNY'
       call mpas_dmpar_abort(domain_ptr % dminfo)
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write(iceLogUnit,*) subname//'Need to set attribute ScalarFieldIdxNextSwCday'
       call mpas_dmpar_abort(domain_ptr % dminfo)
    endif

    call ice_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    use ESMF               , only: ESMF_VMGet
      use mpas_stream_manager, only : MPAS_stream_mgr_init, MPAS_build_stream_filename, MPAS_stream_mgr_validate_streams
      use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
      use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string
      use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_time
      use mpas_bootstrapping, only : mpas_bootstrap_framework_phase1, mpas_bootstrap_framework_phase2
      use mpas_log

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local variables
    integer                 :: lmpicom
    integer                 :: localPet
    integer                 :: npes
    type(ESMF_DistGrid)     :: ice_distgrid
    type(ESMF_Mesh)         :: ice_mesh
    type(ESMF_VM)           :: vm
    type(ESMF_DistGrid)     :: distGrid
    type(ESMF_Mesh)         :: Emesh
    type(ESMF_Time)         :: startTime          ! Start time
    type(ESMF_Time)         :: stopTime           ! Stop time
    type(ESMF_Time)         :: refTime            ! Ref time
    type(ESMF_TimeInterval) :: timeStep           ! Model timestep
    type(ESMF_CalKind_Flag) :: esmf_caltype       ! esmf calendar type
    integer                 :: start_ymd          ! Start date (YYYYMMDD)
    integer                 :: start_tod          ! start time of day (s)
    integer                 :: curr_ymd           ! Current date (YYYYMMDD)
    integer                 :: curr_tod           ! Current time of day (s)
    integer                 :: stop_ymd           ! stop date (YYYYMMDD)
    integer                 :: stop_tod           ! stop time of day (sec)
    integer                 :: ref_ymd            ! Reference date (YYYYMMDD)
    integer                 :: ref_tod            ! reference time of day (s)
    integer                 :: yy,mm,dd           ! Temporaries for time query
    integer                 :: iyear              ! yyyy
    integer                 :: dtime              ! time step
    character(len=cs)       :: starttype          ! infodata start type
    logical                 :: isPresent
    character(CL)           :: cvalue
    character(CL)           :: single_column_lnd_domainfile
    character(len=CL)       :: ice_meshfile
    character(len=CL)       :: ice_maskfile
    type(ESMF_Field)        :: lfield
    character(CL) ,pointer  :: lfieldnamelist(:) => null()
    integer                 :: fieldcount
    real(r8), pointer       :: fldptr1d(:)
    real(r8), pointer       :: fldptr2d(:,:)
    integer                 :: n
    real(r8)                :: scol_lon
    real(r8)                :: scol_lat
    real(r8)                :: scol_spval
    integer                 :: rank
    integer :: icePossibleErrUnit !< unit number to reserve for if a err log needs to be opened
    logical :: readNamelistArg, readStreamsArg
    integer :: shrloglev, shrlogunit, ierr
    type (MPAS_Time_Type) :: alarmStartTime
    type (MPAS_TimeInterval_Type) :: alarmTimeStep
    type(MPAS_Time_Type)         :: currTime           ! Current time
    type (block_type), pointer :: block
    character(len=strKIND) :: &
         calendar_type
    integer :: mesh_iotype, ierr_local
    character(len=strKIND) :: caseid
    character(len=strKIND) :: model_version
    character(len=strKIND) :: username
    character(len=strKIND) :: hostname
    character(len=strKIND) :: institution
    character(len=strKIND) :: curdate
    character(len=strKIND) :: curtime
    character(len=strKIND) :: history
    logical :: streamsExists
    character(kind=c_char), dimension(StrKIND+1) :: c_filename       ! StrKIND+1 for C null-termination character
    integer(kind=c_int) :: c_comm
    integer(kind=c_int) :: c_ierr
    character(len=StrKIND) :: mesh_stream
    character(len=StrKIND) :: mesh_filename
    character(len=StrKIND) :: mesh_filename_temp
    character(len=StrKIND) :: ref_time_temp
    character(len=StrKIND) :: filename_interval_temp
    character(kind=c_char), dimension(StrKIND+1) :: c_mesh_stream
    character(kind=c_char), dimension(StrKIND+1) :: c_mesh_filename_temp
    character(kind=c_char), dimension(StrKIND+1) :: c_ref_time_temp
    character(kind=c_char), dimension(StrKIND+1) :: c_filename_interval_temp
    character(kind=c_char), dimension(StrKIND+1) :: c_iotype
    character(len=8) :: c_inst_index ! instance number
    character(len=8) :: c_npes       ! number of pes
    character(len=StrKIND) :: iotype
    type (c_ptr) :: mgr_p
    type(ESMF_Time)         :: EcurrTime
    type(ESMF_Calendar)     :: Ecalendar
    logical                 :: associated_e
    integer(I8Kind)         :: s_e, sn_e, sd_e
    integer                 :: yy_e
    type(ESMF_CalKind_Flag)   :: type_e
    type (MPAS_Time_type) :: start_time
    type (MPAS_Time_type) :: ref_time
    type (MPAS_TimeInterval_type) :: filename_interval
    type (MPAS_TimeInterval_type) :: denInterval, remInterval, zeroInterval
    real(kind=RKIND), pointer :: &
         dayOfNextShortwaveCalculation ! needed for CESM like coupled simulations
    type (MPAS_Pool_Type), pointer :: shortwave
    character(len=StrKIND)  :: timeStamp
    integer (kind=I8KIND) :: numDivs
    type (mpas_pool_type), pointer :: meshpool
    integer ncells, lsize, num_ice, icell
    integer, dimension(:), pointer :: nCellsArray
    integer, dimension(:), pointer :: indexToCellID
    integer , allocatable   :: gindex_ice(:)
    integer , allocatable   :: gindex(:)
    logical                 :: mastertask
    integer                 :: errorCode       ! error code
    logical, pointer :: tempLogicalConfig
    character(len=StrKIND), pointer :: tempCharConfig
    real (kind=RKIND), pointer :: tempRealConfig
    character(*), parameter     :: F00   = "('(ice_comp_nuopc) ',2a,1x,d21.14)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '

    interface
       subroutine xml_stream_parser(xmlname, mgr_p, comm, ierr) bind(c)
          use iso_c_binding, only : c_char, c_ptr, c_int
          character(kind=c_char), dimension(*), intent(in) :: xmlname
          type (c_ptr), intent(inout) :: mgr_p
          integer(kind=c_int), intent(inout) :: comm
          integer(kind=c_int), intent(out) :: ierr
       end subroutine xml_stream_parser

       subroutine xml_stream_get_attributes(xmlname, streamname, comm, filename, ref_time, filename_interval, io_type, ierr) bind(c)
          use iso_c_binding, only : c_char, c_int
          character(kind=c_char), dimension(*), intent(in) :: xmlname
          character(kind=c_char), dimension(*), intent(in) :: streamname
          integer(kind=c_int), intent(inout) :: comm
          character(kind=c_char), dimension(*), intent(out) :: filename
          character(kind=c_char), dimension(*), intent(out) :: ref_time
          character(kind=c_char), dimension(*), intent(out) :: filename_interval
          character(kind=c_char), dimension(*), intent(out) :: io_type
          integer(kind=c_int), intent(out) :: ierr
       end subroutine xml_stream_get_attributes
    end interface

    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)


    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=iam, PetCount=npes, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, pet=iam, peCount=nthrds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(nthrds==1) then
       call NUOPC_CompAttributeGet(gcomp, "nthreads", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       read(cvalue,*) nthrds
    endif
!$  call omp_set_num_threads(nthrds)

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) caseid

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype


      readNamelistArg = .false.
      readStreamsArg = .false.
!-----------------------------------------------------------------------
!
!   setup mpassi data structures
!
!-----------------------------------------------------------------------
    allocate(corelist)
    nullify(corelist % next)

    allocate(corelist % domainlist)
    nullify(corelist % domainlist % next)

    domain_ptr => corelist % domainlist
    domain_ptr % core => corelist

    call mpas_allocate_domain(domain_ptr)

!-----------------------------------------------------------------------
!
!   first initializaiton phase of mpassi
!   call mpassi initialization routines
!
!-----------------------------------------------------------------------

    call t_startf('mpassi_init_total')
    call t_startf('mpassi_init')

    !io_system => shr_pio_getiosys(iceid)

    !pio_iotype = shr_pio_getiotype(iceid)
!   ! call MPAS_io_set_iotype(domain_ptr % iocontext, pio_iotype)
    !call MPAS_io_unset_iotype(domain_ptr % iocontext)

    ! ----------------
    ! Set up log file information
    ! ----------------
    !inst_suffix = seq_comm_suffix(iceID) ! a suffix to append to log file name
    iceLogUnit = shr_file_getUnit() ! reserve unit number for log unit
    icePossibleErrUnit = shr_file_getUnit() ! reserve unit number for possible error log file

    ! Store shr log unit and level so we can reassign them later when ICE control is complete
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    ! Set the shr log unit to be the iceLogUnit
    !  (I believe this is so that shr routines called from the ICE driver write any messages to the ICE log file.)
    call shr_file_setLogUnit(iceLogUnit)

    ! Write to iceLogUnit here, because the log module is not open yet.
    if (iam==0) write(iceLogUnit,'(a,i6)') '=== Beginning ice_init_nuopc: rank=',iam

    call mpas_framework_init_phase1(domain_ptr % dminfo, lmpicom)
 
    ! Setup function pointers for MPASSI core and domain types
    call seaice_setup_core(corelist)
    call seaice_setup_domain(domain_ptr)

    ! ===========
    ! Initialize log manager
    call mpas_log_init(domain_ptr % logInfo, domain_ptr, unitNumbers=(/iceLogUnit, icePossibleErrUnit/), err=ierr)
    if ( ierr /= 0 ) then
       write(iceLogUnit,*) 'ERROR: log init failed for core ' // trim(domain_ptr % core % coreName)
       call mpas_dmpar_abort(domain_ptr % dminfo)
    end if

    ! Set core specific options here
    ! Disable output from all but the master task for E3SM!
    ! (This overrides the default set by mpas_log_init based on MPAS_DEBUG setting.)
    if (iam /= 0) then
       domain_ptr % logInfo % outputLog % isActive = .false.
    endif

    ! After core has had a chance to modify log defaults, open the output log
    call mpas_log_open(err=ierr)
    if ( ierr /= 0 ) then
       write(iceLogUnit,*) 'ERROR: log open failed for core ' // trim(domain_ptr % core % coreName)
       call mpas_dmpar_abort(domain_ptr % dminfo)
    end if
    ! ===========


    ! ----------
    ! Process namelist and streams files
    ! ----------
    ! Override the names of the stream and namelist files
    domain_ptr % namelist_filename = 'mpassi_in'
    domain_ptr % streams_filename = 'streams.seaice'

    ! Setup namelist variables, and read the namelist
    ierr = domain_ptr % core % setup_namelist(domain_ptr % configs, domain_ptr % namelist_filename, domain_ptr % dminfo)
    if ( ierr /= 0 ) then
       call mpas_log_write('Namelist setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    call mpas_framework_init_phase2(domain_ptr)!, io_system)

    ! Define package variables
    ierr = domain_ptr % core % define_packages(domain_ptr % packages)
    if ( ierr /= 0 ) then
       call mpas_log_write('Package definition failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Setup packages (i.e. determine if they should be on or off)
    ierr = domain_ptr % core % setup_packages(domain_ptr % configs, domain_ptr % packages, domain_ptr % iocontext)
    if ( ierr /= 0 ) then
       call mpas_log_write('Package setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Setup decompositions available for dimensions
    ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
    if ( ierr /= 0 ) then
       call mpas_log_write('Decomposition setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if


    ! Determine runtype and possibly nextsw_cday
    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       read(cvalue,*) starttype
       if (trim(starttype) == trim('startup')) then
          runtype = "initial"
       else if (trim(starttype) == trim('continue') ) then
          runtype = "continue"
       else if (trim(starttype) == trim('branch')) then
          runtype = "continue"
       else
          write(iceLogUnit,*) subname//' ERROR: unknown starttype'
          call mpas_dmpar_abort(domain_ptr % dminfo)
       end if

       ! Note that in the mct version the atm was initialized first so that nextsw_cday could be passed to the other
       ! components - this assumed that cam or datm was ALWAYS initialized first.
       ! In the nuopc version it will be easier to assume that on startup - nextsw_cday is just the current time

       ! TOOD (mvertens, 2019-03-21): need to get the perpetual run working

       if (trim(runtype) == 'initial') then
         ! Turn off restart
         call mpas_pool_get_config(domain_ptr % configs, "config_do_restart", tempLogicalConfig)
         tempLogicalConfig = .false.

         ! Setup start time. Will be over written later when clocks are synchronized
         call mpas_pool_get_config(domain_ptr % configs, "config_start_time", tempCharConfig)
         tempCharConfig = trim(tempCharConfig) // "0:00:00"

         ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
         call mpas_pool_get_config(domain_ptr % configs, "config_run_duration", tempCharConfig)
         tempCharConfig = "0001-00-00_00:00:00"
       else
         ! Turn on restarts
         call mpas_pool_get_config(domain_ptr % configs, "config_do_restart", tempLogicalConfig)
         tempLogicalConfig = .true.
         call mpas_pool_get_config(domain_ptr % configs, "config_do_restart_hbrine", tempLogicalConfig)
         tempLogicalConfig = .true.
         call mpas_pool_get_config(domain_ptr % configs, "config_do_restart_snow_density", tempLogicalConfig)
         tempLogicalConfig = .true.
         call mpas_pool_get_config(domain_ptr % configs, "config_do_restart_snow_grain_radius", tempLogicalConfig)
         tempLogicalConfig = .true.
         call mpas_pool_get_config(domain_ptr % configs, "config_do_restart_zsalinity", tempLogicalConfig)
         tempLogicalConfig = .true.

         ! Set start time to be read from file
         call mpas_pool_get_config(domain_ptr % configs, "config_start_time", tempCharConfig)
         tempCharConfig = "file"

         ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
         call mpas_pool_get_config(domain_ptr % configs, "config_run_duration", tempCharConfig)
         tempCharConfig = "0001-00-00_00:00:00"
          call ESMF_ClockGet( clock, currTime=EcurrTime, rc=rc )
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          !DDcall ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
          !DDif (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    !DD! Determine runid
    !DDcall NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, isPresent=isPresent, rc=rc)
    !DDif (ChkErr(rc,__LINE__,u_FILE_u)) return
    !DDif (isPresent) then
    !DD   read(cvalue,*) runid
    !DDelse
    !DD   runid = 'unknown'  ! read in from the namelist in ice_init.F90 if CESMCOUPLED is not defined
    !DDend if

    ! Get clock information before call to cice_init

    call ESMF_ClockGet( clock, &
         currTime=EcurrTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( EcurrTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    !DDcall ESMF_TimeIntervalGet( timeStep, s=dtime, rc=rc )
    !DDif (ChkErr(rc,__LINE__,u_FILE_u)) return
    !DDdt = real(dtime)

    call ESMF_TimeGet( EcurrTime, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar_type = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar_type = shr_cal_gregorian
    else
       write(iceLogUnit,*) subname//'ERROR:: bad calendar for ESMF'
       call mpas_dmpar_abort(domain_ptr % dminfo)
    end if

    ! Setup MPASSI simulation clock
    ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
    if ( ierr /= 0 ) then
    end if

    call mpas_log_write('Reading streams configuration from file '//trim(domain_ptr % streams_filename))
    inquire(file=trim(domain_ptr % streams_filename), exist=streamsExists)

    if ( .not. streamsExists ) then
       call mpas_log_write('Streams file '//trim(domain_ptr % streams_filename)//' does not exist.', MPAS_LOG_CRIT)
    end if

    !
    ! Using information from the namelist, a graph.info file, and a file containing
    !    mesh fields, build halos and allocate blocks in the domain
    !
    ierr = domain_ptr % core % get_mesh_stream(domain_ptr % configs, mesh_stream)
    if ( ierr /= 0 ) then
       call mpas_log_write('Failed to find mesh stream for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if


    call mpas_f_to_c_string(domain_ptr % streams_filename, c_filename)
    call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
    c_comm = domain_ptr % dminfo % comm
    call xml_stream_get_attributes(c_filename, c_mesh_stream, c_comm, &
                                   c_mesh_filename_temp, c_ref_time_temp, &
                                   c_filename_interval_temp, c_iotype, c_ierr)
    if (c_ierr /= 0) then
       call mpas_log_write('xml_stream_get_attributes failed.', MPAS_LOG_CRIT)
    end if
    call mpas_c_to_f_string(c_mesh_filename_temp, mesh_filename_temp)
    call mpas_c_to_f_string(c_ref_time_temp, ref_time_temp)
    call mpas_c_to_f_string(c_filename_interval_temp, filename_interval_temp)
    call mpas_c_to_f_string(c_iotype, iotype)

    if (trim(iotype) == 'pnetcdf') then
       mesh_iotype = MPAS_IO_PNETCDF
    else if (trim(iotype) == 'pnetcdf,cdf5') then
       mesh_iotype = MPAS_IO_PNETCDF5
    else if (trim(iotype) == 'netcdf') then
       mesh_iotype = MPAS_IO_NETCDF
    else if (trim(iotype) == 'netcdf4') then
       mesh_iotype = MPAS_IO_NETCDF4
    else
       mesh_iotype = MPAS_IO_PNETCDF
    end if

    start_time = mpas_get_clock_time(domain_ptr % clock, MPAS_START_TIME, ierr)
    if ( trim(ref_time_temp) == 'initial_time' ) then
        call mpas_get_time(start_time, dateTimeString=ref_time_temp, ierr=ierr)
    end if

    if ( trim(filename_interval_temp) == 'none' ) then
        call mpas_expand_string(ref_time_temp, -1, mesh_filename_temp, mesh_filename)
    else
        call mpas_set_time(ref_time, dateTimeString=ref_time_temp, ierr=ierr)
        call mpas_set_timeInterval(filename_interval, timeString=filename_interval_temp, ierr=ierr)
        call mpas_build_stream_filename(ref_time, start_time, filename_interval, mesh_filename_temp, -1, mesh_filename, ierr)
    end if

    ! Bootstrap framework (1). Here data structures are setup, but dimensions and arrays are not finalized.
    call mpas_log_write(' ** Attempting to bootstrap MPAS framework using stream: ' // trim(mesh_stream))
    call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype)

    !
    ! Set up run-time streams
    !
    call MPAS_stream_mgr_init(domain_ptr % streamManager, domain_ptr % ioContext, domain_ptr % clock, domain_ptr % blocklist % allFields, &
                              domain_ptr % packages, domain_ptr % blocklist % allStructs)

    call datetime(curdate, curtime)
    history = 'created on ' // trim(curdate) // ' ' // trim(curtime)

    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'case', trim(caseid))
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'source_id', trim(model_version))
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'realm', 'seaIce')
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'product', 'model-output')
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'title', 'MPAS-Seaice History file information')
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'source', 'E3SM Sea Ice Model')
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'institution',trim(institution))
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'institution_id', 'E3SM-Project')
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'contact', 'e3sm-data-support@listserv.llnl.gov')
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'username', trim(username))
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'hostname', trim(hostname))

    ! Other attributes
    call add_stream_attributes(domain_ptr)

    ! overwrite MPAS attributes
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'history', trim(history))
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'Conventions', 'CF-1.7')
    call MPAS_stream_mgr_add_att(domain_ptr % streamManager, 'git_version', trim(model_version))

    ! Setup all immutable streams for the core
    ierr = domain_ptr % core % setup_immutable_streams(domain_ptr % streamManager)
    if ( ierr /= 0 ) then
       call mpas_log_write('Immutable streams setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Parse / read all streams configuration
    mgr_p = c_loc(domain_ptr % streamManager)
    call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
    if (c_ierr /= 0) then
       call mpas_log_write('xml_stream_parser failed.', MPAS_LOG_CRIT)
    end if

    my_task = domain_ptr % dminfo % my_proc_id
    
    !
    ! Finalize the setup of blocks and fields
    !
    call mpas_bootstrap_framework_phase2(domain_ptr)

    !DD! Determine coupling type
    !DD What do I do here?
    !DDcall seq_infodata_GetData(infodata, cpl_seq_option=cpl_seq_option)

    ! Determine time of next atmospheric shortwave calculation
    block => domain_ptr % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_array(shortwave, "dayOfNextShortwaveCalculation", dayOfNextShortwaveCalculation)
       !DD What do I do here?
       !DDcall seq_infodata_GetData(infodata, nextsw_cday=dayOfNextShortwaveCalculation )

       ! Set dayOfNextShortwaveCalculation to -1 for continue and branch runs.
       if (trim(runtype) /= 'initial') dayOfNextShortwaveCalculation = -1
       call mpas_log_write('dayOfNextShortwaveCalculation = $r', realArgs=(/dayOfNextShortwaveCalculation/))

       block => block % next
    end do


    ! Initialize the MPASSI core
    ierr = domain_ptr % core % core_init(domain_ptr, timeStamp)
    if ( ierr /= 0 ) then
       call mpas_log_write('Core init failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

!-----------------------------------------------------------------------
!
!   initialize time-stamp information
!
!-----------------------------------------------------------------------
    call t_stopf ('mpassi_init')

!-----------------------------------------------------------------------
!
!   check for consistency of mpassi and sync clock initial time
!
!-----------------------------------------------------------------------

    call ESMF_ClockGet( clock, currTime=EcurrTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    !DD there hase to be a better way to go from esmf type to MPAS_Time_Type
    call ESMF_TimeGet(Ecurrtime, s_i8=s_e, sn_i8=sn_e, sd_i8=sd_e, yy=yy_e, calendar = ecalendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !currtime%t%basetime%S  = s_e
    !currtime%t%basetime%Sn = sn_e
    !currtime%t%basetime%Sd = sd_e
    !currtime%t%yr = yy_e
    currtime%t = EcurrTime

!    call ESMF_CalendarGet(Ecalendar, calkindflag=type_e, rc=rc)
!    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!    if(type_e == ESMF_CALKIND_NOLEAP) then
!       currTime%t%calendar => noleapCal
!    elseif(type_e == ESMF_CALKIND_GREGORIAN) then
!       currTime%t%calendar => gregorianCal
!    else
!       rc = 1
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    endif

    if (runtype == 'initial') then
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_NOW, ierr)
    else if (runtype == 'continue' .or. runtype == 'branch') then
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_NOW, ierr)
    end if
    ! Doublecheck that clocks are synced here.  
    ! (Note that if this is an initial run, a section below will advance the ocean clock
    ! by one coupling interval, at which point we expect the clocks to be OUT of sync.)
    !DD add the check later
    !?!?call check_clocks_sync(domain_ptr % clock, Eclock, ierr)
    
    !-----------------------------------------------------------------------
    ! initialize necessary coupling info
    !-----------------------------------------------------------------------

    call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet( timeStep, s=ice_cpl_dt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    call mpas_set_timeInterval(alarmTimeStep, S=ice_cpl_dt, ierr=ierr)
    call mpas_get_timeInterval(alarmTimeStep, timeString=coupleTimeStamp, ierr=ierr)

    ! Verify the mpas time step fits into a coupling interval
    call mpas_pool_get_config(domain_ptr % configs, 'config_dt', tempCharConfig)
    call mpas_set_timeInterval(denInterval, timeString=tempCharConfig, ierr=ierr)
    call mpas_set_timeInterval(zeroInterval, S=0, ierr=ierr)
    call mpas_interval_division(start_time, alarmTimeStep, denInterval, numDivs, remInterval)

    ierr = 0

    if ( alarmTimeStep < denInterval ) then
       ierr = 1
    end if
    ierr_local = ierr
    call mpas_dmpar_max_int(domain_ptr % dminfo, ierr_local, ierr)

    if ( ierr == 1 ) then
       call mpas_log_write('Coupling interval is: ' // trim(coupleTimeStamp), MPAS_LOG_ERR)
       call mpas_log_write('        ICE Model time step is: ' // trim(tempCharConfig), MPAS_LOG_ERR)
       call mpas_log_write('        The model time step cannot be longer then the coupling interval', MPAS_LOG_ERR)
       call mpas_log_write('Model is not properly configured for coupling interval.', MPAS_LOG_CRIT)
    end if

    if ( remInterval > zeroInterval ) then
       ierr = 1
    end if

    ierr_local = ierr
    call mpas_dmpar_max_int(domain_ptr % dminfo, ierr_local, ierr)

    if ( ierr == 1 ) then
       call mpas_log_write('Coupling interval is: ' // trim(coupleTimeStamp), MPAS_LOG_ERR)
       call mpas_log_write('        ICE Model time step is: ' // trim(tempCharConfig), MPAS_LOG_ERR)
       call mpas_log_write('        These are not synchronized, so time steps will not match to coupling interval boundaries.', MPAS_LOG_ERR)
       call mpas_log_write('        Please reconfigure either the coupling interval or the time step.', MPAS_LOG_ERR)
       call mpas_log_write('Model is not properly configured for coupling interval.', MPAS_LOG_CRIT)
    end if

    ! set coupling alarm
    alarmStartTime = currTime
    call mpas_add_clock_alarm(domain_ptr % clock, coupleAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
    call mpas_print_alarm(domain_ptr % clock, coupleAlarmID, ierr)
    call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)
    
    !---------------------------------------------------------------------------
    ! Determine the global index space needed for the distgrid
    !---------------------------------------------------------------------------

    n = 0
    block => domain_ptr % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray( 1 )
       n = n + ncells
       block => block % next
    end do
    lsize = n
    allocate(gindex_ice(lsize))

    n = 0
    block => domain_ptr % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
       nCells = nCellsArray( 1 )
       do iCell = 1, nCells
          gindex_ice(n+iCell) =indexToCellID(iCell)
       enddo
       n = n + ncells
       block => block % next
    end do

!       ! No eliminated land blocks
       num_ice = size(gindex_ice)
       allocate(gindex(num_ice))
       do n = 1,num_ice
          gindex(n) = gindex_ice(n)
       end do

    !---------------------------------------------------------------------------
    ! Create distGrid from global index array
    !---------------------------------------------------------------------------
    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Create the MPAS-SI mesh
    !---------------------------------------------------------------------------

    ! read in the mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ice', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    EMesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistgrid=Distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !mastertask = iam == domain_ptr % dminfo % my_proc_id
    mastertask = iam == 0
    if (mastertask) then
       write(stdout,*)'mesh file for mpassi domain is ',trim(cvalue)
    end if


    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call ice_realize_fields(importState, exportState, mesh=Emesh, &
         flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, &
         my_task=iam, mastertask=mastertask, lmpicom=lmpicom,  domain=domain_ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Prescribed ice initialization
    !-----------------------------------------------------------------

    !DD will assume this is false for now and will implement if needed
    !DDcall ice_prescribed_init(clock, ice_mesh, rc)
    !DDif (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------
    ! Create cice export state
    !-----------------------------------------------------------------

    call ice_export (exportState, flds_scalar_name, domain_ptr,    &
                     errorCode, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !DDcall State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, &
    !DD     flds_scalar_name, flds_scalar_num, rc)
    !DDif (ChkErr(rc,__LINE__,u_FILE_u)) return
    !DDcall State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, &
    !DD     flds_scalar_name, flds_scalar_num, rc)
    !DDif (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! TODO (mvertens, 2018-12-21): fill in iceberg_prognostic as .false.

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    call t_stopf ('mpassi_init_total')

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !---------------------------------------------------------------------------
    ! Run MPASSI
    !---------------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State)             :: importState
    type(ESMF_State)             :: exportState
    type(ESMF_Clock)             :: clock
    type(ESMF_Alarm)             :: alarm
    type(ESMF_Time)              :: nextTime
    character(CL)          :: cvalue
    integer                :: shrlogunit ! original log unit
    integer                :: k,n        ! index
    logical                :: stop_now   ! .true. ==> stop at the end of this run phase
    integer                :: ymd        ! Current date (YYYYMMDD)
    integer                :: tod        ! Current time of day (sec)
    integer                :: curr_ymd   ! Current date (YYYYMMDD)
    integer                :: curr_tod   ! Current time of day (s)
    integer                :: yy,mm,dd   ! year, month, day, time of day
    integer                :: ymd_sync   ! Sync date (YYYYMMDD)
    integer                :: yr_sync    ! Sync current year
    integer                :: mon_sync   ! Sync current month
    integer                :: day_sync   ! Sync current day
    integer                :: tod_sync   ! Sync current time of day (sec)
      integer :: ihour, iminute, isecond
      integer :: iyear, imonth, iday
    character(CL)          :: restart_filename
    character(*)   , parameter :: F00   = "('(ice_comp_nuopc) ',2a,i8,d21.14)"
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    logical                      :: debugon
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      integer SHRLOGLEV
    integer                      :: errorCode, ierr ! error flag
    type (block_type), pointer :: block
    type (MPAS_Pool_Type), pointer :: shortwave
    real(kind=RKIND), pointer :: &
         dayOfNextShortwaveCalculation ! needed for CESM like coupled simulations
      type (MPAS_Time_Type) :: currTime
      type (MPAS_timeInterval_type) :: timeStep
      logical, save :: first=.true.
      character(len=StrKIND) :: timeStamp, streamName, WCstring
      real (kind=RKIND) :: dt, current_wallclock_time
      integer  :: streamDirection
      logical :: streamActive
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

      iam = domain_ptr % dminfo % my_proc_id

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Note this logic triggers off of the component clock rather than the internal pop time
    ! The component clock does not get advanced until the end of the loop - not at the beginning

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

      debugOn = .false.
#ifdef MPAS_DEBUG
      debugOn = .true.
#endif
    !--------------------------------
    ! Single column logic if nearest neighbor point has a mask of zero
    !--------------------------------

    !if (single_column .and. .not. scol_valid) then
    !   RETURN
    !end if

!$  call omp_set_num_threads(nthrds)

      ! Set MPAS Log module instance
      mpas_log_info => domain_ptr % logInfo
      if (debugOn) call mpas_log_write("=== Beginning ice_run_nuopc ===")

      call mpas_pool_get_config(domain_ptr % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Setup log information.
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (iceLogUnit)

      ! reinitialize fluxes
      call seaice_column_reinitialize_fluxes(domain_ptr)


    call ice_import(importState, flds_scalar_name, domain_ptr,    &
                         errorCode, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! Post coupling calls
      block => domain_ptr % blocklist
      do while (associated(block))

         ! Determine time of next atmospheric shortwave calculation
         call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
         call MPAS_pool_get_array(shortwave, "dayOfNextShortwaveCalculation", dayOfNextShortwaveCalculation)
         !DD what do I do heere?
         !DDcall seq_infodata_GetData(infodata, nextsw_cday=dayOfNextShortwaveCalculation )

         ! perform post coupling operations
         call post_atmospheric_coupling(block)
         call post_oceanic_coupling(block)

         block => block % next
      end do

      ! reset coupler alarm before we start
      call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)

     ! Get current time
     currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)

     timeStep = mpas_get_clock_timestep(domain_ptr % clock, ierr=ierr)
     call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)

     itimestep = 0  ! We may want to initialize this in init and make it a module, save variable

     ! During integration, time level 1 stores the model state at the beginning of the
     !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
     do while (.not. mpas_is_alarm_ringing(domain_ptr % clock, coupleAlarmID, ierr=ierr))
        call mpas_stream_mgr_read(domain_ptr % streamManager, ierr=ierr)
        call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

        itimestep = itimestep + 1
        call mpas_advance_clock(domain_ptr % clock)

        ! final initialization after clock advance 
        if (first) then
           call seaice_forcing_get(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .true.)

           call seaice_init_post_clock_advance(domain_ptr, domain_ptr % clock)
           first = .false.
        else
           ! forcing
           call seaice_forcing_get(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .false.)
        endif

        currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
        call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
        ! write time stamp on every step
        call mpas_dmpar_get_time(current_wallclock_time)
        write (WCstring,'(F18.3)') current_wallclock_time
        call mpas_log_write(trim(timeStamp) // '  WC time:' // WCstring)

        !DD assumed false for now
        !DD! get prescribed ice coverage
        !DDcall ice_prescribed_run(domain_ptr, currTime)

        ! pre-timestep analysis computation
        if (debugOn) call mpas_log_write('   Starting analysis precompute', masterOnly=.true.)
        call seaice_analysis_precompute(domain_ptr, ierr)
        if (debugOn) call mpas_log_write('   Finished analysis precompute', masterOnly=.true.)

        if (debugOn) call mpas_log_write('   Starting forward update', masterOnly=.true.)
        call mpas_timer_start("time integration", .false.)
        ierr = 0
        call seaice_timestep(domain_ptr, domain_ptr % clock, itimestep)
        call mpas_timer_stop("time integration")
        if (debugOn) call mpas_log_write('   Finished forward update', masterOnly=.true.)

        ! update analysis members
        if (debugOn) call mpas_log_write('   Starting AM compute', masterOnly=.true.)
        call seaice_analysis_compute(domain_ptr, ierr)
        if (debugOn) call mpas_log_write('   Finished AM compute', masterOnly=.true.)
        if (debugOn) call mpas_log_write('   Starting AM restart', masterOnly=.true.)
        call seaice_analysis_restart(domain_ptr, ierr)
        if (debugOn) call mpas_log_write('   Finished AM restart', masterOnly=.true.)
        if (debugOn) call mpas_log_write('   Starting AM write', masterOnly=.true.)
        call seaice_analysis_write(domain_ptr, ierr)
        if (debugOn) call mpas_log_write('   Finished AM write', masterOnly=.true.)

        ! Reset the restart alarm to prevent restart files being written without the coupler requesting it.
        if (debugOn) call mpas_log_write('   Resetting restart stream alarms', masterOnly=.true.)
        call mpas_stream_mgr_begin_iteration(domain_ptr % streamManager)
        do while ( mpas_stream_mgr_get_next_stream(domain_ptr % streamManager, streamID=streamName, &
                   directionProperty=streamDirection, activeProperty=streamActive) )
           if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
              call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, streamID=streamName, ierr=ierr)
           end if
        end do
        if (debugOn) call mpas_log_write('   Resetting restart stream alarms complete', masterOnly=.true.)
        if (debugOn) call mpas_log_write('   Writing output streams', masterOnly=.true.)
        call mpas_stream_mgr_write(domain_ptr % streamManager, streamID='output', ierr=ierr)
        call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, streamID='output', ierr=ierr)

        call mpas_stream_mgr_write(domain_ptr % streamManager, ierr=ierr)
        call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
        if (debugOn) call mpas_log_write('   Finished writing output streams', masterOnly=.true.)
        if (debugOn) call mpas_log_write('Completed timestep ' // trim(timeStamp))
     end do

      ! Check if coupler wants us to write a restart file.
      ! We only write restart files at the end of a coupling interval
      
      if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         call ESMF_AlarmRingerOff( alarm, rc=rc )
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         ! Write a restart file, because the coupler asked for it.
         if (debugOn) call mpas_log_write('Writing restart streams', masterOnly=.true.)
         call seaice_forcing_write_restart_times(domain_ptr)
         call mpas_stream_mgr_begin_iteration(domain_ptr % streamManager)
         do while ( mpas_stream_mgr_get_next_stream(domain_ptr % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive) )
            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               if (debugOn) call mpas_log_write('    Writing stream: ' // trim(streamName), masterOnly=.true.)
               call mpas_stream_mgr_write(domain_ptr % streamManager, forceWriteNow=.true., streamID=streamName, iErr=iErr)
               if (debugOn) call mpas_log_write('    Finished Writing stream: ' // trim(streamName), masterOnly=.true.)
            end if
         end do

         if ( iam == 0 ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) trim(timeStamp)
            close(22)
         end if
         if (debugOn) call mpas_log_write('Finished writing restart streams', masterOnly=.true.)
      end if

    !--------------------------------
    ! Create export state
    !--------------------------------


    if (debugOn) call mpas_log_write('Exporting state', masterOnly=.true.)
    call ice_export(exportState, flds_scalar_name, domain_ptr,    &
                         errorCode, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (debugOn) call mpas_log_write('Finished exporting state', masterOnly=.true.)


      ! Check if clocks are in sync
!TODO?      call check_clocks_sync(domain_ptr % clock, Eclock, ierr)
      ! Check if clocks are in sync
      currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, YYYY=iyear, MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond, ierr=ierr)
      !DDcall seq_timemgr_EClockGetData(EClock, curr_ymd=curr_ymd, curr_tod=curr_tod)

      !DDymd = iyear * 10000 + imonth * 100 + iday
      !DDtod = ihour * 3600 + iminute * 60 + isecond
      !DDif (.not. seq_timemgr_EClockDateInSync( EClock, ymd, tod)) then
      !DD   call mpas_log_write('MPAS ymd=$i MPAS tod=$i', MPAS_LOG_ERR, intArgs=(/ymd,tod/))
      !DD   call mpas_log_write('sync ymd=$i sync tod=$i', MPAS_LOG_ERR, intArgs=(/curr_ymd, curr_tod/))
      !DD   call mpas_log_write('Internal mpas clock not in sync with sync clock', MPAS_LOG_CRIT)
      !DDend if

      call mpas_log_write('=== Completed coupling interval in ice_run_nuopc. ===', flushNow=.true.)

      ! Reset I/O logs
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

  105  format( A, 2i8, A, f10.2, A, f10.2, A)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! intput/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(ice_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(ice_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !--------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (my_task == 0) then
       write(stdout,F91)
       write(stdout,F00) 'MPASSI: end of main integration loop'
       write(stdout,F91)
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================

!***********************************************************************
!BOP
!
! !IROUTINE: datetime
!
! !INTERFACE:
   subroutine datetime(cdate, ctime)
!
! !DESCRIPTION:
! Calculate current date and time for metadata history 
!
! !USES:

! !INPUT PARAMETERS:
      character(len=*), intent(out) :: cdate
      character(len=*), intent(out) :: ctime

! !OUTPUT PARAMETERS:
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
      character(len=strKIND) :: date
      character(len=strKIND) :: time
      character(len=strKIND) :: zone
      integer, dimension(8) :: values

      call date_and_time (date, time, zone, values)

      cdate = ""
      cdate(1:2) = date(5:6)
      cdate(3:3) = '/'
      cdate(4:5) = date(7:8)
      cdate(6:6) = '/'
      cdate(7:8) = date(3:4)

      ctime = ""
      ctime(1:2) = time(1:2)
      ctime(3:3) = ':'
      ctime(4:5) = time(3:4)
      ctime(6:6) = ':'
      ctime(7:8) = time(5:6)

!-----------------------------------------------------------------------
!EOC

   end subroutine datetime!}}}

   subroutine add_stream_attributes(domain)!{{{

      type (domain_type), intent(inout) :: domain

      type (MPAS_Pool_iterator_type) :: itr
      integer, pointer :: intAtt
      logical, pointer :: logAtt
      character (len=StrKIND), pointer :: charAtt
      real (kind=RKIND), pointer :: realAtt
      character (len=StrKIND) :: histAtt

      integer :: local_ierr

      if (domain % dminfo % nProcs < 10) then
          write(histAtt, '(A,I1,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100) then
          write(histAtt, '(A,I2,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 1000) then
          write(histAtt, '(A,I3,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 10000) then
          write(histAtt, '(A,I4,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100000) then
          write(histAtt, '(A,I5,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else
          write(histAtt, '(A,I6,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      end if
     

      call MPAS_stream_mgr_add_att(domain % streamManager, 'on_a_sphere', domain % on_a_sphere)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'model_name', domain % core % modelName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'core_name', domain % core % coreName)
      ! DWJ 10/01/2014: Eventually add the real history attribute, for now (due to length restrictions)
      ! add a shortened version.
!     call MPAS_stream_mgr_add_att(domain % streamManager, 'history', domain % history)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'history', histAtt)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'Conventions', domain % core % Conventions)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'parent_id', domain % parent_id)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'mesh_spec', domain % mesh_spec)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'git_version', domain % core % git_version)

      call mpas_pool_begin_iteration(domain % configs)

      do while (mpas_pool_get_next_member(domain % configs, itr))

         if ( itr % memberType == MPAS_POOL_CONFIG) then

            if ( itr % dataType == MPAS_POOL_REAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, realAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, realAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_INTEGER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, intAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, intAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, charAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, charAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, logAtt)
               if (logAtt) then
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'YES', ierr=local_ierr)
               else
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'NO', ierr=local_ierr)
               end if
            end if

          end if
      end do

   end subroutine add_stream_attributes!}}}

!***********************************************************************
end module ice_comp_nuopc
