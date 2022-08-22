module ice_import_export

  use ESMF
  use NUOPC
  use NUOPC_Model
  use shr_sys_mod        , only : shr_sys_abort, shr_sys_flush
  use shr_kind_mod       , only : cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod,       only : shr_const_spval, radius=>SHR_CONST_REARTH, pi=>shr_const_pi
  use shr_mpi_mod,         only : shr_mpi_min, shr_mpi_max
  use perf_mod           , only : t_startf, t_stopf, t_barrierf
  use nuopc_shr_methods,     only: chkerr

   ! MPAS framework modules
   use mpas_framework
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_stream_manager
   !use mpas_kind_types, only : r8 => r8kind
   use mpas_io_units
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_log

   use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
   use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string

   ! MPASSI modules
   use seaice_analysis_driver
   use seaice_column, only: seaice_column_reinitialize_fluxes, &
                            seaice_column_coupling_prep
   use ice_kinds_mod, only : dbl_kind
   use seaice_constants
   use seaice_core
   use seaice_core_interface
   use seaice_forcing
   use seaice_initialize, only: seaice_init_post_clock_advance
   use seaice_mesh, only: seaice_latlon_vector_rotation_backward
   use seaice_time_integration
   use seaice_error, only: seaice_check_critical_error

   use ice_constants_colpkg, only : Tffresh, ice_ref_salinity, p001
   use ice_colpkg, only: colpkg_sea_freezing_temperature

  implicit none
  public

  public  :: ice_advertise_fields
  public  :: ice_realize_fields
  public  :: ice_import
  public  :: ice_export

  integer, public :: ice_cpl_dt    ! length of coupling interval in seconds - set by coupler/ESMF

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_FldChk

  interface state_getfldptr
     module procedure state_getfldptr_1d
     module procedure state_getfldptr_2d
  end interface state_getfldptr
  private :: state_getfldptr

  ! Private module data

  type fld_list_type
    character(len=128) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  ! area correction factors for fluxes send and received from mediator
  real(r8), allocatable :: mod2med_areacor(:) ! ratios of model areas to input mesh areas
  real(r8), allocatable :: med2mod_areacor(:) ! ratios of input mesh areas to model areas

  integer, parameter       :: fldsMax = 100
  integer                  :: fldsToIce_num = 0
  integer                  :: fldsFrIce_num = 0
  type (fld_list_type)     :: fldsToIce(fldsMax)
  type (fld_list_type)     :: fldsFrIce(fldsMax)
  logical                  :: atm_prognostic

  integer, parameter  :: stdout = 6
  integer     , parameter :: dbug = 1        ! i/o debug messages
  real(R8)    , parameter  :: czero = 0.0_R8
  character(*), parameter  :: u_FILE_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine ice_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    integer       :: n
    character(CS) :: stdname
    character(CS) :: cvalue
    logical       :: flds_i2o_per_cat  ! .true. => select per ice thickness category
    character(len=*), parameter :: subname='(ice_import_export:ice_advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !DDcall NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
    !DDif (ChkErr(rc,__LINE__,u_FILE_u)) return
    !DDread(cvalue,*) send_i2x_per_cat
    !DDcall ESMF_LogWrite('flds_i2o_per_cat = '// trim(cvalue), ESMF_LOGMSG_INFO)

    !-----------------
    ! advertise import fields
    !-----------------

    call fldlist_add(fldsToIce_num, fldsToIce, trim(flds_scalar_name))

    ! from ocean
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_t'    ) !seaSurfaceTemperature
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_s'    ) !seaSurfaceSalinity
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_u'    ) !uOceanVelocity
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_v'    ) !vOceanVelocity
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_dhdx' ) !seaSurfaceTiltU
    call fldlist_add(fldsToIce_num, fldsToIce, 'So_dhdy' ) !seaSurfaceTiltV
    call fldlist_add(fldsToIce_num, fldsToIce, 'Fioo_q'  ) !freezingMeltingPotential
    call fldlist_add(fldsToIce_num, fldsToIce, 'Fioo_frazil') !frazilMassFlux

    ! from atmosphere
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_z'       ) !airLevelHeight
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_ptem'    ) !airPotentialTemperature
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_tbot'    ) !airTemperature
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_shum'    ) !airSpecificHumidity
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_dens'    ) !airDensity
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swvdr' ) !shortwaveVisibleDirectDown
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swvdf' ) !shortwaveVisibleDiffuseDown
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swndr' ) !shortwaveIRDirectDown
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_swndf' ) !shortwaveIRDiffuseDown
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_lwdn'  ) !longwaveDown
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_rain'  ) !rainfallRate
    call fldlist_add(fldsToIce_num, fldsToIce, 'Faxa_snow'  ) !snowfallRate
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_u'       ) !uAirVelocity
    call fldlist_add(fldsToIce_num, fldsToIce, 'Sa_v'       ) !vAirVelocity

    !DD! set aerosols, if configured
    !DDif (config_use_aerosols .or. config_use_column_biogeochemistry.and.config_use_zaerosols) then
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????' )       !bcphodry
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !bcphidry
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !bcphiwet
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !dstwet1
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !dstwet2
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !dstwet3
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !dstwet4
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !dstdry1
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !dstdry2
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !dstdry3
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !dstdry4.
    !DDendif

    !DD! import biogeochemistry fields, if configured
    !DDif (config_use_column_biogeochemistry) then
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_algae1
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_algae2
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_algae3
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_doc1
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_doc2
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_dic1
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_don1
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_no3
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_sio3
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_nh4
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_dms
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_dmsp
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_docr
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_fep1
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_fep2
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_fed1
    !DD   call fldlist_add(fldsToIce_num, fldsToIce, '????'       ) !So_fed2
    !DD   endif
    !DDendif

    do n = 1,fldsToIce_num
       call NUOPC_Advertise(importState, standardName=fldsToIce(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !-----------------
    ! advertise export fields
    !-----------------

    call NUOPC_CompAttributeGet(gcomp, name='ATM_model', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(cvalue) == 'satm' .or. trim(cvalue) == 'datm') then
       atm_prognostic = .false.
    else
       atm_prognostic = .true.
    end if

    !DDcall fldlist_add(fldsFrIce_num, fldsFrIce, trim(flds_scalar_name))

    ! ice states
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_imask'            )  
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_ifrac'            )  !ailohi
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_t') !surfaceTemperatureCell
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_bpress'         ) !basalPressure
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_vsno'        ) !volume of snow per unit area
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_u10'                 ) !atmosReferenceSpeed10m
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_tref'                 ) !atmosReferenceTemperature2m
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_snowh'                ) !Si_snowh
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_qref'                  ) !atmosReferenceHumidity2m
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_swnet'  ) !absorbedShortwaveFlux
    !if (config_use_column_shortwave) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_avsdr' ) !albedoVisibleDirectCell
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_anidr'  ) !albedoIRDirectCell
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_avsdf' ) !albedoVisibleDiffuseCell
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_anidf'  ) !albedoIRDiffuseCell
    !endif
    !DDif (send_i2x_per_cat) then
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, 'Si_ifrac_n', &
    !DD        ungridded_lbound=1, ungridded_ubound=ncat)  !ice_fraction_n
    !DDend if

    !DDif (config_use_data_icebergs) then
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !bergFreshwaterFlux
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???'  ) !bergLatentHeatFlux
    !DDendif

    ! ice/atm fluxes computed by ice
    if (atm_prognostic) then
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_taux') !tauxa
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_tauy') !tauya
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_lat' ) !latentHeatFlux
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_sen' ) !sensibleHeatFlux
       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_lwup') !longwaveUp
!       call fldlist_add(fldsFrIce_num, fldsFrIce, 'Faii_evap') !evaporativeWaterFlux
    end if

    ! ice/ocn fluxes computed by ice
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_melth'  ) !oceanHeatFlux
    call fldlist_add(fldsFrIce_num, fldsFrIce, 'Fioi_swpen'  ) !oceanShortwaveFlux
    !DDif (send_i2x_per_cat) then
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, 'mean_sw_pen_to_ocn_ifrac_n', &
    !DD        ungridded_lbound=1, ungridded_ubound=ncat)
    !DDend if
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_meltw' ) !oceanFreshWaterFlux
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_salt'  ) !oceanSaltFlux
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_taux'  ) !tauxo
    call fldlist_add(fldsFrIce_num , fldsFrIce, 'Fioi_tauy'  ) !tauyo

    !DD! export biogeochemistry fields, if configured
    !DDif (config_use_column_biogeochemistry) then
    !DD   ! convert from mmol N/m^3 to mmol C/m^3
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_algae1
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_algae2
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_algae3
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_doc1
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_doc2
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_doc3
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_dic1
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_don1
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_no3
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_sio3
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_nh4
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_dms
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_dmspp
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_dmspd
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_docr
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_fep1
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_fep2
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_fed1
    !DD   call fldlist_add(fldsFrIce_num, fldsFrIce, '???' ) !Fioi_fed2
    !DDendif

    do n = 1,fldsFrIce_num
       call NUOPC_Advertise(exportState, standardName=fldsFrIce(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine ice_advertise_fields

!==============================================================================

  subroutine ice_realize_fields(importState, exportState,  mesh, grid,          &
                                flds_scalar_name, flds_scalar_num,              &
                                my_task, mastertask, lmpicom,  domain, rc)

    use shr_mpi_mod  , only : shr_mpi_min, shr_mpi_max

    ! input/output variables
    type(ESMF_State)      :: importState
    type(ESMF_State)      :: exportState
    type(ESMF_Mesh) , optional , intent(in)  :: mesh
    type(ESMF_Grid) , optional , intent(in)  :: grid
    character(len=*)           , intent(in)  :: flds_scalar_name
    integer                    , intent(in)  :: flds_scalar_num
    integer          , intent(in)  :: lmpicom  ! the ocean mpi communicator
    logical          , intent(in)  :: mastertask           
    integer          , intent(in)  :: my_task           
    type (domain_type), pointer, intent(in)  :: domain
    integer                    , intent(out) :: rc

    ! local variables
    type(ESMF_Field)      :: lfield
    integer               :: numOwnedElements
    integer               :: spatialDim
    integer               :: i, j, iblk, n
    integer               :: ilo, ihi, jlo, jhi ! beginning and end of physical domain
    type (block_type), pointer :: block
    real(r8), allocatable :: mesh_areas(:)
    real(r8), allocatable :: model_areas(:)
    real(r8), pointer     :: dataptr(:)
    real(r8)              :: max_mod2med_areacor
    real(r8)              :: max_med2mod_areacor
    real(r8)              :: min_mod2med_areacor
    real(r8)              :: min_med2mod_areacor
    real(r8)              :: max_mod2med_areacor_glob
    real(r8)              :: max_med2mod_areacor_glob
    real(r8)              :: min_mod2med_areacor_glob
    real(r8)              :: min_med2mod_areacor_glob
    character(len=*), parameter :: subname='(ice_import_export:realize_fields)'

    real(r8), pointer     :: ownedElemCoords(:)
    real(r8), pointer     :: latModel(:), latMesh(:)
    real(r8), pointer     :: lonModel(:), lonMesh(:)
    real(r8)              :: diff_lon
    real(r8)              :: diff_lat
    integer :: iCell, nCells
    integer, dimension(:), pointer :: nCellsArray
    real (kind=RKIND), dimension(:), pointer :: lonCell,                       &
                                                latCell,                       &
                                                areaCell 

    type (mpas_pool_type), pointer :: meshPool
    integer gcell, cell_offset
    real(r8) :: rad2deg = 180._r8/pi

    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrIce, &
         numflds=fldsFrIce_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':MPASSI_Export',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToIce, &
         numflds=fldsToIce_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':MPASSI_Import',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine mesh lats and lons
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    allocate(lonMesh(numOwnedElements))
    allocate(latMesh(numOwnedElements))
    allocate(lonModel(numOwnedElements))
    allocate(latModel(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do
    deallocate(ownedElemCoords)

    ! Compare mesh lats/lons to model generated lats/lons
    cell_offset = 0
    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray( 1 )
       call mpas_pool_get_array(meshPool, 'lonCell', lonModel)  
       call mpas_pool_get_array(meshPool, 'latCell', latModel)  
       
       do iCell = 1, nCells
          gcell = iCell + cell_offset

          diff_lon = abs(lonMesh(gcell) - lonModel(iCell)*rad2deg)
          if ( (diff_lon > 1.e2  .and. abs(diff_lon - 360.) > 1.e-1) .or.&
               (diff_lon > 1.e-3 .and. diff_lon < 1._r8) ) then
             write(stdout,'(a,i6,2(f21.13,3x),d21.5)') &
                 'ERROR: MPASSI n, lonMesh, lonModel, diff_lon = ',&
                 gcell,lonMesh(gcell),lonModel(icell)*rad2deg, diff_lon
             call shr_sys_abort()
          end if
          if (abs(latMesh(gcell) - latModel(icell)*rad2deg) > 1.e-1) then
             write(stdout,'(a,i6,2(f21.13,3x),d21.5)') &
                  'ERROR: MPASSI n, latMesh, latModel, diff_lat = ', &
                  gcell,latMesh(gcell),latModel(icell)*rad2deg, abs(latMesh(gcell)-latModel(iCell)*rad2deg)
             call shr_sys_abort()
          end if
       end do

       cell_offset = cell_offset + nCells
       
       block => block % next
    end do

    ! Determine mesh areas used in regridding
    lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(mesh_areas(numOwnedElements))
    mesh_areas(:) = dataptr(:)
    call ESMF_FieldDestroy(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine flux correction factors (module variables)
    !DD is this needed?
    allocate(model_areas(numOwnedElements))
    allocate(mod2med_areacor(numOwnedElements))
    allocate(med2mod_areacor(numOwnedElements))
    mod2med_areacor(:) = 1._r8
    med2mod_areacor(:) = 1._r8
    cell_offset = 0
    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray( 1 )
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)  
       
       do iCell = 1, nCells
          gcell = iCell + cell_offset

          model_areas(gcell) = areaCell(icell)/(radius*radius)   
          mod2med_areacor(gcell) = model_areas(gcell) / mesh_areas(gcell)
          med2mod_areacor(gcell) = mesh_areas(gcell) / model_areas(gcell)
       end do

       cell_offset = cell_offset + nCells
       
       block => block % next
    end do

    min_mod2med_areacor = minval(mod2med_areacor)
    max_mod2med_areacor = maxval(mod2med_areacor)
    min_med2mod_areacor = minval(med2mod_areacor)
    max_med2mod_areacor = maxval(med2mod_areacor)
    call shr_mpi_max(max_mod2med_areacor, max_mod2med_areacor_glob, lmpicom)
    call shr_mpi_min(min_mod2med_areacor, min_mod2med_areacor_glob, lmpicom)
    call shr_mpi_max(max_med2mod_areacor, max_med2mod_areacor_glob, lmpicom)
    call shr_mpi_min(min_med2mod_areacor, min_med2mod_areacor_glob, lmpicom)

    if (my_task == mastertask) then
       write(stdout,'(2A,2g23.15,A )') trim(subname),' :  min_mod2med_areacor, max_mod2med_areacor ',&
            min_mod2med_areacor_glob, max_mod2med_areacor_glob, 'MPASSI'
       write(stdout,'(2A,2g23.15,A )') trim(subname),' :  min_med2mod_areacor, max_med2mod_areacor ',&
            min_med2mod_areacor_glob, max_med2mod_areacor_glob, 'MPASSI'
    end if
    
    deallocate(mesh_areas, model_areas)

  end subroutine ice_realize_fields

  !==============================================================================

  subroutine ice_import( importState, flds_scalar_name, domain,    &
                         errorCode, rc )

    ! input/output variables
    type(ESMF_State) , intent(in)  :: importState
    character(len=*)   , intent(in)  :: flds_scalar_name
    type (domain_type), pointer, intent(in) :: domain
    integer            , intent(out) :: errorCode
    integer            , intent(out) :: rc

    ! local variables
    character (cl) :: label,  message
    integer              :: i,j,k,n,ncol,iblock,nfld,nf
    real (r8)            :: m2percm2, gsum
    real (r8), pointer   :: dataptr1d(:)
    real (r8), pointer   :: dataptr2d(:,:)
    character (cl), allocatable :: fieldNameList(:)
    character(len=*),   parameter    :: subname = 'ice_import'

    ! from atm
    real (r8), pointer   :: Sa_z(:)
    real (r8), pointer   :: Sa_ptem(:)
    real (r8), pointer   :: Sa_tbot(:)
    real (r8), pointer   :: Sa_shum(:)
    real (r8), pointer   :: Sa_dens(:)
    real (r8), pointer   :: Sa_u(:)
    real (r8), pointer   :: Sa_v(:)
    real (r8), pointer   :: Faxa_rain(:)
    real (r8), pointer   :: Faxa_snow(:)
    real (r8), pointer   :: Faxa_lwdn(:)
    real (r8), pointer   :: Faxa_swvdr(:)
    real (r8), pointer   :: Faxa_swvdf(:)
    real (r8), pointer   :: Faxa_swndr(:)
    real (r8), pointer   :: Faxa_swndf(:)
    
    ! from ocean
    real (r8), pointer   :: So_t(:)
    real (r8), pointer   :: So_s(:)
    real (r8), pointer   :: So_u(:)
    real (r8), pointer   :: So_v(:)
    real (r8), pointer   :: So_dhdx(:)
    real (r8), pointer   :: So_dhdy(:)
    real (r8), pointer   :: Fioo_q(:)
    real (r8), pointer   :: Fioo_frazil(:)

    integer              :: fieldCount
    character (cl) :: fldname
    type(ESMF_StateItem_Flag) :: itemflag

    ! local variables - mpas names
    integer :: iCell
    integer, pointer :: nCellsSolve
    type (mpas_pool_type), pointer :: meshPool, forcingPool, configs, &
                                      atmosCoupling, oceanCoupling

    type (block_type), pointer :: block
    integer gcell, cell_offset

    real (kind=RKIND), dimension(:), pointer :: airLevelHeight, airPotentialTemperature, &
                                                evaporativeWaterFlux, airTemperature, &
                                                airSpecificHumidity, airDensity, &
                                                shortwaveVisibleDirectDown, shortwaveVisibleDiffuseDown, &
                                                shortwaveIRDirectDown, shortwaveIRDiffuseDown, &
                                                rainfallRate, snowfallRate, &
                                                uAirVelocity, vAirVelocity, &
                                                seaSurfaceTemperature, seaSurfaceSalinity, &
                                                seaFreezingTemperature, freezingMeltingPotential, &
                                                frazilMassAdjust, &
                                                uOceanVelocity, vOceanVelocity, &
                                                seaSurfaceTiltU, seaSurfaceTiltV, &
                                                longwaveDown
    type (field1DReal),         pointer :: airLevelHeightField, airPotentialTemperatureField, &
                                           evaporativeWaterFluxField, airTemperatureField, &
                                           airSpecificHumidityField, airDensityField, &
                                           shortwaveVisibleDirectDownField, shortwaveVisibleDiffuseDownField, &
                                           shortwaveIRDirectDownField, shortwaveIRDiffuseDownField, &
                                           rainfallRateField, snowfallRateField, &
                                           uAirVelocityField, vAirVelocityField, &
                                           seaSurfaceTemperatureField, seaSurfaceSalinityField, &
                                           seaFreezingTemperatureField, freezingMeltingPotentialField, &
                                           frazilMassAdjustField, &
                                           uOceanVelocityField, vOceanVelocityField, &
                                           seaSurfaceTiltUField, seaSurfaceTiltVField
                                                                                        
   logical, pointer :: &
      config_use_aerosols,              &
      config_use_modal_aerosols,        &
      config_use_zaerosols,             &
      config_use_column_biogeochemistry

   character(len=strKIND), pointer ::   &
      config_thermodynamics_type,       &
      config_ocean_surface_type
      
   real(r8) ::  frazilMassFlux, frazilMassFluxRev

    !-----------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call t_startf ('mpassi_imp_atm')

    !-----------------------------------------------------------------------
    ! from atmosphere
    !-----------------------------------------------------------------------
    call state_getfldptr(importState, 'Faxa_lwdn', Faxa_lwdn, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Faxa_rain', Faxa_rain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Faxa_snow', Faxa_snow, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_z', Sa_z, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_ptem', Sa_ptem, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_tbot', Sa_tbot, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_shum', Sa_shum, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_dens', Sa_dens, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Faxa_swvdr', Faxa_swvdr, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Faxa_swvdf', Faxa_swvdf, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Faxa_swndr', Faxa_swndr, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Faxa_swndf', Faxa_swndf, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_u', Sa_u, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Sa_v', Sa_v, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------------
    ! from ocean
    !-----------------------------------------------------------------------
    call state_getfldptr(importState, 'So_t', So_t, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_s', So_s, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_u', So_u, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_v', So_v, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_dhdx', So_dhdx, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'So_dhdy', So_dhdy, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Fioo_q', Fioo_q, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Fioo_frazil', Fioo_frazil, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    !-----------------------------------------------------------------------
    ! copy into component variables
    !-----------------------------------------------------------------------
    cell_offset = 0
    block => domain % blocklist
    do while (associated(block))

       configs => block % configs
       call mpas_pool_get_config(configs, "config_thermodynamics_type", config_thermodynamics_type)
       call mpas_pool_get_config(configs, "config_ocean_surface_type", config_ocean_surface_type)
       call mpas_pool_get_config(configs, "config_use_aerosols", config_use_aerosols)
       call mpas_pool_get_config(configs, "config_use_modal_aerosols", config_use_modal_aerosols)
       call mpas_pool_get_config(configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'ocean_coupling', oceanCoupling)
       call mpas_pool_get_subpool(block % structs, 'atmos_coupling', atmosCoupling)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_array(oceanCoupling, 'seaSurfaceTemperature', seaSurfaceTemperature)
       call mpas_pool_get_array(oceanCoupling, 'seaSurfaceSalinity', seaSurfaceSalinity)
       call mpas_pool_get_array(oceanCoupling, 'seaFreezingTemperature', seaFreezingTemperature)
       call mpas_pool_get_array(oceanCoupling, 'freezingMeltingPotential', freezingMeltingPotential)
       call mpas_pool_get_array(oceanCoupling, 'frazilMassAdjust', frazilMassAdjust)
       call mpas_pool_get_array(oceanCoupling, 'uOceanVelocity', uOceanVelocity)
       call mpas_pool_get_array(oceanCoupling, 'vOceanVelocity', vOceanVelocity)
       call mpas_pool_get_array(oceanCoupling, 'seaSurfaceTiltU', seaSurfaceTiltU)
       call mpas_pool_get_array(oceanCoupling, 'seaSurfaceTiltV', seaSurfaceTiltV)

       call mpas_pool_get_array(atmosCoupling, 'airLevelHeight', airLevelHeight)
       call mpas_pool_get_array(atmosCoupling, 'airPotentialTemperature', airPotentialTemperature)
       call mpas_pool_get_array(atmosCoupling, 'airTemperature', airTemperature)
       call mpas_pool_get_array(atmosCoupling, 'airSpecificHumidity', airSpecificHumidity)
       call mpas_pool_get_array(atmosCoupling, 'airDensity', airDensity)
       call mpas_pool_get_array(atmosCoupling, 'shortwaveVisibleDirectDown', shortwaveVisibleDirectDown)
       call mpas_pool_get_array(atmosCoupling, 'shortwaveVisibleDiffuseDown', shortwaveVisibleDiffuseDown)
       call mpas_pool_get_array(atmosCoupling, 'shortwaveIRDirectDown', shortwaveIRDirectDown)
       call mpas_pool_get_array(atmosCoupling, 'shortwaveIRDiffuseDown', shortwaveIRDiffuseDown)
       call mpas_pool_get_array(atmosCoupling, 'longwaveDown', longwaveDown)
       call mpas_pool_get_array(atmosCoupling, 'rainfallRate', rainfallRate)
       call mpas_pool_get_array(atmosCoupling, 'snowfallRate', snowfallRate)
       call mpas_pool_get_array(atmosCoupling, 'uAirVelocity', uAirVelocity)
       call mpas_pool_get_array(atmosCoupling, 'vAirVelocity', vAirVelocity)


       do iCell = 1, nCellsSolve
       
          gcell = iCell + cell_offset
             
          seaSurfaceTemperature(iCell)       = So_t    (gcell)
          seaSurfaceSalinity(iCell)          = So_s    (gcell)

          seaFreezingTemperature(iCell) = colpkg_sea_freezing_temperature(seaSurfaceSalinity(iCell))

          uOceanVelocity(iCell)              = So_u    (gcell)
          vOceanVelocity(iCell)              = So_v    (gcell)
          seaSurfaceTiltU(iCell)             = So_dhdx (gcell)
          seaSurfaceTiltV(iCell)             = So_dhdy (gcell)

          if (trim(config_ocean_surface_type) == "free") then ! free surface (MPAS-O)
        
             ! freezingMeltingPotential(i) is the ocean energy associated with frazil formation 
             ! when it is positive and frazilMassFlux is positive. Conversely, freezingMeltingPotential(i)
             ! is negative when there is the melting potential in which case frazilMassFlux is zero.

             freezingMeltingPotential(iCell) = Fioo_q  (gcell)

             frazilMassFlux              = Fioo_frazil (gcell)

             ! Now determine the sea ice mass associated with the frazil heat flux given when
             ! freezingMeltingPotential(i) is positive. This produces a revised mass flux, given 
             ! in frazilMassFluxRev for the given sea surface salinity. The resulting difference
             ! is assigned to frazilMassAdjust(i) which is exported to the ocean in the subsequent
             ! coupling step as a freshwater and salt flux. This step is required to balance mass
             ! and heat with the ocean.

             call frazil_mass(freezingMeltingPotential(iCell), frazilMassFluxRev, seaSurfaceSalinity(iCell), &
                              config_thermodynamics_type)

             frazilMassAdjust(iCell) = frazilMassFlux-frazilMassFluxRev 

          else ! non-free surface (SOM)

             freezingMeltingPotential(iCell) = Fioo_q  (gcell)

          endif

          airLevelHeight(iCell)              = Sa_z      (gcell)
          airPotentialTemperature(iCell)     = Sa_ptem   (gcell)
          airTemperature(iCell)              = Sa_tbot   (gcell)
          airSpecificHumidity(iCell)         = Sa_shum   (gcell)
          airDensity(iCell)                  = Sa_dens   (gcell)
          shortwaveVisibleDirectDown(iCell)  = Faxa_swvdr(gcell)
          shortwaveVisibleDiffuseDown(iCell) = Faxa_swvdf(gcell)
          shortwaveIRDirectDown(iCell)       = Faxa_swndr(gcell)
          shortwaveIRDiffuseDown(iCell)      = Faxa_swndf(gcell)
          longwaveDown(iCell)                = Faxa_lwdn (gcell)
          rainfallRate(iCell)                = Faxa_rain (gcell)
          snowfallRate(iCell)                = Faxa_snow (gcell)
          uAirVelocity(iCell)                = Sa_u      (gcell)
          vAirVelocity(iCell)                = Sa_v      (gcell)

         ! set aerosols, if configured
         !DDif (config_use_aerosols) then
         !DD   if (config_use_modal_aerosols) then
         !DD      atmosAerosolFlux(1,iCell) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_bcphidry, n)
         !DD      atmosAerosolFlux(2,iCell) = x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
         !DD      ! combine all the dust into one category
         !DD      atmosAerosolFlux(3,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
         !DD   else
         !DD      atmosAerosolFlux(1,iCell) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n)
         !DD      atmosAerosolFlux(2,iCell) = x2i_i % rAttr(index_x2i_Faxa_bcphidry, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
         !DD      ! combine all the dust into one category
         !DD      atmosAerosolFlux(3,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
         !DD   endif
         !DDendif

         !DD! import biogeochemistry fields, if configured
         !DDif (config_use_column_biogeochemistry) then
         !DD   oceanAlgaeConc(1,iCell)           = x2i_i % rAttr(index_x2i_So_algae1, n)
         !DD   oceanAlgaeConc(2,iCell)           = x2i_i % rAttr(index_x2i_So_algae2, n)
         !DD   oceanAlgaeConc(3,iCell)           = x2i_i % rAttr(index_x2i_So_algae3, n)
         !DD   oceanDOCConc(1,iCell)             = x2i_i % rAttr(index_x2i_So_doc1, n)
         !DD   oceanDOCConc(2,iCell)             = x2i_i % rAttr(index_x2i_So_doc2, n)
         !DD   oceanDOCConc(3,iCell)             = 0.0_RKIND
         !DD   oceanDICConc(1,iCell)             = x2i_i % rAttr(index_x2i_So_dic1, n) !JW not used, set to 0?
         !DD   oceanDONConc(1,iCell)             = x2i_i % rAttr(index_x2i_So_don1, n)
         !DD   oceanNitrateConc(iCell)           = x2i_i % rAttr(index_x2i_So_no3, n)
         !DD   oceanSilicateConc(iCell)          = x2i_i % rAttr(index_x2i_So_sio3, n)
         !DD   oceanAmmoniumConc(iCell)          = x2i_i % rAttr(index_x2i_So_nh4, n)
         !DD   oceanDMSConc(iCell)               = x2i_i % rAttr(index_x2i_So_dms, n)
         !DD   oceanDMSPConc(iCell)              = x2i_i % rAttr(index_x2i_So_dmsp, n)
         !DD   oceanHumicsConc(iCell)            = x2i_i % rAttr(index_x2i_So_docr, n)
         !DD   oceanParticulateIronConc(1,iCell) = x2i_i % rAttr(index_x2i_So_fep1, n)
         !DD   oceanParticulateIronConc(2,iCell) = x2i_i % rAttr(index_x2i_So_fep2, n)
         !DD   oceanDissolvedIronConc(1,iCell)   = x2i_i % rAttr(index_x2i_So_fed1, n)
         !DD   oceanDissolvedIronConc(2,iCell)   = x2i_i % rAttr(index_x2i_So_fed2, n)
         !DD   oceanZAerosolConc(1,iCell)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer1, n)
         !DD   oceanZAerosolConc(2,iCell)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer2, n)
         !DD   oceanZAerosolConc(3,iCell)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer3, n)
         !DD   oceanZAerosolConc(4,iCell)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer4, n)
         !DD   oceanZAerosolConc(5,iCell)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer5, n)
         !DD   oceanZAerosolConc(6,iCell)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer6, n)
         !DD   ! set aerosols, if configured
         !DD   if (config_use_zaerosols) then
         !DD      if (config_use_modal_aerosols) then
         !DD         atmosBlackCarbonFlux(1,iCell) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n) &
         !DD                                   + x2i_i % rAttr(index_x2i_Faxa_bcphidry, n)
         !DD         atmosBlackCarbonFlux(2,iCell) = x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
         !DD         ! combine wet and dry dust
         !DD         atmosDustFlux(1,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n)
         !DD         atmosDustFlux(2,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n)
         !DD         atmosDustFlux(3,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n)
         !DD         atmosDustFlux(4,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
         !DD      else
         !DD         atmosBlackCarbonFlux(1,iCell) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n)
         !DD         atmosBlackCarbonFlux(2,iCell) = x2i_i % rAttr(index_x2i_Faxa_bcphidry, n) &
         !DD                                   + x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
         !DD         ! combine wet and dry dust
         !DD         atmosDustFlux(1,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n)
         !DD         atmosDustFlux(2,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n)
         !DD         atmosDustFlux(3,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n)
         !DD         atmosDustFlux(4,iCell) = x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
         !DD                            + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
         !DD      endif
         !DD   endif
         !DDendif
       end do

       !----------------------------------------------------------------------- 
       !
       !  unit conversions and any manipulation of coupled fields
       !
       !-----------------------------------------------------------------------
       do iCell = 1, nCellsSolve
          seaSurfaceTemperature(iCell)  = seaSurfaceTemperature(iCell) - seaiceFreshWaterFreezingPoint

          !DDif (config_use_column_biogeochemistry) then
          !DD  ! convert from mmol C/m^3 to mmol N/m^3
          !DD  oceanAlgaeConc(1,iCell)           = oceanAlgaeConc(1,iCell) / carbonToNitrogenRatioAlgae(1)
          !DD  oceanAlgaeConc(2,iCell)           = oceanAlgaeConc(2,iCell) / carbonToNitrogenRatioAlgae(2)
          !DD  oceanAlgaeConc(3,iCell)           = oceanAlgaeConc(3,iCell) / carbonToNitrogenRatioAlgae(3)
          !DD  ! convert from mmol Fe/m^3 to umol Fe/m^3
          !DD  oceanParticulateIronConc(1,iCell) = oceanParticulateIronConc(1,iCell) * 1000._RKIND
          !DD  oceanParticulateIronConc(2,iCell) = oceanParticulateIronConc(2,iCell) * 1000._RKIND
          !DD  oceanDissolvedIronConc(1,iCell)   = oceanDissolvedIronConc(1,iCell)   * 1000._RKIND
          !DD  oceanDissolvedIronConc(2,iCell)   = oceanDissolvedIronConc(2,iCell)   * 1000._RKIND
          !DDendif
       end do

       cell_offset = cell_offset + nCellsSolve
       
       block => block % next
    end do

    !-----------------------------------------------------------------------
    ! update ghost cells for fluxes received from the mediator
    !-----------------------------------------------------------------------
    call mpas_pool_get_subpool(domain % blocklist % structs, 'ocean_coupling', oceanCoupling)
    call mpas_pool_get_subpool(domain % blocklist % structs, 'atmos_coupling', atmosCoupling)

    call mpas_pool_get_field(atmosCoupling, 'airLevelHeight', airLevelHeightField)
    call mpas_pool_get_field(atmosCoupling, 'airPotentialTemperature', airPotentialTemperatureField)
!    call mpas_pool_get_field(atmosCoupling, 'evaporativeWaterFlux', evaporativeWaterFluxField)
    call mpas_pool_get_field(atmosCoupling, 'airTemperature', airTemperatureField)
    call mpas_pool_get_field(atmosCoupling, 'airSpecificHumidity', airSpecificHumidityField)
    call mpas_pool_get_field(atmosCoupling, 'airDensity', airDensityField)
    call mpas_pool_get_field(atmosCoupling, 'shortwaveVisibleDirectDown', shortwaveVisibleDirectDownField)
    call mpas_pool_get_field(atmosCoupling, 'shortwaveVisibleDiffuseDown', shortwaveVisibleDiffuseDownField)
    call mpas_pool_get_field(atmosCoupling, 'shortwaveIRDirectDown', shortwaveIRDirectDownField)
    call mpas_pool_get_field(atmosCoupling, 'shortwaveIRDiffuseDown', shortwaveIRDiffuseDownField)
    call mpas_pool_get_field(atmosCoupling, 'rainfallRate', rainfallRateField)
    call mpas_pool_get_field(atmosCoupling, 'snowfallRate', snowfallRateField)
    call mpas_pool_get_field(atmosCoupling, 'uAirVelocity', uAirVelocityField)
    call mpas_pool_get_field(atmosCoupling, 'vAirVelocity', vAirVelocityField)
    
    call mpas_pool_get_field(oceanCoupling, 'seaSurfaceTemperature', seaSurfaceTemperatureField)
    call mpas_pool_get_field(oceanCoupling, 'seaSurfaceSalinity', seaSurfaceSalinityField)
    call mpas_pool_get_field(oceanCoupling, 'seaFreezingTemperature', seaFreezingTemperatureField)
    call mpas_pool_get_field(oceanCoupling, 'freezingMeltingPotential', freezingMeltingPotentialField)
    call mpas_pool_get_field(oceanCoupling, 'frazilMassAdjust', frazilMassAdjustField)
    call mpas_pool_get_field(oceanCoupling, 'uOceanVelocity', uOceanVelocityField)
    call mpas_pool_get_field(oceanCoupling, 'vOceanVelocity', vOceanVelocityField)
    call mpas_pool_get_field(oceanCoupling, 'seaSurfaceTiltU', seaSurfaceTiltUField)
    call mpas_pool_get_field(oceanCoupling, 'seaSurfaceTiltV', seaSurfaceTiltVField)

    !DDif (config_use_aerosols) then
    !DD   call mpas_pool_get_subpool(domain % blocklist % structs, 'aerosols', aerosols)

    !DD   call mpas_pool_get_field(aerosols, "atmosAerosolFlux", atmosAerosolFluxField)
    !DDendif

    !DDif (config_use_column_biogeochemistry) then
    !DD   call mpas_pool_get_subpool(domain % blocklist % structs, 'biogeochemistry', biogeochemistry)

    !DD   call mpas_pool_get_field(biogeochemistry, 'oceanAlgaeConc', oceanAlgaeConcField)
    !DD + more (see ice_comp_mct.F)
    !DD   if (config_use_zaerosols) then
    !DD      call mpas_pool_get_field(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFluxField)
    !DD + more (see ice_comp_mct.F)
    !DD   endif
    !DDendif

    if ( seaSurfaceTemperatureField % isActive ) &
       call mpas_dmpar_exch_halo_field(seaSurfaceTemperatureField)
    if ( seaSurfaceSalinityField % isActive ) &
       call mpas_dmpar_exch_halo_field(seaSurfaceSalinityField)
    if ( seaFreezingTemperatureField % isActive ) &
       call mpas_dmpar_exch_halo_field(seaFreezingTemperatureField)
    if ( freezingMeltingPotentialField % isActive ) &
       call mpas_dmpar_exch_halo_field(freezingMeltingPotentialField)
    if ( frazilMassAdjustField % isActive ) &
       call mpas_dmpar_exch_halo_field(frazilMassAdjustField)
    if ( uOceanVelocityField % isActive ) &
       call mpas_dmpar_exch_halo_field(uOceanVelocityField)
    if ( vOceanVelocityField % isActive ) &
       call mpas_dmpar_exch_halo_field(vOceanVelocityField)
    if ( seaSurfaceTiltUField % isActive ) &
       call mpas_dmpar_exch_halo_field(seaSurfaceTiltUField)
    if ( seaSurfaceTiltVField % isActive ) &
       call mpas_dmpar_exch_halo_field(seaSurfaceTiltVField)

    if ( airLevelHeightField % isActive ) &
       call mpas_dmpar_exch_halo_field(airLevelHeightField)
    if ( airPotentialTemperatureField % isActive ) &
       call mpas_dmpar_exch_halo_field(airPotentialTemperatureField)
!    if ( evaporativeWaterFluxField % isActive ) &
!       call mpas_dmpar_exch_halo_field(evaporativeWaterFluxField)
    if ( airTemperatureField % isActive ) &
       call mpas_dmpar_exch_halo_field(airTemperatureField)
    if ( airSpecificHumidityField % isActive ) &
       call mpas_dmpar_exch_halo_field(airSpecificHumidityField)
    if ( airDensityField % isActive ) &
       call mpas_dmpar_exch_halo_field(airDensityField)
    if ( shortwaveVisibleDirectDownField % isActive ) &
       call mpas_dmpar_exch_halo_field(shortwaveVisibleDirectDownField)
    if ( shortwaveVisibleDiffuseDownField % isActive ) &
       call mpas_dmpar_exch_halo_field(shortwaveVisibleDiffuseDownField)
    if ( shortwaveIRDirectDownField % isActive ) &
       call mpas_dmpar_exch_halo_field(shortwaveIRDirectDownField)
    if ( shortwaveIRDiffuseDownField % isActive ) &
       call mpas_dmpar_exch_halo_field(shortwaveIRDiffuseDownField)
    if ( rainfallRateField % isActive ) &
       call mpas_dmpar_exch_halo_field(rainfallRateField)
    if ( snowfallRateField % isActive ) &
       call mpas_dmpar_exch_halo_field(snowfallRateField)
    if ( uAirVelocityField % isActive ) &
       call mpas_dmpar_exch_halo_field(uAirVelocityField)
    if ( vAirVelocityField % isActive ) &
       call mpas_dmpar_exch_halo_field(vAirVelocityField)

    !DDif (config_use_aerosols) then
    !DDendif

    !DDif (config_use_column_biogeochemistry) then
    !DD   if (config_use_zaerosols) then
    !DD   endif
    !DDendif

    call t_stopf ('mpassi_imp_atm')

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ice_import

  !===============================================================================

  subroutine ice_export( exportState, flds_scalar_name, domain,    &
                         errorCode, rc )

    !-----------------------------------------------------------------------
    ! Create export state
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)                 :: exportState
    character(len=*)   , intent(in)  :: flds_scalar_name
    type (domain_type), intent(in), pointer :: domain
    integer            , intent(out) :: errorCode  ! pop error code
    integer            , intent(out) :: rc         ! returned error code

    ! local variables
    integer              :: n,i,j,k,iblock,nfld,lev
    character (cl) :: label
    real (r8)            :: m2percm2
    real (r8)            :: gsum
    integer              :: fieldCount
    character (cl), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname='(ice_import_export:ice_export)'

    ! local variables - mpas names
    integer :: iCell
    integer, pointer :: nCellsSolve
    real (kind=RKIND) :: ailohi, Tsrf, &
                         tauxa, tauya, tauxo, tauyo, basalPressure,   &
                         snowVolumeToSWE
                
    type (mpas_pool_type), pointer :: &
       configs,          &
       meshPool,         &
       tracersAggregate, &
       velocitySolver,   &
       shortwave,        &
       atmosCoupling,    &
       oceanCoupling,    &
       atmosFluxes,      &
       oceanFluxes,      &
       icebergFluxes,    &
       biogeochemistry
    type (block_type), pointer :: block
    integer gcell, cell_offset

    ! local variables - coupler names
    real (r8), pointer   :: Si_imask  (:) 
    real (r8), pointer   :: Si_ifrac  (:) 
    real (r8), pointer   :: Si_bpress (:) 
    real (r8), pointer   :: Si_t      (:) 
    real (r8), pointer   :: Si_vsno   (:) 
    real (r8), pointer   :: Si_u10    (:) 
    real (r8), pointer   :: Si_tref   (:) 
    real (r8), pointer   :: Si_snowh  (:) 
    real (r8), pointer   :: Si_qref   (:) 
    real (r8), pointer   :: Faii_swnet(:) 
    real (r8), pointer   :: Fioi_melth(:) 
    real (r8), pointer   :: Fioi_swpen(:) 
    real (r8), pointer   :: Fioi_meltw(:) 
    real (r8), pointer   :: Fioi_salt (:) 
    real (r8), pointer   :: Fioi_taux (:) 
    real (r8), pointer   :: Fioi_tauy (:) 
    real (r8), pointer   :: Si_avsdr  (:)
    real (r8), pointer   :: Si_anidr  (:)
    real (r8), pointer   :: Si_avsdf  (:)
    real (r8), pointer   :: Si_anidf  (:)
    real (r8), pointer   :: Faii_taux (:)
    real (r8), pointer   :: Faii_tauy (:)
    real (r8), pointer   :: Faii_lat  (:)
    real (r8), pointer   :: Faii_sen  (:)
    real (r8), pointer   :: Faii_lwup (:)
    real (r8), pointer   :: Faii_evap (:)

   logical, pointer :: &
      config_rotate_cartesian_grid,      &
      config_use_column_biogeochemistry, &
      config_use_column_shortwave,       &
      config_use_topo_meltponds,         &
      config_use_data_icebergs

   real(kind=RKIND), pointer :: &
      sphere_radius

   real (kind=RKIND), dimension(:), pointer :: &
      latCell,                     &
      lonCell,                     &
      xCell,                       &
      yCell,                       &
      zCell,                       &
      iceAreaCell,                 &
      iceVolumeCell,               &
      snowVolumeCell,              &
      pondDepthCell,               &
      pondLidThicknessCell,        &
      pondAreaCell,                &
      surfaceTemperatureCell,      &
      airStressCellU,              &
      airStressCellV,              &
      oceanStressCellU,            &
      oceanStressCellV,            &
      albedoVisibleDirectCell,     &
      albedoIRDirectCell,          &
      albedoVisibleDiffuseCell,    &
      albedoIRDiffuseCell,         &
      atmosReferenceSpeed10m,      &
      atmosReferenceTemperature2m, &
      atmosReferenceHumidity2m,    &
      latentHeatFlux,              &
      sensibleHeatFlux,            &
      longwaveUp,                  &
      evaporativeWaterFlux,        &
      absorbedShortwaveFlux,       &
      oceanHeatFlux,               &
      oceanShortwaveFlux,          &
      oceanFreshWaterFlux,         &
      oceanSaltFlux,               &
      frazilMassAdjust,            &
      bergFreshwaterFlux,          &
      bergLatentHeatFlux,          &
      oceanNitrateFlux,            &
      oceanSilicateFlux,           &
      oceanAmmoniumFlux,           &
      oceanDMSFlux,                &
      oceanDMSPpFlux,              &
      oceanDMSPdFlux,              &
      oceanHumicsFlux,             &
      carbonToNitrogenRatioAlgae,  &
      carbonToNitrogenRatioDON

    !-----------------------------------------------------
    rc = ESMF_SUCCESS
    errorCode = 0

    ! ice mask set one everywhere
    call state_getfldptr(exportState, 'Si_imask ', Si_imask , rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_imask  (:) = 1.0_RKIND ! all points are potentially ice points

    ! zero globally first
    call state_getfldptr(exportState, 'Si_ifrac ', Si_ifrac , rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_ifrac  (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Si_t'     , Si_t     , rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_t      (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Si_bpress', Si_bpress, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_bpress (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Si_vsno'  , Si_vsno  , rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_vsno   (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Si_u10'   , Si_u10   , rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_u10    (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Si_tref'  , Si_tref  , rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_tref   (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Si_snowh' , Si_snowh , rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_snowh  (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Si_qref'  , Si_qref  , rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Si_qref   (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Faii_swnet', Faii_swnet, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Faii_swnet(:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Fioi_melth', Fioi_melth, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Fioi_melth(:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Fioi_swpen', Fioi_swpen, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Fioi_swpen(:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Fioi_meltw', Fioi_meltw, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Fioi_meltw(:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Fioi_salt', Fioi_salt, rc)
      Fioi_salt (:) = 0.0_RKIND
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, 'Fioi_taux', Fioi_taux, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Fioi_taux (:) = 0.0_RKIND
    call state_getfldptr(exportState, 'Fioi_tauy', Fioi_tauy, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      Fioi_tauy (:) = 0.0_RKIND

    !if (config_use_column_shortwave) then
       call state_getfldptr(exportState, 'Si_avsdr', Si_avsdr, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Si_avsdr  (:) = 0.0_RKIND
       call state_getfldptr(exportState, 'Si_anidr', Si_anidr, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Si_anidr  (:) = 0.0_RKIND
       call state_getfldptr(exportState, 'Si_avsdf', Si_avsdf, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Si_avsdf  (:) = 0.0_RKIND
       call state_getfldptr(exportState, 'Si_anidf', Si_anidf, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Si_anidf  (:) = 0.0_RKIND
    !endif

    !DDif (send_i2x_per_cat) then
    !DD  call state_getfldptr(exportState, 'So_omask', So_omask, rc)
    !DD    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !DD      ????  (:) = 0.0_RKIND
    !DDendif
      
    !DDif (config_use_data_icebergs) then
    !DD  call state_getfldptr(exportState, 'So_omask', So_omask, rc)
    !DD    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !DD      ????  (:) = 0.0_RKIND
    !DDendif

    if (atm_prognostic) then
       call state_getfldptr(exportState, 'Faii_taux', Faii_taux, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Faii_taux (:) = 0.0_RKIND
       call state_getfldptr(exportState, 'Faii_tauy', Faii_tauy, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Faii_tauy (:) = 0.0_RKIND
       call state_getfldptr(exportState, 'Faii_lat', Faii_lat, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Faii_lat  (:) = 0.0_RKIND
       call state_getfldptr(exportState, 'Faii_sen', Faii_sen, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Faii_sen  (:) = 0.0_RKIND
       call state_getfldptr(exportState, 'Faii_lwup', Faii_lwup, rc)
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         Faii_lwup (:) = 0.0_RKIND
!       call state_getfldptr(exportState, 'Faii_evap', Faii_evap, rc)
!         if (ChkErr(rc,__LINE__,u_FILE_u)) return
!         Faii_evap (:) = 0.0_RKIND
    endif

    snowVolumeToSWE = seaiceDensitySnow * 0.001 
    
    cell_offset = 0
    block => domain % blocklist
    do while(associated(block))


       configs => block % configs
       call MPAS_pool_get_config(configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
       call MPAS_pool_get_config(configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(configs, "config_use_column_shortwave", config_use_column_shortwave)
       call MPAS_pool_get_config(configs, "config_use_data_icebergs", config_use_data_icebergs)
       call MPAS_pool_get_config(configs, "config_use_topo_meltponds", config_use_topo_meltponds)

       call MPAS_pool_get_subpool(block % structs, 'mesh', meshPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregate)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, 'atmos_coupling', atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, 'ocean_coupling', oceanCoupling)
       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxes)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxes)

       call MPAS_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call MPAS_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call MPAS_pool_get_config(meshPool, "sphere_radius", sphere_radius)
       call MPAS_pool_get_array(meshPool, "latCell", latCell)
       call MPAS_pool_get_array(meshPool, "lonCell", lonCell)
       call MPAS_pool_get_array(meshPool, "xCell", xCell)
       call MPAS_pool_get_array(meshPool, "yCell", yCell)
       call MPAS_pool_get_array(meshPool, "zCell", zCell)

       call MPAS_pool_get_array(tracersAggregate, 'iceAreaCell', iceAreaCell)
       call MPAS_pool_get_array(tracersAggregate, 'iceVolumeCell', iceVolumeCell)
       call MPAS_pool_get_array(tracersAggregate, 'snowVolumeCell', snowVolumeCell)
       call MPAS_pool_get_array(tracersAggregate, 'pondDepthCell', pondDepthCell)
       call MPAS_pool_get_array(tracersAggregate, 'pondLidThicknessCell', pondLidThicknessCell)
       call MPAS_pool_get_array(tracersAggregate, 'pondAreaCell', pondAreaCell)
       call MPAS_pool_get_array(tracersAggregate, 'surfaceTemperatureCell', surfaceTemperatureCell)

       call MPAS_pool_get_array(velocitySolver, 'airStressCellU', airStressCellU)
       call MPAS_pool_get_array(velocitySolver, 'airStressCellV', airStressCellV)
       call MPAS_pool_get_array(velocitySolver, 'oceanStressCellU', oceanStressCellU)
       call MPAS_pool_get_array(velocitySolver, 'oceanStressCellV', oceanStressCellV)

       call MPAS_pool_get_array(shortwave, 'albedoVisibleDirectCell', albedoVisibleDirectCell)
       call MPAS_pool_get_array(shortwave, 'albedoIRDirectCell', albedoIRDirectCell)
       call MPAS_pool_get_array(shortwave, 'albedoVisibleDiffuseCell', albedoVisibleDiffuseCell)
       call MPAS_pool_get_array(shortwave, 'albedoIRDiffuseCell', albedoIRDiffuseCell)
       call MPAS_pool_get_array(shortwave, 'absorbedShortwaveFlux', absorbedShortwaveFlux)

       call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceSpeed10m', atmosReferenceSpeed10m)
       call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceTemperature2m', atmosReferenceTemperature2m)
       call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceHumidity2m', atmosReferenceHumidity2m)

       call MPAS_pool_get_array(oceanCoupling, 'frazilMassAdjust', frazilMassAdjust)

       call MPAS_pool_get_array(atmosFluxes, 'latentHeatFlux', latentHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, 'sensibleHeatFlux', sensibleHeatFlux)
       call MPAS_pool_get_array(atmosFluxes, 'longwaveUp', longwaveUp)
!       call MPAS_pool_get_array(atmosFluxes, 'evaporativeWaterFlux', evaporativeWaterFlux)

       call MPAS_pool_get_array(oceanFluxes, 'oceanHeatFlux', oceanHeatFlux)
       call MPAS_pool_get_array(oceanFluxes, 'oceanShortwaveFlux', oceanShortwaveFlux)
       call MPAS_pool_get_array(oceanFluxes, 'oceanFreshWaterFlux', oceanFreshWaterFlux)
       call MPAS_pool_get_array(oceanFluxes, 'oceanSaltFlux', oceanSaltFlux)

       !DDif (config_use_data_icebergs) then
       !DD  call MPAS_pool_get_subpool(block_ptr % structs, "berg_fluxes", icebergFluxes)

       !DD  call MPAS_pool_get_array(icebergFluxes, "bergFreshwaterFlux", bergFreshwaterFlux)
       !DD  call MPAS_pool_get_array(icebergFluxes, "bergLatentHeatFlux", bergLatentHeatFlux)
       !DDendif

       !DDif (config_use_column_biogeochemistry) then
       !DD   call mpas_pool_get_subpool(block_ptr % structs, 'biogeochemistry', biogeochemistry)

       !DD   call mpas_pool_get_array(biogeochemistry, 'oceanAlgaeFlux', oceanAlgaeFlux)
       !DD      and more
       !DDendif

       do i = 1, nCellsSolve
          gcell = i + cell_offset
 
         ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
         call seaice_latlon_vector_rotation_backward(&
              tauxo,                &
              tauyo,                &
              -oceanStressCellU(i), &
              -oceanStressCellV(i), &
              latCell(i),           &
              lonCell(i),           &
              xCell(i),             &
              yCell(i),             &
              zCell(i),             &
              sphere_radius,        &
              config_rotate_cartesian_grid)

          ! ice fraction
          ailohi = min(iceAreaCell(i), 1.0_RKIND)

          !TODO: CICE has a check for ailohi < 0

          if ( ailohi > 0.0_RKIND ) then
             ! surface temperature
             Tsrf = seaiceFreshWaterFreezingPoint + surfaceTemperatureCell(i)

             ! basal pressure
             call basal_pressure(&
                  basalPressure,           &
                  iceVolumeCell(i),        &
                  snowVolumeCell(i),       &
                  pondDepthCell(i),        &
                  pondLidThicknessCell(i), &
                  pondAreaCell(i),         &
                  config_use_topo_meltponds)

             Si_ifrac  (gcell) = ailohi
             Si_t      (gcell) = Tsrf
             Si_bpress (gcell) = basalPressure
             Si_vsno   (gcell) = snowVolumeCell(i)
             Si_u10    (gcell) = atmosReferenceSpeed10m(i)
             Si_tref   (gcell) = atmosReferenceTemperature2m(i)
             Si_snowh  (gcell) = snowVolumeCell(i) * snowVolumeToSWE / ailohi
             Si_qref   (gcell) = atmosReferenceHumidity2m(i)
             Faii_swnet(gcell) = absorbedShortwaveFlux(i)
             Fioi_melth(gcell) = oceanHeatFlux(i)
             Fioi_swpen(gcell) = oceanShortwaveFlux(i)
             Fioi_meltw(gcell) = oceanFreshWaterFlux(i) + frazilMassAdjust(i)/ailohi
             Fioi_salt (gcell) = oceanSaltFlux(i) + ice_ref_salinity*p001*frazilMassAdjust(i)/ailohi
             Fioi_taux (gcell) = tauxo
             Fioi_tauy (gcell) = tauyo
          endif
       enddo

       if (config_use_column_shortwave) then
          do i = 1, nCellsSolve
             gcell = i + cell_offset

             Si_avsdr(gcell) = albedoVisibleDirectCell(i)
             Si_anidr(gcell) = albedoIRDirectCell(i)
             Si_avsdf(gcell) = albedoVisibleDiffuseCell(i)
             Si_anidf(gcell) = albedoIRDiffuseCell(i)
          enddo
       endif

       !DDif (send_i2x_per_cat) then
       !DD   do i = 1, nCellsSolve
       !DD      gcell = iCell + cell_offset

       !DD      ????  (gcell) = 
       !DD   enddo
       !DDendif
      
       !DDif (config_use_data_icebergs) then
       !DD   do i = 1, nCellsSolve
       !DD      gcell = iCell + cell_offset

       !DD      ????  (gcell) = 
       !DD   enddo
       !DDendif

       ! ------
       ! ice/atm fluxes computed by ice
       ! ------
       if (atm_prognostic) then
          do i = 1, nCellsSolve
             gcell = i + cell_offset

             ! wind stress  (on T-grid:  convert to lat-lon)
             call seaice_latlon_vector_rotation_backward(&
                  tauxa,             &
                  tauya,             &
                  airStressCellU(i), &
                  airStressCellV(i), &
                  latCell(i),        &
                  lonCell(i),        &
                  xCell(i),          &
                  yCell(i),          &
                  zCell(i),          &
                  sphere_radius,     &
                  config_rotate_cartesian_grid)

             Faii_taux (gcell) = tauxa
             Faii_tauy (gcell) = tauya

             Faii_lat  (gcell) = latentHeatFlux(i)
             Faii_sen  (gcell) = sensibleHeatFlux(i)
             Faii_lwup (gcell) = longwaveUp(i)
!             Faii_evap (gcell) = evaporativeWaterFlux(i)
          enddo
       endif

       cell_offset = cell_offset + nCellsSolve

       block => block % next
    enddo

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)


  end subroutine ice_export

  !===============================================================================

  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer             , intent(inout) :: num
    type(fld_list_type) , intent(inout) :: fldlist(:)
    character(len=*)    , intent(in)    :: stdname
    integer, optional   , intent(in)    :: ungridded_lbound
    integer, optional   , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call shr_sys_abort(trim(subname)//": ERROR num > fldsMax "//trim(stdname))
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================

  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, grid, tag, rc)

    use NUOPC, only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU
    use ESMF , only : ESMF_VM

    ! input/output variables
    type(ESMF_State)          , intent(inout) :: state
    type(fld_list_type)       , intent(in)    :: fldList(:)
    integer                   , intent(in)    :: numflds
    character(len=*)          , intent(in)    :: flds_scalar_name
    integer                   , intent(in)    :: flds_scalar_num
    character(len=*)          , intent(in)    :: tag
    type(ESMF_Mesh), optional , intent(in)    :: mesh
    type(ESMF_Grid), optional , intent(in)    :: grid
    integer                   , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(ice_import_export:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)

             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                                         ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                                         ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                                         gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if ! if not scalar field

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(ice_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================
  logical function State_FldChk(State, fldname)
    ! ----------------------------------------------
    ! Determine if field is in state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(in)  :: State
    character(len=*) , intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag) :: itemType
    ! ----------------------------------------------

    call ESMF_StateGet(State, trim(fldname), itemType)
    State_FldChk = (itemType /= ESMF_STATEITEM_NOTFOUND)

  end function State_FldChk

  !===============================================================================
  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)              , intent(in)     :: State
    character(len=*)              , intent(in)     :: fldname
    real(kind=dbl_kind) , pointer , intent(inout)  :: fldptr(:)
    integer, optional             , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ice_import_export:State_GetFldPtr_1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end subroutine State_GetFldPtr_1d

  !===============================================================================
  subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    ,            intent(in)     :: State
    character(len=*)    ,            intent(in)     :: fldname
    real(kind=dbl_kind) , pointer ,  intent(inout)  :: fldptr(:,:)
    integer             , optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ice_import_export:State_GetFldPtr_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
  end subroutine State_GetFldPtr_2d

!***********************************************************************
!BOP
!
! !IROUTINE: frazil_mass
!
! !INTERFACE
   subroutine frazil_mass(freezingPotential, frazilMassFlux, seaSurfaceSalinity, &
                                 config_thermodynamics_type)
!
! !DESCRIPTION:
! Calculate frazil mass based on on the sea surface salinity, and frazil heat flux 
! from the ocean, otherwise referred to as to as the freeze-melt potential. When 
! freezingPotential is positive, it gives the heat flux, according to the ocean model
! associated with the frazil mass passed from the ocean.  This function calculates
! the frazil mass based on the freezingPotential according to sea ice model thermodynamics,
!
! !USES:
      use ice_mushy_physics, only:  &
         liquidus_temperature_mush, &
         enthalpy_mush

      use ice_colpkg_shared, only: &
         dSin0_frazil,             &
         phi_init

! !INPUT PARAMETERS:
      real (kind=RKIND),      intent(in)   :: freezingPotential
      real (kind=RKIND),      intent(in)   :: seaSurfaceSalinity
      character(len=strKIND), intent(in)   :: config_thermodynamics_type

! !OUTPUT PARAMETERS:
      real (kind=RKIND),      intent(out)  :: frazilMassFlux

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real(kind=RKIND) :: &
      Si0new,          &
      Ti,              &
      qi0new,          &
      vi0new


   if (freezingPotential > 0.0_RKIND) then

      if (trim(config_thermodynamics_type) == "mushy") then  ! mushy
         if (seaSurfaceSalinity > 2.0_RKIND * dSin0_frazil) then
             Si0new = seaSurfaceSalinity - dSin0_frazil
         else
             Si0new = seaSurfaceSalinity**2 / (4.0_RKIND*dSin0_frazil)
         endif
         Ti = liquidus_temperature_mush(Si0new/phi_init)
         qi0new = enthalpy_mush(Ti, Si0new)
      else
         qi0new = -seaiceDensityIce*seaiceLatentHeatMelting
      endif    ! ktherm

      frazilMassFlux = -freezingPotential*seaiceDensityIce/qi0new

   else

      frazilMassFlux = 0.0_RKIND

   endif

! REVISION HISTORY:
! Revised Andrew Roberts May 2021
!-----------------------------------------------------------------------
!EOC

   end subroutine frazil_mass!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: basal_pressure
!
! !INTERFACE:
   subroutine basal_pressure(basalPressure,  iceVolumeCell,  snowVolumeCell,  pondDepthCell, &
                             pondLidThicknessCell,  pondAreaCell, config_use_topo_meltponds)!{{{
!
! !DESCRIPTION:
! Calculate basal pressure for a cell
!

! !INPUT PARAMETERS:
      real (kind=RKIND), intent(in)  :: iceVolumeCell
      real (kind=RKIND), intent(in)  :: snowVolumeCell
      real (kind=RKIND), intent(in)  :: pondDepthCell
      real (kind=RKIND), intent(in)  :: pondLidThicknessCell
      real (kind=RKIND), intent(in)  :: pondAreaCell
      logical,           intent(in)  :: config_use_topo_meltponds

! !OUTPUT PARAMETERS:
      real (kind=RKIND), intent(out) :: basalPressure

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real(kind=RKIND) :: &
      seaIceSpecificMass

   ! sea ice and snow specific mass
   seaIceSpecificMass = &
      iceVolumeCell  * seaiceDensityIce + &
      snowVolumeCell * seaiceDensitySnow

   ! only topo ponds have real pond volume
   if (config_use_topo_meltponds) then

      ! add pond specific weight
      seaIceSpecificMass = seaIceSpecificMass + &
         pondDepthCell        * pondAreaCell * seaiceDensitySeaWater + &
         pondLidThicknessCell * pondAreaCell * seaiceDensityIce

   endif ! config_use_topo_meltponds

   ! convert specific mass to pressure at sea ice base
   basalPressure = seaIceSpecificMass * seaiceGravity

!-----------------------------------------------------------------------
!EOC

   end subroutine basal_pressure!}}}

end module ice_import_export
