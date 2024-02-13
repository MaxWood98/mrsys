!MRsys Subdivision Main Program
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version : 0.2.4
!Updated : 25-01-2024

!Program
program mrsys
use mrsys_io_mod
use mrsys_oqc_mod
use mrsys_catmull_clark_mod
use mrsys_linear_subdivide_mod
use mrsys_emc_multiresolution_mod
use mrsys_smooth_reverse_catmull_clark_mod
implicit none 

!Variables
integer(in64) :: rr
real(dp) :: mr_parerror_log
type(mesh_data) :: mesh_base
type(options_data) :: options 
type(mressystem_data) :: mres_system
type(subdsystem_data), dimension(:), allocatable :: subd_system

!Set option config defaults 
call set_option_defaults(options)

!Get command arguments 
call get_process_arguments(options)
if (options%mode == 'exit') then 
    write(*,'(A)') '    {error exit}'
    stop 
end if 

!Import options 
call import_options(options)

!Display splash
if (options%console_disp == 'yes') then
    write(*,'(A)') ' '
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)')'|                   MRsys                    |'
    write(*,'(A)')'|     Multi-Resolution Subdivision Code      |'
    write(*,'(A)')'|        Version 0.2.6 || 13/02/2024         |'
    write(*,'(A)')'|                 Max Wood                   |'
    write(*,'(A)')'|           University of Bristol            |'
    write(*,'(A)')'|    Department of Aerospace Engineering     |'
    write(*,'(A)')'+--------------------------------------------+'
    write(*,'(A)') ' '
end if 

!Construct data folder
call construct_data_folder(options)

!Mode switch
if ((options%mode == 'sd') .OR. (options%mode == 'sdm'))then !Subdivide supplied geometry mesh (and export subdivision matricies if requested)

    !Display mode
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> mode -> subdivision'
    end if 

    !Import the base geometry mesh 
    call read_mesh(mesh_base,options,options%mesh_path)

    !Construct the subdivision matricies and propagate the subdivision
    if (options%subd_method == 'catmull_clark') then
        call propagate_catmull_clark(subd_system,mesh_base,options)
    elseif (options%subd_method == 'interpolating_catmull_clark') then

    elseif (options%subd_method == 'linear') then
        call propagate_linear_sd(subd_system,mesh_base,options)
    end if 

    !Export the final surface 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> exporting surface geometry'
    end if 
    call export_surface(options%data_path//'/mrsys_surface',subd_system(options%n_sd_levels+1)%mesh)

    !Export the subdivision matrices
    if (options%mode == 'sdm') then
        if (options%console_disp == 'yes') then
            write(*,'(A)') '--> exporting subdivision matricies'
        end if 
        do rr=1,options%n_sd_levels
            call export_csr_matrix(options%data_path//'/Ru/Ru'//int2str(rr),subd_system(rr)%Ru)
        end do 
    end if 
elseif (options%mode == 'bmrs') then !Construct multi-resolution system 

    !Display mode
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> mode -> multi-resolution system construct'
    end if 

    !Import the target geometry mesh 
    call read_mesh(mesh_base,options,options%tgtmesh_path)

    !Coarsen surface mesh 
    call ordered_quad_coarsen(mres_system%nrefine,mres_system%subd_system,mesh_base,options)
    mres_system%nlevel = mres_system%nrefine + 1

    !Evaluate the connectivity for each mesh level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> evaluating mesh connectivity at each level'
    end if
    do rr=1,mres_system%nlevel 
        call get_valence(mres_system%subd_system(rr)%mesh)
        call get_connectivity(mres_system%subd_system(rr)%mesh%connectivity,mres_system%subd_system(rr)%mesh)
    end do 

    !Store the target vertex positions at each level 
    do rr=1,mres_system%nlevel 
        allocate(mres_system%subd_system(rr)%mesh%vertices_tgt(&
        mres_system%subd_system(rr)%mesh%nvertex,mres_system%subd_system(rr)%mesh%ndim))
        mres_system%subd_system(rr)%mesh%vertices_tgt(:,:) = mres_system%subd_system(rr)%mesh%vertices(:,:)
    end do 

    !Construct the refinement and coarsening matrices for each level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> constructing refinement matricies'
    end if
    if (options%subd_method == 'catmull_clark') then
        call evaluate_mrsys_catmull_clark_refinement_matrices(mres_system,options)
    elseif (options%subd_method == 'interpolating_catmull_clark') then

    end if 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> constructing coarsening matricies'
    end if
    if (options%subd_method == 'catmull_clark') then
        call evaluate_mrsys_catmull_clark_coarsening_matrices(mres_system,options)
    elseif (options%subd_method == 'interpolating_catmull_clark') then

    end if 

    !Perform smooth reverse subdivision from the finest level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> reverse subdividing'
    end if
    if (options%subd_method == 'catmull_clark') then
        call propagate_mrsys_sr_catmull_clark(mres_system,options)
    elseif (options%subd_method == 'interpolating_catmull_clark') then

    end if 

    !Construct the detail matricies for each level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> constructing detail matrices'
    end if
    if (options%subd_method == 'catmull_clark') then
        call construct_emc_detail_matricies(mres_system,options)
    elseif (options%subd_method == 'interpolating_catmull_clark') then

    end if 

    !Construct the full refinement matricies for each level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> constructing multi-resolution refinement matrices'
    end if
    do rr=1,mres_system%nrefine
        call csrhrzcat(mres_system%subd_system(rr)%RuQ,mres_system%subd_system(rr)%Ru,mres_system%subd_system(rr)%Q)
    end do 
    
    !Construct the detail vectors for each level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> constructing multi-resolution detail vectors'
    end if
    if (options%subd_method == 'catmull_clark') then
        call construct_emc_detail_vectors(mres_system,options)
    elseif (options%subd_method == 'interpolating_catmull_clark') then

    end if 

    !Evaluate the final parameterisation error 
    mr_parerror_log = norm2(mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,1) - &
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,1))
    mr_parerror_log = mr_parerror_log + norm2(mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,2) - &
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,2))
    if (mesh_base%ndim == 3) then 
        mr_parerror_log = mr_parerror_log + norm2(mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,3) - &
        mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,3))
    end if 
    if (mr_parerror_log .NE. 0.0d0) then 
        mr_parerror_log = log10(mr_parerror_log)
    end if 
    if (options%console_disp == 'yes') then
        write(*,'(A,A,A)') '    < final surface error lognorm = ',real2F0_Xstring(mr_parerror_log,8_in64),' >'
    end if

    !Export the multiresolution system 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> exporting multi-resolution data'
    end if
    call export_mrsys_data(options%data_path//'/'//options%mrsystem_name,mres_system)

    !Export surfaces if requested 
    if (options%mrcon_export_meshes == 'yes') then 
        call export_surface(options%data_path//'/mrsys_surface',mres_system%subd_system(mres_system%nlevel)%mesh)
        do rr=1,mres_system%nlevel 
            call export_surface(options%data_path//'/mrsys_surface_level_'//int2str(rr),mres_system%subd_system(rr)%mesh)
        end do 
    end if 
elseif (options%mode == 'pmrs') then !Propagate deformation though multi-resolution system 

    !Display mode
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> mode -> multi-resolution system deform'
    end if 

    !Import multi-resolution system 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> importing multi-resolution data'
    end if
    call import_mrsys_data(options%mrsystemf,mres_system)

    !Import deformation file 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> importing deformations'
    end if
    call import_deformation_file(options%deformationf,mres_system%deformation)

    !Display
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> propagating deformations'
    end if

    !Apply the deformation 
    call emc_multiresolution_deform(mres_system,mres_system%deformation,options%mrdeflevel)

    !Propagate deformation to the final surface 
    call emc_multiresolution_refine_2_end(mres_system,options,options%mrdeflevel)

    !Display
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> exporting deformed geometry'
    end if

    !Export the final deformed surface 
    call export_surface(options%data_path//'/mrsys_surface',mres_system%subd_system(mres_system%nlevel)%mesh)
elseif (options%mode == 'rmrs') then !Reset multi-resolution system control nets to match the current surface geometry 

    !Display mode
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> mode -> multi-resolution system control nets reset'
    end if 

    !Import the target geometry mesh 
    call read_mrsys_surface_mesh(mesh_base,options,options%tgtmesh_path)

    !Import multi-resolution system 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> importing multi-resolution data'
    end if
    call import_mrsys_data(options%mrsystemf,mres_system)

    !Evaluate the connectivity for each mesh level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> evaluating mesh connectivity at each level'
    end if
    do rr=1,mres_system%nlevel 
        call get_valence(mres_system%subd_system(rr)%mesh)
        call get_connectivity(mres_system%subd_system(rr)%mesh%connectivity,mres_system%subd_system(rr)%mesh)
    end do 

    !Set the finest level of the system to the new target surface 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> inserting new target surface vertex positions'
    end if
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,:) = mesh_base%vertices(:,:)
    allocate(mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(mesh_base%nvertex,mesh_base%ndim))
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,:) = mesh_base%vertices(:,:)

    !Perform smooth reverse subdivision from the new finest level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> reverse subdividing'
    end if
    if (options%subd_method == 'catmull_clark') then
        call propagate_mrsys_sr_catmull_clark(mres_system,options)
    elseif (options%subd_method == 'interpolating_catmull_clark') then

    end if 

    !Construct the detail vectors for each level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> updating multi-resolution detail vectors'
    end if
    if (options%subd_method == 'catmull_clark') then
        call construct_emc_detail_vectors(mres_system,options)
    elseif (options%subd_method == 'interpolating_catmull_clark') then

    end if 

    !Evaluate the final parameterisation error 
    mr_parerror_log = norm2(mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,1) - &
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,1))
    mr_parerror_log = mr_parerror_log + norm2(mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,2) - &
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,2))
    if (mesh_base%ndim == 3) then 
        mr_parerror_log = mr_parerror_log + norm2(mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,3) - &
        mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,3))
    end if 
    if (mr_parerror_log .NE. 0.0d0) then 
        mr_parerror_log = log10(mr_parerror_log)
    end if 
    if (options%console_disp == 'yes') then
        write(*,'(A,A,A)') '    < final surface error lognorm = ',real2F0_Xstring(mr_parerror_log,8_in64),' >'
    end if

    !Export the multiresolution system 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> exporting multi-resolution data'
    end if
    call export_mrsys_data(options%data_path//'/'//options%mrsystem_name,mres_system)

    !Export surfaces if requested 
    if (options%mrcon_export_meshes == 'yes') then 
        do rr=1,mres_system%nlevel 
            call export_surface(options%data_path//'/mrsys_surface_level_'//int2str(rr),mres_system%subd_system(rr)%mesh)
        end do 
    end if 
elseif (options%mode == 'gradients') then 

    !Display mode
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> mode -> surface response gradient evaluation'
    end if

    !Import multi-resolution system 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> importing multi-resolution data'
    end if
    call import_mrsys_data(options%mrsystemf,mres_system)

    !Evaluate gradient propagation matrix from the target level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> constructing surface sensitivity matrix'
    end if
    call emc_build_gradient_matrix(mres_system,options,options%mrdeflevel)

    !Evaluate gradients WRT each control point in the control net at this level 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> evaluating surface position gradients'
    end if
    call emc_evaluate_surface_gradients(mres_system,options%mrdeflevel)

    !Export gradients 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '--> exporting surface position gradients'
    end if
    call export_surface_gradients(options%data_path//'/'//options%gradient_fname,&
    mres_system%subd_system(options%mrdeflevel)%surface_grad)
end if 

!Display complete
if (options%console_disp == 'yes') then
    write(*,'(A)') '    {complete}'
end if 
stop 
end program mrsys