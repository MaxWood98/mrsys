!MRsys io module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.7
!Updated 14-12-2023

!Module
module mrsys_io_mod
use io_utilities
use mrsys_connectivity_mod
contains 


!Read command arguments subroutine ===========================
subroutine get_process_arguments(options)
implicit none

!Variables - Import
type(options_data) :: options

!Variables - Local  
integer(in32) :: arglen,argstat,nargs
character(len=:), allocatable :: arg_nm1

!Check and process supplied command arguments
nargs = command_argument_count()
if (nargs == 0) then 
    write(*,'(A)') '** one operation must be requested :' 
    write(*,'(A)') '[sd "surface" "nlevel" -> subdivide "surface" "nlevel" times]' 
    write(*,'(A)') '[sdm "surface" "nlevel" -> subdivide "surface" "nlevel" times and export subdivision matricies]' 
    write(*,'(A)') '[bmrs "surface" "nclevel" -> build multi-resolution system from mesh "surface" with at most "nclevel" levels]' 
    write(*,'(A)') '[pmrs "system" "levdef" "deformation" -> propagate multi-resolution system "system" from control net level &
    &"levdef" with deformations file "deformation"]'
    write(*,'(A)') '[rmrs "system" "surface" -> reset all control nets in the multi-resolution system "system" to match the &
    &surface geometry "surface" with no applied deformation]' 
    write(*,'(A)') '[gradients "system" "levdef" -> evaluate surface position with control point movement gradients for the&
    & multi-resolution system "system" when deformed from level "levdef"]' 
    write(*,'(A)') '[optionally the path to an options file can be specified with "-o path" as the final arguments &
    &(the default is "io/MRsys_options.dat")]'
    allocate(character(len=4) :: options%mode)
    options%mode = 'exit'
else

    !Read mode
    call get_command_argument(number=1, length=arglen)
    allocate(character(len=arglen) :: options%mode)
    call get_command_argument(number=1, value=options%mode, status=argstat)

    !Display error mode 
    if ((options%mode .NE. 'sd') .AND. (options%mode .NE. 'sdm') .AND. &
    (options%mode .NE. 'bmrs') .AND. (options%mode .NE. 'pmrs') .AND. &
    (options%mode .NE. 'gradients') .AND. (options%mode .NE. 'rmrs')) then 
        write(*,'(A)') '** unknown mode request "'//options%mode//'" specified'
        deallocate(options%mode)
        allocate(character(len=4) :: options%mode)
        options%mode = 'exit'
        return 
    end if 

    !Read additional arguments
    if ((options%mode == 'sd') .OR. (options%mode == 'sdm')) then 
        if (nargs .LT. 3) then 
            write(*,'(A)') '** missing arguments' 
            write(*,'(A)') '[mrsys sd "surface" "nlevel"]' 
            write(*,'(A)') '[mrsys sdm "surface" "nlevel"]' 
            options%mode = 'exit'
        else
            call get_command_argument(number=2, length=arglen)
            allocate(character(len=arglen) :: options%mesh_path)
            call get_command_argument(number=2, value=options%mesh_path, status=argstat)
            options%n_sd_levels = get_command_argument_n_int(3_in32) 
        end if 
    end if 
    if (options%mode == 'bmrs') then 
        if (nargs .LT. 3) then 
            write(*,'(A)') '** missing arguments' 
            write(*,'(A)') '[mrsys bmrs "surface" "nclevel"]' 
            options%mode = 'exit'
        else
            call get_command_argument(number=2, length=arglen)
            allocate(character(len=arglen) :: options%tgtmesh_path)
            call get_command_argument(number=2, value=options%tgtmesh_path, status=argstat)
            options%n_crs_levels_max = get_command_argument_n_int(3_in32) 
        end if 
    end if 
    if (options%mode == 'pmrs') then 
        if (nargs .LT. 4) then 
            write(*,'(A)') '** missing arguments' 
            write(*,'(A)') '[mrsys pmrs "system" "levdef" "deformation"]'
            options%mode = 'exit'
        else
            call get_command_argument(number=2, length=arglen)
            allocate(character(len=arglen) :: options%mrsystemf)
            call get_command_argument(number=2, value=options%mrsystemf, status=argstat)
            options%mrdeflevel = get_command_argument_n_int(3_in32) 
            call get_command_argument(number=4, length=arglen)
            allocate(character(len=arglen) :: options%deformationf)
            call get_command_argument(number=4, value=options%deformationf, status=argstat)
        end if 
    end if 
    if (options%mode == 'rmrs') then 
        if (nargs .LT. 3) then 
            write(*,'(A)') '** missing arguments' 
            write(*,'(A)') '[rmrs "system" "surface"]'
            options%mode = 'exit'
        else
            call get_command_argument(number=2, length=arglen)
            allocate(character(len=arglen) :: options%mrsystemf)
            call get_command_argument(number=2, value=options%mrsystemf, status=argstat)
            call get_command_argument(number=3, length=arglen)
            allocate(character(len=arglen) :: options%tgtmesh_path)
            call get_command_argument(number=3, value=options%tgtmesh_path, status=argstat)
        end if 
    end if
    if (options%mode == 'gradients') then 
        if (nargs .LT. 3) then 
            write(*,'(A)') '** missing arguments' 
            write(*,'(A)') '[mrsys gradients "system" "levdef"]' 
            options%mode = 'exit'
        else
            call get_command_argument(number=2, length=arglen)
            allocate(character(len=arglen) :: options%mrsystemf)
            call get_command_argument(number=2, value=options%mrsystemf, status=argstat)
            options%mrdeflevel = get_command_argument_n_int(3_in32) 
        end if 
    end if 

    !Check if options path is specified then read and set if so 
    call get_command_argument(number=nargs-1, length=arglen)
    allocate(character(len=arglen) :: arg_nm1)
    call get_command_argument(number=nargs-1, value=arg_nm1, status=argstat)
    if (arg_nm1 == '-o') then 
        if (allocated(options%options_path)) then 
            deallocate(options%options_path)
        end if 
        call get_command_argument(number=nargs, length=arglen)
        allocate(character(len=arglen) :: options%options_path)
        call get_command_argument(number=nargs, value=options%options_path, status=argstat)
    end if 
end if 
return 
end subroutine get_process_arguments




!Set option defaults subroutine =========================
subroutine set_option_defaults(options)
implicit none 

!Variables - Import
type(options_data) :: options 

!Set option defaults
options%options_path = 'MRsys_options.dat'
options%console_disp = 'yes'
options%data_path = 'mrsys_data'
options%mrsystem_name = 'mrsystem_data'
options%subd_method = 'catmull_clark'
options%crs_base_vlnc = 'auto'
options%srcc_omega = 0.5d0 
options%mrcon_export_meshes = 'yes'
options%gradient_fname = 'mr_surfgrad'
return 
end subroutine set_option_defaults




!Import options subroutine =========================
subroutine import_options(options)
implicit none 

!Variables - Import
type(options_data) :: options 

!Variables - Local 
integer(in32) :: fh
logical :: ofexists

!Check if base file exists 
inquire(file=options%options_path,exist=ofexists)

!Return if options file does not exist
if (.NOT.ofexists) then 
    write(*,'(A)') '** options file "'//options%options_path//'" does not exist -> using default options'
    return 
end if 

!Open file 
open(newunit=fh,file=options%options_path)

!Scan for and set options that are present 
call set_str_opt(options%console_disp,fh,'con_disp')
call set_str_opt(options%data_path,fh,'data_path')
call set_str_opt(options%mrsystem_name,fh,'mrsystem_name')
call set_str_opt(options%subd_method,fh,'sd_method')
call set_str_opt(options%crs_base_vlnc,fh,'cb_vlnc')
call set_real_opt(options%srcc_omega,fh,'srcc_omega')
call set_str_opt(options%mrcon_export_meshes,fh,'exp_mr_surfs')
call set_str_opt(options%gradient_fname,fh,'gradient_fname')

!Close file 
close(fh)
return 
end subroutine import_options




!Import mesh subroutine =========================
subroutine read_mesh(mesh,options,filename)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(mesh_data) :: mesh
type(options_data) :: options 

!Variables - Local
integer(in64) :: ii,jj
integer(in64) :: iostatus,idxsp1
real(dp) :: edgeline(3)
real(dp), dimension(:), allocatable :: vertexline
character(len=1000) :: rtemp 

!Display
if (options%console_disp == 'yes') then
    write(*,'(A,I0,A)') '--> importing mesh'
end if 

!Open mesh file 
open(11,file=filename)

!Read number of dimensions
mesh%ndim = 0 
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:4) == 'ndim') then 
        read(rtemp(6:len_trim(rtemp)),*) mesh%ndim
        exit 
    end if 
end do
rewind(11)

!Read vertices 
mesh%nvertex = 0 
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:7) == 'nvertex') then 
        read(rtemp(9:len_trim(rtemp)),*) mesh%nvertex
        allocate(mesh%vertices(mesh%nvertex,mesh%ndim))
        allocate(mesh%vertex_sharp(mesh%nvertex))
        allocate(vertexline(mesh%ndim+1))
        do ii=1,mesh%nvertex
            read(11,*) vertexline 
            mesh%vertices(ii,:) = vertexline(1:mesh%ndim)
            mesh%vertex_sharp(ii) = int(vertexline(mesh%ndim+1),in64)
        end do 
        exit 
    end if 
end do
rewind(11)

!Read edges 
mesh%nedge = 0 
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:5) == 'nedge') then 
        read(rtemp(7:len_trim(rtemp)),*) mesh%nedge
        if (mesh%nedge .GT. 0) then 
            allocate(mesh%edges(mesh%nedge,2))
            allocate(mesh%edge_sharp(mesh%nedge))
            do ii=1,mesh%nedge
                read(11,*) edgeline 
                mesh%edges(ii,:) = int(edgeline(1:2),in64)
                mesh%edge_sharp(ii) = edgeline(3)
            end do 
        end if 
        exit 
    end if 
end do
rewind(11)

!Read faces 
mesh%nface = 0 
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:5) == 'nface') then 
        read(rtemp(7:len_trim(rtemp)),*) mesh%nface
        if (mesh%nface .GT. 0) then 
            allocate(mesh%faces(mesh%nface))
            do ii=1,mesh%nface
                read(11,'(A)') rtemp 
                idxsp1 = 0 
                do jj=1,1000
                    if (rtemp(jj:jj) == ' ') then 
                        read(rtemp(1:jj-1),*) mesh%faces(ii)%nvertex 
                        idxsp1 = jj
                        exit 
                    end if 
                end do 
                allocate(mesh%faces(ii)%vertices(mesh%faces(ii)%nvertex))
                read(rtemp(idxsp1+1:len_trim(rtemp)),*) mesh%faces(ii)%vertices(:)
            end do 
        end if 
        exit 
    end if 
end do

!Close file 
close(11)

!Display mesh metrics 
if (options%console_disp == 'yes') then
    write(*,'(A,I0,A)') '    {mesh contains ',mesh%ndim,' dimensions}'
    write(*,'(A,I0,A)') '    {imported = ',mesh%nvertex,' vertices}'
    write(*,'(A,I0,A)') '    {imported = ',mesh%nedge,' edges}'
    write(*,'(A,I0,A)') '    {imported = ',mesh%nface,' faces}'
end if

!If no edges then construct edges and assume none are sharp 
if (mesh%nedge == 0) then 
    if (options%console_disp == 'yes') then
        write(*,'(A)') '    {constructing edges}'
    end if 
    call get_valence(mesh)
    call build_edges_from_faces(mesh,mesh%connectivity%max_valence)
    if (options%console_disp == 'yes') then
        write(*,'(A,I0,A)') '    {constructed = ',mesh%nedge,' edges}'
    end if 
end if 
return 
end subroutine read_mesh




!Import mrsys surface mesh subroutine =========================
subroutine read_mrsys_surface_mesh(mesh,options,filename)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(mesh_data) :: mesh
type(options_data) :: options 

!Variables - Local
integer(in64) :: ii,jj
integer(in64) :: nedge_face,nvertex,ndim,lenreadface

!Display
if (options%console_disp == 'yes') then
    write(*,'(A,I0,A)') '--> importing mesh'
end if 

!Open mesh file 
open(11,file=filename)

!Read mesh data 
read(11,*) nvertex,nedge_face,ndim
mesh%ndim = ndim
mesh%nvertex = nvertex 
if (ndim == 2) then 
    mesh%nface = 0 
    mesh%nedge = nedge_face
elseif (ndim == 3) then 
    mesh%nface = nedge_face
    mesh%nedge = 0
end if 
read(11,*) lenreadface

!Read mesh vertices 
allocate(mesh%vertices(mesh%nvertex,mesh%ndim))
do ii=1,mesh%nvertex
    read(11,*) mesh%vertices(ii,:)
end do 

!Read mesh edges/faces
allocate(mesh%faces(mesh%nface))
allocate(mesh%edges(mesh%nedge,2))
if (ndim == 2) then !Read edges
    do ii=1,mesh%nedge
        read(11,*) lenreadface
        do jj=1,2
            read(11,*) mesh%edges(ii,jj)
        end do 
    end do 
elseif (ndim == 3) then !Read faces
    do ii=1,mesh%nface
        read(11,*) lenreadface
        mesh%faces(ii)%nvertex = lenreadface
        allocate(mesh%faces(ii)%vertices(mesh%faces(ii)%nvertex))
        do jj=1,lenreadface
            read(11,*) mesh%faces(ii)%vertices(jj)
        end do 
    end do 
end if 

!Close file 
close(11)

!Display mesh metrics 
if (options%console_disp == 'yes') then
    write(*,'(A,I0,A)') '    {mesh contains ',mesh%ndim,' dimensions}'
    write(*,'(A,I0,A)') '    {imported = ',mesh%nvertex,' vertices}'
    write(*,'(A,I0,A)') '    {imported = ',mesh%nedge,' edges}'
    write(*,'(A,I0,A)') '    {imported = ',mesh%nface,' faces}'
end if
return 
end subroutine read_mrsys_surface_mesh




!Export surface subroutine =========================
subroutine export_surface(filename,mesh)
implicit none 

!Variables - Import 
character(*), intent(in) :: filename
type(mesh_data) :: mesh

!Variables - Local 
integer(in64) :: ii,jj
integer(in64) :: lenreadface
character(len=50) :: vertices_str(mesh%nvertex)

!Find face read length 
lenreadface = 0
if (mesh%ndim == 2) then !2D -> edges
    do ii=1,mesh%nedge
        lenreadface = lenreadface + 1
        lenreadface = lenreadface + 2
    end do 
elseif (mesh%ndim == 3) then !3D -> faces
    do ii=1,mesh%nface
        lenreadface = lenreadface + 1
        lenreadface = lenreadface + mesh%faces(ii)%nvertex
    end do 
end if 

!Open file
open(11,file=filename)

!Item counts
if (mesh%ndim == 2) then !2D -> edges
    write(11,'(I0,A,I0,A,I0)') mesh%nvertex,' ',mesh%nedge,' ',mesh%ndim
elseif (mesh%ndim == 3) then !3D -> faces
    write(11,'(I0,A,I0,A,I0)') mesh%nvertex,' ',mesh%nface,' ',mesh%ndim
end if 
write(11,'(I0)') lenreadface

!Convert vertices to strings 
call array_real2F0_Xstring(vertices_str,mesh%vertices,8_in64) 

!Surface vertices 
do ii=1,mesh%nvertex
    write(11,'(A)') trim(vertices_str(ii))
end do

!Surface faces or edges 
if (mesh%ndim == 2) then !2D -> edges
    do ii=1,mesh%nedge
        write(11,'(I0)') 2
        do jj=1,2
            write(11,'(I0)') mesh%edges(ii,jj)
        end do 
    end do 
elseif (mesh%ndim == 3) then !3D -> faces
    do ii=1,mesh%nface
        write(11,'(I0)') mesh%faces(ii)%nvertex
        do jj=1,mesh%faces(ii)%nvertex
            write(11,'(I0)') mesh%faces(ii)%vertices(jj)
        end do 
    end do 
end if 

!Close file
close(11)
return 
end subroutine export_surface




!Export csr matrix subroutine =========================
subroutine export_csr_matrix(filename,mat)
implicit none 

!Variables - Import 
character(*), intent(in) :: filename
type(csrmatrix) :: mat

!Open file
open(11,file=filename,status='replace',access='stream')

!Export 
call dump_csr_mat(mat,11)

!Close file 
close(11)
return 
end subroutine export_csr_matrix




!Export multi-resolution system data subroutine =========================
subroutine export_mrsys_data(filename,mres_system)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(mressystem_data) :: mres_system

!Variables - Local 
integer(in64) :: ii,jj,rr 

!Open file
open(11,file=filename,status='replace',access='stream')

!Write system data 
write(11) mres_system%nlevel
write(11) mres_system%nrefine

!Export each sd refinement matrix Ru
do rr=1,mres_system%nrefine
    call dump_csr_mat(mres_system%subd_system(rr)%Ru,11)
end do 

!Export each mr refinement matrix RuQ
do rr=1,mres_system%nrefine
    call dump_csr_mat(mres_system%subd_system(rr)%RuQ,11)
end do 

!Export each mr coarsening matrix Ri
do rr=2,mres_system%nrefine+1
    call dump_csr_mat(mres_system%subd_system(rr)%Ri,11)
end do 

!Export each detail vector 
do rr=1,mres_system%nrefine+1
    call dump_realarr_bin(mres_system%subd_system(rr)%d,11)
end do 

!Export all surface vertex positions as constructed from the base geometry at each level
do rr=1,mres_system%nrefine+1
    call dump_realarr_bin(mres_system%subd_system(rr)%mesh%vertices,11)
end do 

!Export all surface vertex sharpness data from the base geometry at each level
do rr=1,mres_system%nrefine+1
    call dump_intvec_bin(mres_system%subd_system(rr)%mesh%vertex_sharp,11)
end do 

!Export the mesh faces (and edges) of each surface 
if (mres_system%subd_system(1)%mesh%ndim == 2) then !2D -> edges
    do rr=1,mres_system%nlevel
        write(11) mres_system%subd_system(rr)%mesh%nedge
        do ii=1,mres_system%subd_system(rr)%mesh%nedge
            write(11) mres_system%subd_system(rr)%mesh%edges(ii,:)
            write(11) mres_system%subd_system(rr)%mesh%edge_sharp(ii)
        end do 
    end do 
elseif (mres_system%subd_system(1)%mesh%ndim == 3) then !3D -> faces and edges
    do rr=1,mres_system%nlevel !edges
        write(11) mres_system%subd_system(rr)%mesh%nedge
        do ii=1,mres_system%subd_system(rr)%mesh%nedge
            write(11) mres_system%subd_system(rr)%mesh%edges(ii,:)
            write(11) mres_system%subd_system(rr)%mesh%edge_sharp(ii)
        end do 
    end do 
    do rr=1,mres_system%nlevel !faces
        write(11) mres_system%subd_system(rr)%mesh%nface
        do ii=1,mres_system%subd_system(rr)%mesh%nface
            write(11) mres_system%subd_system(rr)%mesh%faces(ii)%nvertex
            do jj=1,mres_system%subd_system(rr)%mesh%faces(ii)%nvertex
                write(11) mres_system%subd_system(rr)%mesh%faces(ii)%vertices(jj)
            end do 
        end do 
    end do    
end if 

!Close file 
close(11)
return 
end subroutine export_mrsys_data




!Import multi-resolution system data subroutine =========================
subroutine import_mrsys_data(filename,mres_system)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
type(mressystem_data) :: mres_system

!Variables - Local 
integer(in64) :: ii,jj,rr 

!Open file
open(11,file=filename,access='stream')

!Read system data 
read(11) mres_system%nlevel
read(11) mres_system%nrefine

!Allocate system 
allocate(mres_system%subd_system(mres_system%nlevel))

!Read each sd refinement matrix Ru
do rr=1,mres_system%nrefine
    call read_csr_mat(mres_system%subd_system(rr)%Ru,11)
end do 

!Read each mr refinement matrix RuQ
do rr=1,mres_system%nrefine
    call read_csr_mat(mres_system%subd_system(rr)%RuQ,11)
end do 

!Read each mr coarsening matrix Ri
do rr=2,mres_system%nrefine+1
    call read_csr_mat(mres_system%subd_system(rr)%Ri,11)
end do 

!Read each detail vector 
do rr=1,mres_system%nrefine+1
    call read_realarr_bin(mres_system%subd_system(rr)%d,11)
    mres_system%subd_system(rr)%n_d = size(mres_system%subd_system(rr)%d(:,1))
end do 

!Read all surface vertex positions as constructed from the base geometry at each level
do rr=1,mres_system%nrefine+1
    call read_realarr_bin(mres_system%subd_system(rr)%mesh%vertices,11)
    mres_system%subd_system(rr)%mesh%ndim = size(mres_system%subd_system(rr)%mesh%vertices(1,:))
    mres_system%subd_system(rr)%mesh%nvertex = size(mres_system%subd_system(rr)%mesh%vertices(:,1))
end do 

!Read all surface vertex sharpness data from the base geometry at each level
do rr=1,mres_system%nrefine+1
    call read_intvec_bin(mres_system%subd_system(rr)%mesh%vertex_sharp,11)
end do 

!Import the mesh faces (and edges) of each surface 
mres_system%subd_system(:)%mesh%nedge = 0 
mres_system%subd_system(:)%mesh%nface = 0 
if (mres_system%subd_system(1)%mesh%ndim == 2) then !2D -> edges
    do rr=1,mres_system%nlevel
        read(11) mres_system%subd_system(rr)%mesh%nedge
        allocate(mres_system%subd_system(rr)%mesh%edges(&
        mres_system%subd_system(rr)%mesh%nedge,2))
        allocate(mres_system%subd_system(rr)%mesh%edge_sharp(&
        mres_system%subd_system(rr)%mesh%nedge))
        do ii=1,mres_system%subd_system(rr)%mesh%nedge
            read(11) mres_system%subd_system(rr)%mesh%edges(ii,:)
            read(11) mres_system%subd_system(rr)%mesh%edge_sharp(ii)
        end do
    end do 
elseif (mres_system%subd_system(mres_system%nlevel)%mesh%ndim == 3) then !3D -> faces and edges
    do rr=1,mres_system%nlevel !edges
        read(11) mres_system%subd_system(rr)%mesh%nedge
        allocate(mres_system%subd_system(rr)%mesh%edges(&
        mres_system%subd_system(rr)%mesh%nedge,2))
        allocate(mres_system%subd_system(rr)%mesh%edge_sharp(&
        mres_system%subd_system(rr)%mesh%nedge))
        do ii=1,mres_system%subd_system(rr)%mesh%nedge
            read(11) mres_system%subd_system(rr)%mesh%edges(ii,:)
            read(11) mres_system%subd_system(rr)%mesh%edge_sharp(ii)
        end do
    end do 
    do rr=1,mres_system%nlevel !faces
        read(11) mres_system%subd_system(rr)%mesh%nface
        allocate(mres_system%subd_system(rr)%mesh%faces(&
        mres_system%subd_system(rr)%mesh%nface))
        do ii=1,mres_system%subd_system(rr)%mesh%nface
            read(11) mres_system%subd_system(rr)%mesh%faces(ii)%nvertex
            allocate(mres_system%subd_system(rr)%mesh%faces(ii)%&
            vertices(mres_system%subd_system(rr)%mesh%faces(ii)%nvertex))
            do jj=1,mres_system%subd_system(rr)%mesh%faces(ii)%nvertex
                read(11) mres_system%subd_system(rr)%mesh%faces(ii)%vertices(jj)
            end do 
        end do 
    end do 
end if

!Close file 
close(11)
return 
end subroutine import_mrsys_data




!Import deformation file =========================
subroutine import_deformation_file(filename,deformation)
implicit none 

!Variables - Import
character(*), intent(in) :: filename
real(dp), dimension(:,:), allocatable :: deformation

!Variables - Local 
integer(in64) :: ii 
integer(in64) :: ndim,nvertex,iostatus
character(len=1000) :: rtemp 

!Open file 
open(11,file=filename)

!Read number of dimensions
ndim = 0 
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:4) == 'ndim') then 
        read(rtemp(6:len_trim(rtemp)),*) ndim
        exit 
    end if 
end do
rewind(11)

!Read vertices 
nvertex = 0 
iostatus = 0 
do while (iostatus == 0)
    read(11,'(A)',iostat=iostatus) rtemp
    if (rtemp(1:7) == 'nvertex') then 
        read(rtemp(9:len_trim(rtemp)),*) nvertex
        allocate(deformation(nvertex,ndim))
        deformation(:,:) = 0.0d0 
        do ii=1,nvertex
            read(11,*) deformation(ii,:) 
        end do 
        exit 
    end if 
end do

!Close file 
close(11)
return 
end subroutine import_deformation_file




!Export surface gradient file =========================
subroutine export_surface_gradients(filename,gradients)
implicit none 

!Variables - Import
real(dp), dimension(:,:) :: gradients
character(*), intent(in) :: filename

!Variables - Local 
integer(in64) :: ii 
integer(in64) :: nrow
character(len=size(gradients,2)*15), dimension(:), allocatable :: gradients_str

!Open file 
open(11,file=filename)

!Convert to correct format
nrow = size(gradients,1)
allocate(gradients_str(nrow))
call array_real2F0_Xstring(gradients_str,gradients,8_in64) 

!Export gradients 
do ii=1,nrow
    write(11,'(A)') trim(gradients_str(ii))
end do 

!Close file 
close(11)
return 
end subroutine export_surface_gradients




!Dump csr matrix subroutine =========================
subroutine dump_csr_mat(mat,fh)
implicit none 

!Variables - Import 
integer(in32) :: fh
type(csrmatrix) :: mat

!Variables - Local 
integer(in64) :: ii 

!Dump data 
write(fh) mat%nrow
write(fh) mat%ncol
write(fh) mat%nnz
do ii=1,mat%nrow+1
    write(fh) mat%i(ii)
end do 
do ii=1,mat%nnz
    write(fh) mat%j(ii)
end do 
do ii=1,mat%nnz
    write(fh) mat%entry(ii)
end do 
return 
end subroutine dump_csr_mat




!Read csr matrix subroutine =========================
subroutine read_csr_mat(mat,fh)
implicit none 

!Variables - Import 
integer(in32) :: fh
type(csrmatrix) :: mat

!Variables - Local 
integer(in64) :: ii 

!Read data 
read(fh) mat%nrow
read(fh) mat%ncol
read(fh) mat%nnz
allocate(mat%i(mat%nrow+1))
do ii=1,mat%nrow+1
    read(fh) mat%i(ii)
end do 
allocate(mat%j(mat%nnz))
do ii=1,mat%nnz
    read(fh) mat%j(ii)
end do 
allocate(mat%entry(mat%nnz))
do ii=1,mat%nnz
    read(fh) mat%entry(ii)
end do 
return 
end subroutine read_csr_mat



!Dump integer vector subroutine =========================
subroutine dump_intvec_bin(intvec,fh)
implicit none 

!Variables - Import 
integer(in32) :: fh
integer(in64), dimension(:) :: intvec

!Variables - Local 
integer(in64) :: ii 
integer(in64) :: lenv

!Vector size
lenv = size(intvec)

!Write details 
write(fh) lenv
do ii=1,lenv
    write(fh) intvec(ii)
end do 
return 
end subroutine dump_intvec_bin




!Dump real array subroutine =========================
subroutine dump_realarr_bin(realvec,fh)
implicit none 

!Variables - Import 
integer(in32) :: fh
real(dp), dimension(:,:) :: realvec

!Variables - Local 
integer(in64) :: ii 
integer(in64) :: lenv,widv

!Vector size
lenv = size(realvec(:,1))
widv = size(realvec(1,:))

!Write details 
write(fh) lenv
write(fh) widv
do ii=1,lenv
    write(fh) realvec(ii,:)
end do 
return 
end subroutine dump_realarr_bin




!Read integer vector binary subroutine =========================
subroutine read_intvec_bin(intvec,fh)
implicit none 

!Variables - Import 
integer(in32) :: fh
integer(in64), dimension(:),allocatable :: intvec

!Variables - Local 
integer(in64) :: ii 
integer(in64) :: lenv

!Vector size
lenv = size(intvec)

!Write details 
read(fh) lenv
allocate(intvec(lenv))
do ii=1,lenv
    read(fh) intvec(ii)
end do 
return 
end subroutine read_intvec_bin




!Read read vector binary subroutine =========================
subroutine read_realarr_bin(realvec,fh)
implicit none 

!Variables - Import 
integer(in32) :: fh
real(dp), dimension(:,:), allocatable :: realvec

!Variables - Local 
integer(in64) :: ii 
integer(in64) :: lenv,widv

!Read details 
read(fh) lenv
read(fh) widv
allocate(realvec(lenv,widv))
do ii=1,lenv
    read(fh) realvec(ii,:)
end do 
return 
end subroutine read_realarr_bin




!Read read vector text subroutine =========================
subroutine read_realvec_txt(realvec,fh)
implicit none 

!Variables - Import 
integer(in32) :: fh
real(dp), dimension(:,:), allocatable :: realvec

!Variables - Local 
integer(in64) :: ii 
integer(in64) :: lenv,widv

!Read details 
read(fh,*) lenv
read(fh,*) widv
allocate(realvec(lenv,widv))
do ii=1,lenv
    read(fh,*) realvec(ii,:)
end do 
return 
end subroutine read_realvec_txt




!Construct data folder subroutine =========================
subroutine construct_data_folder(options)
implicit none 

!Variables - Import
type(options_data) :: options 

!Variables - Local 
logical :: dfexists

!Check if base file exists 
inquire(file=options%data_path,exist=dfexists)

!If it does not exist create the file 
if (.NOT.dfexists) then 
    call system('mkdir '//options%data_path)
end if 

!Construct surface files if they do not exist
! inquire(file=options%data_path//'\surfaces_base',exist=dfexists)
! if (.NOT.dfexists) then 
!     call system('mkdir '//options%data_path//'\surfaces_base')
! end if 
! inquire(file=options%data_path//'\surfaces_deformed',exist=dfexists)
! if (.NOT.dfexists) then 
!     call system('mkdir '//options%data_path//'\surfaces_deformed')
! end if 

!Construct matrix files if they do not exist
inquire(file=options%data_path//'/Ru',exist=dfexists)
if (.NOT.dfexists) then 
    call system('mkdir '//options%data_path//'/Ru')
end if 
! inquire(file=options%data_path//'\RuQ',exist=dfexists)
! if (.NOT.dfexists) then 
!     call system('mkdir '//options%data_path//'\RuQ')
! end if 
! inquire(file=options%data_path//'\Ri',exist=dfexists)
! if (.NOT.dfexists) then 
!     call system('mkdir '//options%data_path//'\Ri')
! end if 
! inquire(file=options%data_path//'\d',exist=dfexists)
! if (.NOT.dfexists) then 
!     call system('mkdir '//options%data_path//'\d')
! end if 
return 
end subroutine construct_data_folder


end module mrsys_io_mod