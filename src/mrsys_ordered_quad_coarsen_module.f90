!MRsys ordered quad coarsen module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.5
!Updated 13-12-2023

!Module
module mrsys_oqc_mod
use mrsys_connectivity_mod
contains 


!Ordered quad coarsen subroutine =========================
subroutine ordered_quad_coarsen(nrefine_sd,subd_system,mesh_base,options)
implicit none 

!Variables - Import
integer(in64) :: nrefine_sd
type(mesh_data) :: mesh_base
type(options_data) :: options 
type(subdsystem_data), dimension(:), allocatable :: subd_system

!Variables - Local 
integer(in64) :: cc,ff 
integer(in64) :: vtype_vtx,vtype_emid,vtype_fmid,coarsen_levels,meshlevel!,levexp
integer(in64), dimension(:), allocatable :: vbase,vtype
type(subdsystem_data), dimension(:), allocatable :: subd_system_coarsen

!Set vertex types for coarse mesh construction 
vtype_vtx = 1
vtype_emid = 2
vtype_fmid = 3

!Display
if (options%console_disp == 'yes') then
    write(*,'(A)') '--> coarsening the target surface mesh'
end if 

!Allocate the coarsening subdivision system 
allocate(subd_system_coarsen(options%n_crs_levels_max+1))

!Store base mesh in the finest level 
subd_system_coarsen(1)%mesh%ndim = mesh_base%ndim
subd_system_coarsen(1)%mesh%nedge = mesh_base%nedge
subd_system_coarsen(1)%mesh%nface = mesh_base%nface
subd_system_coarsen(1)%mesh%nvertex = mesh_base%nvertex
allocate(subd_system_coarsen(1)%mesh%vertices(mesh_base%nvertex,mesh_base%ndim))
allocate(subd_system_coarsen(1)%mesh%vertex_sharp(mesh_base%nvertex))
subd_system_coarsen(1)%mesh%vertices(:,:) = mesh_base%vertices(:,:)
subd_system_coarsen(1)%mesh%vertex_sharp(:) = mesh_base%vertex_sharp(:)
allocate(subd_system_coarsen(1)%mesh%edges(mesh_base%nedge,2))
allocate(subd_system_coarsen(1)%mesh%edge_sharp(mesh_base%nedge))
subd_system_coarsen(1)%mesh%edges(:,:) = mesh_base%edges(:,:)
subd_system_coarsen(1)%mesh%edge_sharp(:) = mesh_base%edge_sharp(:)
allocate(subd_system_coarsen(1)%mesh%faces(mesh_base%nface))
do ff=1,mesh_base%nface
    subd_system_coarsen(1)%mesh%faces(ff)%nvertex = mesh_base%faces(ff)%nvertex
    allocate(subd_system_coarsen(1)%mesh%faces(ff)%vertices(mesh_base%faces(ff)%nvertex))
    subd_system_coarsen(1)%mesh%faces(ff)%vertices(:) = mesh_base%faces(ff)%vertices(:)
end do 
do cc=2,options%n_crs_levels_max+1
    subd_system_coarsen(cc)%mesh%ndim = mesh_base%ndim
    subd_system_coarsen(cc)%mesh%nedge = 0
    subd_system_coarsen(cc)%mesh%nface = 0
    subd_system_coarsen(cc)%mesh%nvertex = 0
end do 

!Evaluate connectivity of the finest level 
call get_valence(subd_system_coarsen(1)%mesh)
call get_connectivity(subd_system_coarsen(1)%mesh%connectivity,subd_system_coarsen(1)%mesh,1.0d0)

!Coarsen the mesh 
coarsen_levels = 0 
do cc=1,options%n_crs_levels_max

    !Identify base vertices for coarsening 
    if (mesh_base%ndim == 2) then 
        call identify_coarsening_base_vertices_2D(vbase,subd_system_coarsen(cc)%mesh,options) 
    elseif (mesh_base%ndim == 3) then 
        call identify_coarsening_base_vertices_3D(vbase,subd_system_coarsen(cc)%mesh,options) 
    end if 

    !Flood vertex types through the mesh 
    if (mesh_base%ndim == 2) then 
        call flood_vertex_types_2D(vtype,vbase,subd_system_coarsen(cc)%mesh,vtype_vtx,vtype_emid)
    elseif (mesh_base%ndim == 3) then 
        call flood_vertex_types_3D(vtype,vbase,subd_system_coarsen(cc)%mesh,vtype_vtx,vtype_emid,vtype_fmid)
    end if 

    !Check if mesh can be coarsened one level 
    if (.NOT.can_be_coarsend(subd_system_coarsen(cc)%mesh,vtype,vtype_vtx,vtype_emid,vtype_fmid)) then 
        if (options%console_disp == 'yes') then
            write(*,'(A)') '    <** mesh coarsening limit reached **> '
        end if
        exit 
    end if 
    
    !Construct coarse level mesh and its connectivity if mesh can be coarsened 
    if (mesh_base%ndim == 2) then 
        call build_coarse_mesh_2D(subd_system_coarsen(cc+1)%mesh,subd_system_coarsen(cc)%mesh,vtype,vtype_vtx,vtype_emid)
    elseif (mesh_base%ndim == 3) then 
        call build_coarse_mesh_3D(subd_system_coarsen(cc+1)%mesh,subd_system_coarsen(cc)%mesh,vtype,vtype_vtx,vtype_emid,vtype_fmid)
    end if 

    !Increment number of coarsening levels 
    coarsen_levels = coarsen_levels + 1

    !Display
    if (options%console_disp == 'yes') then
        write(*,'(A,I0,A,I0,A,I0,A)') '    < coarsening = ',cc,' - vertices = ',subd_system_coarsen(cc)%mesh%nvertex,&
        ' -> ',subd_system_coarsen(cc+1)%mesh%nvertex,' >'
    end if
end do 

!Re-order vertices at each level to the catmull-clark structure 
call reorder_vertices_2cc(subd_system_coarsen,coarsen_levels,vtype_vtx,vtype_emid,vtype_fmid)

!Allocate the extracted subdivision system structure 
allocate(subd_system(coarsen_levels+1))

!Extract meshes to the return subdivision structure 
do cc=1,coarsen_levels+1

    !Level of mesh to extract (coarsest first)
    meshlevel = coarsen_levels - cc + 2

    !Extract mesh 
    subd_system(cc)%mesh%ndim = subd_system_coarsen(meshlevel)%mesh%ndim 
    subd_system(cc)%mesh%nedge = subd_system_coarsen(meshlevel)%mesh%nedge
    subd_system(cc)%mesh%nface = subd_system_coarsen(meshlevel)%mesh%nface
    subd_system(cc)%mesh%nvertex = subd_system_coarsen(meshlevel)%mesh%nvertex 
    allocate(subd_system(cc)%mesh%vertices(subd_system(cc)%mesh%nvertex,subd_system(cc)%mesh%ndim))
    allocate(subd_system(cc)%mesh%vertex_sharp(subd_system(cc)%mesh%nvertex))
    subd_system(cc)%mesh%vertices(:,:) = subd_system_coarsen(meshlevel)%mesh%vertices(:,:)
    subd_system(cc)%mesh%vertex_sharp(:) = subd_system_coarsen(meshlevel)%mesh%vertex_sharp(:)
    allocate(subd_system(cc)%mesh%edges(subd_system(cc)%mesh%nedge,2))
    allocate(subd_system(cc)%mesh%edge_sharp(subd_system(cc)%mesh%nedge))
    subd_system(cc)%mesh%edges(:,:) = subd_system_coarsen(meshlevel)%mesh%edges(:,:)
    subd_system(cc)%mesh%edge_sharp(:) = subd_system_coarsen(meshlevel)%mesh%edge_sharp(:)
    allocate(subd_system(cc)%mesh%faces(subd_system(cc)%mesh%nface))
    do ff=1,subd_system(cc)%mesh%nface
        subd_system(cc)%mesh%faces(ff)%nvertex = subd_system_coarsen(meshlevel)%mesh%faces(ff)%nvertex
        allocate(subd_system(cc)%mesh%faces(ff)%vertices(subd_system(cc)%mesh%faces(ff)%nvertex))
        subd_system(cc)%mesh%faces(ff)%vertices(:) = subd_system_coarsen(meshlevel)%mesh%faces(ff)%vertices(:)
    end do 
end do 

!Set the number of refinement levels in the resulting subdivision system 
nrefine_sd = coarsen_levels

! !DEBUG =======================
! !Test export vertex types 
! levexp = 1
! open(11,file='io\vtype1')
! do cc=1,subd_system_coarsen(levexp)%mesh%nvertex
!     if (vtype(cc) == 1) then 
!         write(11,*) subd_system_coarsen(levexp)%mesh%vertices(cc,:)
!     end if 
! end do 
! close(11)
! open(11,file='io\vtype2')
! do cc=1,subd_system_coarsen(levexp)%mesh%nvertex
!     if (vtype(cc) == 2) then 
!         write(11,*) subd_system_coarsen(levexp)%mesh%vertices(cc,:)
!     end if 
! end do 
! close(11)
! open(11,file='io\vtype3')
! do cc=1,subd_system_coarsen(levexp)%mesh%nvertex
!     if (vtype(cc) == 3) then 
!         write(11,*) subd_system_coarsen(levexp)%mesh%vertices(cc,:)
!     end if 
! end do 
! close(11)
return 
end subroutine ordered_quad_coarsen




!Re-order vertices subroutine =========================
subroutine reorder_vertices_2cc(subd_system_coarsen,coarsen_levels,vtype_vtx,vtype_emid,vtype_fmid)
implicit none

!Variables - Import
integer(in64) :: coarsen_levels
integer(in64) :: vtype_vtx,vtype_emid,vtype_fmid
type(subdsystem_data), dimension(:) :: subd_system_coarsen

!Variables - Local 
integer(in64) :: rr,vv,ee,ff,aa
integer(in64) :: level_coarse,level_fine,vtx_idx_crs,nvertex_crs,nedge_crs
integer(in64), dimension(:), allocatable :: vertex_sharp_temp
real(dp), dimension(:,:), allocatable :: vertices_temp
type(vtxodrdata) :: vtx_order(coarsen_levels+1)

!Reorder vertices at the finer level of each coarsening to correspond to the base vertex, face and edge indexing of the coarser level 
do rr=1,coarsen_levels

    !Index of the coarse and fine level 
    level_coarse = coarsen_levels - rr + 2
    level_fine = level_coarse - 1

    !Initialise the vertex mapping on the coarsest level to show no change 
    if (rr == 1) then 
        allocate(vtx_order(level_coarse)%vtx_idx_odr(subd_system_coarsen(level_coarse)%mesh%nvertex))
        do vv=1,subd_system_coarsen(level_coarse)%mesh%nvertex
            vtx_order(level_coarse)%vtx_idx_odr(vv) = vv
        end do 
    end if 

    !Number of vertices and edges in the coarse level 
    nvertex_crs = subd_system_coarsen(level_coarse)%mesh%nvertex
    nedge_crs = subd_system_coarsen(level_coarse)%mesh%nedge

    !Find new poisition for each vertex on the fine level (orignal vertices -> new edge vertices -> new face vertices) in order of the parent edge/face index
    allocate(vtx_order(level_fine)%vtx_idx_odr(subd_system_coarsen(level_fine)%mesh%nvertex))
    vtx_order(level_fine)%vtx_idx_odr(:) = 0 
    do vv=1,subd_system_coarsen(level_fine)%mesh%nvertex !Original vertices from the coarse level -> set to the same position as in the coarse level
        if (subd_system_coarsen(level_fine)%mesh%vertex_felink2crs(vv,2) == vtype_vtx) then 
            vtx_idx_crs = subd_system_coarsen(level_fine)%mesh%vtxmap_fin2crs(vv)
            vtx_order(level_fine)%vtx_idx_odr(vv) = vtx_order(level_coarse)%vtx_idx_odr(vtx_idx_crs)
        end if 
    end do 
    do vv=1,subd_system_coarsen(level_fine)%mesh%nvertex !Edge midpoint vertices-> set to the number of vertices in the coarse level + the index of the parent edge 
        if (subd_system_coarsen(level_fine)%mesh%vertex_felink2crs(vv,2) == vtype_emid) then 
            vtx_order(level_fine)%vtx_idx_odr(vv) = nvertex_crs + subd_system_coarsen(level_fine)%mesh%vertex_felink2crs(vv,1) 
        end if 
    end do 
    do vv=1,subd_system_coarsen(level_fine)%mesh%nvertex !Face midpoint vertices-> set to the number of vertices in the coarse level + to the number of edges in the coarse level + the index of the parent face 
        if (subd_system_coarsen(level_fine)%mesh%vertex_felink2crs(vv,2) == vtype_fmid) then 
            vtx_order(level_fine)%vtx_idx_odr(vv) = nvertex_crs + nedge_crs + &
            subd_system_coarsen(level_fine)%mesh%vertex_felink2crs(vv,1) 
        end if 
    end do 

    !Update the refined level vertex positions 
    allocate(vertices_temp(subd_system_coarsen(level_fine)%mesh%nvertex,subd_system_coarsen(level_fine)%mesh%ndim))
    vertices_temp(:,:) = subd_system_coarsen(level_fine)%mesh%vertices(:,:)
    allocate(vertex_sharp_temp(subd_system_coarsen(level_fine)%mesh%nvertex))
    vertex_sharp_temp(:) = subd_system_coarsen(level_fine)%mesh%vertex_sharp(:)
    do vv=1,subd_system_coarsen(level_fine)%mesh%nvertex 
        subd_system_coarsen(level_fine)%mesh%vertices(vtx_order(level_fine)%vtx_idx_odr(vv),:) = vertices_temp(vv,:)
        subd_system_coarsen(level_fine)%mesh%vertex_sharp(vtx_order(level_fine)%vtx_idx_odr(vv)) = vertex_sharp_temp(vv)
    end do 
    deallocate(vertices_temp)
    deallocate(vertex_sharp_temp)

    !Re-map the refined level edges and faces 
    do ee=1,subd_system_coarsen(level_fine)%mesh%nedge
        subd_system_coarsen(level_fine)%mesh%edges(ee,1) = &
        vtx_order(level_fine)%vtx_idx_odr(subd_system_coarsen(level_fine)%mesh%edges(ee,1))
        subd_system_coarsen(level_fine)%mesh%edges(ee,2) = &
        vtx_order(level_fine)%vtx_idx_odr(subd_system_coarsen(level_fine)%mesh%edges(ee,2))
    end do 
    do ff=1,subd_system_coarsen(level_fine)%mesh%nface
        do aa=1,subd_system_coarsen(level_fine)%mesh%faces(ff)%nvertex
            subd_system_coarsen(level_fine)%mesh%faces(ff)%vertices(aa) = &
            vtx_order(level_fine)%vtx_idx_odr(subd_system_coarsen(level_fine)%mesh%faces(ff)%vertices(aa))
        end do 
    end do 
end do 
return 
end subroutine reorder_vertices_2cc




!Check if mesh can be coarsened =========================
function can_be_coarsend(mesh,vtype,vtype_vtx,vtype_emid,vtype_fmid) result(can_coarsen)
implicit none 

!Variables - Import
logical :: can_coarsen
integer(in64) :: vtype_vtx,vtype_emid,vtype_fmid
integer(in64), dimension(:) :: vtype
type(mesh_data) :: mesh

!Variables - Local 
integer(in64) :: ff,ee,vv
integer(in64) :: v1,v2,has_vfmid,vtgt,fm_minvalence,ncrs_edge
integer(in64) :: fcv_valence(mesh%nvertex)

!Set base state
can_coarsen = .true.

!Check for faces that do not contain any face centre vertices 
if (mesh%ndim == 3) then 
    do ff=1,mesh%nface
        has_vfmid = 0 
        do vv=1,mesh%faces(ff)%nvertex
            if (vtype(mesh%faces(ff)%vertices(vv)) == vtype_fmid) then 
                has_vfmid = 1
                exit 
            end if 
        end do 
        if (has_vfmid == 0) then !Cannot be coarsened 
            can_coarsen = .false.
            return 
        end if 
    end do 
end if 

!Check for edges that connect two vertices that are retained on the coarse level 
do ee=1,mesh%nedge
    v1 = mesh%edges(ee,1)
    v2 = mesh%edges(ee,2)
    if ((vtype(v1) == vtype_vtx) .AND. (vtype(v2) == vtype_vtx)) then !Cannot be coarsened 
        can_coarsen = .false.
        return 
    end if 
end do 

!Check for edges that connect two edge mid vertices
do ee=1,mesh%nedge
    v1 = mesh%edges(ee,1)
    v2 = mesh%edges(ee,2)
    if ((vtype(v1) == vtype_emid) .AND. (vtype(v2) == vtype_emid)) then !Cannot be coarsened 
        can_coarsen = .false.
        return 
    end if 
end do 

!Check for edges that connect two face centre vertices
if (mesh%ndim == 3) then 
    do ee=1,mesh%nedge
        v1 = mesh%edges(ee,1)
        v2 = mesh%edges(ee,2)
        if ((vtype(v1) == vtype_fmid) .AND. (vtype(v2) == vtype_fmid)) then !Cannot be coarsened 
            can_coarsen = .false.
            return 
        end if 
    end do 
end if 

!Check face valence of all face centre vertices for any less than three
if (mesh%ndim == 3) then 
    fcv_valence(:) = 0 
    do ff=1,mesh%nface
        do vv=1,mesh%faces(ff)%nvertex
            vtgt = mesh%faces(ff)%vertices(vv)
            if (vtype(vtgt) == vtype_fmid) then 
                fcv_valence(vtgt) = fcv_valence(vtgt) + 1
            end if 
        end do 
    end do 
    fm_minvalence = 10*mesh%nvertex
    do vv=1,mesh%nvertex
        if (vtype(vv) == vtype_fmid) then 
            if (fcv_valence(vv) .LT. fm_minvalence) then 
                fm_minvalence = fcv_valence(vv)
            end if 
        end if 
    end do 
    if (fm_minvalence .LT. 3) then
        can_coarsen = .false.
        return 
    end if 
end if 

!Check if any vertices have not been tagged with a type 
if (minval(vtype) == 0) then 
    can_coarsen = .false.
    return 
end if 

!Check if the resulting number of mesh edges is less than three
if (mesh%ndim == 2) then 
    ncrs_edge = 0 
    do vv=1,mesh%nvertex
        if (vtype(vv) == vtype_emid) then 
            ncrs_edge = ncrs_edge + 1
        end if 
    end do 
    if (ncrs_edge .LT. 3) then 
        can_coarsen = .false.
        return 
    end if 
end if 
return 
end function can_be_coarsend




!Construct coarse level mesh subroutine 2D =========================
subroutine build_coarse_mesh_2D(mesh_coarse,mesh_fine,vtype,vtype_vtx,vtype_emid)
implicit none 

!Variables - Import
integer(in64), dimension(:) :: vtype
integer(in64) :: vtype_vtx,vtype_emid
type(mesh_data) :: mesh_coarse,mesh_fine

!Variables - Local 
integer(in64) :: vv,ee,aa
integer(in64) :: eins,vins,etgt,edge_1,edge_2,vtx_1,vtx_2

!Find number of coarse level vertices 
mesh_coarse%nvertex = 0 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_vtx) then 
        mesh_coarse%nvertex = mesh_coarse%nvertex + 1
    end if 
end do 
allocate(mesh_coarse%vertices(mesh_coarse%nvertex,mesh_fine%ndim))
allocate(mesh_coarse%vertex_sharp(mesh_coarse%nvertex))

!Find number of coarse level edges 
mesh_coarse%nedge = 0 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_emid) then 
        mesh_coarse%nedge = mesh_coarse%nedge + 1
    end if 
end do 
allocate(mesh_coarse%edges(mesh_coarse%nedge,2))
allocate(mesh_coarse%edge_sharp(mesh_coarse%nedge))

!Set number of coarse level faces
mesh_coarse%nface = 0 

!Extract and map coarse level vertices from the fine level mesh 
allocate(mesh_coarse%vtxmap_crs2fin(mesh_coarse%nvertex))
allocate(mesh_fine%vtxmap_fin2crs(mesh_fine%nvertex))
mesh_coarse%vtxmap_crs2fin(:) = 0 
mesh_fine%vtxmap_fin2crs(:) = 0 
vins = 0 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_vtx) then 
        vins = vins + 1
        mesh_coarse%vertices(vins,:) = mesh_fine%vertices(vv,:)
        mesh_coarse%vertex_sharp(vins) = mesh_fine%vertex_sharp(vv)
        mesh_coarse%vtxmap_crs2fin(vins) = vv 
        mesh_fine%vtxmap_fin2crs(vv) = vins 
    end if 
end do 

!Allocate face and edge links from the fine mesh to the corse mesh 
allocate(mesh_fine%vertex_felink2crs(mesh_fine%nvertex,2))
mesh_fine%vertex_felink2crs(:,:) = 0 
mesh_fine%vertex_felink2crs(:,2) = vtype(:)

!Construct coarse level edges 
eins = 0 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_emid) then 

        !Find start vertex of this coarse edge (fine edge with vv as end 2)
        edge_1 = 0 
        do aa=1,mesh_fine%connectivity%valence(vv)
            etgt = mesh_fine%connectivity%v2e(vv,aa)
            if (mesh_fine%edges(etgt,2) == vv) then 
                edge_1 = etgt 
                exit 
            end if 
        end do 
        vtx_1 = mesh_fine%edges(edge_1,1)

        !Find end vertex of this edge (fine edge with vv as end 1)
        edge_2 = 0 
        do aa=1,mesh_fine%connectivity%valence(vv)
            etgt = mesh_fine%connectivity%v2e(vv,aa)
            if (mesh_fine%edges(etgt,1) == vv) then 
                edge_2 = etgt 
                exit 
            end if 
        end do 
        vtx_2 = mesh_fine%edges(edge_2,2)

        !Construct edge and assign sharpness 
        eins = eins + 1
        mesh_coarse%edges(eins,1) = vtx_1
        mesh_coarse%edges(eins,2) = vtx_2
        mesh_coarse%edge_sharp(eins) = max(mesh_fine%edge_sharp(edge_1),mesh_fine%edge_sharp(edge_2))

        !Assign coarse level edge index to the fine level vertex at this coarse level edge midpoint 
        mesh_fine%vertex_felink2crs(vv,1) = eins
    end if 
end do 

!Map edges to coarse mesh indecies 
do ee=1,mesh_coarse%nedge
    mesh_coarse%edges(ee,1) = mesh_fine%vtxmap_fin2crs(mesh_coarse%edges(ee,1))
    mesh_coarse%edges(ee,2) = mesh_fine%vtxmap_fin2crs(mesh_coarse%edges(ee,2))
end do 

!Build coarse level connectivity 
call get_valence(mesh_coarse)
call get_connectivity(mesh_coarse%connectivity,mesh_coarse,1.0d0)
return 
end subroutine build_coarse_mesh_2D




!Construct coarse level mesh subroutine 3D =========================
subroutine build_coarse_mesh_3D(mesh_coarse,mesh_fine,vtype,vtype_vtx,vtype_emid,vtype_fmid)
implicit none 

!Variables - Import
integer(in64), dimension(:) :: vtype
integer(in64) :: vtype_vtx,vtype_emid,vtype_fmid
type(mesh_data) :: mesh_coarse,mesh_fine

!Variables - Local 
integer(in64) :: ii,vv,ee,ff,aa
integer(in64) :: fins,vins,nvface,fpi,fvtx,etgt,v1crs,v2crs,eparent
integer(in64) :: vfcentre,vstart,vstart_fi,fstart,fcurrent,ecurrent,vcurrent,vcurrent_fi,vf,vtgt
integer(in64) :: face_temp(mesh_fine%connectivity%max_valence),vretatt(2)
real(dp) :: eattsharp(2)

!Find number of coarse level vertices 
mesh_coarse%nvertex = 0 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_vtx) then 
        mesh_coarse%nvertex = mesh_coarse%nvertex + 1
    end if 
end do 
allocate(mesh_coarse%vertices(mesh_coarse%nvertex,mesh_fine%ndim))
allocate(mesh_coarse%vertex_sharp(mesh_coarse%nvertex))

!Find number of coarse level faces 
mesh_coarse%nface = 0 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_fmid) then 
        mesh_coarse%nface = mesh_coarse%nface + 1
    end if 
end do 
allocate(mesh_coarse%faces(mesh_coarse%nface))

!Extract and map coarse level vertices from the fine level mesh 
allocate(mesh_coarse%vtxmap_crs2fin(mesh_coarse%nvertex))
allocate(mesh_fine%vtxmap_fin2crs(mesh_fine%nvertex))
mesh_coarse%vtxmap_crs2fin(:) = 0 
mesh_fine%vtxmap_fin2crs(:) = 0 
vins = 0 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_vtx) then 
        vins = vins + 1
        mesh_coarse%vertices(vins,:) = mesh_fine%vertices(vv,:)
        mesh_coarse%vertex_sharp(vins) = mesh_fine%vertex_sharp(vv)
        mesh_coarse%vtxmap_crs2fin(vins) = vv 
        mesh_fine%vtxmap_fin2crs(vv) = vins 
    end if 
end do 

!Allocate face and edge links from the fine mesh to the corse mesh 
allocate(mesh_fine%vertex_felink2crs(mesh_fine%nvertex,2))
mesh_fine%vertex_felink2crs(:,:) = 0 
mesh_fine%vertex_felink2crs(:,2) = vtype(:)

!Construct coarse level faces
fins = 0 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_fmid) then 

        !Face centre vertex
        vfcentre = vv

        !Select starting vertex for this face
        vstart = 0
        vstart_fi = 0
        fstart = mesh_fine%connectivity%v2f(vfcentre,1)
        do ii=1,mesh_fine%faces(fstart)%nvertex 
            vf = mesh_fine%faces(fstart)%vertices(ii)
            if (vtype(vf) == vtype_vtx) then 
                vstart_fi = ii
                vstart = vf
                exit
            end if 
        end do 

        !Construct face
        nvface = 1 
        face_temp(:) = 0 
        face_temp(1) = vstart
        fcurrent = fstart
        vcurrent_fi = vstart_fi
        do ee=1,mesh_fine%connectivity%valence(vfcentre)

            !Find edge in fcurrent subseqent to vcurrent that contains vfcentre    
            ecurrent = 0
            fpi = vcurrent_fi
            do ii=1,mesh_fine%faces(fcurrent)%nvertex

                !Index of current position in face
                fvtx = mod(fpi-1,mesh_fine%faces(fcurrent)%nvertex) + 1

                !Edge
                etgt = mesh_fine%faces(fcurrent)%edges(fvtx) 
                if ((mesh_fine%edges(etgt,1) == vfcentre) .OR. (mesh_fine%edges(etgt,2) == vfcentre)) then 
                    ecurrent = etgt
                    exit
                end if 
                
                !Increment face vertex index position
                fpi = fpi + 1
            end do 

            !Find next face on this edge and update this to the current face
            if (mesh_fine%connectivity%e2f(ecurrent,1) == fcurrent) then 
                fcurrent = mesh_fine%connectivity%e2f(ecurrent,2)
            else
                fcurrent = mesh_fine%connectivity%e2f(ecurrent,1)
            end if 

            !Find type vtype_vtx vertex in this face as the next vertex for this face
            vcurrent = 0
            vcurrent_fi = 0
            do ii=1,mesh_fine%faces(fcurrent)%nvertex
                vf = mesh_fine%faces(fcurrent)%vertices(ii)
                if (vtype(vf) == vtype_vtx) then 
                    vcurrent_fi = ii
                    vcurrent = vf
                    exit
                end if 
            end do

            !Exit if vcurrent == vstart
            if (vcurrent == vstart) then 
                exit 
            end if 

            !Add this vertex to this face
            nvface = nvface + 1
            face_temp(nvface) = vcurrent
        end do 

        !Store face 
        fins = fins + 1
        mesh_coarse%faces(fins)%nvertex = nvface
        allocate(mesh_coarse%faces(fins)%vertices(nvface))
        mesh_coarse%faces(fins)%vertices(:) = face_temp(1:nvface)
        ! print *, mesh_coarse%faces(fins)%nvertex
        ! print *, mesh_coarse%faces(fins)%vertices

        !Assign coarse level face index to the fine level vertex at this coarse level faces centre 
        mesh_fine%vertex_felink2crs(vfcentre,1) = fins
    end if 
end do 

!Map faces to coarse mesh indecies 
do ff=1,mesh_coarse%nface
    do ii=1,mesh_coarse%faces(ff)%nvertex
        mesh_coarse%faces(ff)%vertices(ii) = mesh_fine%vtxmap_fin2crs(mesh_coarse%faces(ff)%vertices(ii))
    end do 
end do 

!Construct coarse level edges 
call get_valence(mesh_coarse)
call build_edges_from_faces(mesh_coarse,mesh_coarse%connectivity%max_valence)

!Build coarse level connectivity 
call get_connectivity(mesh_coarse%connectivity,mesh_coarse,1.0d0)

!Assign coarse level edge index to the fine level vertex at each coarse level edge midpoint 
do vv=1,mesh_fine%nvertex
    if (vtype(vv) == vtype_emid) then 

        !Find the two attached vtype_vtx vertices
        vins = 0
        vretatt(:) = 0
        do aa=1,mesh_fine%connectivity%valence(vv)
            if (vtype(mesh_fine%connectivity%v2v(vv,aa)) == vtype_vtx) then 
                vins = vins + 1
                vretatt(vins) = mesh_fine%connectivity%v2v(vv,aa)
                if (vins == 2) then 
                    exit
                end if 
            end if 
        end do 
        
        !Map these vertices to the coarser level 
        v1crs = mesh_fine%vtxmap_fin2crs(vretatt(1))
        v2crs = mesh_fine%vtxmap_fin2crs(vretatt(2))

        !Identify the edge that connects these two vertices on the coarse level 
        eparent = 0
        do aa=1,mesh_coarse%connectivity%valence(v1crs)
            etgt = mesh_coarse%connectivity%v2e(v1crs,aa)
            if ((mesh_coarse%edges(etgt,1) == v2crs) .OR. (mesh_coarse%edges(etgt,2) == v2crs)) then 
                eparent = etgt
                exit 
            end if 
        end do 

        !Store
        mesh_fine%vertex_felink2crs(vv,1) = eparent

        !Assign sharpness to the coarse level edge  mesh%edge_sharp
        vins = 0
        eattsharp(:) = 0.0d0 
        do aa=1,mesh_fine%connectivity%valence(vv)
            etgt = mesh_fine%connectivity%v2e(vv,aa)
            if (mesh_fine%edges(etgt,1) == vv) then 
                vtgt = mesh_fine%edges(etgt,2)
            else
                vtgt = mesh_fine%edges(etgt,1)
            end if 
            if ((vtgt == vretatt(1)) .OR. (vtgt == vretatt(2))) then 
                vins = vins + 1
                eattsharp(vins) = mesh_fine%edge_sharp(etgt)
                if (vins == 2) then 
                    exit 
                end if 
            end if 
        end do 
        mesh_coarse%edge_sharp(eparent) = maxval(eattsharp)
    end if 
end do 
return 
end subroutine build_coarse_mesh_3D




!Flood vertex types subroutine 2D =========================
subroutine flood_vertex_types_2D(vtype,vbase,mesh,vtype_vtx,vtype_emid)
implicit none 

!Variables - Import
integer(in64) :: vtype_vtx,vtype_emid
integer(in64), dimension(:) :: vbase
integer(in64), dimension(:), allocatable :: vtype
type(mesh_data) :: mesh

!Variables - Local 
integer(in64) :: rr,ee
integer(in64) :: v1,v2,Nupdate

!Reset
if (allocated(vtype)) then 
    deallocate(vtype)
end if 
allocate(vtype(mesh%nvertex))
vtype(:) = 0 

!Set state of base vertices 
vtype(vbase) = vtype_vtx

!Flood states
do rr=1,2*mesh%nvertex

    !Set number of updated vertices to zero 
    Nupdate = 0 

    !Set all vtype 0 vertices connected to a vtype_vtx vertex by an edge as vtype_emid
    do ee=1,mesh%nedge
        
        !Edge ends 
        v1 = mesh%edges(ee,1)
        v2 = mesh%edges(ee,2)

        !Flood 
        if ((vtype(v1) == 0) .AND. (vtype(v2) == vtype_vtx)) then 
            vtype(v1) = vtype_emid
            Nupdate = Nupdate + 1
            cycle 
        elseif ((vtype(v2) == 0) .AND. (vtype(v1) == vtype_vtx)) then 
            vtype(v2) = vtype_emid
            Nupdate = Nupdate + 1
            cycle 
        end if 
    end do 

    !Set all vtype 0 vertices connected to a vtype_emid vertex by an edge as vtype_vtx
    do ee=1,mesh%nedge
        
        !Edge ends 
        v1 = mesh%edges(ee,1)
        v2 = mesh%edges(ee,2)

        !Flood 
        if ((vtype(v1) == 0) .AND. (vtype(v2) == vtype_emid)) then 
            vtype(v1) = vtype_vtx
            Nupdate = Nupdate + 1
            cycle 
        elseif ((vtype(v2) == 0) .AND. (vtype(v1) == vtype_emid)) then 
            vtype(v2) = vtype_vtx
            Nupdate = Nupdate + 1
            cycle 
        end if 
    end do 

    !Exit at completed flood 
    if (Nupdate == 0) then 
        exit 
    end if 
    !print *, Nupdate
end do 
return 
end subroutine flood_vertex_types_2D




!Flood vertex types subroutine 3D =========================
subroutine flood_vertex_types_3D(vtype,vbase,mesh,vtype_vtx,vtype_emid,vtype_fmid)
implicit none 

!Variables - Import
integer(in64) :: vtype_vtx,vtype_emid,vtype_fmid
integer(in64), dimension(:) :: vbase
integer(in64), dimension(:), allocatable :: vtype
type(mesh_data) :: mesh

!Variables - Local 
integer(in64) :: rr,ff,vv,aa
integer(in64) :: vtgt,has_vtype_vtx,connected2vtype_vtx,has_vtype_fmid,connected2vtype_fmid,Nupdate
integer(in64) :: face_front(mesh%nface),face_front2(mesh%nface)

!Reset
if (allocated(vtype)) then 
    deallocate(vtype)
end if 
allocate(vtype(mesh%nvertex))
vtype(:) = 0 

!Set state of base vertices 
vtype(vbase) = vtype_vtx

!Flood states
face_front(:) = 0 
face_front2(:) = 0 
do rr=1,2*mesh%nvertex

    !Set number of updated vertices to zero 
    Nupdate = 0 

    !Set all vtype 0 vertices connected to a vtype_vtx vertex by an edge as vtype_emid
    do vv=1,mesh%nvertex
        if (vtype(vv) == vtype_vtx) then 
            do aa=1,mesh%connectivity%valence(vv) 
                if (vtype(mesh%connectivity%v2v(vv,aa)) == 0) then !Tag vertices
                    vtype(mesh%connectivity%v2v(vv,aa)) = vtype_emid
                    Nupdate = Nupdate + 1
                end if 
                if (mesh%connectivity%v2f(vv,aa) .GT. 0) then !Build search front of faces
                    face_front(mesh%connectivity%v2f(vv,aa)) = 1 
                end if 
            end do 
        end if 
    end do 

    !Set all vtype 0 vertices connected to a vtype_vtx vertex by a face and not an edge as vtype_fmid
    do ff=1,mesh%nface
        if (face_front(ff) == 1) then !If face on front boundary 

            !Reset this face front tag
            face_front(ff) = 0 
            
            !Check if this face contains a vtype_vtx vertex
            has_vtype_vtx = 0 
            do vv=1,mesh%faces(ff)%nvertex
                if (vtype(mesh%faces(ff)%vertices(vv)) == vtype_vtx) then 
                    has_vtype_vtx = 1
                    exit 
                end if 
            end do 
            if (has_vtype_vtx == 0) then 
                cycle
            end if 

            !If it contains a vtype_vtx vertex then set vtype 0 vertices connected to a vtype_vtx vertex by this face and not an edge as vtype_fmid
            do vv=1,mesh%faces(ff)%nvertex
                if (vtype(mesh%faces(ff)%vertices(vv)) == 0) then 

                    !Target vertex 
                    vtgt = mesh%faces(ff)%vertices(vv)

                    !Check if connected to a vtype_vtx by an edge 
                    connected2vtype_vtx = 0 
                    do aa=1,mesh%connectivity%valence(vtgt) 
                        if (vtype(mesh%connectivity%v2v(vtgt,aa)) == vtype_vtx) then  
                            connected2vtype_vtx = 1
                        end if 
                    end do 
                    if (connected2vtype_vtx == 1) then 
                        cycle 
                    end if 

                    !If not connected to a vtype_vtx by an edge then set as vtype_fmid vertex
                    vtype(vtgt) = vtype_fmid
                    Nupdate = Nupdate + 1

                    !Add faces on this vertex to the second face front
                    do aa=1,mesh%connectivity%valence(vtgt) 
                        if (mesh%connectivity%v2f(vtgt,aa) .GT. 0) then 
                            face_front2(mesh%connectivity%v2f(vtgt,aa)) = 1 
                        end if 
                    end do 
                end if 
            end do 
        end if 
    end do 

    !Set any vtype 0 vertex connected to a vtype_fmid vertex by a face and not an edge as vtype_vtx
    do ff=1,mesh%nface
        if (face_front2(ff) == 1) then !If face on front boundary 

            !Reset this face front tag
            face_front2(ff) = 0 

            !Check if this face contains a vtype_fmid vertex
            has_vtype_fmid = 0 
            do vv=1,mesh%faces(ff)%nvertex
                if (vtype(mesh%faces(ff)%vertices(vv)) == vtype_fmid) then 
                    has_vtype_fmid = 1
                    exit 
                end if 
            end do 
            if (has_vtype_fmid == 0) then 
                cycle
            end if 

            !If it contains a vtype_vtx vertex then set vtype 0 vertices connected to a vtype_fmid vertex by this face and not an edge as vtype_vtx
            do vv=1,mesh%faces(ff)%nvertex
                if (vtype(mesh%faces(ff)%vertices(vv)) == 0) then 

                    !Target vertex 
                    vtgt = mesh%faces(ff)%vertices(vv)

                    !Check if connected to a vtype_fmid by an edge 
                    connected2vtype_fmid = 0 
                    do aa=1,mesh%connectivity%valence(vtgt) 
                        if (vtype(mesh%connectivity%v2v(vtgt,aa)) == vtype_fmid) then  
                            connected2vtype_fmid = 1
                        end if 
                    end do 
                    if (connected2vtype_fmid == 1) then 
                        cycle 
                    end if 

                    !If not connected to a vtype_fmid by an edge then set as a vtype_vtx
                    vtype(vtgt) = vtype_vtx
                    Nupdate = Nupdate + 1
                end if 
            end do 
        end if 
    end do 

    !Exit at completed flood 
    if (Nupdate == 0) then 
        exit 
    end if 
    !print *, Nupdate
end do 
return 
end subroutine flood_vertex_types_3D




!Identify base vertex for coarsening 2D =========================
subroutine identify_coarsening_base_vertices_2D(vbase,mesh,options) 
implicit none 

!Variables - Import
integer(in64), dimension(:), allocatable :: vbase
type(mesh_data) :: mesh
type(options_data) :: options 

!Variables - Local 
integer(in64) :: vv
integer(in64) :: nexrtvtx
integer(in64) :: maxvalence,minvalence

!Initialise
if (allocated(vbase)) then 
    deallocate(vbase)
end if 

!Find valence bounds
maxvalence = maxval(mesh%connectivity%valence)
minvalence = minval(mesh%connectivity%valence)

!Return if maxvalence < 2 as mesh cannot be coarsened further
if (maxvalence .LT. 2) then 
    allocate(vbase(1))
    vbase(1) = -1
    return 
end if 

!Force auto setting as 2D
if (options%crs_base_vlnc == 'gt4') then 
    options%crs_base_vlnc = 'auto'
    if (options%console_disp == 'yes') then
        write(*,'(A)') '    {selecting auto coarsening base vertex selection as mesh is two dimensional}'
    end if  
elseif (options%crs_base_vlnc == 'lt4') then 
    options%crs_base_vlnc = 'auto'
    if (options%console_disp == 'yes') then
        write(*,'(A)') '    {selecting auto coarsening base vertex selection as mesh is two dimensional}'
    end if 
elseif (options%crs_base_vlnc == 'ne4') then 
    options%crs_base_vlnc = 'auto'
    if (options%console_disp == 'yes') then
        write(*,'(A)') '    {selecting auto coarsening base vertex selection as mesh is two dimensional}'
    end if 
end if 

!Select base vertices 
if (minvalence .LT. 2) then !Select any lower valence vertices 
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if (mesh%connectivity%valence(vv) .LT. 2) then 
            nexrtvtx = nexrtvtx + 1
        end if 
    end do 
    allocate(vbase(nexrtvtx))
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if (mesh%connectivity%valence(vv) .LT. 2) then 
            nexrtvtx = nexrtvtx + 1
            vbase(nexrtvtx) = vv 
        end if 
    end do 
elseif (maxvalence .GT. 2) then !Select any higher valence vertices 
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if (mesh%connectivity%valence(vv) .GT. 2) then 
            nexrtvtx = nexrtvtx + 1
        end if 
    end do 
    allocate(vbase(nexrtvtx))
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if (mesh%connectivity%valence(vv) .GT. 2) then 
            nexrtvtx = nexrtvtx + 1
            vbase(nexrtvtx) = vv 
        end if 
    end do 
else !Select first vertex
    allocate(vbase(1))
    vbase(1) = 1
end if 
return 
end subroutine identify_coarsening_base_vertices_2D




!Identify base vertex for coarsening 3D =========================
subroutine identify_coarsening_base_vertices_3D(vbase,mesh,options) 
implicit none 

!Variables - Import
integer(in64), dimension(:), allocatable :: vbase
type(mesh_data) :: mesh
type(options_data) :: options 

!Variables - Local 
integer(in64) :: vv,aa,ee,pp 
integer(in64) :: nexrtvtx,vftypeC,ftgt,vf,vbatt
integer(in64) :: maxvalence,minvalence,tgtvalence

!Initialise
if (allocated(vbase)) then 
    deallocate(vbase)
end if 

!Find valence bounds
maxvalence = maxval(mesh%connectivity%valence)
minvalence = minval(mesh%connectivity%valence)

!Return if maxvalence < 4 as mesh cannot be coarsened further
if (maxvalence .LT. 4) then 
    allocate(vbase(1))
    vbase(1) = -1
    return 
end if 

!Select base vertices (>4/<4/!=4/auto)
if (options%crs_base_vlnc == 'gt4') then 
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if (mesh%connectivity%valence(vv) .GT. 4) then 
            nexrtvtx = nexrtvtx + 1
        end if 
    end do 
    allocate(vbase(nexrtvtx))
    vbase(:) = 0 
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if (mesh%connectivity%valence(vv) .GT. 4) then 
            nexrtvtx = nexrtvtx + 1
            vbase(nexrtvtx) = vv 
        end if 
    end do 
elseif (options%crs_base_vlnc == 'lt4') then 
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if (mesh%connectivity%valence(vv) .LT. 4) then 
            nexrtvtx = nexrtvtx + 1
        end if 
    end do 
    allocate(vbase(nexrtvtx))
    vbase(:) = 0 
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if (mesh%connectivity%valence(vv) .LT. 4) then 
            nexrtvtx = nexrtvtx + 1
            vbase(nexrtvtx) = vv 
        end if 
    end do 
elseif (options%crs_base_vlnc == 'ne4') then 
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if ((mesh%connectivity%valence(vv) .GT. 4) .OR. (mesh%connectivity%valence(vv) .LT. 4)) then 
            nexrtvtx = nexrtvtx + 1
        end if 
    end do 
    allocate(vbase(nexrtvtx))
    vbase(:) = 0 
    nexrtvtx = 0 
    do vv=1,mesh%nvertex
        if ((mesh%connectivity%valence(vv) .GT. 4) .OR. (mesh%connectivity%valence(vv) .LT. 4)) then 
            nexrtvtx = nexrtvtx + 1
            vbase(nexrtvtx) = vv 
        end if 
    end do 
elseif (options%crs_base_vlnc == 'auto') then 
    allocate(vbase(1))
    vbase(1) = 0 
    tgtvalence = minvalence
    do vv=1,mesh%nvertex 
        if (mesh%connectivity%valence(vv) == tgtvalence) then 
            vftypeC = 0
            do aa=1,mesh%connectivity%valence(vv)
                ftgt = mesh%connectivity%v2f(vv,aa)
                if (ftgt .GT. 0) then  
                    do ee=1,mesh%faces(ftgt)%nvertex 
        
                        !Vertex in face
                        vf = mesh%faces(ftgt)%vertices(ee)
        
                        !Check this is not attached to vv
                        vbatt = 0
                        do pp=1,mesh%connectivity%valence(vf)
                            if (mesh%connectivity%v2v(vf,pp) == vv) then 
                                vbatt = 1
                                exit
                            end if 
                        end do 
                        
                        !Check if valence 4
                        if (vbatt == 0) then 
                            ! if (mesh%connectivity%valence(vf) == TGTValence) then 
                            if (mesh%connectivity%valence(vf) == 4) then 
                                vftypeC = 1
                            end if 
                        end if 
                    end do 
                end if 
            end do 
            if (vftypeC == 1) then 
                vbase(1) = vv
                exit
            end if 
        end if 
    end do 
else
    write(*,'(A)') '** incorrect base vertex type for coarsening specified: '//options%crs_base_vlnc
    stop
end if 
return 
end subroutine identify_coarsening_base_vertices_3D


end module mrsys_oqc_mod