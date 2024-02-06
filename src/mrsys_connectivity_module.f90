!MRsys mesh connectivity module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.4
!Updated 12-12-2023

!Module
module mrsys_connectivity_mod
use mrsys_data_mod
contains 


!Get valence subroutine =========================
subroutine get_valence(mesh)
implicit none 

!Variables - Import
type(mesh_data) :: mesh

!Switch on dimensions 
if (mesh%ndim == 2) then 
    call get_valence_2D(mesh%connectivity%valence,mesh%connectivity%max_valence,mesh%edges,mesh%nedge,mesh%nvertex)
elseif (mesh%ndim == 3) then 
    call get_valence_3D(mesh%connectivity%valence,mesh%connectivity%max_valence,mesh%faces,mesh%nface,mesh%nvertex)
end if 
return 
end subroutine get_valence




!Get valence 2D subroutine =========================
subroutine get_valence_2D(valence,maxValence,edges,nedge,nvtx)
implicit none 

!Variables - Import
integer(in64) :: maxValence,nedge,nvtx
integer(in64), dimension(:), allocatable :: valence
integer(in64), dimension(:,:) :: edges

!Variables - Local
integer(in64) :: ee

!Initialse valence array
if (allocated(valence)) then 
    deallocate(valence)
end if 
allocate(valence(nvtx))
valence(:) = 0 

!Evaluate valence 
do ee=1,nedge
    valence(edges(ee,1)) = valence(edges(ee,1)) + 1
    valence(edges(ee,2)) = valence(edges(ee,2)) + 1
end do 

!Set maximum valence
maxValence = maxval(valence(:))
return 
end subroutine get_valence_2D




!Get valence 3D subroutine =========================
subroutine get_valence_3D(valence,maxValence,faces,nface,nvtx)
implicit none 

!Variables - Import
integer(in64) :: maxValence,nface,nvtx
integer(in64), dimension(:), allocatable :: valence
type(face_data), dimension(:) :: faces

!Variables - Local
integer(in64) :: ff,ee,vv
integer(in64) :: ev1,ev2,ubValence,evalid
integer(in64), dimension(:,:), allocatable :: vconnect

!Initialse valence array
if (allocated(valence)) then 
    deallocate(valence)
end if 
allocate(valence(nvtx))
valence(:) = 0 

!Upper bound of maximum valence 
valence(:) = 0 
do ff=1,nface
    do ee=1,faces(ff)%nvertex

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,faces(ff)%nvertex) + 1
        ev1 = faces(ff)%vertices(ev1) 
        ev2 = faces(ff)%vertices(ev2) 

        !Accumulate valence 
        valence(ev1) = valence(ev1) + 1
        valence(ev2) = valence(ev2) + 1
    end do 
end do 
ubValence = 2*maxval(valence)

!Construct actual valence of each vertex
allocate(vconnect(nvtx,ubValence))
vconnect(:,:) = 0 
valence(:) = 0
do ff=1,nface
    do ee=1,faces(ff)%nvertex

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,faces(ff)%nvertex) + 1
        ev1 = faces(ff)%vertices(ev1) 
        ev2 = faces(ff)%vertices(ev2)  

        !Check against vconnect 
        evalid = 1
        do vv=1,ubValence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !Add valence if new edge
        if (evalid == 1) then 
            
            !Increment valence on each vertex
            valence(ev1) = valence(ev1) + 1
            valence(ev2) = valence(ev2) + 1

            !Update vconnect
            do vv=1,ubValence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    exit 
                end if 
            end do 
            do vv=1,ubValence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Set maximum valence
maxValence = maxval(valence(:))
return 
end subroutine get_valence_3D




!Get connectivity subroutine =========================
subroutine get_connectivity(connectivity,mesh,set_sharp_shell_edges)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(connectivity_data) :: connectivity
real(dp), optional :: set_sharp_shell_edges

!Build v2v and v2e
call build_v2v_v2e(mesh,connectivity)

!Build v2f
call build_v2f(mesh,connectivity)

!Build f2e
call build_f2e(mesh,connectivity)

!Build e2f 
call build_e2f(mesh,connectivity)

!Set shell edges as sharp with the value set_sharp_shell_edges if present and requested 
if (present(set_sharp_shell_edges)) then 
    call set_shell_edges_as_sharp(mesh,connectivity,set_sharp_shell_edges)
end if 
return 
end subroutine get_connectivity




!Build v2v and v2e =========================
subroutine build_v2v_v2e(mesh,connectivity)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(connectivity_data) :: connectivity

!Variables - Local 
integer(in64) :: ee
integer(in64) :: v1,v2
integer(in64) :: vins(mesh%nvertex)

!Build v2v and v2e
vins(:) = 0 
allocate(connectivity%v2v(mesh%nvertex,connectivity%max_valence))
allocate(connectivity%v2e(mesh%nvertex,connectivity%max_valence))
connectivity%v2v(:,:) = 0 
connectivity%v2e(:,:) = 0 
do ee=1,mesh%nedge

    !Edge vertices 
    v1 = mesh%edges(ee,1)
    v2 = mesh%edges(ee,2)

    !Increment counts
    vins(v1) = vins(v1) + 1
    vins(v2) = vins(v2) + 1

    !Add to v2v
    connectivity%v2v(v1,vins(v1)) = v2
    connectivity%v2v(v2,vins(v2)) = v1

    !Add to v2e
    connectivity%v2e(v1,vins(v1)) = ee
    connectivity%v2e(v2,vins(v2)) = ee
end do
return 
end subroutine build_v2v_v2e




!Build v2f subroutine =========================
subroutine build_v2f(mesh,connectivity)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(connectivity_data) :: connectivity

!Variables - Local 
integer(in64) :: ff,vv
integer(in64) :: v1
integer(in64) :: vins(mesh%nvertex)

!Build v2f
vins(:) = 0 
allocate(connectivity%v2f(mesh%nvertex,connectivity%max_valence))
connectivity%v2f(:,:) = 0 
do ff=1,mesh%nface
    do vv=1,mesh%faces(ff)%nvertex
        v1 = mesh%faces(ff)%vertices(vv)
        vins(v1) = vins(v1) + 1
        connectivity%v2f(v1,vins(v1)) = ff
    end do 
end do 
return
end subroutine build_v2f




!Build f2e subroutine =========================
subroutine build_f2e(mesh,connectivity)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(connectivity_data) :: connectivity

!Variables - Local 
integer(in64) :: ff,vv,aa
integer(in64) :: v1,v2,etgt

!Build f2e
do ff=1,mesh%nface
    allocate(mesh%faces(ff)%edges(mesh%faces(ff)%nvertex))
    mesh%faces(ff)%edges(:) = 0 
    do vv=1,mesh%faces(ff)%nvertex

        !End vertices of this edge 
        v1 = vv 
        v2 = mod(vv,mesh%faces(ff)%nvertex) + 1
        v1 = mesh%faces(ff)%vertices(v1)
        v2 = mesh%faces(ff)%vertices(v2)

        !Search v2e for connection 
        do aa=1,connectivity%valence(v1)
            etgt = connectivity%v2e(v1,aa)
            if ((mesh%edges(etgt,1) == v2) .OR. (mesh%edges(etgt,2) == v2)) then
                mesh%faces(ff)%edges(vv) = etgt 
                exit 
            end if 
        end do 
    end do 
end do 
return 
end subroutine build_f2e




!Build e2f subroutine =========================
subroutine build_e2f(mesh,connectivity)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(connectivity_data) :: connectivity

!Variables - Local 
integer(in64) :: ff,vv
integer(in64) :: etgt
integer(in64) :: eins(mesh%nedge)

!Build e2f 
eins(:) = 0 
allocate(connectivity%e2f(mesh%nedge,2))
connectivity%e2f(:,:) = 0 
do ff=1,mesh%nface
    do vv=1,mesh%faces(ff)%nvertex
        etgt = mesh%faces(ff)%edges(vv)
        eins(etgt) = eins(etgt) + 1
        connectivity%e2f(etgt,eins(etgt)) = ff
    end do 
end do 
return 
end subroutine build_e2f




!Tag shell edges as sharp subroutine =========================
subroutine set_shell_edges_as_sharp(mesh,connectivity,sharpval)
implicit none 

!Variables - Import
real(dp) :: sharpval
type(mesh_data) :: mesh
type(connectivity_data) :: connectivity

!Variables - Local 
integer(in64) :: ee 

!Set shell edges as sharp 
do ee=1,mesh%nedge
    if ((connectivity%e2f(ee,1) == 0) .OR. (connectivity%e2f(ee,2) == 0)) then 
        mesh%edge_sharp(ee) = sharpval
    end if 
end do 
return 
end subroutine set_shell_edges_as_sharp




!Build edges from faces subroutine =========================
subroutine build_edges_from_faces(mesh,maxvalence)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
integer(in64) :: maxvalence

!Variables - Local 
integer(in64) :: ff,ee,vv 
integer(in64) :: ev1,ev2,evalid,edge_idx
integer(in64), dimension(:,:), allocatable :: vconnect,edgeidx 

!Initialise
allocate(vconnect(mesh%nvertex,maxvalence))
allocate(edgeidx(mesh%nvertex,maxvalence))
vconnect(:,:) = 0
edgeidx(:,:) = 0

!Index edges 
mesh%nedge = 0 
do ff=1,mesh%nface
    do ee=1,mesh%faces(ff)%nvertex

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,mesh%faces(ff)%nvertex) + 1
        ev1 = mesh%faces(ff)%vertices(ev1)
        ev2 = mesh%faces(ff)%vertices(ev2)

        !Check against vconnect 
        evalid = 1 
        do vv=1,maxvalence
            if (vconnect(ev1,vv) == ev2) then 
                evalid = 0
                exit 
            end if 
            if (vconnect(ev2,vv) == ev1) then 
                evalid = 0
                exit 
            end if 
        end do 

        !Add if valid 
        if (evalid == 1) then 

            !Increment edge count 
            mesh%nedge = mesh%nedge + 1

            !Add edge to connection structure     
            do vv=1,maxvalence
                if (vconnect(ev1,vv) == 0) then 
                    vconnect(ev1,vv) = ev2 
                    edgeidx(ev1,vv) = mesh%nedge
                    exit 
                end if 
            end do 
            do vv=1,maxvalence
                if (vconnect(ev2,vv) == 0) then 
                    vconnect(ev2,vv) = ev1 
                    edgeidx(ev2,vv) = mesh%nedge
                    exit 
                end if 
            end do 
        end if 
    end do 
end do 

!Build edges 
if (allocated(mesh%edges)) then 
    deallocate(mesh%edges)
end if 
if (allocated(mesh%edge_sharp)) then 
    deallocate(mesh%edge_sharp)
end if 
allocate(mesh%edges(mesh%nedge,2))
allocate(mesh%edge_sharp(mesh%nedge))
mesh%edges(:,:) = 0 
mesh%edge_sharp(:) = 0.0d0 
do ff=1,mesh%nface
    do ee=1,mesh%faces(ff)%nvertex

        !Edge end vertices
        ev1 = ee
        ev2 = mod(ee,mesh%faces(ff)%nvertex) + 1
        ev1 = mesh%faces(ff)%vertices(ev1)
        ev2 = mesh%faces(ff)%vertices(ev2)

        !Check index of this edge 
        edge_idx = 0 
        do vv=1,maxvalence
            if (vconnect(ev1,vv) == ev2) then 
                edge_idx = edgeidx(ev1,vv)
                exit 
            end if 
        end do 

        !Construct
        if (edge_idx .GT. 0) then !If valid index then remove from edgeidx and construct 

            !Remove edge index 
            do vv=1,maxvalence
                if (edgeidx(ev1,vv) == edge_idx) then 
                    edgeidx(ev1,vv) = -1*edgeidx(ev1,vv)
                    exit 
                end if 
            end do 
            do vv=1,maxvalence
                if (edgeidx(ev2,vv) == edge_idx) then 
                    edgeidx(ev2,vv) = -1*edgeidx(ev2,vv)
                    exit 
                end if 
            end do 

            !Build edge 
            mesh%edges(edge_idx,1) = ev1
            mesh%edges(edge_idx,2) = ev2
        end if 
    end do 
end do 
return 
end subroutine build_edges_from_faces




!Subroutine to construct refined level mesh =========================
subroutine construct_refined_mesh(mesh_refined,mesh_base)
implicit none 

!Variables - Import 
type(mesh_data) :: mesh_base
type(mesh_data) :: mesh_refined 

!Switch on dimensions 
if (mesh_base%ndim == 2) then 
    call construct_refined_mesh_2D(mesh_refined,mesh_base)
elseif (mesh_base%ndim == 3) then 
    call construct_refined_mesh_3D(mesh_refined,mesh_base)
end if 
return 
end subroutine construct_refined_mesh




!Subroutine to construct refined level mesh 2D =========================
subroutine construct_refined_mesh_2D(mesh_refined,mesh_base)
implicit none 

!Variables - Import 
type(mesh_data) :: mesh_base
type(mesh_data) :: mesh_refined 

!Variables - Local 
integer(in64) :: ee
integer(in64) :: eins
integer(in64) :: nedgeN

!Number of new edges 
nedgeN = 2*mesh_base%nedge

!Allocate required new variables 
allocate(mesh_refined%edges(nedgeN,2))
allocate(mesh_refined%edge_sharp(nedgeN))
allocate(mesh_refined%faces(0))

!Mesh new surface 
eins = 0
do ee=1,mesh_base%nedge

    !Build first edge (v1->vmid)
    eins = eins + 1
    mesh_refined%edges(eins,1) = mesh_base%edges(ee,1)
    mesh_refined%edges(eins,2) = mesh_base%nvertex + ee
    mesh_refined%edge_sharp(eins) = mesh_base%edge_sharp(ee)

    !Build second edge (vmid->v2)
    eins = eins + 1
    mesh_refined%edges(eins,1) = mesh_base%nvertex + ee
    mesh_refined%edges(eins,2) = mesh_base%edges(ee,2) 
    mesh_refined%edge_sharp(eins) = mesh_base%edge_sharp(ee)
end do 

!Set properties in the new mesh 
mesh_refined%ndim = mesh_base%ndim
mesh_refined%nedge = nedgeN
mesh_refined%nface = 0
return 
end subroutine construct_refined_mesh_2D




!Subroutine to construct refined level mesh 3D =========================
subroutine construct_refined_mesh_3D(mesh_refined,mesh_base)
implicit none 

!Variables - Import 
type(mesh_data) :: mesh_base
type(mesh_data) :: mesh_refined 

!Variables - Local  
integer(in64) :: ff,vv,ee,aa
integer(in64) :: eins,fins,e1,e2,v1,v2
integer(in64) :: nedgeN,nfaceN,maxValenceN,evalid
integer(in64) :: eNcand(4,2)
integer(in64), dimension(:,:), allocatable :: vconnect
real(dp) :: ecShrp
real(dp) :: eNcand_sharp(4)

!Number of new edges and faces 
nfaceN = sum(mesh_base%faces(:)%nvertex)
nedgeN = 2*mesh_base%nedge + sum(mesh_base%faces(:)%nvertex)

!Estimate maximum valence of the new mesh 
maxValenceN = max(4,mesh_base%connectivity%max_valence,maxval(mesh_base%faces(:)%nvertex))

!Allocate required new variables 
allocate(mesh_refined%edges(nedgeN,2))
allocate(mesh_refined%edge_sharp(nedgeN))
allocate(mesh_refined%faces(nfaceN))
allocate(vconnect(mesh_refined%nvertex,maxValenceN))
vconnect(:,:) = 0 

!Mesh new surface 
eins = 0
fins = 0 
do ff=1,mesh_base%nface
    do vv=1,mesh_base%faces(ff)%nvertex

        !Increment face count 
        fins = fins + 1

        !Allocate new face
        mesh_refined%faces(fins)%nvertex = 4
        allocate(mesh_refined%faces(fins)%vertices(4))

        !Edges in face adjacent to vertex vv
        e1 = modulo(vv-2,mesh_base%faces(ff)%nvertex) + 1
        e2 = vv 

        !Populate new face
        mesh_refined%faces(fins)%vertices(1) = mesh_base%faces(ff)%vertices(vv)
        mesh_refined%faces(fins)%vertices(2) = mesh_base%nvertex + mesh_base%faces(ff)%edges(e2)
        mesh_refined%faces(fins)%vertices(3) = mesh_base%nvertex + mesh_base%nedge + ff
        mesh_refined%faces(fins)%vertices(4) = mesh_base%nvertex + mesh_base%faces(ff)%edges(e1)

        !Candidate new edges and their sharpness
        eNcand(1,1) = mesh_refined%faces(fins)%vertices(1)
        eNcand(1,2) = mesh_refined%faces(fins)%vertices(2)
        eNcand_sharp(1) = mesh_base%edge_sharp(mesh_base%faces(ff)%edges(e2)) - 1.0d0 !Decriment

        eNcand(2,1) = mesh_refined%faces(fins)%vertices(2)
        eNcand(2,2) = mesh_refined%faces(fins)%vertices(3)
        eNcand_sharp(2) = 0.0d0 

        eNcand(3,1) = mesh_refined%faces(fins)%vertices(3)
        eNcand(3,2) = mesh_refined%faces(fins)%vertices(4)
        eNcand_sharp(3) = 0.0d0 

        eNcand(4,1) = mesh_refined%faces(fins)%vertices(4)
        eNcand(4,2) = mesh_refined%faces(fins)%vertices(1)
        eNcand_sharp(4) = mesh_base%edge_sharp(mesh_base%faces(ff)%edges(e1)) - 1.0d0 !Decriment

        !Add candidate edges (if new) to Vconnect and eN  
        do ee=1,4
            
            !Vertices on edge ends of the new candidate edge 
            v1 = eNcand(ee,1)
            v2 = eNcand(ee,2)
            ecShrp = eNcand_sharp(ee)

            !Check against vconnect 
            evalid = 1
            do aa=1,MaxValenceN
                if (vconnect(v1,aa) == v2) then 
                    evalid = 0
                    exit 
                end if 
                if (vconnect(v2,aa) == v1) then 
                    evalid = 0
                    exit 
                end if 
            end do 
            
            !Add edge if valid
            if (evalid == 1) then 

                !Create edge
                mesh_refined%edges(eins+1,1) = v1
                mesh_refined%edges(eins+1,2) = v2
                mesh_refined%edge_sharp(eins+1) = ecShrp

                !Update vconnect
                do aa=1,MaxValenceN
                    if (vconnect(v1,aa) == 0) then 
                        vconnect(v1,aa) = v2
                        exit
                    end if 
                end do 
                do aa=1,MaxValenceN
                    if (vconnect(v2,aa) == 0) then 
                        vconnect(v2,aa) = v1
                        exit 
                    end if 
                end do 

                !Increment edge count 
                eins = eins + 1
            end if 
        end do 
    end do 
end do 

!Set properties in the new mesh 
mesh_refined%ndim = mesh_base%ndim
mesh_refined%nedge = nedgeN
mesh_refined%nface = nfaceN
return 
end subroutine construct_refined_mesh_3D


end module mrsys_connectivity_mod