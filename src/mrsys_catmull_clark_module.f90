!MRsys catmull_clark module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.5
!Updated 15-12-2023

!Module
module mrsys_catmull_clark_mod
use mrsys_io_mod
use sparse_matrix_mod
use mrsys_connectivity_mod
contains 


!Subroutine to evaluate the Catmull-Clark coarsening matricies for a multiresolution system =========================
subroutine evaluate_mrsys_catmull_clark_coarsening_matrices(mres_system,options)
implicit none 

!Variables - Import
type(mressystem_data) :: mres_system
type(options_data) :: options 

!Variables - Local 
integer(in64) :: rr
real(dp) :: sparsity

!For each coarsening 
do rr=2,mres_system%nlevel

    !Construct 
    call construct_catmull_clark_crsmat(mres_system%subd_system(rr)%Ri,mres_system%subd_system(rr-1)%mesh,&
    mres_system%subd_system(rr)%mesh)

    !Display
    if (options%console_disp == 'yes') then
        sparsity = (real(mres_system%subd_system(rr)%Ri%nnz,dp)/&
        (real(mres_system%subd_system(rr)%Ri%nrow*mres_system%subd_system(rr)%Ri%ncol,dp)))*100.0d0
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,A,A)') '    < coarsening = ',rr,' -> ',rr-1,' - vertices = ',&
        mres_system%subd_system(rr)%mesh%nvertex,' -> ',mres_system%subd_system(rr-1)%mesh%nvertex,&
        ' - nnz fraction = ',real2F0_Xstring(sparsity,8_in64)//'%',' >'
    end if 
end do 
return 
end subroutine evaluate_mrsys_catmull_clark_coarsening_matrices




!Subroutine to evaluate the Catmull-Clark refinement matricies for a multiresolution system =========================
subroutine evaluate_mrsys_catmull_clark_refinement_matrices(mres_system,options)
implicit none 

!Variables - Import
type(mressystem_data) :: mres_system
type(options_data) :: options 

!Variables - Local 
integer(in64) :: rr
real(dp) :: sparsity

!For each refinement 
do rr=1,mres_system%nrefine

    !Construct 
    call construct_catmull_clark_refmat(mres_system%subd_system(rr)%Ru,mres_system%subd_system(rr)%mesh)

    !Display
    if (options%console_disp == 'yes') then
        sparsity = (real(mres_system%subd_system(rr)%Ru%nnz,dp)/&
        (real(mres_system%subd_system(rr)%Ru%nrow*mres_system%subd_system(rr)%Ru%ncol,dp)))*100.0d0
        write(*,'(A,I0,A,I0,A,I0,A,A,A)') '    < refinement = ',rr,' - vertices = ',mres_system%subd_system(rr)%mesh%nvertex,&
        ' -> ',mres_system%subd_system(rr+1)%mesh%nvertex,' - nnz fraction = ',real2F0_Xstring(sparsity,8_in64)//'%',' >'
    end if 
end do 
return 
end subroutine evaluate_mrsys_catmull_clark_refinement_matrices




!Subroutine to propagate Catmull-Clark subdivision =========================
subroutine propagate_catmull_clark(subd_system,mesh_base,options)
implicit none 

!Variables - Import
type(mesh_data) :: mesh_base
type(options_data) :: options 
type(subdsystem_data), dimension(:), allocatable :: subd_system

!Variables - Local 
integer(in64) :: rr,ff
real(dp) :: sparsity

!Display
if (options%console_disp == 'yes') then
    write(*,'(A)') '--> propagating catmull-clark refinement'
end if 

!Allocate the subdivision system 
allocate(subd_system(options%n_sd_levels+1))

!Copy the base mesh into the base level of the subdivision system
subd_system(1)%mesh%ndim = mesh_base%ndim
subd_system(1)%mesh%nedge = mesh_base%nedge
subd_system(1)%mesh%nface = mesh_base%nface
subd_system(1)%mesh%nvertex = mesh_base%nvertex
allocate(subd_system(1)%mesh%vertices(mesh_base%nvertex,mesh_base%ndim))
allocate(subd_system(1)%mesh%vertex_sharp(mesh_base%nvertex))
subd_system(1)%mesh%vertices(:,:) = mesh_base%vertices(:,:)
subd_system(1)%mesh%vertex_sharp(:) = mesh_base%vertex_sharp(:)
allocate(subd_system(1)%mesh%edges(mesh_base%nedge,2))
allocate(subd_system(1)%mesh%edge_sharp(mesh_base%nedge))
subd_system(1)%mesh%edges(:,:) = mesh_base%edges(:,:)
subd_system(1)%mesh%edge_sharp(:) = mesh_base%edge_sharp(:)*real(options%n_sd_levels+1,dp)
allocate(subd_system(1)%mesh%faces(mesh_base%nface))
do ff=1,mesh_base%nface
    subd_system(1)%mesh%faces(ff)%nvertex = mesh_base%faces(ff)%nvertex
    allocate(subd_system(1)%mesh%faces(ff)%vertices(mesh_base%faces(ff)%nvertex))
    subd_system(1)%mesh%faces(ff)%vertices(:) = mesh_base%faces(ff)%vertices(:)
end do 

!Subdivision refinement loop 
do rr=1,options%n_sd_levels

    !Evaluate connectivity of the current level 
    call get_valence(subd_system(rr)%mesh)
    if (rr == 1) then !set to tag any shell edges as fully sharp on the first iteration 
        call get_connectivity(subd_system(rr)%mesh%connectivity,subd_system(rr)%mesh,real(options%n_sd_levels+1,dp))
    else
        call get_connectivity(subd_system(rr)%mesh%connectivity,subd_system(rr)%mesh)
    end if 

    !Build the refinement matrix for level rr
    call construct_catmull_clark_refmat(subd_system(rr)%Ru,subd_system(rr)%mesh)

    !Construct refined level vertices 
    subd_system(rr+1)%mesh%nvertex = subd_system(rr)%Ru%nrow
    allocate(subd_system(rr+1)%mesh%vertices(subd_system(rr+1)%mesh%nvertex,subd_system(rr)%mesh%ndim))
    allocate(subd_system(rr+1)%mesh%vertex_sharp(subd_system(rr+1)%mesh%nvertex))
    subd_system(rr+1)%mesh%vertex_sharp(1:subd_system(rr)%mesh%nvertex) = subd_system(rr)%mesh%vertex_sharp(:)
    subd_system(rr+1)%mesh%vertex_sharp(subd_system(rr)%mesh%nvertex+1:subd_system(rr+1)%mesh%nvertex) = 0 
    call matmul_csr_dvec(subd_system(rr+1)%mesh%vertices(:,1),subd_system(rr)%Ru,subd_system(rr)%mesh%vertices(:,1))
    call matmul_csr_dvec(subd_system(rr+1)%mesh%vertices(:,2),subd_system(rr)%Ru,subd_system(rr)%mesh%vertices(:,2))
    if (subd_system(rr)%mesh%ndim == 3) then 
        call matmul_csr_dvec(subd_system(rr+1)%mesh%vertices(:,3),subd_system(rr)%Ru,subd_system(rr)%mesh%vertices(:,3))
    end if 

    !Construct refined level mesh 
    call construct_refined_mesh(subd_system(rr+1)%mesh,subd_system(rr)%mesh)

    !Display
    if (options%console_disp == 'yes') then
        sparsity = (real(subd_system(rr)%Ru%nnz,dp)/(real(subd_system(rr)%Ru%nrow*subd_system(rr)%Ru%ncol,dp)))*100.0d0
        write(*,'(A,I0,A,I0,A,I0,A,A,A)') '    < refinement = ',rr,' - vertices = ',subd_system(rr)%mesh%nvertex,&
        ' -> ',subd_system(rr+1)%mesh%nvertex,' - nnz fraction = ',real2F0_Xstring(sparsity,8_in64)//'%',' >'
    end if 
end do 
return 
end subroutine propagate_catmull_clark




!Construct catmull-clark refinement matrix subroutine =========================
subroutine construct_catmull_clark_refmat(Ru,mesh)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(csrmatrix) :: Ru

!Variables - Local 
integer(in64) :: ii,vv,ee,ff
integer(in64) :: Rs_entry,Ra_entry
integer(in64) :: nentry_rs,nentry_ra
integer(in64) :: n_sharp,nvertex_n,rins,esins,Eatt,Vatt
integer(in64) :: EsharpI(2)
integer(in64), dimension(:,:), allocatable :: Rs_ij,Ra_ij
real(dp) :: Npfinv,k,insSharp
real(dp), dimension(:), allocatable :: Rs_V,Ra_V
type(csrmatrix) :: Rs,Ra

!Find number of seperate entries in refinement matrix components ---
!Rs ------
!Original vertex contribution
nentry_rs = mesh%nvertex

!Edge vertex creation contribution
nentry_rs = nentry_rs + 2*mesh%nedge

!Face vertex creation contribution
nentry_rs = nentry_rs + sum(mesh%faces(:)%nvertex)

!Ra ------
!Original vertex contribution
nentry_ra = 0
do vv=1,mesh%nvertex

    !Sharp incident edges
    n_sharp = 0
    do ee=1,mesh%connectivity%valence(vv)
        if (mesh%edge_sharp(mesh%connectivity%v2e(vv,ee)) .GT. 0.0d0) then 
            n_sharp = n_sharp + 1
        end if 
    end do 

    !Update
    if ((n_sharp .LT. 2) .AND. (mesh%vertex_sharp(vv) == 0)) then !Normal
        if (mesh%connectivity%valence(vv) == 2) then
            nentry_ra = nentry_ra + 4
        else
            nentry_ra = nentry_ra + 2*mesh%connectivity%valence(vv) + 2
        end if 
    elseif ((n_sharp == 2) .AND. (mesh%vertex_sharp(vv) == 0)) then !Sharp edge
        nentry_ra = nentry_ra + 3
    elseif ((n_sharp .GT. 2) .OR. (mesh%vertex_sharp(vv) == 1)) then !Sharp vertex
        nentry_ra = nentry_ra + 1
    end if 
end do 

!Edge vertex averaging contribution
do ee=1,mesh%nedge
    if (mesh%edge_sharp(ee) .GT. 1.0d0) then 
        nentry_ra = nentry_ra + 1
    else
        nentry_ra = nentry_ra + 3
    end if 
end do 

!Face vertex fixing contribution
nentry_ra = nentry_ra + mesh%nface

!Populate Catmull-Clark refinement matricies Ra Rs ====================
!New vertex quantity 
nvertex_n = mesh%nvertex + mesh%nedge + mesh%nface

!Rs sparse structure
allocate(Rs_ij(nentry_rs,2))
allocate(Rs_V(nentry_rs))
Rs_ij(:,:) = 0
Rs_V(:) = 0.0d0 

!Ra sparse structure
allocate(Ra_ij(nentry_ra,2))
allocate(Ra_V(nentry_ra))
Ra_ij(:,:) = 0
Ra_V(:) = 0.0d0 

!Populate splitting matrix Rs -----------------------------------------
Rs_entry = 1
rins = mesh%nvertex + 1

!Current verticies ------
do vv=1,mesh%nvertex
    Rs_V(Rs_entry) = 1
    Rs_ij(Rs_entry,:) = vv
    Rs_entry = Rs_entry + 1
end do 

!New edge verticies ------
do ee=1,mesh%nedge

    !End 1
    Rs_V(Rs_entry) = 0.5d0
    Rs_ij(Rs_entry,1) = rins
    Rs_ij(Rs_entry,2) = mesh%edges(ee,1)
    Rs_entry = Rs_entry + 1

    !End 2
    Rs_V(Rs_entry) = 0.5d0
    Rs_ij(Rs_entry,1) = rins
    Rs_ij(Rs_entry,2) = mesh%edges(ee,2)
    Rs_entry = Rs_entry + 1

    !Increment row 
    rins = rins + 1
end do 

!New face verticies ------
do ff=1,mesh%nface

    !Face points
    Npfinv = 1.0d0/real(mesh%faces(ff)%nvertex,dp)
    do vv=1,mesh%faces(ff)%nvertex
        Rs_V(Rs_entry) = Npfinv
        Rs_ij(Rs_entry,1) = rins
        Rs_ij(Rs_entry,2) = mesh%faces(ff)%vertices(vv)
        Rs_entry = Rs_entry + 1
    end do 

    !Increment
    rins = rins + 1
end do 

!Populate averaging matrix Ra -----------------------------------------
Ra_entry = 1

!Original verticies
do vv=1,mesh%nvertex
    
    !Identify any incident sharp edges 
    n_sharp = 0
    do ee=1,mesh%connectivity%valence(vv)
        if (mesh%edge_sharp(mesh%connectivity%v2e(vv,ee)) .GT. 0.0d0) then  
            n_sharp = n_sharp + 1
        end if 
    end do 

    !If this vertex in on a crease -> identify attached sharp edges
    if (n_sharp == 2) then  
        esins = 0
        EsharpI(:) = 0
        do ee=1,mesh%connectivity%valence(vv)
            if (mesh%edge_sharp(mesh%connectivity%v2e(vv,ee)) .GT. 0.0d0) then  
                esins = esins + 1
                EsharpI(esins) = mesh%connectivity%v2e(vv,ee)   
            end if
            if (esins .GE. 2) then 
                exit 
            end if 
        end do 
    end if 

    !Valence of this vertex (real valued)
    k = real(mesh%connectivity%valence(vv),dp)

    !Cases
    if ((n_sharp .LT. 2) .AND. (mesh%vertex_sharp(vv) == 0)) then !Normal 
        if (mesh%connectivity%valence(vv) == 2) then !Valence 2 approximation

            !Midpoint verticies on attached edges
            Ra_V(Ra_entry) = 0.25d0
            Ra_ij(Ra_entry,1) = vv
            Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%connectivity%v2e(vv,1)
            Ra_entry = Ra_entry + 1
            Ra_V(Ra_entry) = 0.25d0
            Ra_ij(Ra_entry,1) = vv
            Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%connectivity%v2e(vv,2)
            Ra_entry = Ra_entry + 1

            !Midpoint verticies on attached faces
            Ra_V(Ra_entry) = 0.25d0
            Ra_ij(Ra_entry,1) = vv
            Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%v2f(vv,1)
            Ra_entry = Ra_entry + 1
            Ra_V(Ra_entry) = 0.25d0
            Ra_ij(Ra_entry,1) = vv
            Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%v2f(vv,2)
            Ra_entry = Ra_entry + 1

        else !Standard Catmull-Clark

            !Base vertex
            Ra_V(Ra_entry) = (k - 2.0d0)/k
            Ra_ij(Ra_entry,1) = vv
            Ra_ij(Ra_entry,2) = vv
            Ra_entry = Ra_entry + 1
    
            !Attached vertices
            do ii=1,mesh%connectivity%valence(vv)
                Ra_V(Ra_entry) = 1.0d0/(k*k)
                Ra_ij(Ra_entry,1) = vv
                Ra_ij(Ra_entry,2) = mesh%connectivity%v2v(vv,ii)
                Ra_entry = Ra_entry + 1
            end do
    
            !Attached face midpoints
            do ii=1,mesh%connectivity%valence(vv)
                Ra_V(Ra_entry) = 1.0d0/(k*k)
                Ra_ij(Ra_entry,1) = vv
                Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%v2f(vv,ii)
                Ra_entry = Ra_entry + 1
            end do 
        end if 
    elseif ((n_sharp == 2) .AND. (mesh%vertex_sharp(vv) == 0)) then !On a sharp edge/crease
        
        !Average incident sharpness
        insSharp = 0.5d0*(mesh%edge_sharp(EsharpI(1)) + mesh%edge_sharp(EsharpI(2)))

        !Cases
        if (insSharp .GE. 1.0d0) then !Fully sharp

            !Base vertex
            Ra_V(Ra_entry) = 0.5d0
            Ra_ij(Ra_entry,1) = vv
            Ra_ij(Ra_entry,2) = vv
            Ra_entry = Ra_entry + 1

            !Attached midpoint verticies on sharp edges
            Ra_V(Ra_entry) = 0.25d0
            Ra_ij(Ra_entry,1) = vv
            Ra_ij(Ra_entry,2) = mesh%nvertex + EsharpI(1)
            Ra_entry = Ra_entry + 1
            Ra_V(Ra_entry) = 0.25d0
            Ra_ij(Ra_entry,1) = vv
            Ra_ij(Ra_entry,2) = mesh%nvertex + EsharpI(2)
            Ra_entry = Ra_entry + 1

        elseif ((insSharp .LT. 1.0d0) .AND. (insSharp .GT. 0.0d0)) then !Semi sharp
            if (mesh%connectivity%valence(vv) == 2) then !Valence 2 approximation   
                
                !Base vertex (linearly interpolated)   
                Ra_V(Ra_entry) = 0.5d0*insSharp
                Ra_ij(Ra_entry,1) = vv
                Ra_ij(Ra_entry,2) = vv
                Ra_entry = Ra_entry + 1

                !Midpoint verticies on attached edges (constant between both cases)
                Ra_V(Ra_entry) = 0.25d0
                Ra_ij(Ra_entry,1) = vv
                Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%connectivity%v2e(vv,1)
                Ra_entry = Ra_entry + 1
                Ra_V(Ra_entry) = 0.25d0
                Ra_ij(Ra_entry,1) = vv
                Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%connectivity%v2e(vv,2)
                Ra_entry = Ra_entry + 1

                !Midpoint verticies on attached faces (linearly interpolated) 
                Ra_V(Ra_entry) = 0.25d0*(1.0d0 - insSharp)
                Ra_ij(Ra_entry,1) = vv
                Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%v2f(vv,1)
                Ra_entry = Ra_entry + 1
                Ra_V(Ra_entry) = 0.25d0*(1.0d0 - insSharp)
                Ra_ij(Ra_entry,1) = vv
                Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%v2f(vv,2)
                Ra_entry = Ra_entry + 1

            else !Standard Catmull-Clark

                !Base vertex (linearly interpolated)       
                Ra_V(Ra_entry) = ((k - 2.0d0)/k)*(1.0d0 - insSharp) + 0.5d0*insSharp
                Ra_ij(Ra_entry,1) = vv
                Ra_ij(Ra_entry,2) = vv
                Ra_entry = Ra_entry + 1

                !Attached vertices (linearly interpolated)
                do ii=1,mesh%connectivity%valence(vv)
    
                    !Attached edge and vertex
                    Eatt = mesh%connectivity%v2e(vv,ii)
                    if (mesh%edges(Eatt,1) == vv) then 
                        Vatt = mesh%edges(Eatt,2)
                    else
                        Vatt = mesh%edges(Eatt,1)
                    end if 
                    
                    !Attached vertex
                    Ra_V(Ra_entry) = (1/(k*k))*(1 - insSharp)
                    Ra_ij(Ra_entry,1) = vv
                    Ra_ij(Ra_entry,2) = Vatt
                    Ra_entry = Ra_entry + 1
    
                    !Attached edge midpoint vertex if edge is sharp
                    if (mesh%edge_sharp(Eatt) .GT. 0.0d0) then  
                        Ra_V(Ra_entry) = 0.25d0*insSharp
                        Ra_ij(Ra_entry,1) = vv
                        Ra_ij(Ra_entry,2) = mesh%nvertex + Eatt
                        Ra_entry = Ra_entry + 1
                    end if 
    
                    !Attached face midpoint
                    Ra_V(Ra_entry) = (1.0d0/(k*k))*(1.0d0 - insSharp)
                    Ra_ij(Ra_entry,1) = vv
                    Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%v2f(vv,ii)
                    Ra_entry = Ra_entry + 1
                end do 
            end if 
        end if 
    elseif ((n_sharp .GT. 2) .OR. (mesh%vertex_sharp(vv) == 1)) then !Sharp vertex (implied || specified)
        Ra_V(Ra_entry) = 1.0d0 
        Ra_ij(Ra_entry,1) = vv
        Ra_ij(Ra_entry,2) = vv
        Ra_entry = Ra_entry + 1
    end if 
end do 

!Edge midpoints
rins = mesh%nvertex + 1
do ee=1,mesh%nedge
    if (mesh%edge_sharp(ee) .LE. 0.0d0) then !Non sharp standard
        
        !Edge midpoint
        Ra_V(Ra_entry) = 0.5d0
        Ra_ij(Ra_entry,1) = rins
        Ra_ij(Ra_entry,2) = rins
        Ra_entry = Ra_entry + 1

        !Attached face points
        Ra_V(Ra_entry) = 0.25d0
        Ra_ij(Ra_entry,1) = rins
        Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%e2f(ee,1)
        Ra_entry = Ra_entry + 1
        Ra_V(Ra_entry) = 0.25d0
        Ra_ij(Ra_entry,1) = rins
        Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%e2f(ee,2)
        Ra_entry = Ra_entry + 1

    else !Sharp/Semi sharp
        if (mesh%edge_sharp(ee) .GE. 1.0d0) then !Fully sharp
            Ra_V(Ra_entry) = 1.0d0 
            Ra_ij(Ra_entry,1) = rins
            Ra_ij(Ra_entry,2) = rins
            Ra_entry = Ra_entry + 1
        elseif ((mesh%edge_sharp(ee) .LT. 1.0d0) .AND. (mesh%edge_sharp(ee) .GT. 0)) then !Semi sharp
            
            !Edge midpoint (linearly interpolated)
            Ra_V(Ra_entry) = mesh%edge_sharp(ee) + (1.0d0 - mesh%edge_sharp(ee))*0.5d0
            Ra_ij(Ra_entry,1) = rins
            Ra_ij(Ra_entry,2) = rins
            Ra_entry = Ra_entry + 1

            !Attached face points (linearly interpolated)
            Ra_V(Ra_entry) = 0.25d0*(1.0d0 - mesh%edge_sharp(ee))
            Ra_ij(Ra_entry,1) = rins
            Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%e2f(ee,1)
            Ra_entry = Ra_entry + 1
            Ra_V(Ra_entry) = 0.25d0*(1.0d0 - mesh%edge_sharp(ee))
            Ra_ij(Ra_entry,1) = rins
            Ra_ij(Ra_entry,2) = mesh%nvertex + mesh%nedge + mesh%connectivity%e2f(ee,2)
            Ra_entry = Ra_entry + 1
        end if 
    end if 
    rins = rins + 1
end do 

!Face midpoints
do ff=1,mesh%nface

    !Entry
    Ra_V(Ra_entry) = 1.0d0 
    Ra_ij(Ra_entry,1) = rins
    Ra_ij(Ra_entry,2) = rins
    Ra_entry = Ra_entry + 1

    !Increment
    rins = rins + 1
end do 

!Set final entry counts 
Rs_entry = Rs_entry - 1
Ra_entry = Ra_entry - 1

!Construct sparse Rs 
call build_csr_matrix(Rs,Rs_V(1:Rs_entry),Rs_ij(1:Rs_entry,1),Rs_ij(1:Rs_entry,2),nvertex_n,mesh%nvertex)

!Construct sparse Ra
call build_csr_matrix(Ra,Ra_V(1:Ra_entry),Ra_ij(1:Ra_entry,1),Ra_ij(1:Ra_entry,2),nvertex_n,nvertex_n)

!Construct refinement matrix 
call matmul_csr(Ru,Ra,Rs)
return 
end subroutine construct_catmull_clark_refmat




!Construct catmull-clark coarsening matrix subroutine =========================
subroutine construct_catmull_clark_crsmat(Ri,mesh_coarse,mesh_fine)
implicit none 

!Variables - Import
type(mesh_data) :: mesh_fine,mesh_coarse
type(csrmatrix) :: Ri

!Variables - Local 
integer(in64) :: vv,ee
integer(in64) :: Ri_entry,nentry_ri,fvidx,etgt,vtgt,ftgt,fa,fb
integer(in64) :: n_sharp,esins,vertex_state,edge_v3x,vtx_v3x
integer(in64) :: EsharpI(2),VsharpI(2)
integer(in64), dimension(:,:), allocatable :: Ri_ij
real(dp) :: k,ka
real(dp), dimension(:), allocatable :: Ri_V

!Find upper bound of non zero elements in Ri ------
!Original vertex contribution
nentry_ri = mesh_coarse%nvertex

!Vertex contribution 
nentry_ri = nentry_ri + 2*sum(mesh_fine%connectivity%valence(:))

!Ri sparse structure
allocate(Ri_ij(nentry_ri,2))
allocate(Ri_V(nentry_ri))
Ri_ij(:,:) = 0
Ri_V(:) = 0.0d0 

!Construct direct reverse Catmull-Clark subdivision matrix
Ri_entry = 1
do vv=1,mesh_coarse%nvertex

    !Identify any incident sharp edges 
    n_sharp = 0
    do ee=1,mesh_fine%connectivity%valence(vv)
        if (mesh_fine%edge_sharp(mesh_fine%connectivity%v2e(vv,ee)) .GT. 0.0d0) then  
            n_sharp = n_sharp + 1
        end if 
    end do 

    !If this vertex in on a crease -> identify attached sharp edges
    if (n_sharp == 2) then  
        esins = 0
        EsharpI(:) = 0
        do ee=1,mesh_fine%connectivity%valence(vv)
            if (mesh_fine%edge_sharp(mesh_fine%connectivity%v2e(vv,ee)) .GT. 0.0d0) then  
                esins = esins + 1
                EsharpI(esins) = mesh_fine%connectivity%v2e(vv,ee)   
            end if
            if (esins .GE. 2) then 
                exit 
            end if 
        end do 
        if (mesh_fine%edges(EsharpI(1),1) == vv) then 
            VsharpI(1) = mesh_fine%edges(EsharpI(1),2)
        else
            VsharpI(1) = mesh_fine%edges(EsharpI(1),1)
        end if 
        if (mesh_fine%edges(EsharpI(2),1) == vv) then 
            VsharpI(2) = mesh_fine%edges(EsharpI(2),2)
        else
            VsharpI(2) = mesh_fine%edges(EsharpI(2),1)
        end if 
    end if 

    !Valence of this vertex (real valued)
    k = real(mesh_fine%connectivity%valence(vv),dp)

    !Identify vertex state (0 = normal | 1 = on sharp edge | 2 = sharp vertex) 
    if ((n_sharp .LT. 2) .AND. (mesh_fine%vertex_sharp(vv) == 0)) then !Normal
        vertex_state = 0 
    elseif ((n_sharp == 2) .AND. (mesh_fine%vertex_sharp(vv) == 0)) then !On a sharp edge/crease
        vertex_state = 1 
    elseif ((n_sharp .GT. 2) .OR. (mesh_fine%vertex_sharp(vv) == 1)) then !Sharp vertex (implied || specified)   
        vertex_state = 2 
    end if 

    !Populate reverse subdivision matrix
    if (vertex_state == 0) then !Normal 
        if (mesh_fine%connectivity%valence(vv) == 3) then !Valence 3 -> find approximation 

            !Check if there is a non valence 3 vetex attached on the coarse level and if so select this vertex and the edge joining the two 
            vtx_v3x = 0
            edge_v3x = 0 
            fa = 0 
            fb = 0 
            do ee=1,mesh_coarse%connectivity%valence(vv)
                etgt = mesh_coarse%connectivity%v2e(vv,ee)
                if (mesh_coarse%edges(etgt,1) == vv) then 
                    vtgt = mesh_coarse%edges(etgt,2)
                else
                    vtgt = mesh_coarse%edges(etgt,1)
                end if 
                if (mesh_coarse%connectivity%valence(vtgt) .NE. 3) then 
                    vtx_v3x = vtgt
                    edge_v3x = etgt
                    fa = mesh_coarse%connectivity%e2f(edge_v3x,1)
                    fb = mesh_coarse%connectivity%e2f(edge_v3x,2)
                    exit 
                end if 
            end do 

            !Approximation 
            if (vtx_v3x == 0) then !No attached non valence 3 -> keep fixed
                Ri_ij(Ri_entry,1) = vv
                Ri_ij(Ri_entry,2) = vv
                Ri_V(Ri_entry) = 1.0d0 
                Ri_entry = Ri_entry + 1
            else !Attached non valence 3 -> approximate 

                !Extract adjacent vertex valence 
                ka = real(mesh_fine%connectivity%valence(vtx_v3x),dp)

                !Attached vertex vtx_v3x
                Ri_ij(Ri_entry,1) = vv
                Ri_ij(Ri_entry,2) = vtx_v3x
                Ri_V(Ri_entry) = -ka/(ka - 3.0d0)
                Ri_entry = Ri_entry + 1

                !Edge verticies on vtx_v3x
                do ee=1,mesh_coarse%connectivity%valence(vtx_v3x)
                    etgt = mesh_coarse%connectivity%v2e(vtx_v3x,ee)
                    vtgt = mesh_coarse%nvertex + etgt
                    if (etgt == edge_v3x) then !perturb with extra entry
                        Ri_ij(Ri_entry,1) = vv
                        Ri_ij(Ri_entry,2) = vtgt
                        Ri_V(Ri_entry) = (4.0d0/(ka*(ka - 3.0d0))) + 4.0d0
                        Ri_entry = Ri_entry + 1
                    else !use normal entry
                        Ri_ij(Ri_entry,1) = vv
                        Ri_ij(Ri_entry,2) = vtgt
                        Ri_V(Ri_entry) = 4.0d0/(ka*(ka - 3.0d0))
                        Ri_entry = Ri_entry + 1
                    end if 
                end do 

                !Face verticies on vtx_v3x
                do ee=1,mesh_coarse%connectivity%valence(vtx_v3x)
                    ftgt = mesh_coarse%connectivity%v2f(vtx_v3x,ee) 
                    if (ftgt .GT. 0) then 
                        fvidx = mesh_coarse%nvertex + mesh_coarse%nedge + ftgt
                        if ((ftgt == fa) .OR. (ftgt == fb)) then !perturb with extra entry
                            Ri_ij(Ri_entry,1) = vv
                            Ri_ij(Ri_entry,2) = fvidx
                            Ri_V(Ri_entry) = (-1.0d0/(ka*(ka - 3.0d0))) - 1.0d0 
                            Ri_entry = Ri_entry + 1
                        else !use normal entry
                            Ri_ij(Ri_entry,1) = vv
                            Ri_ij(Ri_entry,2) = fvidx
                            Ri_V(Ri_entry) = -1.0d0/(ka*(ka - 3.0d0))
                            Ri_entry = Ri_entry + 1
                        end if 
                    end if 
                end do
            end if 
        else !Normal

            !Base vertex
            Ri_ij(Ri_entry,1) = vv
            Ri_ij(Ri_entry,2) = vv
            Ri_V(Ri_entry) = k/(k - 3.0d0)
            Ri_entry = Ri_entry + 1

            !Edge verticies
            do ee=1,mesh_fine%connectivity%valence(vv)
                Ri_ij(Ri_entry,1) = vv
                Ri_ij(Ri_entry,2) = mesh_fine%connectivity%v2v(vv,ee)
                Ri_V(Ri_entry) = -4.0d0/(k*(k - 3.0d0))
                Ri_entry = Ri_entry + 1
            end do 

            !Face verticies
            do ee=1,mesh_fine%connectivity%valence(vv)
                if (mesh_coarse%connectivity%v2f(vv,ee) .GT. 0) then 
                    fvidx = mesh_coarse%nvertex + mesh_coarse%nedge + mesh_coarse%connectivity%v2f(vv,ee) 
                    Ri_ij(Ri_entry,1) = vv
                    Ri_ij(Ri_entry,2) = fvidx
                    Ri_V(Ri_entry) = 1.0d0/(k*(k - 3.0d0))
                    Ri_entry = Ri_entry + 1
                end if 
            end do 
        end if 
    elseif (vertex_state == 1) then !On a sharp edge/crease
        
        !Base vertex
        Ri_ij(Ri_entry,1) = vv
        Ri_ij(Ri_entry,2) = vv
        Ri_V(Ri_entry) = 2.0d0 
        Ri_entry = Ri_entry + 1

        !Attached edge verticies
        Ri_ij(Ri_entry,1) = vv
        Ri_ij(Ri_entry,2) = VsharpI(1)
        Ri_V(Ri_entry) = -0.5d0
        Ri_entry = Ri_entry + 1
        Ri_ij(Ri_entry,1) = vv
        Ri_ij(Ri_entry,2) = VsharpI(2)
        Ri_V(Ri_entry) = -0.5d0
        Ri_entry = Ri_entry + 1
    elseif (vertex_state == 2) then !Sharp vertex (implied || specified)     
        Ri_ij(Ri_entry,1) = vv
        Ri_ij(Ri_entry,2) = vv
        Ri_V(Ri_entry) = 1.0d0 
        Ri_entry = Ri_entry + 1
    end if 
end do 

!Set final entry count
Ri_entry = Ri_entry - 1

!Construct sparse Ri 
call build_csr_matrix(Ri,Ri_V(1:Ri_entry),Ri_ij(1:Ri_entry,1),Ri_ij(1:Ri_entry,2),mesh_coarse%nvertex,mesh_fine%nvertex)
return 
end subroutine construct_catmull_clark_crsmat


end module mrsys_catmull_clark_mod