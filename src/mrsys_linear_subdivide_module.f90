!MRsys linear_subdivide module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.1
!Updated 13-02-2024

!Module
module mrsys_linear_subdivide_mod
use mrsys_io_mod
use sparse_matrix_mod
use mrsys_connectivity_mod
contains 


!Subroutine to propagate linear subdivision =========================
subroutine propagate_linear_sd(subd_system,mesh_base,options)
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
    write(*,'(A)') '--> propagating linear refinement'
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
    call construct_linear_refmat(subd_system(rr)%Ru,subd_system(rr)%mesh)

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
end subroutine propagate_linear_sd




!Construct linear refinement matrix subroutine =========================
subroutine construct_linear_refmat(Ru,mesh)
implicit none 

!Variables - Import
type(mesh_data) :: mesh
type(csrmatrix) :: Ru

!Variables - Local 
integer(in64) :: vv,ee,ff
integer(in64) :: nvertex_n,rins,nentry_rs,Rs_entry
integer(in64), dimension(:,:), allocatable :: Rs_ij
real(dp) :: Npfinv
real(dp), dimension(:), allocatable :: Rs_V

!Find number of seperate entries in refinement matrix ---
!Rs ------
!Original vertex contribution
nentry_rs = mesh%nvertex

!Edge vertex creation contribution
nentry_rs = nentry_rs + 2*mesh%nedge

!Face vertex creation contribution
nentry_rs = nentry_rs + sum(mesh%faces(:)%nvertex)

!New vertex quantity 
nvertex_n = mesh%nvertex + mesh%nedge + mesh%nface

!Rs sparse structure
allocate(Rs_ij(nentry_rs,2))
allocate(Rs_V(nentry_rs))
Rs_ij(:,:) = 0
Rs_V(:) = 0.0d0 

!Populate splitting matrix Rs (linear refinement matrix) -----------------------------------------
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

!Set final entry counts
Rs_entry = Rs_entry - 1

!Construct refinement matrix 
call build_csr_matrix(Ru,Rs_V(1:Rs_entry),Rs_ij(1:Rs_entry,1),Rs_ij(1:Rs_entry,2),nvertex_n,mesh%nvertex)
return 
end subroutine construct_linear_refmat


end module mrsys_linear_subdivide_mod