!MRsys energy minimised coarsening multiresolution module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.5
!Updated 05-02-2024

!Module
module mrsys_emc_multiresolution_mod
use mrsys_io_mod
use sparse_matrix_mod
contains 


!Evaluate surface position gradients subroutine =========================
subroutine emc_evaluate_surface_gradients(mres_system,level_2def)
implicit none 

!Variables - Import
integer(in64) :: level_2def
type(mressystem_data) :: mres_system
    
!Variables - Local 
integer(in64) :: vv 
real(dp) :: vec_response(mres_system%subd_system(level_2def)%mesh%nvertex)

!Allocate surface_grad structure for this level 
allocate(mres_system%subd_system(level_2def)%surface_grad(mres_system%subd_system(mres_system%nlevel)%mesh%nvertex,&
mres_system%subd_system(level_2def)%mesh%nvertex))

!Evaluate surface gradients for each vertex in the control net of this level 
vec_response(:) = 0.0d0 
do vv=1,mres_system%subd_system(level_2def)%mesh%nvertex
    vec_response(vv) = 1.0d0 
    call matmul_csr_dvec(mres_system%subd_system(level_2def)%surface_grad(:,vv),mres_system%gradient_mat(level_2def),vec_response)
    vec_response(vv) = 0.0d0 
end do 
return 
end subroutine emc_evaluate_surface_gradients




!Construct gradient matrix from level level_2def subroutine =========================
subroutine emc_build_gradient_matrix(mres_system,options,level_2def)
implicit none 

!Variables - Import
integer(in64) :: level_2def
type(mressystem_data) :: mres_system
type(options_data) :: options 

!Variables - Local 
integer(in64) :: rr 
real(dp) :: sparsity
type(csrmatrix) :: gradm_current,gradm_new

!Allocate gradient matrix structure if required 
if (.NOT.allocated(mres_system%gradient_mat)) then 
    allocate(mres_system%gradient_mat(mres_system%nrefine))
end if 

!Set current gradient matrix 
call csr_copy(gradm_current,mres_system%subd_system(level_2def)%Ru)

!Display
if (options%console_disp == 'yes') then
    sparsity = (real(gradm_current%nnz,dp)/&
    (real(gradm_current%nrow*gradm_current%ncol,dp)))*100.0d0 
    write(*,'(A,I0,A,I0,A,I0,A,A,A)') '    < refinement = ',level_2def,' - vertices = ',&
    mres_system%subd_system(level_2def)%mesh%nvertex,' -> ',mres_system%subd_system(level_2def+1)%mesh%nvertex,&
    ' - nnz fraction = ',real2F0_Xstring(sparsity,8_in64)//'%',' >'
end if 

!Project matrix though from level level_2def to the final refinement 
if (level_2def+1 .LE. mres_system%nrefine) then 
    do rr=level_2def+1,mres_system%nrefine

        !Multiply current gradient matrix by the refinement matrix at this level 
        call matmul_csr(gradm_new,mres_system%subd_system(rr)%Ru,gradm_current)

        !Free and update current gradient matrix 
        call csr_deallocate(gradm_current)
        call csr_copy(gradm_current,gradm_new)

        !Free new matrix 
        call csr_deallocate(gradm_new)

        !Display 
        if (options%console_disp == 'yes') then
            sparsity = (real(gradm_current%nnz,dp)/&
            (real(gradm_current%nrow*gradm_current%ncol,dp)))*100.0d0 
            write(*,'(A,I0,A,I0,A,I0,A,A,A)') '    < refinement = ',rr,' - vertices = ',mres_system%subd_system(rr)%mesh%nvertex,&
            ' -> ',mres_system%subd_system(rr+1)%mesh%nvertex,' - nnz fraction = ',real2F0_Xstring(sparsity,8_in64)//'%',' >'
        end if 
    end do 
end if 

!Store gradient matrix from level_2def
call csr_deallocate(mres_system%gradient_mat(level_2def))
call csr_copy(mres_system%gradient_mat(level_2def),gradm_current)
return 
end subroutine emc_build_gradient_matrix




!EMC multi-resolution apply deformation at level level_2def =========================
subroutine emc_multiresolution_deform(mres_system,deformation,level_2def)
implicit none 

!Variables - Import
integer(in64) :: level_2def
real(dp), dimension(:,:) :: deformation
type(mressystem_data) :: mres_system
    
!Apply deformation to target level 
mres_system%subd_system(level_2def)%mesh%vertices(:,:) = &
mres_system%subd_system(level_2def)%mesh%vertices(:,:) + deformation(:,:)
return 
end subroutine emc_multiresolution_deform




!EMC multi-resolution refinement from from level_2ref to the end surface subroutine =========================
subroutine emc_multiresolution_refine_2_end(mres_system,options,level_2initref)
implicit none 

!Variables - Import
integer(in64) :: level_2initref
type(mressystem_data) :: mres_system
type(options_data) :: options 
    
!Variables - Local 
integer(in64) :: rr 

!Refine 
do rr=level_2initref,mres_system%nrefine
    call emc_multiresolution_refine(mres_system,options,rr,1_in64)
end do 

!Apply final level detail correction 
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,1) = &
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,1) + mres_system%subd_system(mres_system%nlevel)%d(:,1)
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,2) = &
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,2) + mres_system%subd_system(mres_system%nlevel)%d(:,2)
if (mres_system%subd_system(mres_system%nlevel)%mesh%ndim == 3) then
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,3) = &
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,3) + mres_system%subd_system(mres_system%nlevel)%d(:,3)
end if 
return 
end subroutine emc_multiresolution_refine_2_end




!EMC multi-resolution refinement one level from level_2ref subroutine =========================
subroutine emc_multiresolution_refine(mres_system,options,level_2ref,disptoggle)
implicit none 

!Variables - Import
integer(in64) :: level_2ref,disptoggle
type(mressystem_data) :: mres_system
type(options_data) :: options 

!Variables - Local 
integer(in64) :: vdlen
real(dp), dimension(:,:), allocatable :: vertices_d

!Construct concatonated vertices and details vector
vdlen = mres_system%subd_system(level_2ref)%mesh%nvertex + mres_system%subd_system(level_2ref)%n_d 
allocate(vertices_d(vdlen,mres_system%subd_system(level_2ref)%mesh%ndim))
vertices_d(1:mres_system%subd_system(level_2ref)%mesh%nvertex,:) = mres_system%subd_system(level_2ref)%mesh%vertices(:,:)
vertices_d(mres_system%subd_system(level_2ref)%mesh%nvertex+1:vdlen,:) = mres_system%subd_system(level_2ref)%d(:,:)

!Perform refinement to the new level 
call matmul_csr_dvec(mres_system%subd_system(level_2ref+1)%mesh%vertices(:,1),&
mres_system%subd_system(level_2ref)%RuQ,vertices_d(:,1))
call matmul_csr_dvec(mres_system%subd_system(level_2ref+1)%mesh%vertices(:,2),&
mres_system%subd_system(level_2ref)%RuQ,vertices_d(:,2))
if (mres_system%subd_system(level_2ref)%mesh%ndim == 3) then 
    call matmul_csr_dvec(mres_system%subd_system(level_2ref+1)%mesh%vertices(:,3),&
    mres_system%subd_system(level_2ref)%RuQ,vertices_d(:,3))
end if 

!Display 
if ((options%console_disp == 'yes') .AND. (disptoggle == 1)) then
    write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '    < refinement = ',level_2ref,' -> ',level_2ref+1,' - vertices = ',&
    mres_system%subd_system(level_2ref)%mesh%nvertex,' -> ',mres_system%subd_system(level_2ref+1)%mesh%nvertex,' >'
end if 
return 
end subroutine emc_multiresolution_refine




!Construct emc detail vectors subroutine =========================
subroutine construct_emc_detail_vectors(mres_system,options)
implicit none 

!Variables - Import
type(mressystem_data) :: mres_system
type(options_data) :: options 

!Variables - Local 
integer(in64) :: rr
integer(in64) :: lend
real(dp) :: detail_norm
real(dp), dimension(:,:), allocatable :: vertices_ccr,delta

!For each refinement 
do rr=1,mres_system%nrefine

    !Perform normal catmull-clark refinement from the coarse level 
    allocate(vertices_ccr(mres_system%subd_system(rr+1)%mesh%nvertex,mres_system%subd_system(rr+1)%mesh%ndim))
    call matmul_csr_dvec(vertices_ccr(:,1),mres_system%subd_system(rr)%Ru,mres_system%subd_system(rr)%mesh%vertices(:,1))
    call matmul_csr_dvec(vertices_ccr(:,2),mres_system%subd_system(rr)%Ru,mres_system%subd_system(rr)%mesh%vertices(:,2))
    if (mres_system%subd_system(rr)%mesh%ndim == 3) then 
        call matmul_csr_dvec(vertices_ccr(:,3),mres_system%subd_system(rr)%Ru,mres_system%subd_system(rr)%mesh%vertices(:,3))
    end if

    !Calculate errors 
    allocate(delta(mres_system%subd_system(rr+1)%mesh%nvertex,mres_system%subd_system(rr+1)%mesh%ndim))
    delta(:,1) = mres_system%subd_system(rr+1)%mesh%vertices(:,1) - vertices_ccr(:,1)
    delta(:,2) = mres_system%subd_system(rr+1)%mesh%vertices(:,2) - vertices_ccr(:,2)
    if (mres_system%subd_system(rr)%mesh%ndim == 3) then
        delta(:,3) = mres_system%subd_system(rr+1)%mesh%vertices(:,3) - vertices_ccr(:,3)
    end if 

    !Construct detail vector for level rr
    lend = mres_system%subd_system(rr+1)%mesh%nvertex - mres_system%subd_system(rr)%mesh%nvertex
    if (.NOT.allocated(mres_system%subd_system(rr)%d)) then 
        mres_system%subd_system(rr)%n_d = lend
        allocate(mres_system%subd_system(rr)%d(lend,mres_system%subd_system(rr)%mesh%ndim))
    end if 
    mres_system%subd_system(rr)%d(:,1) = delta(mres_system%subd_system(rr)%mesh%nvertex+1:&
    mres_system%subd_system(rr+1)%mesh%nvertex,1)
    mres_system%subd_system(rr)%d(:,2) = delta(mres_system%subd_system(rr)%mesh%nvertex+1:&
    mres_system%subd_system(rr+1)%mesh%nvertex,2)
    if (mres_system%subd_system(rr)%mesh%ndim == 3) then
        mres_system%subd_system(rr)%d(:,3) = delta(mres_system%subd_system(rr)%mesh%nvertex+1:&
        mres_system%subd_system(rr+1)%mesh%nvertex,3)
    end if 
    if (mres_system%subd_system(rr)%mesh%ndim == 2) then
        detail_norm = sqrt(norm2(mres_system%subd_system(rr)%d(:,1))**2 + norm2(mres_system%subd_system(rr)%d(:,2))**2)
    elseif (mres_system%subd_system(rr)%mesh%ndim == 3) then
        detail_norm = sqrt(norm2(mres_system%subd_system(rr)%d(:,1))**2 + &
        norm2(mres_system%subd_system(rr)%d(:,2))**2 + norm2(mres_system%subd_system(rr)%d(:,3))**2)
    end if 
    if (abs(detail_norm) .GE. 1e-40) then 
        detail_norm = log10(abs(detail_norm))
    else 
        detail_norm = abs(detail_norm)
    end if 

    !Consrtuct actual refined level rr+1 vertices 
    call emc_multiresolution_refine(mres_system,options,rr,0_in64)

    !Deallocate locals for this level 
    deallocate(delta)
    deallocate(vertices_ccr)

    !Display
    if (options%console_disp == 'yes') then
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,A,A)') '    < refinement = ',rr,' -> ',rr+1,' - vertices = ',&
        mres_system%subd_system(rr)%mesh%nvertex,' -> ',mres_system%subd_system(rr+1)%mesh%nvertex,&
        ' - lognorm(details) = ',real2F0_Xstring(detail_norm,8_in64),' >'
    end if 
end do 

!Construct final level detail correction 
lend = mres_system%subd_system(mres_system%nlevel)%mesh%nvertex
if (.NOT.allocated(mres_system%subd_system(mres_system%nlevel)%d)) then 
    mres_system%subd_system(mres_system%nlevel)%n_d = lend
    allocate(mres_system%subd_system(mres_system%nlevel)%d(lend,mres_system%subd_system(mres_system%nlevel)%mesh%ndim))
end if 
mres_system%subd_system(mres_system%nlevel)%d(:,1) = mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,1) - &
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,1)
mres_system%subd_system(mres_system%nlevel)%d(:,2) = mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,2) - &
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,2)
if (mres_system%subd_system(mres_system%nlevel)%mesh%ndim == 3) then
    mres_system%subd_system(mres_system%nlevel)%d(:,3) = mres_system%subd_system(mres_system%nlevel)%mesh%vertices_tgt(:,3) - &
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,3)
end if 
if (mres_system%subd_system(mres_system%nlevel)%mesh%ndim == 2) then
    detail_norm = sqrt(norm2(mres_system%subd_system(mres_system%nlevel)%d(:,1))**2 + &
    norm2(mres_system%subd_system(mres_system%nlevel)%d(:,2))**2)
elseif (mres_system%subd_system(mres_system%nlevel)%mesh%ndim == 3) then
    detail_norm = sqrt(norm2(mres_system%subd_system(mres_system%nlevel)%d(:,1))**2 + &
    norm2(mres_system%subd_system(mres_system%nlevel)%d(:,2))**2 + norm2(mres_system%subd_system(mres_system%nlevel)%d(:,3))**2)
end if
if (abs(detail_norm) .GE. 1e-40) then 
    detail_norm = log10(abs(detail_norm))
else 
    detail_norm = abs(detail_norm)
end if 

!Apply final level correction 
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,1) = &
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,1) + mres_system%subd_system(mres_system%nlevel)%d(:,1)
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,2) = &
mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,2) + mres_system%subd_system(mres_system%nlevel)%d(:,2)
if (mres_system%subd_system(mres_system%nlevel)%mesh%ndim == 3) then
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,3) = &
    mres_system%subd_system(mres_system%nlevel)%mesh%vertices(:,3) + mres_system%subd_system(mres_system%nlevel)%d(:,3)
end if 

!Display
if (options%console_disp == 'yes') then
    write(*,'(A,A,A)') '    < final error correction lognorm = ',real2F0_Xstring(detail_norm,8_in64),'>'
end if 
return 
end subroutine construct_emc_detail_vectors




!Construct emc detail matricies subroutine =========================
subroutine construct_emc_detail_matricies(mres_system,options)
implicit none 

!Variables - Import
type(mressystem_data) :: mres_system
type(options_data) :: options 

!Variables - Local 
integer(in64) :: rr
real(dp) :: sparsity

!For each refinement 
do rr=1,mres_system%nrefine

    !Construct detail matrix for refinement rr 
    call build_emc_detail_matrix(mres_system%subd_system(rr)%Q,mres_system%subd_system(rr)%mesh,&
    mres_system%subd_system(rr+1)%mesh,options)

    !Display
    if (options%console_disp == 'yes') then
        sparsity = (real(mres_system%subd_system(rr)%Q%nnz,dp)/&
        (real(mres_system%subd_system(rr)%Q%nrow*mres_system%subd_system(rr)%Q%ncol,dp)))*100.0d0
        write(*,'(A,I0,A,I0,A,I0,A,A,A)') '    < refinement = ',rr,' - vertices = ',mres_system%subd_system(rr)%mesh%nvertex,&
        ' -> ',mres_system%subd_system(rr+1)%mesh%nvertex,' - nnz fraction = ',real2F0_Xstring(sparsity,8_in64)//'%',' >'
    end if 
end do 
return 
end subroutine construct_emc_detail_matricies




!Construct emc detail matrix subroutine =========================
subroutine build_emc_detail_matrix(Q,mesh_coarse,mesh_fine,options)
implicit none 

!Variables - Import
type(mesh_data) :: mesh_coarse,mesh_fine
type(csrmatrix) :: Q
type(options_data) :: options

!Variables - Local 
integer(in64) :: vv,ee
integer(in64) :: Q_entry,nentry_q,ncol,nrow,fvidx,ftgt,etgt,n_sharp,sinsert,v1sr,v2sr
integer(in64) :: EsharpI(2)
integer(in64), dimension(:,:), allocatable :: Q_ij
real(dp) :: k,alphae,alphaf,omega
real(dp), dimension(:), allocatable :: Q_V

!Extract omega smoothing parameter 
omega = options%srcc_omega

!Upper bound of non zero elements
nrow = mesh_fine%nvertex 
ncol = mesh_fine%nvertex - mesh_coarse%nvertex
nentry_q = ncol + 3*sum(mesh_fine%connectivity%valence(1:mesh_coarse%nvertex))

!Q sparse structure
allocate(Q_ij(nentry_q,2))
allocate(Q_V(nentry_q))
Q_ij(:,:) = 0
Q_V(:) = 0.0d0 

!Populate detail matrix Q for level rr-1 vertex details
Q_entry = 1
do vv=1,mesh_coarse%nvertex

    !Extract valence
    k = mesh_coarse%connectivity%valence(vv) 

    !Find number of sharp incident edges to identify vertex type on refined level
    n_sharp = 0
    do ee=1,mesh_fine%connectivity%valence(vv)
        etgt = mesh_fine%connectivity%v2e(vv,ee)
        if (mesh_fine%edge_sharp(etgt) .GT. 0.0d0) then
            n_sharp = n_sharp + 1
        end if
    end do 

    !If 2 sharp edges the locate these edges and their repsective edge vertices on the fine level 
    if (n_sharp == 2) then 
        EsharpI(:) = 0
        sinsert = 0
        do ee=1,mesh_fine%connectivity%valence(vv)
            etgt = mesh_fine%connectivity%v2e(vv,ee)
            if (mesh_fine%edge_sharp(etgt) .GT. 0.0d0) then
                sinsert = sinsert + 1
                EsharpI(sinsert) = etgt
            end if 
            if (sinsert == 2) then 
                exit 
            end if 
        end do 
        if (mesh_fine%edges(EsharpI(1),1) == vv) then 
            v1sr = mesh_fine%edges(EsharpI(1),2)
        else
            v1sr = mesh_fine%edges(EsharpI(1),1)
        end if 
        if (mesh_fine%edges(EsharpI(2),1) == vv) then 
            v2sr = mesh_fine%edges(EsharpI(2),2)
        else
            v2sr = mesh_fine%edges(EsharpI(2),1)
        end if 
    end if 

    !Cases
    if ((n_sharp .LT. 2) .AND. (mesh_coarse%vertex_sharp(vv) == 0)) then !Normal 

        !alpha e and f
        alphae = 4.0d0/(k*k)
        alphaf = -1.0d0/(k*k)

        !Construct entries of the detail matrix ---
        !Edge vertices 
        do ee=1,mesh_coarse%connectivity%valence(vv) 
            Q_ij(Q_entry,1) = vv
            Q_ij(Q_entry,2) = mesh_coarse%connectivity%v2e(vv,ee)  !Exclude mesh_coarse%nvertex from column position as this is added through concatonation with Ru
            Q_V(Q_entry) = alphae + ((1.0d0 - omega)/(k*k)) !as omega decreases the coarsening adds dependance on each ek of (1-omega)/k, so here we add this dependance averaged over all k ek
            Q_entry = Q_entry + 1
        end do 

        !Face vertices 
        do ee=1,mesh_coarse%connectivity%valence(vv) 
            ftgt = mesh_coarse%connectivity%v2f(vv,ee)
            if (ftgt .GT. 0) then 
                fvidx = mesh_coarse%nedge + ftgt !Exclude mesh_coarse%nvertex from column position as this is added through concatonation with Ru
                Q_ij(Q_entry,1) = vv
                Q_ij(Q_entry,2) = fvidx
                Q_V(Q_entry) = alphaf 
                Q_entry = Q_entry + 1
            end if 
        end do 
    elseif ((n_sharp == 2) .AND. (mesh_coarse%vertex_sharp(vv) == 0)) then !On a sharp edge/crease

        !Vertex 1
        Q_ij(Q_entry,1) = vv
        Q_ij(Q_entry,2) = v1sr - mesh_coarse%nvertex
        Q_V(Q_entry) = 0.5d0 
        Q_entry = Q_entry + 1

        !Vertex 2
        Q_ij(Q_entry,1) = vv
        Q_ij(Q_entry,2) = v2sr - mesh_coarse%nvertex
        Q_V(Q_entry) = 0.5d0 
        Q_entry = Q_entry + 1
    elseif ((n_sharp .GT. 2) .OR. (mesh_coarse%vertex_sharp(vv) == 1)) then !Sharp vertex (implied || specified)  
        Q_ij(Q_entry,1) = vv
        Q_ij(Q_entry,2) = vv - mesh_coarse%nvertex
        Q_V(Q_entry) = 0.0d0 !dummy entry to force correct row structure in the csr matrix format
        Q_entry = Q_entry + 1
    end if 
end do 

!Populate matrix for level rr vertex details
do vv=mesh_coarse%nvertex+1,mesh_fine%nvertex 
    Q_ij(Q_entry,1) = vv
    Q_ij(Q_entry,2) = vv - mesh_coarse%nvertex !Subtract mesh_coarse%nvertex from column position as this is added through concatonation with Ru
    Q_V(Q_entry) = 1.0d0 
    Q_entry = Q_entry + 1
end do 

!Set final entry count
Q_entry = Q_entry - 1

!Construct sparse Q 
call build_csr_matrix(Q,Q_V(1:Q_entry),Q_ij(1:Q_entry,1),Q_ij(1:Q_entry,2),nrow,ncol)
return 
end subroutine build_emc_detail_matrix



! !Construct emc detail matrix subroutine =========================               UPDATE for sharp vertices ===================================
! subroutine build_emc_detail_matrix(Q,mesh_coarse,mesh_fine)
! implicit none 

! !Variables - Import
! type(mesh_data) :: mesh_coarse,mesh_fine
! type(csrmatrix) :: Q

! !Variables - Local 
! integer(in64) :: vv,ee
! integer(in64) :: Q_entry,nentry_q,ncol,nrow,fvidx,ftgt
! integer(in64), dimension(:,:), allocatable :: Q_ij
! real(dp) :: k
! real(dp), dimension(:), allocatable :: Q_V

! !Upper bound of non zero elements
! nrow = mesh_fine%nvertex 
! ncol = mesh_fine%nvertex - mesh_coarse%nvertex
! nentry_q = ncol + 2*sum(mesh_fine%connectivity%valence(1:mesh_coarse%nvertex))

! !Q sparse structure
! allocate(Q_ij(nentry_q,2))
! allocate(Q_V(nentry_q))
! Q_ij(:,:) = 0
! Q_V(:) = 0.0d0 

! !Populate detail matrix Q for level rr-1 vertex details
! Q_entry = 1
! do vv=1,mesh_coarse%nvertex

!     !Extract valence
!     k = mesh_coarse%connectivity%valence(vv) 

!     !Refined level edge attached verticies
!     do ee=1,mesh_coarse%connectivity%valence(vv) 
!         Q_ij(Q_entry,1) = vv
!         Q_ij(Q_entry,2) = mesh_fine%connectivity%v2v(vv,ee) - mesh_coarse%nvertex !Subtract mesh_coarse%nvertex from column position as this is added through concatonation with Ru
!         Q_V(Q_entry) = 4.0d0/(k*k)
!         Q_entry = Q_entry + 1
!     end do 

!     !Refined level face attached verticies
!     do ee=1,mesh_coarse%connectivity%valence(vv) 
!         ftgt = mesh_coarse%connectivity%v2f(vv,ee)
!         if (ftgt .GT. 0) then 
!             fvidx = mesh_coarse%nedge + ftgt !Exclude mesh_coarse%nvertex from column position as this is added through concatonation with Ru
!             Q_ij(Q_entry,1) = vv
!             Q_ij(Q_entry,2) = fvidx
!             Q_V(Q_entry) = -1.0d0/(k*k)
!             Q_entry = Q_entry + 1
!         end if 
!     end do 
! end do 

! !Populate matrix for level rr vertex details
! do vv=mesh_coarse%nvertex+1,mesh_fine%nvertex 
!     Q_ij(Q_entry,1) = vv
!     Q_ij(Q_entry,2) = vv - mesh_coarse%nvertex !Subtract mesh_coarse%nvertex from column position as this is added through concatonation with Ru
!     Q_V(Q_entry) = 1.0d0 
!     Q_entry = Q_entry + 1
! end do 

! !Set final entry count
! Q_entry = Q_entry - 1

! !Construct sparse Q 
! call build_csr_matrix(Q,Q_V(1:Q_entry),Q_ij(1:Q_entry,1),Q_ij(1:Q_entry,2),nrow,ncol)
! return 
! end subroutine build_emc_detail_matrix


end module mrsys_emc_multiresolution_mod