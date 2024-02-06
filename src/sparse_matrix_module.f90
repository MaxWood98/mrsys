!Sparse Matrix Operation Subroutines Module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 4.0
!Updated 14-12-2023

!Sparse matrix opteration module
module sparse_matrix_mod

!Use Modules 
use mrsys_data_mod
use sparse_idx_hash_mod

!Data types 
! type csrmatrix 
!     integer(in64) :: nrow,ncol
!     integer(in64), dimension(:), allocatable :: i,j
!     real(dp), dimension(:), allocatable :: entry
! end type csrmatrix 

!Routines
contains 


!Subroutine to deallocate a csr matrix ===================================================
subroutine csr_deallocate(mat)
implicit none 

!Variables - Import
type(csrmatrix) :: mat

!Reset values
mat%nnz = 0 
mat%nrow = 0
mat%ncol = 0

!Dellaocte
if (allocated(mat%entry)) then 
    deallocate(mat%entry)
end if 
if (allocated(mat%j)) then 
    deallocate(mat%j)
end if 
if (allocated(mat%i)) then 
    deallocate(mat%i)
end if 
return 
end subroutine csr_deallocate




!Subroutine to copy a csr matrix ===================================================
subroutine csr_copy(mat_cpyto,mat_cpyfrom)
implicit none 

!Variables - Import
type(csrmatrix) :: mat_cpyto,mat_cpyfrom

!Assign matrix sizes
mat_cpyto%nnz = mat_cpyfrom%nnz
mat_cpyto%nrow = mat_cpyfrom%nrow 
mat_cpyto%ncol = mat_cpyfrom%ncol 
allocate(mat_cpyto%entry(mat_cpyfrom%nnz))
allocate(mat_cpyto%j(mat_cpyfrom%nnz))
allocate(mat_cpyto%i(mat_cpyfrom%nrow+1))

!Copy entry data 
mat_cpyto%entry = mat_cpyfrom%entry
mat_cpyto%j = mat_cpyfrom%j
mat_cpyto%i = mat_cpyfrom%i
return 
end subroutine csr_copy




!Subroutine to construct CSR matrix from list of entries ===================================================
subroutine build_csr_matrix(mat,entry_val,i,j,nrow,ncol,collate_entries)
implicit none 

!Variables - Import
integer(in64) :: nrow,ncol
integer(in64), optional :: collate_entries
integer(in64), dimension(:) :: i,j
real(dp), dimension(:) :: entry_val
type(csrmatrix) :: mat

!Variables - Local 
integer(in64) :: ii
integer(in64) :: nentry,rowc

!Evaluate number of entries nentry and collate and order row wise if requested
if (present(collate_entries)) then 
    if (collate_entries == 1) then 

    else 
        nentry = size(entry_val,1)
    end if 
else
    nentry = size(entry_val,1)
end if 

!Assign matrix sizes
mat%nnz = nentry
mat%nrow = nrow 
mat%ncol = ncol 
allocate(mat%entry(nentry))
allocate(mat%j(nentry))
allocate(mat%i(nrow+1))

!Assign entries and columns
mat%entry(:) = entry_val(:)
mat%j(:) = j(:)

!Assign row pointers
rowc = i(1)
mat%i(:) = 0
mat%i(1) = 1
do ii=2,nentry
    if ((i(ii) .NE. rowc) .AND. (i(ii-1) == rowc)) then !update row
        mat%i(rowc+1) = ii
        rowc = rowc + 1
    end if 
end do 
mat%i(rowc+1) = nentry + 1
return 
end subroutine build_csr_matrix




!Wrapper subroutine to multiply two CSR matricies (C = A*B) ===================================================
subroutine matmul_csr(mat_C,mat_A,mat_B)
implicit none 

!Variables - Import 
type(csrmatrix) :: mat_A,mat_B,mat_C

!Variables - Local 
integer(in64) :: nrowC,ncolC,ierr,nzmax,nnz_c
integer(in64) :: iw(mat_A%ncol)
integer(in64), dimension(:), allocatable :: ic,jc
real(dp), dimension(:), allocatable :: c_entry

!Size of matrix c
nrowC = mat_A%nrow
ncolC = mat_B%ncol

!Set maximum size of C
nzmax = 2*(mat_A%nnz + mat_B%nnz) 

!Intialise C arrays 
allocate(c_entry(nzmax))
allocate(ic(nzmax))
allocate(jc(nzmax))

!Call multiplication subroutine 
call csrmucsr(nrowC,ncolC,mat_A%entry,mat_A%j,mat_A%i,mat_B%entry,mat_B%j,mat_B%i,c_entry,jc,ic,nnz_c,nzmax,iw,ierr)

!Assign output matrix sizes
mat_C%nnz = nnz_c
mat_C%nrow = nrowC 
mat_C%ncol = ncolC 
allocate(mat_C%entry(nnz_c))
allocate(mat_C%j(nnz_c))
allocate(mat_C%i(nrowC+1))

!Extract matrix components 
mat_C%entry(:) = c_entry(1:nnz_c)
mat_C%j(:) = jc(1:nnz_c)
mat_C%i(:) = ic(1:nrowC+1)
return 
end subroutine matmul_csr




!Subroutine to multiply a sparse csr matrix by a dense vector ===================================================
subroutine matmul_csr_dvec(vec_out,mat,vec_in)
implicit none 

!Variables - Import 
real(dp), dimension(:) :: vec_in,vec_out
type(csrmatrix) :: mat

!Variables - Local 
integer(in64) :: ii,jj

!Multiply
vec_out(:) = 0.0d0 
do ii=1,mat%nrow
    do jj=mat%i(ii),mat%i(ii+1)-1
        vec_out(ii) = vec_out(ii) + mat%entry(jj)*vec_in(mat%j(jj))
    end do 
end do 
return 
end subroutine matmul_csr_dvec




!Subroutine to perform the multiplication of two CSR matricies ===================================================
!(based on the method in sparsekit.f90)
subroutine csrmucsr(nrow,ncol,a,ja,ia,b,jb,ib,c,jc,ic,nnz_c,nzmax,work,ierr)
implicit none 

!Variables - Import 
integer(in64) :: nnz_c,ncol,nrow,nzmax,ierr
integer(in64) :: ia(nrow+1),ib(nrow+1),ic(nrow+1),jc(nzmax),work(ncol)
integer(in64), dimension(:) :: ja,jb
real(dp) :: c(nzmax)
real(dp), dimension(:) :: a,b

!Variables - Local 
integer(in64) :: ii,jj,kk
integer(in64) :: jcol,jpos,ka,kb,entry_ins
real(dp) :: val_a

!Initialise entries of matrix c
ic(:) = 0
jc(:) = 0
ic(1) = 1
entry_ins = 0

!Initialise error state and work array 
ierr = 0
work(:) = 0

!Multiply 
do ii=1,nrow !each row of a
    do ka=ia(ii),ia(ii+1)-1 !each entry in this row of a

        !Entry value in a 
        val_a = a(ka)

        !Column in a
        jj = ja(ka) 

        !Each entry in this row of b
        do kb=ib(jj),ib(jj+1)-1 

            !Column in b
            jcol = jb(kb) 

            !Entry location in c
            jpos = work(jcol) 

            !Accumulate to c
            if (jpos == 0) then !if no entry for this location in c (in c -> row = ii col = jj, entry in c arrays jpos)
                entry_ins = entry_ins + 1 !increment length 
                if (nzmax .LT. entry_ins) then !fail if nnz max exceeded
                    ierr = ii
                    return
                end if
                jc(entry_ins) = jcol !store column in c
                work(jcol)= entry_ins !store entry position 
                c(entry_ins) = val_a * b(kb) !store value 
            else !if entry already exists add value to this entry 
                c(jpos) = c(jpos) + val_a * b(kb)
            end if
        end do
    end do

    !Reset work array 
    do kk=ic(ii),entry_ins
        work(jc(kk)) = 0
    end do

    !Add pointer to the start of the next row in c
    ic(ii+1) = entry_ins + 1
end do

!Assign number of non zero elements in the result matrix c
nnz_c = entry_ins
return
end subroutine csrmucsr




!Subroutine to horizontally concatonate two sparse matricies in the form C = [A | B] ===================================================
subroutine csrhrzcat(mat_C,mat_A,mat_B)
implicit none 

!Variables - Import 
type(csrmatrix) :: mat_A,mat_B,mat_C

!Variables - Local
integer(in64) :: rr,cc
integer(in64) :: nrowC,ncolC,nentryC,cins
integer(in64), dimension(:,:), allocatable :: C_ij
real(dp), dimension(:), allocatable :: C_V

!Set the size of C
nrowC = mat_A%nrow
ncolC = mat_A%ncol + mat_B%ncol
nentryC = mat_A%nnz + mat_B%nnz 

!Allocate components of C
allocate(C_ij(nentryC,2))
allocate(C_V(nentryC))

!Populate C row-wise 
cins = 0 
do rr=1,mat_A%nrow

    !From A
    do cc=mat_A%i(rr),mat_A%i(rr+1)-1 !each entry in this row of A
        cins = cins + 1
        C_V(cins) = mat_A%entry(cc)
        C_ij(cins,1) = rr 
        C_ij(cins,2) = mat_A%j(cc)
    end do 

    !From B
    do cc=mat_B%i(rr),mat_B%i(rr+1)-1 !each entry in this row of A
        cins = cins + 1
        C_V(cins) = mat_B%entry(cc)
        C_ij(cins,1) = rr 
        C_ij(cins,2) = mat_B%j(cc) + mat_A%ncol !shift across by the number of columns in A
    end do 
end do 

!Construct matrix c
call build_csr_matrix(mat_C,C_V,C_ij(:,1),C_ij(:,2),nrowC,ncolC)
return 
end subroutine csrhrzcat




! !Subroutine to collate a list of entries into a sparse matrix ===================================================
! subroutine collate_sparse_entires(NNZ,Nentry,entry_idx,entry_val,Ncol)
! implicit none 

! !Variables - Import
! integer(in64) :: Nentry,Ncol,NNZ
! integer(in64), dimension(:,:), allocatable :: entry_idx
! real(dp), dimension(:), allocatable :: entry_val

! !Variables - Local 
! integer(in64) :: ii
! integer(in64) :: sparse_idx
! integer(in64) :: entry_map(Nentry)
! integer(in64) :: entry_idx_col(Nentry,2)
! real(dp) :: entry_val_col(Nentry)
! type(sparse_mat_data) :: mat

! !Allocate temporary entry hash structure 
! mat%entry_hash%table_length = Nentry + 1
! mat%entry_hash%NNZ = Nentry
! mat%entry_hash%col_total = Ncol
! allocate(mat%entry_hash%table(Nentry+1))
! mat%entry_hash%table(:)%Nentry = 0
! mat%entry_hash%table(:)%size = 0

! !Count non-zero entries and map multiple entries to the same spare index where required
! NNZ = 0 
! entry_map(:) = 0 
! do ii=1,Nentry
!     call hash_table_retrive(sparse_idx,entry_idx(ii,1),entry_idx(ii,2),mat%entry_hash)
!     if (sparse_idx .LT. 0) then !no entry here -> increment non-zero count 
!         NNZ = NNZ + 1
!         call hash_table_add(NNZ,entry_idx(ii,1),entry_idx(ii,2),mat%entry_hash)
!         entry_map(ii) = NNZ
!     else !entry here -> map this entry to this sparse index 
!         entry_map(ii) = sparse_idx
!     end if 
! end do 

! !Collapse list by combining entries at the same indecies 
! entry_val_col(:) = 0 
! entry_idx_col(:,:) = 0 
! do ii=1,Nentry
!     entry_val_col(entry_map(ii)) = entry_val_col(entry_map(ii)) + entry_val(ii)
!     entry_idx_col(entry_map(ii),:) = entry_idx(ii,:)
! end do 

! !Rebuild collated entry lists 
! deallocate(entry_val)
! deallocate(entry_idx)
! allocate(entry_val(NNZ))
! allocate(entry_idx(NNZ,2))
! entry_val(:) = entry_val_col(1:NNZ)
! entry_idx(:,:) = entry_idx_col(1:NNZ,:)
! return
! end subroutine collate_sparse_entires

end module sparse_matrix_mod