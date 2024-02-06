!Hash table for sparse matrix indexing module
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 1.1
!Updated 28-01-2022

!Module
module sparse_idx_hash_mod

!Define data types
use mrsys_data_mod

!Subroutines 
contains 


!Subroutine to retrive item from table =====================
subroutine hash_table_retrive(sparse_idx,row,col,hash_table)
implicit none 

!Variables - Import
integer(in64) :: row,col,sparse_idx
type(hash_table_data) :: hash_table

!Variables - Local
integer(in64) :: hashval,linear_idx,NNZ,col_total,ii 

!Find item in table with same hash and linear index 
sparse_idx = -1 !Default to no entry here 
NNZ = hash_table%NNZ
col_total = hash_table%col_total
hashval = hash_val(row,col,NNZ,col_total)
linear_idx = col + (row - 1)*col_total
if(hash_table%table(hashval)%Nentry .GT. 0) then !Find entry
    if (hash_table%table(hashval)%Nentry .GT. 1) then !Multiple entries for this hash
        do ii=1,hash_table%table(hashval)%Nentry
            if (hash_table%table(hashval)%linear_idx(ii) == linear_idx) then 
                sparse_idx = hash_table%table(hashval)%sparse_idx(ii)
                exit
            end if
            if (ii == hash_table%table(hashval)%Nentry) then 
                sparse_idx = -1
                !print *, '** no entries at this hash for given linear index **'
            end if 
        end do 
    else !1 entry for this hash 
        if (hash_table%table(hashval)%linear_idx(1) == linear_idx) then !If referencing the same matrix entry
            sparse_idx = hash_table%table(hashval)%sparse_idx(1)
        else
            sparse_idx = -1
            !print *, '** no entries at this hash for given linear index **'
        end if
    end if
else !No entries for this input
    sparse_idx = -1
    !print *, '** no entries for this hash **'
end if
return     
end subroutine hash_table_retrive




!Interface to add item to table =====================
subroutine hash_table_add(sparse_idx,row,col,hash_table)
implicit none 

!Variables - Import
integer(in64) :: row,col,sparse_idx
type(hash_table_data) :: hash_table

!Variables - Local
integer(in64) :: hashval,linear_idx,NNZ,col_total

!Add to table
NNZ = hash_table%NNZ
col_total = hash_table%col_total
hashval = hash_val(row,col,NNZ,col_total)
linear_idx = col + (row - 1)*col_total
call hash_insert(sparse_idx,linear_idx,hashval,hash_table)
return 
end subroutine hash_table_add




!Subroutine to add item to table =====================
subroutine hash_insert(sparse_idx,linear_idx,hashval,hash_table)
implicit none 

!Variables - Import
integer(in64) :: sparse_idx,linear_idx,hashval
type(hash_table_data) :: hash_table

!Variables - Local 
integer(in64) :: Nentry,Nentry_N,ii,exist,sizeN 
integer(in64), dimension(:), allocatable :: entry_temp

!Check if any entries at this hash value to select addition case 
if (hash_table%table(hashval)%Nentry .NE. 0) then !Existing entries here 

    !Check if this entry already exists 
    exist = 0
    do ii=1,hash_table%table(hashval)%Nentry
        if (hash_table%table(hashval)%sparse_idx(ii) == sparse_idx) then
            exist = 1
        end if
    end do 

    !Add to table 
    if (exist == 0) then 

        !Current and new number of entries 
        Nentry = hash_table%table(hashval)%Nentry
        Nentry_N = Nentry + 1

        !If new entry exceeds the size at this entry then reallocate else just add
        if (Nentry_N .GT. hash_table%table(hashval)%size) then 

            !Set new size 
            sizeN = 2*hash_table%table(hashval)%size 
            ! print *, hash_table%table(hashval)%size,'->',sizeN
            
            !Copy new data 
            allocate(entry_temp(hash_table%table(hashval)%Nentry))
            entry_temp(:) = hash_table%table(hashval)%sparse_idx(:)
            deallocate(hash_table%table(hashval)%sparse_idx)
            allocate(hash_table%table(hashval)%sparse_idx(sizeN))
            hash_table%table(hashval)%sparse_idx(1:Nentry) = entry_temp(:)
            entry_temp(:) = hash_table%table(hashval)%linear_idx(:)
            deallocate(hash_table%table(hashval)%linear_idx)
            allocate(hash_table%table(hashval)%linear_idx(sizeN))
            hash_table%table(hashval)%linear_idx(1:Nentry) = entry_temp(:)
            hash_table%table(hashval)%sparse_idx(Nentry_N) = sparse_idx
            hash_table%table(hashval)%linear_idx(Nentry_N) = linear_idx
            hash_table%table(hashval)%Nentry = Nentry_N
            hash_table%table(hashval)%size = sizeN
            deallocate(entry_temp)

            !Initialise the new region 
            hash_table%table(hashval)%sparse_idx(Nentry_N+1:sizeN) = 0
            hash_table%table(hashval)%linear_idx(Nentry_N+1:sizeN) = 0
        else
            hash_table%table(hashval)%sparse_idx(Nentry_N) = sparse_idx
            hash_table%table(hashval)%linear_idx(Nentry_N) = linear_idx
            hash_table%table(hashval)%Nentry = Nentry_N
        end if 
    end if 
else !No entries at this hash 
    allocate(hash_table%table(hashval)%sparse_idx(hashtab_bsize))
    allocate(hash_table%table(hashval)%linear_idx(hashtab_bsize))
    hash_table%table(hashval)%sparse_idx(1) = sparse_idx
    hash_table%table(hashval)%linear_idx(1) = linear_idx
    hash_table%table(hashval)%sparse_idx(2:hashtab_bsize) = 0
    hash_table%table(hashval)%linear_idx(2:hashtab_bsize) = 0
    hash_table%table(hashval)%Nentry = 1
    hash_table%table(hashval)%size = hashtab_bsize
end if
return 
end subroutine hash_insert




!Hash function =====================
function hash_val(row,col,NNZ,col_total) result(hashval)
implicit none 

!Variables - Import
integer(in64) :: row,col,NNZ,col_total

!Variables - Local 
integer(in64) :: il,hashval

!Hash value 
il = col + (row - 1)*col_total
hashval = mod(il,NNZ) + 1
return     
end function hash_val 

end module sparse_idx_hash_mod