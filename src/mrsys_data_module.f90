!MRsys data types module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.4
!Updated 14-12-2023

!Module
module mrsys_data_mod

!Integer data type 
use ISO_FORTRAN_ENV, only: in32=>int32
use ISO_FORTRAN_ENV, only: in64=>int64

!Real data types 
use ISO_FORTRAN_ENV, only: sp=>real32     
use ISO_FORTRAN_ENV, only: dp=>real64     

!Sparse matrix and hash table types =======
!Hash table base allocation size
integer(in64) :: hashtab_bsize = 4

!Hash table entry
type hash_entry_data
    integer(in64) :: Nentry,size
    integer(in64), dimension(:), allocatable :: sparse_idx
    integer(in64), dimension(:), allocatable :: linear_idx
end type hash_entry_data

!Hash table type
type hash_table_data
    integer(in64) :: table_length
    integer(in64) :: NNZ,col_total
    type(hash_entry_data), dimension(:), allocatable :: table 
end type hash_table_data

!CSR matrix data type
type csrmatrix 
    integer(in64) :: nrow,ncol,nnz
    integer(in64), dimension(:), allocatable :: i,j
    real(dp), dimension(:), allocatable :: entry
end type csrmatrix 
!Sparse matrix and hash table types =======

!Options data type
type options_data 
    character(len=:), allocatable :: console_disp,mode,options_path,mesh_path,tgtmesh_path,data_path,subd_method,crs_base_vlnc
    character(len=:), allocatable :: mrsystemf,deformationf,mrcon_export_meshes,mrsystem_name,gradient_fname
    integer(in64) :: n_sd_levels,n_crs_levels_max,mrdeflevel
    real(dp) :: srcc_omega
end type options_data

!Face data type 
type face_data
    integer(in64) :: nvertex 
    integer(in64), dimension(:), allocatable :: vertices,edges
end type face_data

!Connectivity data type 
type connectivity_data
    integer(in64) :: max_valence
    integer(in64), dimension(:), allocatable :: valence
    integer(in64), dimension(:,:), allocatable :: v2v,v2e,v2f,e2f
end type connectivity_data

!Mesh data type
type mesh_data 
    integer(in64) :: ndim,nvertex,nedge,nface
    integer(in64), dimension(:), allocatable :: vertex_sharp
    integer(in64), dimension(:), allocatable :: vtxmap_fin2crs,vtxmap_crs2fin
    integer(in64), dimension(:,:), allocatable :: edges,vertex_felink2crs
    real(dp), dimension(:), allocatable :: edge_sharp
    real(dp), dimension(:,:), allocatable :: vertices,vertices_tgt
    type(face_data), dimension(:), allocatable :: faces
    type(connectivity_data) :: connectivity
end type mesh_data 

!Vertex ordering data type
type vtxodrdata
    integer(in64), dimension(:), allocatable :: vtx_idx_odr
end type vtxodrdata

!Subdivision system data type
type subdsystem_data
    integer(in64) :: n_d
    real(dp), dimension(:,:), allocatable :: d,surface_grad
    type(mesh_data) :: mesh 
    type(csrmatrix) :: Ru,Ri,Q,RuQ
end type subdsystem_data

!Multi-resolution system data type
type mressystem_data
    integer(in64) :: nlevel,nrefine
    real(dp), dimension(:,:), allocatable :: deformation
    type(csrmatrix), dimension(:), allocatable :: gradient_mat
    type(subdsystem_data), dimension(:), allocatable :: subd_system
end type mressystem_data

end module mrsys_data_mod