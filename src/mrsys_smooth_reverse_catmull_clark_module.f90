!MRsys smooth reverse catmull_clark module 
!Max Wood - mw16116@bristol.ac.uk
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 0.2
!Updated 12-12-2023

!Module
module mrsys_smooth_reverse_catmull_clark_mod
use mrsys_io_mod
use sparse_matrix_mod
contains 


!Subroutine to perform smooth reverse Catmull-Clark subdivision of a multiresolution system =========================
subroutine propagate_mrsys_sr_catmull_clark(mres_system,options)
implicit none 

!Variables - Import
type(mressystem_data) :: mres_system
type(options_data) :: options 

!Variables - Local 
integer(in64) :: rr 
integer(in64) :: level_coarse,level_fine
real(dp) :: deltanorm 
real(dp), dimension(:,:), allocatable :: vertices_apprx,def,delta

!For each coarsening 
do rr=1,mres_system%nrefine

    !Index of the coarse and fine level 
    level_coarse = mres_system%nlevel - rr 
    level_fine = level_coarse + 1

    !Perform trial coarening with reverse catmull-clark
    call matmul_csr_dvec(mres_system%subd_system(level_coarse)%mesh%vertices(:,1),&
    mres_system%subd_system(level_fine)%Ri,mres_system%subd_system(level_fine)%mesh%vertices(:,1))
    call matmul_csr_dvec(mres_system%subd_system(level_coarse)%mesh%vertices(:,2),&
    mres_system%subd_system(level_fine)%Ri,mres_system%subd_system(level_fine)%mesh%vertices(:,2))
    if (mres_system%subd_system(level_coarse)%mesh%ndim == 3) then 
        call matmul_csr_dvec(mres_system%subd_system(level_coarse)%mesh%vertices(:,3),&
        mres_system%subd_system(level_fine)%Ri,mres_system%subd_system(level_fine)%mesh%vertices(:,3))
    end if 

    !Foreward subdivide to approximate the fine level 
    allocate(vertices_apprx(mres_system%subd_system(level_fine)%mesh%nvertex,mres_system%subd_system(level_fine)%mesh%ndim))
    call matmul_csr_dvec(vertices_apprx(:,1),mres_system%subd_system(level_coarse)%Ru,&
    mres_system%subd_system(level_coarse)%mesh%vertices(:,1))
    call matmul_csr_dvec(vertices_apprx(:,2),mres_system%subd_system(level_coarse)%Ru,&
    mres_system%subd_system(level_coarse)%mesh%vertices(:,2))
    if (mres_system%subd_system(level_coarse)%mesh%ndim == 3) then 
        call matmul_csr_dvec(vertices_apprx(:,3),mres_system%subd_system(level_coarse)%Ru,&
        mres_system%subd_system(level_coarse)%mesh%vertices(:,3))
    end if 

    !Find def for the approximated refined level vertices 
    allocate(def(mres_system%subd_system(level_fine)%mesh%nvertex,mres_system%subd_system(level_fine)%mesh%ndim))
    def(:,1) = mres_system%subd_system(level_fine)%mesh%vertices(:,1) - vertices_apprx(:,1)
    def(:,2) = mres_system%subd_system(level_fine)%mesh%vertices(:,2) - vertices_apprx(:,2)
    if (mres_system%subd_system(level_coarse)%mesh%ndim == 3) then 
        def(:,3) = mres_system%subd_system(level_fine)%mesh%vertices(:,3) - vertices_apprx(:,3)
    end if 

    !Find deltas for the coarse level vertices 
    allocate(delta(mres_system%subd_system(level_coarse)%mesh%nvertex,mres_system%subd_system(level_coarse)%mesh%ndim))
    call srcc_evaluate_coarse_deltas(delta,def,level_coarse,level_fine,mres_system%subd_system,options)
    if (mres_system%subd_system(level_coarse)%mesh%ndim == 3) then 
        deltanorm = sqrt(norm2(delta(:,1))**2 + norm2(delta(:,2))**2)!/real(mres_system%subd_system(level_coarse)%mesh%nvertex,dp)
    elseif (mres_system%subd_system(level_coarse)%mesh%ndim == 3) then 
        deltanorm = sqrt(norm2(delta(:,1))**2 + norm2(delta(:,2))**2 + norm2(delta(:,3))**2)!/real(mres_system%subd_system(level_coarse)%mesh%nvertex,dp)
    end if 
    if (abs(deltanorm) .GE. 1e-40) then 
        deltanorm = log10(abs(deltanorm))
    else 
        deltanorm = abs(deltanorm)
    end if 

    !Update coarse level vertices to their final positions using delta
    mres_system%subd_system(level_coarse)%mesh%vertices(:,:) = mres_system%subd_system(level_coarse)%mesh%vertices(:,:) + delta(:,:)

    !Deallocate locals for this level 
    deallocate(def)
    deallocate(delta)
    deallocate(vertices_apprx)

    !Display
    if (options%console_disp == 'yes') then
        write(*,'(A,I0,A,I0,A,I0,A,I0,A,A,A)') '    < coarsening = ',level_fine,' -> ',level_coarse,' - vertices = ',&
        mres_system%subd_system(level_fine)%mesh%nvertex,' -> ',mres_system%subd_system(level_coarse)%mesh%nvertex,&
        ' - lognorm(delta) = ',real2F0_Xstring(deltanorm,8_in64),' >'
    end if 
end do 
return 
end subroutine propagate_mrsys_sr_catmull_clark




!Evaluate coarse level deltas =========================
subroutine srcc_evaluate_coarse_deltas(delta_coarse,def,level_coarse,level_fine,subd_system,options)
implicit none 

!Weighting parameter omega = [0,1]
!Suggested omega in the range [0.25,0.5]
!omega = 1 --> fully weighted to minimising the subdivision residuals
!omega = 0 --> fully weighted to minimising the coarsened surface energy

!Variables - Import 
integer(in64) :: level_coarse,level_fine
real(dp), dimension(:,:) :: delta_coarse,def
type(subdsystem_data), dimension(:) :: subd_system
type(options_data) :: options 

!Variables - Local 
integer(in64) :: vv,ee
integer(in64) :: n_sharp,sinsert,k_i
integer(in64) :: v1sr,v2sr,v1sc,v2sc,ftgt,etgt,vtgt,fvidx
integer(in64) :: EsharpI(2)
real(dp) :: k,omega,rf,r,s_si2,si,s_ti2,npfr,ti
real(dp) :: deltae(1,subd_system(level_coarse)%mesh%ndim),deltaf(1,subd_system(level_coarse)%mesh%ndim)
real(dp) :: wv(1,subd_system(level_coarse)%mesh%ndim),we(1,subd_system(level_coarse)%mesh%ndim)
real(dp) :: coarse_vtxp(subd_system(level_coarse)%mesh%nvertex,subd_system(level_coarse)%mesh%ndim)

!Extract omega smoothing parameter 
omega = options%srcc_omega

!Extract coarse level approximated vertices 
coarse_vtxp(:,:) = subd_system(level_coarse)%mesh%vertices(:,:)

!Construct deltas
delta_coarse(:,:) = 0.0d0 
do vv=1,subd_system(level_coarse)%mesh%nvertex

    !Extract vertex valence
    k_i = subd_system(level_coarse)%mesh%connectivity%valence(vv)
    k = real(subd_system(level_coarse)%mesh%connectivity%valence(vv),dp)

    !Find number of sharp incident edges to identify vertex type on refined level
    n_sharp = 0
    do ee=1,subd_system(level_fine)%mesh%connectivity%valence(vv)
        etgt = subd_system(level_fine)%mesh%connectivity%v2e(vv,ee)
        if (subd_system(level_fine)%mesh%edge_sharp(etgt) .GT. 0.0d0) then
            n_sharp = n_sharp + 1
        end if
    end do 

    !If 2 sharp edges the locate these edges and adjacent vertices on the coarse and refined levels 
    if (n_sharp == 2) then 

        !On fine level
        EsharpI(:) = 0
        sinsert = 0
        do ee=1,subd_system(level_fine)%mesh%connectivity%valence(vv)
            etgt = subd_system(level_fine)%mesh%connectivity%v2e(vv,ee)
            if (subd_system(level_fine)%mesh%edge_sharp(etgt) .GT. 0.0d0) then
                sinsert = sinsert + 1
                EsharpI(sinsert) = etgt
            end if 
            if (sinsert == 2) then 
                exit 
            end if 
        end do 
        if (subd_system(level_fine)%mesh%edges(EsharpI(1),1) == vv) then 
            v1sr = subd_system(level_fine)%mesh%edges(EsharpI(1),2)
        else
            v1sr = subd_system(level_fine)%mesh%edges(EsharpI(1),1)
        end if 
        if (subd_system(level_fine)%mesh%edges(EsharpI(2),1) == vv) then 
            v2sr = subd_system(level_fine)%mesh%edges(EsharpI(2),2)
        else
            v2sr = subd_system(level_fine)%mesh%edges(EsharpI(2),1)
        end if 
        
        !On coarse level
        EsharpI(:) = 0
        sinsert = 0
        do ee=1,subd_system(level_coarse)%mesh%connectivity%valence(vv)
            etgt = subd_system(level_coarse)%mesh%connectivity%v2e(vv,ee)
            if (subd_system(level_coarse)%mesh%edge_sharp(etgt) .GT. 0.0d0) then
                sinsert = sinsert + 1
                EsharpI(sinsert) = etgt
            end if 
            if (sinsert == 2) then 
                exit 
            end if 
        end do 
        if (subd_system(level_coarse)%mesh%edges(EsharpI(1),1) == vv) then 
            v1sc = subd_system(level_coarse)%mesh%edges(EsharpI(1),2)
        else
            v1sc = subd_system(level_coarse)%mesh%edges(EsharpI(1),1)
        end if 
        if (subd_system(level_coarse)%mesh%edges(EsharpI(2),1) == vv) then 
            v2sc = subd_system(level_coarse)%mesh%edges(EsharpI(2),2)
        else
            v2sc = subd_system(level_coarse)%mesh%edges(EsharpI(2),1)
        end if 
    end if 

    !Valence cases
    if (k_i == 2) then ! -> Valence 2 
        if (subd_system(level_coarse)%mesh%vertex_sharp(vv) == 1) then !Sharp vertex
            delta_coarse(vv,:) = 0.0d0 
        else !Normal vertex

            !Attached vertices on coarse and refined levels 
            v1sc = subd_system(level_coarse)%mesh%connectivity%v2v(vv,1)
            v2sc = subd_system(level_coarse)%mesh%connectivity%v2v(vv,2)
            v1sr = subd_system(level_fine)%mesh%connectivity%v2v(vv,1) 
            v2sr = subd_system(level_fine)%mesh%connectivity%v2v(vv,2) 

            !Construct delta
            delta_coarse(vv,:) = (11.0d0*omega/(64.0d0 - 47.0d0*omega))*(def(v1sr,:) + def(v2sr,:)) + &
            ((32.0d0 - 32.0d0*omega)/(64.0d0 - 47.0d0*omega))*(coarse_vtxp(v1sc,:) - 2.0d0*coarse_vtxp(vv,:) + coarse_vtxp(v2sc,:))
        end if 
    elseif (k_i == 3) then ! -> Valence 3 case
        if ((n_sharp .LT. 2) .AND. (subd_system(level_coarse)%mesh%vertex_sharp(vv) == 0)) then !Normal 

            !Valence 3 deltae
            deltae(1,:) = 0.0d0 
            do ee=1,subd_system(level_coarse)%mesh%connectivity%valence(vv)
                deltae(1,:) = deltae(1,:) + (1.0d0/k)*coarse_vtxp(subd_system(level_coarse)%mesh%connectivity%v2v(vv,ee),:)
            end do 

            !Construct delta
            delta_coarse(vv,:) = deltae(1,:) - coarse_vtxp(vv,:)

        elseif ((n_sharp == 2) .AND. (subd_system(level_coarse)%mesh%vertex_sharp(vv) == 0)) then !On a sharp edge/crease
            delta_coarse(vv,:) = (11.0d0*omega/(64.0d0 - 47.0d0*omega))*(def(v1sr,:) + def(v2sr,:)) + &
            ((32.0d0 - 32.0d0*omega)/(64.0d0 - 47.0d0*omega))*(coarse_vtxp(v1sc,:) - 2.0d0*coarse_vtxp(vv,:) + coarse_vtxp(v2sc,:))
        elseif ((n_sharp .GT. 2) .OR. (subd_system(level_coarse)%mesh%vertex_sharp(vv) == 1)) then !Sharp vertex (implied || specified)  
            delta_coarse(vv,:) = 0.0d0 
        end if 
    else !k > 3 -> Normal case
        if ((n_sharp .LT. 2) .AND. (subd_system(level_coarse)%mesh%vertex_sharp(vv) == 0)) then !Normal 

            !r
            rf = 0.0d0 
            do ee=1,subd_system(level_coarse)%mesh%connectivity%valence(vv)
                ftgt = subd_system(level_coarse)%mesh%connectivity%v2f(vv,ee)
                if (ftgt .GT. 0) then
                    rf = rf + 1.0d0/real(subd_system(level_coarse)%mesh%faces(ftgt)%nvertex,dp)
                end if 
            end do 
            r = ((k - 2.0d0)/k) + (1.0d0/(k*k))*rf

            !deltae
            s_si2 = 0
            deltae(1,:) = 0.0d0 
            do ee=1,subd_system(level_coarse)%mesh%connectivity%valence(vv) 
                
                !Attached edge on the coarse level and edge vertex on the fine level
                etgt = subd_system(level_coarse)%mesh%connectivity%v2e(vv,ee) 
                vtgt = etgt + subd_system(level_coarse)%mesh%nvertex

                !si for this edge
                si = 1.0d0 
                ftgt = subd_system(level_coarse)%mesh%connectivity%e2f(etgt,1)
                if (ftgt .GT. 0) then 
                    npfr = real(subd_system(level_coarse)%mesh%faces(ftgt)%nvertex,dp)
                    si = si + (1.0d0/npfr)
                end if 
                ftgt = subd_system(level_coarse)%mesh%connectivity%e2f(etgt,2)
                if (ftgt .GT. 0) then 
                    npfr = real(subd_system(level_coarse)%mesh%faces(ftgt)%nvertex,dp)
                    si = si + (1.0d0/npfr)
                end if 
                si = 0.25d0*si

                !Add to sum of squares
                s_si2 = s_si2 + si*si

                !Add to deltae
                deltae(1,:) = deltae(1,:) + ((4.0d0*omega/(k*k))*r + omega*si)*def(vtgt,:)
            end do 

            !deltaf
            s_ti2 = 0.0d0
            deltaf(1,:) = 0.0d0 
            do ee=1,subd_system(level_coarse)%mesh%connectivity%valence(vv) 
                ftgt = subd_system(level_coarse)%mesh%connectivity%v2f(vv,ee)
                if (ftgt .GT. 0) then

                    !ti for this face
                    npfr = real(subd_system(level_coarse)%mesh%faces(ftgt)%nvertex,dp)
                    ti = 1.0d0/npfr

                    !Add to sum of squares
                    s_ti2 = s_ti2 + ti*ti
                    
                    !Add to deltaf
                    fvidx = subd_system(level_coarse)%mesh%nvertex + subd_system(level_coarse)%mesh%nedge + ftgt
                    deltaf(1,:) = deltaf(1,:) + ((-omega/(k*k))*r + omega*ti)*def(fvidx,:)
                end if 
            end do 

            !We
            We(1,:) = 0.0d0
            do ee=1,subd_system(level_coarse)%mesh%connectivity%valence(vv) 
                vtgt = subd_system(level_coarse)%mesh%connectivity%v2v(vv,ee)
                We(1,:) = We(1,:) + coarse_vtxp(vtgt,:)
            end do 
            We(1,:) = ((1.0d0 - omega)/k)*We(1,:) 

            !Wv
            Wv(1,:) = (1.0d0 - omega)*coarse_vtxp(vv,:)

            !Construct delta
            delta_coarse(vv,:) = (deltae(1,:) + deltaf(1,:) + We(1,:) - Wv(1,:))/(omega*(r*r + s_si2 + s_ti2) + (1.0d0 - omega))
        elseif ((n_sharp == 2) .AND. (subd_system(level_coarse)%mesh%vertex_sharp(vv) == 0)) then !On a sharp edge/crease
            delta_coarse(vv,:) = (11.0d0*omega/(64.0d0 - 47.0d0*omega))*(def(v1sr,:) + def(v2sr,:)) + &
            ((32.0d0 - 32.0d0*omega)/(64.0d0 - 47.0d0*omega))*(coarse_vtxp(v1sc,:) - 2.0d0*coarse_vtxp(vv,:) + coarse_vtxp(v2sc,:))
        elseif ((n_sharp .GT. 2) .OR. (subd_system(level_coarse)%mesh%vertex_sharp(vv) == 1)) then !Sharp vertex (implied || specified)  
            delta_coarse(vv,:) = 0.0d0 
        end if 
    end if 
end do 
return 
end subroutine srcc_evaluate_coarse_deltas


end module mrsys_smooth_reverse_catmull_clark_mod