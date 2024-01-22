#define pgf
module gravity_waves_sources
  use derivative_mod, only: derivative_t
  use dimensions_mod, only: np,nlev
  use edgetype_mod,   only: EdgeBuffer_t
  use element_mod,    only: element_t
  use hybrid_mod,     only: hybrid_t
  use shr_kind_mod,   only: r8 => shr_kind_r8

  implicit none
  private
  save

  !! gravity_waves_sources created by S Santos, 10 Aug 2011
  !!
  !! gws_src_fnct starts parallel environment and computes frontogenesis
  !!   for use by WACCM (via dp_coupling)

  public  :: gws_src_fnct
  public  :: gws_init
  private :: compute_frontogenesis
#ifdef pgf
  public  :: pgf_init
  public  :: pgf_src
  type (EdgeBuffer_t) :: edge_pgf
#endif
  type (EdgeBuffer_t) :: edge3
  type (derivative_t)   :: deriv
  real(r8) :: psurf_ref

!----------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------

  subroutine gws_init(elem)
    use parallel_mod, only   : par
    use edge_mod, only       : initEdgeBuffer
    use hycoef, only         : hypi
    use pmgrid, only         : plev
    use thread_mod, only     : horz_num_threads
    implicit none

    ! Elem will be needed for future updates to edge code
    type(element_t), pointer :: elem(:)

    ! Set up variables similar to dyn_comp and prim_driver_mod initializations
    call initEdgeBuffer(par, edge3, elem, 3*nlev,nthreads=1)

    psurf_ref = hypi(plev+1)

  end subroutine gws_init

  subroutine gws_src_fnct(elem, tl, tlq, frontgf, frontga,nphys)
    use derivative_mod, only  : derivinit
    use dimensions_mod, only  : npsq, nelemd
    use dof_mod, only         : UniquePoints
    use hybrid_mod, only      : config_thread_region, get_loop_ranges
    use parallel_mod, only    : par
    use ppgrid, only          : pver
    use thread_mod, only      : horz_num_threads
    use dimensions_mod, only  : fv_nphys
    implicit none
    type (element_t), intent(inout), dimension(:) :: elem
    integer, intent(in)          :: tl, nphys, tlq
    real (kind=r8), intent(out) :: frontgf(nphys*nphys,pver,nelemd)
    real (kind=r8), intent(out) :: frontga(nphys*nphys,pver,nelemd)

    ! Local variables
    type (hybrid_t) :: hybrid
    integer :: nets, nete, ithr, ncols, ie
    real(kind=r8), allocatable  ::  frontgf_thr(:,:,:,:)
    real(kind=r8), allocatable  ::  frontga_thr(:,:,:,:)

    ! This does not need to be a thread private data-structure
    call derivinit(deriv)
    !!$OMP PARALLEL NUM_THREADS(horz_num_threads),  DEFAULT(SHARED), PRIVATE(nets,nete,hybrid,ie,ncols,frontgf_thr,frontga_thr)
!    hybrid = config_thread_region(par,'horizontal')
    hybrid = config_thread_region(par,'serial')
    call get_loop_ranges(hybrid,ibeg=nets,iend=nete)

    allocate(frontgf_thr(nphys,nphys,nlev,nets:nete))
    allocate(frontga_thr(nphys,nphys,nlev,nets:nete))    
    call compute_frontogenesis(frontgf_thr,frontga_thr,tl,tlq,elem,deriv,hybrid,nets,nete,nphys)
    if (fv_nphys>0) then
      do ie=nets,nete
        frontgf(:,:,ie) = RESHAPE(frontgf_thr(:,:,:,ie),(/nphys*nphys,nlev/))
        frontga(:,:,ie) = RESHAPE(frontga_thr(:,:,:,ie),(/nphys*nphys,nlev/))
      end do
    else
      do ie=nets,nete
        ncols = elem(ie)%idxP%NumUniquePts
        call UniquePoints(elem(ie)%idxP, nlev, frontgf_thr(:,:,:,ie), frontgf(1:ncols,:,ie))
        call UniquePoints(elem(ie)%idxP, nlev, frontga_thr(:,:,:,ie), frontga(1:ncols,:,ie))
      end do
    end if
    deallocate(frontga_thr)
    deallocate(frontgf_thr)
    !!$OMP END PARALLEL

  end subroutine gws_src_fnct

  subroutine compute_frontogenesis(frontgf,frontga,tl,tlq,elem,ederiv,hybrid,nets,nete,nphys)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute frontogenesis function F
  !   F =  -gradth dot C
  ! with:
  !   theta  = potential temperature
  !   gradth = grad(theta)
  !   C = ( gradth dot grad ) U
  !
  ! Original by Mark Taylor, July 2011
  ! Change by Santos, 10 Aug 2011:
  ! Integrated into gravity_waves_sources module, several arguments made global
  !  to prevent repeated allocation/initialization
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use physconst,      only: cappa
    use air_composition,only: dry_air_species_num, thermodynamic_active_species_num
    use air_composition,only: thermodynamic_active_species_idx_dycore    
    use derivative_mod, only: gradient_sphere, ugradv_sphere
    use edge_mod,       only: edgevpack, edgevunpack
    use bndry_mod,      only: bndry_exchange
    use dyn_grid,       only: hvcoord
    use dimensions_mod, only: fv_nphys,ntrac
    use fvm_mapping,    only: dyn2phys_vector,dyn2phys
    
    type(hybrid_t),     intent(in)            :: hybrid
    type(element_t),    intent(inout), target :: elem(:)
    type(derivative_t), intent(in)            :: ederiv
    integer,            intent(in)            :: nets,nete,nphys
    integer,            intent(in)            :: tl,tlq
    real(r8),           intent(out)           :: frontgf(nphys,nphys,nlev,nets:nete)
    real(r8),           intent(out)           :: frontga(nphys,nphys,nlev,nets:nete)

    ! local
    real(r8) :: area_inv(fv_nphys,fv_nphys), tmp(np,np)
    real(r8) :: uv_tmp(fv_nphys*fv_nphys,2,nlev)
    real(r8) :: frontgf_gll(np,np,nlev,nets:nete)
    real(r8) :: frontga_gll(np,np,nlev,nets:nete)
    integer  :: k,kptr,i,j,ie,component,h,nq,m_cnst
    real(r8) :: gradth(np,np,2,nlev,nets:nete) ! grad(theta)
    real(r8) :: p(np,np)                       ! pressure at mid points
    real(r8) :: pint(np,np)                    ! pressure at interface points
    real(r8) :: theta(np,np)                   ! potential temperature at mid points
    real(r8) :: C(np,np,2), sum_water(np,np)

    do ie=nets,nete
      ! pressure at model top
      pint(:,:) = hvcoord%hyai(1) 
      do k=1,nlev
        ! moist pressure at mid points
        sum_water(:,:) = 1.0_r8
        do nq=dry_air_species_num+1,thermodynamic_active_species_num
          m_cnst = thermodynamic_active_species_idx_dycore(nq)
          !
          ! make sure Q is updated
          !
          sum_water(:,:) = sum_water(:,:) + elem(ie)%state%Qdp(:,:,k,m_cnst,tlq)/elem(ie)%state%dp3d(:,:,k,tl)
        end do
        p(:,:) = pint(:,:) + 0.5_r8*sum_water(:,:)*elem(ie)%state%dp3d(:,:,k,tl)
        ! moist pressure at interface for next iteration
        pint(:,:) = pint(:,:)+elem(ie)%state%dp3d(:,:,k,tl)
        !
        theta(:,:) = elem(ie)%state%T(:,:,k,tl)*(psurf_ref / p(:,:))**cappa
        ! gradth(:,:,:,k,ie) = gradient_sphere(theta,ederiv,elem(ie)%Dinv)        
        call gradient_sphere(theta,ederiv,elem(ie)%Dinv,gradth(:,:,:,k,ie))        
        ! compute C = (grad(theta) dot grad ) u
        C(:,:,:) = ugradv_sphere(gradth(:,:,:,k,ie), elem(ie)%state%v(:,:,:,k,tl),ederiv,elem(ie))        
        ! gradth dot C
        frontgf_gll(:,:,k,ie) = -( C(:,:,1)*gradth(:,:,1,k,ie) +  C(:,:,2)*gradth(:,:,2,k,ie)  )        
        ! apply mass matrix
        gradth(:,:,1,k,ie)=gradth(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
        gradth(:,:,2,k,ie)=gradth(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
        frontgf_gll(:,:,k,ie)=frontgf_gll(:,:,k,ie)*elem(ie)%spheremp(:,:)        
      enddo
      ! pack
      call edgeVpack(edge3, frontgf_gll(:,:,:,ie),nlev,0,ie)
      call edgeVpack(edge3, gradth(:,:,:,:,ie),2*nlev,nlev,ie)
    enddo
    call bndry_exchange(hybrid,edge3,location='compute_frontogenesis')
    do ie=nets,nete
      call edgeVunpack(edge3, frontgf_gll(:,:,:,ie),nlev,0,ie)
      call edgeVunpack(edge3, gradth(:,:,:,:,ie),2*nlev,nlev,ie)
      ! apply inverse mass matrix,
      do k=1,nlev
        gradth(:,:,1,k,ie)=gradth(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
        gradth(:,:,2,k,ie)=gradth(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
        frontgf_gll(:,:,k,ie)=frontgf_gll(:,:,k,ie)*elem(ie)%rspheremp(:,:)        
      end do
      if (fv_nphys>0) then
        uv_tmp(:,:,:) = dyn2phys_vector(gradth(:,:,:,:,ie),elem(ie))
        do k=1,nlev
          h=0
          do j=1,fv_nphys
            do i=1,fv_nphys
              h=h+1
              frontga(i,j,k,ie) = atan2 ( uv_tmp(h,2,k) , uv_tmp(h,1,k) + 1.e-10_r8 )
            end do
          end do
        end do
        !
        ! compute inverse physgrid area for mapping of scaler
        !
        tmp = 1.0_r8
        area_inv = dyn2phys(tmp,elem(ie)%metdet)
        area_inv = 1.0_r8/area_inv
        do k=1,nlev
          frontgf(:,:,k,ie) = dyn2phys(frontgf_gll(:,:,k,ie),elem(ie)%metdet,area_inv)
        end do        
      else
        do k=1,nlev
          frontgf(:,:,k,ie)=frontgf_gll(:,:,k,ie)
          ! Frontogenesis angle
          frontga(:,:,k,ie) = atan2 ( gradth(:,:,2,k,ie) , gradth(:,:,1,k,ie) + 1.e-10_r8 )
        end do
      end if
    enddo
  end subroutine compute_frontogenesis

#ifdef pgf
  subroutine pgf_init(elem)
    use parallel_mod, only   : par
    use edge_mod, only       : initEdgeBuffer
    use hycoef, only         : hypi
    use pmgrid, only         : plev
    use thread_mod, only     : horz_num_threads
    implicit none

    ! Elem will be needed for future updates to edge code
    type(element_t), pointer :: elem(:)

    ! Set up variables similar to dyn_comp and prim_driver_mod initializations
    call initEdgeBuffer(par, edge_pgf, elem, 2*nlev,nthreads=1)
  end subroutine pgf_init

  subroutine pgf_src(elem, tl, tlq, pgf_u, pgf_v,nphys)
    use derivative_mod, only  : derivinit
    use dimensions_mod, only  : npsq, nelemd
    use dof_mod, only         : UniquePoints
    use hybrid_mod, only      : config_thread_region, get_loop_ranges
    use parallel_mod, only    : par
    use ppgrid, only          : pver
    use thread_mod, only      : horz_num_threads
    use dimensions_mod, only  : fv_nphys
    implicit none
    type (element_t), intent(inout), dimension(:) :: elem
    integer, intent(in)          :: tl, nphys, tlq
    real (kind=r8), intent(out) :: pgf_u(nphys*nphys,pver,nelemd)
    real (kind=r8), intent(out) :: pgf_v(nphys*nphys,pver,nelemd)

    ! Local variables
    type (hybrid_t) :: hybrid
    integer :: nets, nete, ithr, ncols, ie
    real(kind=r8), allocatable  ::  pgfu_thr(:,:,:,:)
    real(kind=r8), allocatable  ::  pgfv_thr(:,:,:,:)

    ! This does not need to be a thread private data-structure
    call derivinit(deriv)
    !!$OMP PARALLEL NUM_THREADS(horz_num_threads),  DEFAULT(SHARED), PRIVATE(nets,nete,hybrid,ie,ncols,frontgf_thr,frontga_thr)
!    hybrid = config_thread_region(par,'horizontal')
    hybrid = config_thread_region(par,'serial')
    call get_loop_ranges(hybrid,ibeg=nets,iend=nete)

    allocate(pgfu_thr(nphys,nphys,nlev,nets:nete))
    allocate(pgfv_thr(nphys,nphys,nlev,nets:nete))    
    call compute_pgf(pgfu_thr,pgfv_thr,tl,tlq,elem,deriv,hybrid,nets,nete,nphys)
    if (fv_nphys>0) then
      do ie=nets,nete
        pgf_u(:,:,ie) = RESHAPE(pgfu_thr(:,:,:,ie),(/nphys*nphys,nlev/))
        pgf_v(:,:,ie) = RESHAPE(pgfv_thr(:,:,:,ie),(/nphys*nphys,nlev/))
      end do
    else
      do ie=nets,nete
        ncols = elem(ie)%idxP%NumUniquePts
        call UniquePoints(elem(ie)%idxP, nlev, pgfu_thr(:,:,:,ie), pgf_u(1:ncols,:,ie))
        call UniquePoints(elem(ie)%idxP, nlev, pgfv_thr(:,:,:,ie), pgf_v(1:ncols,:,ie))
      end do
    end if
    deallocate(pgfu_thr)
    deallocate(pgfv_thr)
    !!$OMP END PARALLEL

  end subroutine pgf_src
  
  subroutine compute_pgf(pgf_u,pgf_v,n0,qn0,elem,ederiv,hybrid,nets,nete,nphys)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute frontogenesis function F
  !   F =  -gradth dot C
  ! with:
  !   theta  = potential temperature
  !   gradth = grad(theta)
  !   C = ( gradth dot grad ) U
  !
  ! Original by Mark Taylor, July 2011
  ! Change by Santos, 10 Aug 2011:
  ! Integrated into gravity_waves_sources module, several arguments made global
  !  to prevent repeated allocation/initialization
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use physconst,      only: cappa, tref,cpair,rga,lapse_rate
    use air_composition,only: dry_air_species_num, thermodynamic_active_species_num
    use air_composition,only: thermodynamic_active_species_idx_dycore
    use air_composition, only: get_cp_dry, get_R_dry
    use derivative_mod, only: derivative_t, gradient_sphere
    use edge_mod,       only: edgevpack, edgevunpack
    use bndry_mod,      only: bndry_exchange
    use dyn_grid,       only: hvcoord
    use dimensions_mod, only: fv_nphys,ntrac
    use fvm_mapping,    only: dyn2phys_vector,dyn2phys
    use control_mod,    only: pgf_formulation
    use hybrid_mod,     only: hybrid_t
    use cam_thermo,     only: get_gz, get_virtual_temp
    use cam_thermo,     only: get_kappa_dry
    use cam_abortutils, only: endrun
 
    type(hybrid_t),     intent(in)            :: hybrid
    type(element_t),    intent(inout), target :: elem(:)
    type(derivative_t), intent(in)            :: ederiv
    integer,            intent(in)            :: nets,nete,nphys
    integer,            intent(in)            :: n0,qn0
    real(r8),           intent(out)           :: pgf_u(nphys,nphys,nlev,nets:nete)
    real(r8),           intent(out)           :: pgf_v(nphys,nphys,nlev,nets:nete)

    ! local
    real(r8) :: uv_tmp(fv_nphys*fv_nphys,2,nlev)
    real(r8) :: pgf_gll(np,np,2,nlev,nets:nete)
    integer  :: k,kptr,i,j,ie,h,nq,m_cnst

    real (kind=r8) :: qwater(np,np,nlev,thermodynamic_active_species_num,nets:nete)
    integer        :: qidx(thermodynamic_active_species_num)
    real (kind=r8) :: kappa(np,np,nlev,nets:nete)

    real (kind=r8), dimension(np,np,nlev)                         :: phi
    real (kind=r8), dimension(np,np,2)                            :: vtemp
    real (kind=r8), dimension(np,np,2)                            :: grad_kappa_term
    real (kind=r8), dimension(np,np,2,nlev)                       :: grad_p_full
    real (kind=r8), dimension(np,np,nlev)                         :: dp_dry       ! delta pressure dry
    real (kind=r8), dimension(np,np,nlev)                         :: R_dry, cp_dry!
    real (kind=r8), dimension(np,np,nlev)                         :: p_full       ! pressure
    real (kind=r8), dimension(np,np,nlev)                         :: dp_full
    real (kind=r8), dimension(np,np)                              :: exner
    real (kind=r8), dimension(np,np,nlev)                         :: T_v(np,np,nlev)
    real (kind=r8), dimension(np,np)                              :: suml
    real (kind=r8), dimension(np,np,2) :: pgf_term
    real (kind=r8), dimension(np,np,2) :: grad_exner,grad_logexner
    real (kind=r8) :: T0,T1
    real (kind=r8), dimension(np,np)   :: theta_v
    
    real (kind=r8) :: sum_water(np,np,nlev), density_inv(np,np)
    real (kind=r8) :: ptop
    
    do nq=1,thermodynamic_active_species_num
       qidx(nq) = nq
    end do
    do ie=nets,nete
       do nq=1,thermodynamic_active_species_num
          m_cnst = thermodynamic_active_species_idx_dycore(nq)
          !
          ! make sure Q is updated
          !
          qwater(:,:,:,nq,ie) = elem(ie)%state%Qdp(:,:,:,m_cnst,qn0)/elem(ie)%state%dp3d(:,:,:,n0)
       end do
    end do
    !
    ! compute Cp and kappa=Rdry/cpdry here and not in RK-stages since Q stays constant
    !
    do ie=nets,nete
       call get_kappa_dry(qwater(:,:,:,:,ie), qidx, kappa(:,:,:,ie))
    end do

    ptop = hvcoord%hyai(1)*hvcoord%ps0
    do ie=nets,nete
       !
       ! compute virtual temperature and sum_water
       !
       call get_virtual_temp(qwater(:,:,:,:,ie), t_v(:,:,:),temp=elem(ie)%state%T(:,:,:,n0),&
            sum_q =sum_water(:,:,:), active_species_idx_dycore=qidx)
       call get_R_dry(qwater(:,:,:,:,ie), qidx, R_dry)
       call get_cp_dry(qwater(:,:,:,:,ie), qidx, cp_dry)
       
       do k=1,nlev
          dp_dry(:,:,k)  = elem(ie)%state%dp3d(:,:,k,n0)
          dp_full(:,:,k) = sum_water(:,:,k)*dp_dry(:,:,k)
       end do
       call get_gz(dp_full, T_v, R_dry, elem(ie)%state%phis, ptop, phi, pmid=p_full)
       do k=1,nlev
          call gradient_sphere(p_full(:,:,k),deriv,elem(ie)%Dinv,grad_p_full(:,:,:,k))
          call gradient_sphere(phi(:,:,k),deriv,elem(ie)%Dinv,vtemp)
          density_inv(:,:) = R_dry(:,:,k)*T_v(:,:,k)/p_full(:,:,k)
          
          if (pgf_formulation==1) then
             if (dry_air_species_num==0) then
                exner(:,:)=(p_full(:,:,k)/hvcoord%ps0)**kappa(:,:,k,ie)
                theta_v(:,:)=T_v(:,:,k)/exner(:,:)
                call gradient_sphere(exner(:,:),deriv,elem(ie)%Dinv,grad_exner)
                pgf_term(:,:,1) = cp_dry(:,:,k)*theta_v(:,:)*grad_exner(:,:,1)
                pgf_term(:,:,2) = cp_dry(:,:,k)*theta_v(:,:)*grad_exner(:,:,2)
             else
                exner(:,:)=(p_full(:,:,k)/hvcoord%ps0)**kappa(:,:,k,ie)
                theta_v(:,:)=T_v(:,:,k)/exner(:,:)
                call gradient_sphere(exner(:,:),deriv,elem(ie)%Dinv,grad_exner)
                call gradient_sphere(kappa(:,:,k,ie),deriv,elem(ie)%Dinv,grad_kappa_term)
                suml = exner(:,:)*LOG(p_full(:,:,k)/hvcoord%ps0)
                grad_kappa_term(:,:,1)=-suml*grad_kappa_term(:,:,1)
                grad_kappa_term(:,:,2)=-suml*grad_kappa_term(:,:,2)
                pgf_term(:,:,1) = cp_dry(:,:,k)*theta_v(:,:)*(grad_exner(:,:,1)+grad_kappa_term(:,:,1))
                pgf_term(:,:,2) = cp_dry(:,:,k)*theta_v(:,:)*(grad_exner(:,:,2)+grad_kappa_term(:,:,2))
             end if
             ! balanced ref profile correction:
             ! reference temperature profile (Simmons and Jiabin, 1991, QJRMS, Section 2a)
             !
             !  Tref = T0+T1*Exner
             !  T1 = .0065*Tref*Cp/g ! = ~191
             !  T0 = Tref-T1         ! = ~97
             !
             T1 = lapse_rate*Tref*cpair*rga
             T0 = Tref-T1
             if (hvcoord%hybm(k)>0) then
                !only apply away from constant pressure levels
                call gradient_sphere(log(exner(:,:)),deriv,elem(ie)%Dinv,grad_logexner)
                pgf_term(:,:,1)=pgf_term(:,:,1) + &
                     cpair*T0*(grad_logexner(:,:,1)-grad_exner(:,:,1)/exner(:,:))
                pgf_term(:,:,2)=pgf_term(:,:,2) + &
                     cpair*T0*(grad_logexner(:,:,2)-grad_exner(:,:,2)/exner(:,:))
             end if
          elseif (pgf_formulation==2) then
             pgf_term(:,:,1)  = density_inv(:,:)*grad_p_full(:,:,1,k)
             pgf_term(:,:,2)  = density_inv(:,:)*grad_p_full(:,:,2,k)
          else
             call endrun('ERROR: bad choice of pgf_formulation (must be 1 or 2)')
          end if
          
          do j=1,np
             do i=1,np
                pgf_gll(i,j,1,k,ie) = - vtemp(i,j,1) - pgf_term(i,j,1)
                
                pgf_gll(i,j,2,k,ie)    = - vtemp(i,j,2) - pgf_term(i,j,2)
             end do
          end do
          pgf_gll(:,:,1,k,ie) = pgf_gll(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
          pgf_gll(:,:,2,k,ie) = pgf_gll(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
       end do
       ! pack
       call edgeVpack(edge3, pgf_gll(:,:,:,:,ie),2*nlev,0,ie)
    enddo
    call bndry_exchange(hybrid,edge3,location='compute_frontogenesis')
    do ie=nets,nete
       call edgeVunpack(edge3, pgf_gll(:,:,:,:,ie),2*nlev,0,ie)
       ! apply inverse mass matrix,
       do k=1,nlev
          pgf_gll(:,:,1,k,ie) = pgf_gll(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
          pgf_gll(:,:,2,k,ie) = pgf_gll(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
       end do
              
       if (fv_nphys>0) then
          uv_tmp(:,:,:) = dyn2phys_vector(pgf_gll(:,:,:,:,ie),elem(ie))
          do k=1,nlev
             h=0
             do j=1,fv_nphys
                do i=1,fv_nphys
                   h=h+1
                   pgf_u(i,j,k,ie) = uv_tmp(h,1,k)
                   pgf_v(i,j,k,ie) = uv_tmp(h,2,k)
                end do
             end do
          end do
       else
          do k=1,nlev
             pgf_u(:,:,k,ie) = pgf_gll(:,:,1,k,ie)
             pgf_v(:,:,k,ie) = pgf_gll(:,:,2,k,ie)
          end do
       end if
    enddo
  end subroutine compute_pgf

#endif

end module gravity_waves_sources
