module dimensions_mod
  use shr_kind_mod, only: r8=>shr_kind_r8
#ifdef FVM_TRACERS
  use constituents, only: ntrac_d=>pcnst ! _EXTERNAL
#else
  use constituents, only: qsize_d=>pcnst ! _EXTERNAL
#endif

  implicit none
  private

! set MAX number of tracers.  actual number of tracers is a run time argument  
#ifdef FVM_TRACERS
  integer, parameter         :: qsize_d =10 ! SE tracers (currently SE supports 10 condensate loading tracers)
#else
  integer, parameter         :: ntrac_d = 0 ! No fvm tracers if CSLAM is off
#endif
  !
  ! The variables below hold indices of water vapor and condensate loading tracers as well as
  ! associated heat capacities (initialized in dyn_init):
  !
  !   qsize_condensate_loading_idx     = index of water tracers included in condensate loading according to CAM physics
  !   qsize_condensate_loading_idx_gll = index of water tracers included in condensate loading terms for SE tracers
  !
  ! Note that when running without CSLAM then
  !
  !   qsize_condensate_loading_idx_gll = qsize_condensate_loading_idx
  !
  ! but when running with CSLAM then SE tracers are only the water tracers included in the condensate loading
  !
  character(len=16),  allocatable, public :: cnst_name_gll(:)     ! constituent names for SE tracers
  character(len=128), allocatable, public :: cnst_longname_gll(:) ! long name of SE tracers

  integer, parameter, public :: np = NP
  integer, parameter, public :: nc = 3       !cslam resolution
  integer           , public :: fv_nphys !physics-grid resolution - the "MAX" is so that the code compiles with NC=0

  integer         :: ntrac = 0           !ntrac is set in dyn_comp
  logical, public :: use_cslam = .false. !logical for CSLAM
  integer         :: qsize = 0           !qsize is set in dyn_comp
  integer, public :: wv_idx_dycore = -1  !water vapor index in dycore Qdp; set in dyn_comp
  logical, public :: del4_cslam_qgll = .false. !2026-07-24 (Peter): REPLACED by del2_qgll_diff -- the del4 acted
                                           !on the full q field right after cslam2gll, i.e. it diffused the
                                           !freshly-mapped CSLAM-consistent field itself ("don't apply del2/4
                                           !to cslam itself"). History: del4-on-target probe 2026-07-23 (job
                                           !6871447) worked (stripe wiggle 13x down); was .true. in works1-4.
  logical, public :: del2_qgll_diff = .true. !2026-07-24 (Peter): "del2 on difference between qdp on GLL and
                                           !cslam. don't apply del2 to cslam itself." At each anchor the
                                           !freshly cslam2gll-mapped q(wv) is STORED as reference (qref,
                                           !prim_advance hypervis_Qdp_diff_setref) -- CSLAM continues to rule
                                           !via the anchor overwrite and the mapped field is never damped.
                                           !Between anchors (per physics step, prim_driver nsubstep==1 site)
                                           !del2 damping acts ONLY on the deviation d = q_gll - qref
                                           !(physics increments/vertical remap drift), tendency
                                           !+nu_q_diff*del2(d) weak-form + DSS (hypervis_Qdp_diff).
                                           !.false. = no deviation damping (works4 state has del4 instead).
  real(r8), public :: nu_q_diff = 1.0e5_r8 !del2 coefficient (m^2/s) for del2_qgll_diff deviation damping;
                                           !grid-scale e-folding ~20 min at ne30, stable (dt*nu*k2max ~ 0.1).
  logical, public :: wvdp_nudge_scalesel = .true. !restored after rev6 (2026-07-24): pointwise anchor did NOT
                                           !change the vertex deviation growth (dev is regenerated per
                                           !supercycle by CSLAM Q noise, not accumulated by wvdp) -- rev7
                                           !tests the del4 Q filter instead.
                                           !2026-07-23 (Peter "do 2"): scale-selective wvdp nudging --
                                           !project the nudge increment (dp3d+Qdp(wv)-wvdp) onto the element-
                                           !bilinear (p=1) subspace (GLL-weighted LS, kills ALL intra-element
                                           !grid-scale modes, preserves element water mass) + DSS, THEN apply
                                           !alpha (wvdp_nudge_p1 in prim_advance). Syncs the scales that carry
                                           !the secular PS/PS_gll divergence without slaving wvdp onto the
                                           !static grid-scale texture of the CSLAM->GLL map (stripes).
                                           !.false. = original pointwise nudge (b4b).
  integer, public :: wvdp_nudge_porder = 0 !projection order for wvdp_nudge_scalesel: 1 = element-bilinear
                                           !(p=1 run: tropical texture 13x down but residual alternating
                                           !pattern from the LINEAR moments of the static mapping error);
                                           !0 = element-MEAN only (physically essential water-mass sync;
                                           !mean-error field is smooth across elements => no texture expected).
                                           !2026-07-24: p=1 RETRIED post-cellmean2gll (job 6882740) => still
                                           !UNSTABLE: vtx p50 0.11(d7)->3.96(d10)->13.7(d11), ~doubling/day,
                                           !exponential blow-up (works5 p=0 bounded 0.289 @d15). cellmean2gll
                                           !did NOT tame it -- the p=1 feedback (nudging wvdp's intra-element
                                           !linear moments toward Qdp every anchor) is a genuine instability,
                                           !not just the old static stripe texture. p=0 stays. Do not retry
                                           !p=1 without damping the fed-back linear moments.
  logical, public :: wvdp_nudge_spheremp = .true. !2026-07-24 (Peter "go"): porder=0 element mean taken with
                                           !spheremp weights (true SPHERICAL mass mean) instead of GLL
                                           !reference-space weights => projected increment carries the
                                           !element's exact spherical water-mass discrepancy; DSS already
                                           !conserves the global spheremp integral exactly, so global
                                           !water-mass sync becomes exact at alpha=1 (per-element exact
                                           !pre-DSS; smooth neighbor redistribution post-DSS). Only the
                                           !projection weights change; still constant-per-element =>
                                           !noise-proof property untouched. porder=0 only (p=1 falls back
                                           !to reference weights). .false. = b4b works2 (reference-space
                                           !mean, standing gmean(PS-PSWET) bias -0.005..-0.045 Pa).
  logical, public :: se_del4_dp_damp = .true. !EXPERIMENT 2026-07-09: del4 damping of dp3d in advance_hypervis_dp
                                           !(the nu_p pathway) AND its two-part mass-flux injection into
                                           !sub_elem_mass_flux (interior dpflux + DSS-increment fluxes).
                                           !.false. = dp3d del4 OFF -- tests whether the hypervis-dp flux
                                           !injected into CSLAM swept areas at element/panel edges seeds the
                                           !panel-edge Q noise.  Does NOT touch nu_t/nu/nu_q_cslam, the
                                           !sponge del2 on dp (nu_top, top levels only), or the CSLAM Q filter
                                           !(which reuses nu_p_lev).  Restore .true. after the experiment.
  logical, public :: cslam2gll_hoc = .true.!2026-07-23 (Peter): high-order cell-mean-CONSISTENT CSLAM->GLL map
                                           !(cellmean2gll in fvm_mapping): 4th-order symmetric edge values (C0
                                           !across elements by construction, no limiter, no DSS) + degree-5
                                           !primitive-function interior nodes. Replaces point-value tensor
                                           !Lagrange interp whose static element-periodic artifact the wvdp
                                           !nudging accumulates into element-edge PS noise (tau=1h run 6870585:
                                           !PS lap p50 5x by day 14). .false. = original fvm2dyn path (b4b).
  real(r8), public :: tau_wvdp = 0._r8     !NUDGE-OFF A/B run 6873720 (Peter "remove p=0 stuff"): tau=-1 brought
                                           !back secular divergence (gmean PS-PSWET -0.02 -> -10.1 Pa d15, local
                                           !354 Pa) with noise metrics unchanged => p=0 nudge necessary for sync,
                                           !injects no noise; its te cost ~3.7e-5 (col4 5.1e-5 off vs 8.8e-5 on).
                                           !TAU=0 TEST 2026-07-23 (Peter): hard anchor -- element-mean wvdp
                                           !overwritten to the CSLAM-consistent value every anchor step
                                           !(alpha=1). Safe only WITH wvdp_nudge_scalesel p=0 (projection
                                           !blocks grid-scale injection). Validated 3600 s config = job
                                           !6871925 (works1 + dycore_snapshots/ps_sync_validated_2026-07-23).
                                           !Requires cslam2gll ON so the target dp3d+Qdp(wv) is fresh.
                                           !Free-run (-1) drifted 350 Pa local by day 14.
                                           !SE_WV_CONTINUITY: relaxation timescale (s) for joining SE-evolved
                                           !MOIST dp (wvdp; = dry dp3d + wv mass since 2026-07-22) to the
                                           !CSLAM-consistent target dp3d+Qdp(wv) on GLL after cslam2gll;
                                           !< 0 nudging OFF (wvdp free-runs); = 0 hard overwrite; > 0 relax
  real(r8), public :: cslam_q_filter_mult = 3.0_r8 !multiplier on nu_p_lev(k) for the CSLAM-grid del4 Q filter
                                           !(apply_cslam_q_filter_del4); also used in print_cfl to auto-set
                                           !cslam_q_filter_nsub (control_mod) from the 2D del4 stability bound
  logical, public :: cslam_q_filter_gradmask = .false. !gradient mask on the del4 Q filter: damp only where a
                                           !ABANDONED after run12 2026-07-09 (symmetric two-pass mask): noise
                                           !no longer pumped catastrophically (run11 fixed) but small-scale Q
                                           !noise still GROWS ~20%/day where the detector fires late -- the
                                           !mask concept fails, not just the one-sided implementation.
                                           !grid-scale oscillation is detected (opposing 1D slopes in x or y);
                                           !monotone stretches (resolved fronts) get NO damping at all.
                                           !Mask applied in BOTH del2 passes (operator stays symmetric
                                           !positive semi-definite squared -> provably dissipative; one-sided
                                           !masking pumps grid-scale noise -- run11 2026-07-09). Composable
                                           !with the limiter. Costs 1 extra ghost exchange per subcycle.
  logical, public :: cslam_q_filter_limiter = .false. !Zalesak-type monotone flux limiter on the del4 Q filter:
                                           !limits face fluxes so no cell exits the pre-filter min/max of its
                                           !5-point neighborhood (kills del4 ringing/blockiness at sharp Q
                                           !fronts and guarantees positivity); .false. recovers the unlimited
                                           !filter (run8 behavior). Costs 1 extra ghost exchange per subcycle.
  logical, public :: cslam_q_filter_xdiff = .true. !cross-diffusion (non-orthogonality) correction on the del4
                                           !Q filter face fluxes: adds the tangential-gradient contribution
                                           !that the two-point flux drops on the skewed gnomonic mesh, making
                                           !each del2 pass a consistent spherical Laplacian (largest effect at
                                           !panel edges/corners).  Uses diagonal ring-1 halo cells (already
                                           !exchanged; no extra communication).  .false. recovers the pure
                                           !two-point (TPFA) fluxes bit-for-bit.  Not wired into the abandoned
                                           !gradmask path (aborts if both on).  NOTE: operator is no longer
                                           !provably dissipative (matrix not symmetric); correction is
                                           !O(mesh skew) and the FCT limiter bounds any local overshoot.
  !
  ! fvm dimensions:
  logical, public :: lprint!for debugging
  integer, parameter, public :: ngpc=3          !number of Gausspoints for the fvm integral approximation   !phl change from 4
  integer, parameter, public :: irecons_tracer=6!=1 is PCoM, =3 is PLM, =6 is PPM for tracer reconstruction
  integer,            public :: irecons_tracer_lev(PLEV)
  integer, parameter, public :: nhe=1           !Max. Courant number
  integer, parameter, public :: nhr=2           !halo width needed for reconstruction - phl
  integer, parameter, public :: nht=nhe+nhr     !total halo width where reconstruction is needed (nht<=nc) - phl
  integer, parameter, public :: ns=3!quadratic halo interpolation - recommended setting for nc=3
  !nhc determines width of halo exchanged with neighboring elements
  integer, parameter, public :: nhc = nhr+(nhe-1)+(ns-MOD(ns,2))/2
                                                !(different from halo needed for elements on edges and corners
  integer, parameter, public :: lbc = 1-nhc
  integer, parameter, public :: ubc = nc+nhc
  logical, public            :: large_Courant_incr

  integer, public :: kmin_jet,kmax_jet !min and max level index for the jet
  integer, public :: fvm_supercycling    
  integer, public :: fvm_supercycling_jet

  integer, allocatable, public :: kord_tr(:), kord_tr_cslam(:)
  
  real(r8), public :: nu_scale_top(PLEV)! scaling of del2 viscosity in sopnge layer (initialized in dyn_comp)
  real(r8), public :: nu_lev(PLEV)      ! level dependent del4 (u,v) damping
  real(r8), public :: nu_t_lev(PLEV)    ! level depedendet del4 T damping
  real(r8), public :: nu_p_lev(PLEV)    ! level depedendet del4 dp damping
  real(r8), public :: nu_omega_del2_lev(PLEV) = 0.0_r8 ! del2 omega damping; ramps from 0 to nu_top
  logical,  public, parameter :: del2omega = .false.    ! apply del2 omega sponge layer damping
  integer,  public :: ksponge_end       ! sponge is active k=1,ksponge_end
  real(r8), public :: nu_div_lev(PLEV) = 1.0_r8 ! scaling of viscosity in sponge layer
                                                      ! (set in prim_state; if applicable)
  real(r8), public :: kmvis_ref(PLEV)        !reference profiles for molecular diffusion 
  real(r8), public :: kmcnd_ref(PLEV)        !reference profiles for molecular diffusion  
  real(r8), public :: rho_ref(PLEV)          !reference profiles for rho
  real(r8), public :: km_sponge_factor(PLEV) !scaling for molecular diffusion (when used as sponge)

  integer,  public :: nhc_phys 
  integer,  public :: nhe_phys 
  integer,  public :: nhr_phys 
  integer,  public :: ns_phys  

  integer, public :: npdg = 0  ! dg degree for hybrid cg/dg element  0=disabled 

  integer, parameter, public :: npsq = np*np
  integer, parameter, public :: nlev=PLEV
  integer, parameter, public :: nlevp=nlev+1


!  params for a mesh 
!  integer, public, parameter :: max_elements_attached_to_node = 7
!  integer, public, parameter :: s_nv = 2*max_elements_attached_to_node 

  !default for non-refined mesh (note that these are *not* parameters now)
  integer, public  :: max_elements_attached_to_node = 4
  integer, public  :: s_nv = 6
  integer, public  :: max_corner_elem               = 1 !max_elements_attached_to_node-3
  integer, public  :: max_neigh_edges               = 8 !4 + 4*max_corner_elem

  public :: qsize,qsize_d,ntrac_d,ntrac

  integer, public :: ne
  integer, public :: nelem       ! total number of elements
  integer, public :: nelemd      ! number of elements per MPI task
  integer, public :: nelemdmax   ! max number of elements on any MPI task
  integer, public :: nPhysProc                          ! This is the number of physics processors/ per dynamics processor
  integer, public :: nnodes,npart,nmpi_per_node
  integer, public :: GlobalUniqueCols

  public :: set_mesh_dimensions

contains

  subroutine set_mesh_dimensions()

    ! new "params"
    max_elements_attached_to_node = 7  ! variable resolution
    s_nv = 2*max_elements_attached_to_node 

    !recalculate these
    max_corner_elem               = max_elements_attached_to_node-3
    max_neigh_edges               = 4 + 4*max_corner_elem


  end subroutine set_mesh_dimensions


end module dimensions_mod

