! $Id: initproblem.F90 7959 2014-08-19 09:38:08Z wolt $
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "user_macro.h"

module initproblem

! Initial condition for Galactic disk
! Written by: D. Woltanski, June 2007
! Restrictions: using only cor bnd_cond, disk radius must be smaller than xmax and ymax

  use constants, only: cbuff_len

   implicit none

   private
   public :: problem_initial_conditions, read_problem_par, problem_pointers, problem_initial_nbody

   real                   :: d0, r_max, r_min, rhoa, alpha, beta_cr, stheight, stfact, dmulti
   character(len=cbuff_len) :: bgfile              !< buildgal file name
   real, parameter        :: met = 1.36
   integer(kind=4)        :: mtrout, mtrin
   integer, parameter     :: mfo_len = 32
   character(len=mfo_len) :: mf_orient
   logical                :: oldpiernik, pres_cor

   namelist /PROBLEM_CONTROL/ oldpiernik, mf_orient, rhoa, d0, r_max, r_min, mtrout, mtrin, pres_cor, alpha, beta_cr, stheight, stfact, dmulti

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user,           only: user_vars_hdf5, user_tsl, user_attrs_wr, user_attrs_rd, user_attrs_pre
      use dodges,                only: interpose_after_restart
      use iosupport,             only: galdisk_vars_hdf5, galdisk_tsl, read_initial_fld_from_restart, write_global_to_restart, galdisk_attrs_pre, galdisk_redostep
      use user_hooks,            only: problem_customize_solution, problem_post_restart, user_vars_arr_in_restart, user_reaction_to_redo_step
#ifdef OVLP_ON_BNDS
      use dodges,                only: overlap_on_bnds
      use fluidboundaries_funcs, only: user_fluidbnd
#endif /* OVLP_ON_BNDS */

      implicit none

      user_vars_hdf5             => galdisk_vars_hdf5
      user_tsl                   => galdisk_tsl
      user_attrs_pre             => galdisk_attrs_pre
      user_attrs_wr              => write_global_to_restart
      user_attrs_rd              => read_initial_fld_from_restart
      user_vars_arr_in_restart   => galdisk_init_supplement
      problem_customize_solution => galdisk_problem_customize_solution
      user_reaction_to_redo_step => galdisk_redostep
      problem_post_restart       => interpose_after_restart
#ifdef OVLP_ON_BNDS
      user_fluidbnd              => overlap_on_bnds
#endif /* OVLP_ON_BNDS */

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,   only: cbuff_len
      use dataio_pub,  only: warn, nh, restarted_sim  ! QA_WARN required for diff_nml
      use funcgaldisk, only: g0, v0, aktugalpos, aktugalvel
      use mpisetup,    only: cbuff, ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast
      use units,       only: pc

      implicit none

      g0             = [0.0, 0.0, 0.0] ! initial position of galaxy
      v0             = [0.0, 0.0, 0.0] ! initial velocity of galaxy
      aktugalpos     = [0.0, 0.0, 0.0] ! current position of galaxy
      aktugalvel     = [0.0, 0.0, 0.0] ! current velocity of galaxy

      oldpiernik     = .false.
      pres_cor       = .false.
      rhoa           = 1.0e-4
      d0             = 1.0
      r_max          = 0.8
      r_min          = 0.1
      mtrout         = 10
      mtrin          = 10
      mf_orient      = 'null' ! 'toroidal', 'vertical'
      alpha          = 0.01          ! alpha parameter of Shakura&Sunyaev
      beta_cr        = 0.01
      stheight       = 100.0
      stfact         = 4.e5
      dmulti         = 1.0
      
      bgfile      = 'SPIRAL'

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         lbuff(1) = oldpiernik
         lbuff(2) = pres_cor

         cbuff(2) = bgfile
         cbuff(3) = mf_orient

         rbuff(1) = rhoa
         rbuff(2) = d0
         rbuff(3) = r_max
         rbuff(4) = r_min
         rbuff(5) = alpha
         rbuff(6) = beta_cr
         rbuff(7) = stheight
         rbuff(8) = stfact
         rbuff(9) = dmulti

         ibuff(1) = mtrout
         ibuff(2) = mtrin

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         oldpiernik   = lbuff(1)
         pres_cor     = lbuff(2)

         bgfile      = cbuff(2)
         mf_orient    = cbuff(3)

         rhoa         = rbuff(1)
         d0           = rbuff(2)
         r_max        = rbuff(3)
         r_min        = rbuff(4)
         alpha        = rbuff(5)
         beta_cr      = rbuff(6)
         stheight     = rbuff(7)
         stfact       = rbuff(8)
         dmulti       = rbuff(9)

         mtrout       = ibuff(1)
         mtrin        = ibuff(2)

      endif
      stheight = stheight * pc
#ifndef COSM_RAYS
      beta_cr = 0.0
#endif /* !COSM_RAYS */

      if (master .and. (mf_orient == 'toroidal')) call warn('Toroidal scheme of initial magnetic field has been choosen. This produces really high valu of divB. Do not be surprised ;)')

      call print_essential_units

      call init_user_modules

      ! Perhaps a more general approach should be implemeted once we decide
      ! to drop most of precompiler conditionals on NBODY.
      ! Right now it is used in just one setup.
      ! restarted_sim is set up early in the initpiernik, in
      ! init_dataio_parameters, much earlier than init_dataio
      if (.not. restarted_sim) call problem_initial_nbody

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine print_essential_units

      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master
      use units,      only: units_set, s_len_u, s_time_u, s_mass_u, miu0, km, au, lyr, mH, gmu, erg, eV, clight, Gs, mGs, cm, ppcm2, ppcm3

      implicit none

      if (master) then

         call printinfo('VINE to PIERNIK convertion of units done:')
         write(msg,'(a,a28       )') units_set, ' is used in this simulation.'                         ; call printinfo(msg)
         write(msg,'(a16,3a      )') 'It is based on: ', trim(s_len_u), trim(s_time_u), trim(s_mass_u) ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                       miu0 = ', miu0                            ; call printinfo(msg)
         call printinfo('Some essential units for this set:')
         write(msg,'(a30,e22.15,a)') '                         km = ', km, s_len_u                     ; call printinfo(msg)
         write(msg,'(a30,e22.15,a)') '                         AU = ', au, s_len_u                     ; call printinfo(msg)
         write(msg,'(a30,e22.15,a)') '                        lyr = ', lyr, s_len_u                    ; call printinfo(msg)
         write(msg,'(a30,e22.15,a)') 'hydrogen atom mass:      mH = ', mH, s_mass_u                    ; call printinfo(msg)
         write(msg,'(a30,e22.15,a)') 'galactic mass unit:     gmu = ', gmu, s_mass_u                   ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '    hydrogen atoms / cm2    = ', mH/cm**2                        ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '    hydrogen atoms / cm3    = ', mH/cm**3                        ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') 'averaged particles / cm2    = ', ppcm2                           ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') 'averaged particles / cm3    = ', ppcm3                           ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                        erg = ', erg                             ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                         eV = ', eV                              ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                  erg / cm3 = ', erg/cm**3                       ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                   eV / cm3 = ', eV/cm**3                        ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') 'speed of light in vacuum: c = ', clight                          ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') 'magnetic induction unit: Gs = ', Gs                              ; call printinfo(msg)
         write(msg,'(a30,e22.15  )') '                    microGs = ', mGs                             ; call printinfo(msg)
      endif

   end subroutine print_essential_units

!-----------------------------------------------------------------------------

   subroutine galdisk_init_supplement

#ifdef ANY_LIMITS
      use dodges,  only: register_mass_vars
#endif /* ANY_LIMITS */
#ifdef SNE_DISTR
      use sndistr, only: register_user_var_sndistr
#endif /* SNE_DISTR */

      implicit none

#ifdef SNE_DISTR
      call register_user_var_sndistr
#endif /* SNE_DISTR */
#ifdef ANY_LIMITS
      call register_mass_vars
#endif /* ANY_LIMITS */

   end subroutine galdisk_init_supplement

!-----------------------------------------------------------------------------

   subroutine init_user_modules

      use gravity_user, only: init_galdisk_grav
      use iosupport,    only: init_iosupport
#ifdef ANY_LIMITS
      use dodges,       only: init_dodges
#endif /* ANY_LIMITS */
#ifdef SNE_DISTR
      use sndistr,      only: init_supernovae
#endif /* SNE_DISTR */

      implicit none

      call init_galdisk_grav
#ifdef ANY_LIMITS
      call init_dodges
#endif /* ANY_LIMITS */
#ifdef SNE_DISTR
      call init_supernovae
#endif /* SNE_DISTR */
      call init_iosupport

   end subroutine init_user_modules

!-----------------------------------------------------------------------------
   subroutine set_galdisk_distributions(rc, coldens)

      use units, only: cm, kpc, mp, r_gc_sun

      implicit none

      real, intent(in)  :: rc
      real, intent(out) :: coldens
      real              :: dcmol         !< column density of molecular H_2 medium
      real              :: dcneut        !< column density of cold and warm neutral medium
      real              :: dcion         !< column density of cold and warm ionized medium
      real              :: dchot         !< column density of hot ionized medium
      real              :: rgalcent      !< Galactocentric radius
      real              :: rstrongm      !< parameter radius (btw radius of smoothly increasing magnetic field ~4.2 miuGs)
      real              :: r2relcent     !< expression value with rgalcent
      real              :: r2relstrm     !< expression value with rstrongm

      rgalcent = 4.5*kpc
      rstrongm = 4.0*kpc
      r2relcent = -1.0 * ((rc - rgalcent)**2 - (r_gc_sun - rgalcent)**2)
      r2relstrm = -1.0 * ((rc - rstrongm)**2 - (r_gc_sun - rstrongm)**2)

      ! neutral   [1.e20/cm**2] nearly constant between ~3.5 and 20 kpc (sum of CNM and WNM)
      dcneut = 6.2 &
                   &      /cosh(min((rc/(20.*kpc))**mtrout,100.0)) &
                   & *(-1./cosh(min((rc/(3.5*kpc))**mtrin ,100.0))+1.0)
      ! molecular [1.e20/cm**2]
      dcmol  = 2.6            * exp(r2relcent / ( 2.9*kpc)**2)
      ! ionized   [1.e20/cm**2]
      dcion  = 1.20e-2        * exp(r2relstrm / ( 2.0*kpc)**2) + 1.46 * exp(-( rc**2 - r_gc_sun**2) / (37.0*kpc)**2)
      ! hot ionized    [1.e20/cm**2]
      dchot  = 4.4e-2 * (0.88 * exp(r2relcent / ( 2.9*kpc)**2) + 0.12 * exp(-( rc    - r_gc_sun)    / ( 4.9*kpc)))
      ! sum * mean number of nuclei * unit * proton mass
#ifdef APJL706
      dcneut = 6.2
      coldens = (dcmol+dcneut+dcion+dchot) * met * 0.8
#else /* !APJL706 */
      coldens = (dcmol+dcneut+dcion+dchot) * met * 1.e20/cm**2 * mp
#endif /* !APJL706 */

   end subroutine set_galdisk_distributions

   subroutine tune_disk(rr,dio)

      use global,         only: smalld

      implicit none

      real, intent(in)    :: rr
      real, intent(inout) :: dio

      if (r_max > 0.0) dio = dio     /cosh(min((rr/r_max)**mtrout,100.0))
      if (r_min > 0.0) dio = dio*(-1./cosh(min((rr/r_min)**mtrin ,100.0))+1.0)
      dio = dio + rhoa
      dio = max(dio, smalld)

   end subroutine tune_disk
!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, LO, HI, cwdlen
      use dataio_pub,       only: die
      use domain,           only: dom, is_multicg
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: operator(.notequals.)
      use funcgaldisk,      only: take_care_about_bnd_indx, g0, v0
      use gravity,          only: r_smooth
      use grid_cont,        only: grid_container
      use hydrostatic,      only: hydrostatic_zeq_coldens, set_default_hsparams, dprof
      use gravity_user,     only: grav_hdf5
#ifndef ISO
      use func,             only: ekin, emag
      use global,           only: smallei
#endif /* !ISO */
#ifdef THERM
      use thermal,          only: itemp, thermal_active
      use units,            only: mH, kboltz
#endif /* THERM */
#ifdef MAGNETIC
      use constants,        only: half, pi, MAXL
      use funcgaldisk,      only: update_b_with_A
      use func,             only: sq_sum3
      use global,           only: smalld
      use named_array_list, only: qna
      use types,            only: value
#ifdef SNE_DISTR
      use constants,        only: mag_n
      use sndistr,          only: A_n
#endif /* SNE_DISTR */
#endif /* MAGNETIC */
#ifdef COSM_RAYS
      use constants,        only: small
      use initcosmicrays,   only: iarr_crn, gamma_crn
#endif /* COSM_RAYS */
#ifdef MASS_COMPENS
      use dodges,           only: save_init_mass
#endif /* MASS_COMPENS */
#ifdef INIT_STATE_REFERENCE
      use dodges,           only: ndn0, nom0, nox0, noy0
      use fluidindex,       only: iarr_all_dn
#endif /* INIT_STATE_REFERENCE */
#if defined(INIT_STATE_REFERENCE) || defined(MAGNETIC)
      use named_array_list, only: wna
#endif /* INIT_STATE_REFERENCE || MAGNETIC */
#ifdef VERBOSE
      use dataio_pub,       only: msg, printinfo
      use units,            only: cm, kpc, mp, units_set
#endif /* VERBOSE */

      implicit none

      class(component_fluid), pointer               :: fl
      integer                                       :: i, j, k, iu, id, ju, jd
      real                                          :: xi, yj, zk, rc, irc, iOmega, csim2, densio
      real                                          :: xgradp, ygradp, xgradgp, ygradgp, sfqx, sfqy
      character(len=cwdlen)                         :: filename
      type(grid_container), pointer                 :: cg
#ifdef MAGNETIC
      type(value)                                   :: b_max, dens_max
#endif /* MAGNETIC */
#ifdef VERBOSE
      real                                          :: cdmax, rdmax, realcdmax
      logical                                       :: changerealcdmax

      call printinfo("[initproblem:problem_initial_conditions] Commencing problem module initialization")

      cdmax     = 0.0
      rdmax     = 0.0
      realcdmax = 0.0
      changerealcdmax = .false.
#endif /* VERBOSE */

      fl => flind%ion
      cg => leaves%first%cg
      if (is_multicg) call die("[initproblem:problem_initial_conditions] multiple grid pieces per procesor not implemented yet") !nontrivial

      call set_default_hsparams(cg)

      do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
         yj = cg%y(j) - g0(ydim)
         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            xi = cg%x(i) - g0(xdim)
            rc = sqrt(xi**2+yj**2+r_smooth)
            irc = 0.0
            if (rc .notequals. 0.0) irc = 1./rc

            call set_galdisk_distributions(rc, d0)
#ifdef VERBOSE
            if (cdmax < d0) then
               cdmax = d0
               rdmax = rc
               changerealcdmax = .true.
            endif
#endif /* VERBOSE */
            csim2 = fl%cs2*(1.0 + alpha + beta_cr)
            if (dom%n_t(zdim) == 1) then
               dprof = d0
            else
               call hydrostatic_zeq_coldens(i, j, d0, csim2)
            endif

            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               zk = cg%z(k) - g0(zdim)
               densio = dmulti*dprof(k)
               call tune_disk(rc,densio)
               cg%u(fl%idn,i,j,k) = densio

               call take_care_about_bnd_indx(i, cg%lhn(xdim,:), iu, id, sfqx)
               call take_care_about_bnd_indx(j, cg%lhn(ydim,:), ju, jd, sfqy)
               xgradgp=(cg%gp(iu,j,k)-cg%gp(id,j,k))*sfqx*cg%idl(xdim)
               ygradgp=(cg%gp(i,ju,k)-cg%gp(i,jd,k))*sfqy*cg%idl(ydim)

               if (pres_cor) then
                  xgradp =-sfqx*fl%cs2/fl%gam/cg%u(fl%idn,i,j,k)*(cg%u(fl%idn,iu,j,k)-cg%u(fl%idn,id,j,k))*cg%idl(xdim)
                  ygradp =-sfqy*fl%cs2/fl%gam/cg%u(fl%idn,i,j,k)*(cg%u(fl%idn,i,ju,k)-cg%u(fl%idn,i,jd,k))*cg%idl(ydim)
                  iOmega=sqrt(abs(sqrt((xgradgp+xgradp)**2+(ygradgp+ygradp)**2))*irc)
               else
                  iOmega=sqrt(abs(sqrt(xgradgp**2+ygradgp**2))*irc)
               endif
               !print *, xgradp, ygradp, irc, iOmega
               cg%u(fl%imx,i,j,k) =-iOmega*yj*cg%u(fl%idn,i,j,k)
               cg%u(fl%imy,i,j,k) = iOmega*xi*cg%u(fl%idn,i,j,k)
               cg%u(fl%imz,i,j,k) = 0.0
               cg%u(fl%imx:fl%imz,i,j,k) = cg%u(fl%imx:fl%imz,i,j,k) + v0(:)*cg%u(fl%idn,i,j,k)
#ifdef INIT_STATE_REFERENCE
               cg%w(wna%ind(ndn0))%arr(:,i,j,k) = cg%u(iarr_all_dn,i,j,k)
               cg%w(wna%ind(nox0))%arr(:,i,j,k) =-iOmega*yj
               cg%w(wna%ind(noy0))%arr(:,i,j,k) = iOmega*xi
               cg%w(wna%ind(nom0))%arr(:,i,j,k) = iOmega
#endif /* INIT_STATE_REFERENCE */
#ifndef ISO
               cg%u(fl%ien,i,j,k) = fl%cs2/fl%gam_1*cg%u(fl%idn,i,j,k)
               cg%u(fl%ien,i,j,k) = max(cg%u(fl%ien,i,j,k), smallei)
#ifdef THERM
               if (thermal_active) cg%q(itemp)%arr(i,j,k) = fl%gam_1 * mH / kboltz * cg%u(fl%ien,i,j,k) / cg%u(fl%idn,i,j,k)
#endif /* THERM */
               cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
#endif /* !ISO */

#ifdef COSM_RAYS
               cg%u(iarr_crn(:),i,j,k)   =  beta_cr*flind%ion%cs2 * (cg%u(flind%ion%idn,i,j,k)-max(rhoa,smalld)+small)/(gamma_crn(1)-1.0)
#endif /* COSM_RAYS */

#ifdef MAGNETIC
               select case (mf_orient)
                  case ('null', 'none')
                     cg%b(xdim,i,j,k) = 0.0
                     cg%b(ydim,i,j,k) = 0.0
                     cg%b(zdim,i,j,k) = 0.0
                  case ('vertical')
                     cg%b(xdim,i,j,k) = 0.0
                     cg%b(ydim,i,j,k) = 0.0
                     cg%b(zdim,i,j,k) = sqrt(2.*alpha*d0*fl%cs2)
                  case ('toroidal')       ! this produces really bad value of divB !!! do not be surprised ;)
                     cg%b(xdim,i,j,k) =-sqrt(2.*alpha*fl%cs2*(abs(cg%u(fl%idn,i,j,k)-max(rhoa,smalld)) ))*yj*irc
                     cg%b(ydim,i,j,k) = sqrt(2.*alpha*fl%cs2*(abs(cg%u(fl%idn,i,j,k)-max(rhoa,smalld)) ))*xi*irc
                     cg%b(zdim,i,j,k) = 0.0
                  case ('vectoroidal')
#ifdef SNE_DISTR
                     cg%w(wna%ind(A_n))%arr(xdim,i,j,k) = 0.0
                     cg%w(wna%ind(A_n))%arr(ydim,i,j,k) = 0.0
                     cg%w(wna%ind(A_n))%arr(zdim,i,j,k) = cos(min(rc/r_max*pi,pi))*exp(-((zk/stheight)**2)*half)/stheight/sqrt(2.*pi)
!                     cg%w(wna%ind(A_n))%arr(zdim,i,j,k) = sqrt(2.*alpha*fl%cs2 * abs(cg%u(fl%idn,i,j,k) * cos(min(rc/r_max*pi,pi))))
#else /* !SNE_DISTR */
#error To use vectoroidal you need code compiled with SNE_DISTR
#endif /* !SNE_DISTR */
               end select
#ifndef ISO
               cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
#endif /* MAGNETIC */
               
            enddo
#ifdef VERBOSE
            if (changerealcdmax) then
               realcdmax = sum(cg%u(fl%idn,i,j,:))*cg%dz
               changerealcdmax = .false.
            endif
#endif /* VERBOSE */
         enddo
      enddo
#if defined(MAGNETIC) && defined(SNE_DISTR)
      if (mf_orient == 'vectoroidal') then
         call update_b_with_A(cg, wna%ind(A_n), wna%ind(mag_n), .true.)

         !> is_multicg is false for whole the routine. If switched to multicg then the following should be taken for all cg lists.
         cg%wa(:,:,:) = sqrt(sq_sum3(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:)))
         call leaves%get_extremum(qna%wai, MAXL, b_max)

         cg%wa = cg%u(fl%idn,:,:,:)
         call leaves%get_extremum(qna%wai, MAXL, dens_max)

         cg%b(xdim,:,:,:)  = sqrt(2.*alpha*fl%cs2*abs(dens_max%val)-max(rhoa,smalld))/b_max%val * cg%b(xdim,:,:,:)
         cg%b(ydim,:,:,:)  = sqrt(2.*alpha*fl%cs2*abs(dens_max%val)-max(rhoa,smalld))/b_max%val * cg%b(ydim,:,:,:)
      endif
#endif /* MAGNETIC && SNE_DISTR */

#ifdef VERBOSE
      write(msg,*) 'Max column density is: ',cdmax/met/mp*cm**2,' [mp/cm^2] at radius ',rdmax/kpc,' [kpc]'
      call printinfo(msg)
      write(msg,*) 'real max coldens.  is: ',realcdmax/met/mp*cm**2,' [mp/cm^2]'
      call printinfo(msg)
      write(msg,*) 'real max density   is: ',maxval(cg%u(fl%idn,:,:,:))/met/mp*cm**3,' [mp/cm^3] (', maxval(cg%u(fl%idn,:,:,:)),') [',trim(units_set),']'
      call printinfo(msg)
#endif /* VERBOSE */

#ifdef MASS_COMPENS
      call save_init_mass
#endif /* MASS_COMPENS */

#ifdef VERBOSE
      call printinfo("[initproblem:problem_initial_conditions] finished. \o/")
#endif /* VERBOSE */

      cg%gp = 0.0

   end subroutine problem_initial_conditions

!=============================================================================

   subroutine galdisk_problem_customize_solution(forward)

#ifdef ANY_LIMITS
      use dodges,         only: dump_mass_update
#endif /* ANY_LIMITS */
#ifdef VZ_LIMITS
      use dodges,         only: vzlimiter
#endif /* VZ_LIMITS */
#ifdef MASS_COMPENS
      use dodges,         only: mass_loss_compensate
#endif /* MASS_COMPENS */

#ifdef SNE_DISTR
#ifdef SN_DISTRIBUTION
      use sndistr,        only: distribute_sn_products
#endif /* SN_DISTRIBUTION */
      use star_formation, only: SF
#ifdef DEBUG
      use piernikiodebug, only: force_dumps
#endif /* DEBUG */
#endif /* SNE_DISTR */

      implicit none

      logical, intent(in) :: forward

#ifdef VZ_LIMITS
      call vzlimiter
#endif /* VZ_LIMITS */

#ifdef SNE_DISTR
#ifdef SN_DISTRIBUTION
      call SF(forward)
      !call distribute_sn_products(forward)
#endif /* SN_DISTRIBUTION */

#ifdef DEBUG
      call force_dumps
#endif /* DEBUG */
#endif /* SNE_DISTR */

#ifdef MASS_COMPENS
      call mass_loss_compensate(forward)
#endif /* MASS_COMPENS */
#ifdef ANY_LIMITS
      call dump_mass_update
#endif /* ANY_LIMITS */

      return
      if (.false. .and. forward) return ! suppress compiler warnings on unused arguments

   end subroutine galdisk_problem_customize_solution

   subroutine problem_initial_nbody

      implicit none
      
      call read_buildgal
      
    end subroutine problem_initial_nbody

    
   subroutine read_buildgal

      use constants,      only: ndims
      use dataio_pub,     only: msg, printio, warn
      use particle_utils, only: add_part_in_proper_cg
      use mpisetup,       only: master, proc, nproc
      use star_formation, only: initialize_id, attribute_id

      implicit none

      integer                           :: i, j
      integer(kind=4)                   :: nbodies, pid
      integer, parameter                :: galfile = 1
      real, dimension(:,:), allocatable :: pos, vel
      real, dimension(:),   allocatable :: mass
      real                              :: tbirth

      call initialize_id()
      open(unit=galfile, file=bgfile, action='read', status='old')
      read(galfile,*) nbodies
         if (master) then
            write(msg,'(3a,i8,a)') 'Reading ', trim(bgfile), ' file with ', nbodies, ' particles'
            call printio(msg)
         endif

         allocate(mass(nbodies),pos(nbodies,ndims),vel(nbodies,ndims))

         read(galfile,*) (mass(i),i=1,nbodies), ((pos(i,j),j=1,ndims),i=1,nbodies), ((vel(i,j),j=1,ndims),i=1,nbodies)

      close(galfile)

      i = 0
      do j = 1, nbodies
         i = i + 1
         if (i > nbodies) exit
#ifdef VERBOSE
         if (modulo(i, 10000) .eq. 0) then
            write(msg,'(i8,a)') i, ' particles read' ; call printio(msg)
         endif
#endif /* VERBOSE */
         if (i .le. nbodies / 2) then
            tbirth=-1000.0
         else
            tbirth=-2000.0
         endif
         call attribute_id(pid)
         call add_part_in_proper_cg(pid, mass(i), pos(i,:), vel(i,:),[0.0, 0.0, 0.0], 0.0, tbirth, 0.0)
      enddo
      deallocate(mass,pos,vel)

    end subroutine read_buildgal
    
end module initproblem
