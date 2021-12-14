! $Id: sndistr.F90 7877 2013-06-03 13:09:27Z wolt $
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

!>
!! \brief [DW] Module containing subroutines and functions that govern supernovae distribution
!<
module sndistr
! pulled by SNE_DISTR
   use constants,  only: dsetnamelen
   implicit none
   private
   public :: init_supernovae, register_user_var_sndistr
   public :: A_n, SFR_n, SFRh_n, SFRp_n, DIP_n, DIPh_n, sfr2plt
#ifdef SN_DISTRIBUTION
   public :: distribute_sn_products, ald_dms, rec_dmass_stars, sum_dmass_stars, tot_dmass_stars
#endif /* SN_DISTRIBUTION */
! Supernovae distribution and drawing
! Written by: D. Woltanski, December 2007

   real                                  :: EexplSN, EcrSN, snenerg
   real                                  :: eps_sf, dens_thresh_sf, mass_per_sn
   logical                               :: add_encr
   character(len=dsetnamelen), parameter :: A_n = "av_n", SFR_N = "SFR_n", SFRh_n = "SFRh_n", SFRp_n = "SFRp_n", DIP_n = "DIP_n", DIPh_n = "DIPh_n"
   logical, save                         :: sfr_dump, sfr2plt, sfr2hdf
#ifdef SN_DISTRIBUTION
   real                                  :: rec_dmass_stars, sum_dmass_stars, tot_dmass_stars
   logical                               :: ald_dms
#endif /* SN_DISTRIBUTION */

   contains

!>
!! \brief Routine that sets the initial values of %supernovae parameters from namelist @c SUPERNOVAE_PARAMS
!!
!! \n \n
!! @b SUPERNOVAE_PARAMS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>g_z         </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::g_z          </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_supernovae

      use dataio_pub,     only: die, nh ! QA_WARN required for diff_nml
      use mpisetup,       only: lbuff, rbuff, master, slave, piernik_MPI_Bcast
      use units,          only: erg
#ifdef COSM_RAYS
      use initcosmicrays, only: cr_eff
#endif /* COSM_RAYS */
#ifdef VERBOSE
      use dataio_pub,     only: printinfo
#endif /* VERBOSE */

      implicit none

      namelist /SUPERNOVAE_PARAMS/ add_encr, snenerg, eps_sf, dens_thresh_sf, mass_per_sn

#ifdef VERBOSE
      call printinfo("[sndistr:init_supernovae] Commencing sndistr module initialization")
#endif /* VERBOSE */

      snenerg           = 1.e51        !<  typical energy of supernova explosion [erg]
      add_encr          = .true.       !<  permission for inserting CR energy inside randomly selected areas
      eps_sf            = 0.05         !<  star formation efficiency: fraction of gas mass that falls into stars within one free fall time.
      dens_thresh_sf    = 0.035        !<  threshold density for star formation (0.025 M_sun/pc**3 corresponds to 1H atom/cm^3, 0.035 includes 1 He atomper 10 H atoms)
      mass_per_sn       = 105.         !<  star forming gas mass corresponding to one SN (the quantity depending on IMF)

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=SUPERNOVAE_PARAMS)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=SUPERNOVAE_PARAMS, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "SUPERNOVAE_PARAMS")
         read(nh%cmdl_nml,nml=SUPERNOVAE_PARAMS, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "SUPERNOVAE_PARAMS", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=SUPERNOVAE_PARAMS)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(180) = snenerg
         rbuff(191) = eps_sf
         rbuff(192) = dens_thresh_sf
         rbuff(193) = mass_per_sn

         lbuff(182) = add_encr

      endif

      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         snenerg            = rbuff(180)
         eps_sf             = rbuff(191)
         dens_thresh_sf     = rbuff(192)
         mass_per_sn        = rbuff(193)

         add_encr           = lbuff(182)

      endif

      EexplSN  = snenerg*erg
#ifdef COSM_RAYS
      EcrSN    = EexplSN*cr_eff
#endif /* COSM_RAYS */

#ifdef SN_DISTRIBUTION
      tot_dmass_stars= 0.0 ; sum_dmass_stars= 0.0
#endif /* SN_DISTRIBUTION */

#ifdef VERBOSE
      call printinfo("[sndistr:init_supernovae] finished. \o/")
#endif /* VERBOSE */

   end subroutine init_supernovae

!-----------------------------------------------------------------------------

   subroutine register_user_var_sndistr

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use common_hdf5,      only: hdf_vars
      use constants,        only: AT_NO_B, ndims
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      sfr_dump = add_encr .and. any(hdf_vars == 'sfrl')
      !sfr_dump = .true.
      sfr2plt  = add_encr .and. any(hdf_vars == 'sfrh')
      sfr2hdf  = add_encr .and. any(hdf_vars == 'sfrh')

                    call all_cg%reg_var(A_n,    dim4 = ndims)
      if (sfr_dump) call all_cg%reg_var(SFR_n,                restart_mode = AT_NO_B)
      if (sfr2plt)  call all_cg%reg_var(SFRp_n,               restart_mode = AT_NO_B)
      if (sfr2hdf)  call all_cg%reg_var(SFRh_n,               restart_mode = AT_NO_B)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (sfr_dump) cg%q(qna%ind(SFR_n ))%arr = 0.0
         if (sfr2plt)  cg%q(qna%ind(SFRp_n))%arr = 0.0
         if (sfr2hdf)  cg%q(qna%ind(SFRh_n))%arr = 0.0
         cgl => cgl%nxt
      enddo

   end subroutine register_user_var_sndistr

!-----------------------------------------------------------------------------

#ifdef SN_DISTRIBUTION
!>
!! \brief Routine to add all explosion products directly from SN distribution
!<
   subroutine distribute_sn_products(forward)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
#ifdef COSM_RAYS
      use constants,        only: pi, LO, HI, xdim, ydim, zdim
      use cr_data,          only: icr_H1, cr_table
      use dataio_pub,       only: last_hdf_time
      use fluidindex,       only: flind
      use global,           only: dt, t
      use units,            only: kpc
      use initcosmicrays,   only: iarr_crn
      use named_array,      only: named_array3d
      use named_array_list, only: qna, wna
      use units,            only: newtong
#endif /* COSM_RAYS */

      implicit none

      logical, intent(in)             :: forward
      real, dimension(:,:,:), pointer :: psfr, pdens, pcrH1
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer  :: i,j,k
#ifdef COSM_RAYS
      type(named_array3d)             :: sfr
      real                            :: dmass_stars, c_tau_ff !, decr_sn
#ifdef INSERTSNONCEASTEP
      real, parameter                 :: dtf = 2.0 !< factor depending on number of calls during double (2*dt) step
#else /* !INSERTSNONCEASTEP */
      real, parameter                 :: dtf = 1.0 !< factor depending on number of calls during double (2*dt) step
#endif /* !INSERTSNONCEASTEP */
#endif /* COSM_RAYS */

#ifdef INSERTSNONCEASTEP
      if (.not.forward) return
#endif /* INSERTSNONCEASTEP */

      if (forward) rec_dmass_stars = sum_dmass_stars
      ald_dms = .false.

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

#ifdef COSM_RAYS
         if (add_encr) then

            dmass_stars = 0.0
            call sfr%init (cg%lhn(:,LO), cg%lhn(:,HI))
            psfr  => sfr%span(cg%ijkse)
            pdens => cg%w(wna%fi)%span(flind%ion%idn,cg%ijkse)
            pcrH1 => cg%w(wna%fi)%span(iarr_crn(cr_table(icr_H1)),cg%ijkse)
!>
!! Number of SN exploding in every cell in a double timestep,
!! resulting from the assumption that gas collapses to stars, within a time period
!! which is related to free-fall time, with a relatively small (~5%) efficiency
!! as soon as the system becomes gravitationally unstable.
!! d rho_stars/dt = sfe*rho_gas/tau_ff, where sfe is the star formation efficiency parameter
!! and tau_ff= sqrt(3 \pi/(32 G \rho_gas)) is the free-fall time.
!! this leads to the  relation : d M_stars/dt \propto (rho_gas)^1.5
!! (see Mac Low & Klessen: eqn. (1) in REVIEWS OF MODERN PHYSICS, VOLUME 76, 2004, p.126)
!<
!            tau_ff = c_tau_ff/sqrt(dens), where
!            dens = mH*(14./10.)/cm**3, assuming 1 atom He per 10 H atoms

            c_tau_ff = sqrt(3.*pi/(32.*newtong))

!            M_cell = dens * dvol
!            M_stars = eps_sf * M_cell * twodt / tau_ff
!            dnsn = M_stars/(mass_per_sn*M_sun)      ! mass_per_sn = 105 M_sun of gas forming stars resulting in 1 SN

! CR energy density increment from one SN in cell
!            decr_sn = EcrSN/cg%dvol

! CR energy density increment from all SN in cell
            sfr%arr(:,:,:)  = 0.0
            where (pdens > dens_thresh_sf)
               psfr  = eps_sf / c_tau_ff * pdens**(3./2.)
            endwhere
            do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
               do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
                  if ((cg%x(i)**2 + cg%y(j)**2 > 625 * kpc * kpc) .and. (t < 1500.0)) then             
                     sfr%arr(i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI)) = 0.0
                     !if ((cg%x(i) .gt. 31.0 * kpc ) .and. (cg%y(j) .gt. 28.0 * kpc) .and. (sum(cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI))) .gt. 0.00000000000000000000001)) then
                     !        print *, '!!!!!!!!!!!!!!!!!!! COSMIC RAYS IN STREAM', sum(cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI))), '!!!!!!!!!!!!'
                     !endif
                     !cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)),i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI)) = 0.0
                     !cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)), i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI)) = 0.0
                     !psfr(i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI)) = 0.0
                  else 
                     cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)), i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI)) = cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)), i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI)) + EcrSN * sfr%arr(i,j,cg%lhn(zdim,LO):cg%lhn(zdim,HI))/mass_per_sn * dtf * dt
                  endif
                  !do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
                  !      if ((cg%x(i)> 30000.0) .and. (cg%y(i)>28000.0) .and. (cg%z(k)>10000.0) .and. (t>290) .and. (t<320)) then
                  !              cg%w(wna%fi)%arr(iarr_crn(cr_table(icr_H1)), i,j,k) = 0.0
                  !      endif
                  !enddo
               end do
            end do
            !do i = cg%lhn(zdim, LO), cg%lhn(zdim, HI) 
            !    if (abs(cg%z(i)) >5) then
            !         sfr%arr(cg%lhn(xdim,LO):cg%lhn(xdim,HI), cg%lhn(ydim,LO):cg%lhn(ydim,HI),i)  = 0.0
            !    endif
            !enddo
            
            pdens = pdens - psfr * dtf * dt
            !pcrH1 = pcrH1 + EcrSN * psfr/mass_per_sn * dtf * dt
            dmass_stars = sum(psfr) * cg%dvol * dtf * dt

! We dump the SFR to .h5 files and dmass_stars (per double timestep) to .tsl files

            if (sfr_dump) cg%q(qna%ind(SFR_n ))%arr(:,:,:) = sfr%arr(:,:,:)

            if (sfr2plt)  cg%q(qna%ind(SFRp_n))%arr(:,:,:) = cg%q(qna%ind(SFRp_n))%arr(:,:,:) + sfr%arr(:,:,:) * cg%dvol * dtf * dt
            if (sfr2hdf)  cg%q(qna%ind(SFRh_n))%arr(:,:,:) = cg%q(qna%ind(SFRh_n))%arr(:,:,:) + sfr%arr(:,:,:) * cg%dvol * dtf * dt

            call sfr%clean
         endif

         sum_dmass_stars = sum_dmass_stars + dmass_stars
#endif /* COSM_RAYS */

         cgl => cgl%nxt
      enddo

   end subroutine  distribute_sn_products
#endif /* SN_DISTRIBUTION */

end module sndistr
