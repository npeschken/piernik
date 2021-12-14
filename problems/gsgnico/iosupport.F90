! $Id: iosupport.F90 7934 2013-06-11 07:59:50Z wolt $
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

module iosupport

   implicit none

   private
   public :: galdisk_vars_hdf5, galdisk_tsl, read_initial_fld_from_restart, write_global_to_restart, galdisk_attrs_pre
   public :: init_iosupport, galdisk_redostep
   real   :: emir, emor, emoh !< emag inner radius, outer radius, outer height
   integer :: zflx_layer      !< index of cells layer to count mass flux through z boundaries
   logical :: t_mcmp, t_mcmp_tot, t_emag, t_emag_tot, t_encr, t_encr_tot, t_dmass_stars, t_dmass_stars_tot, &
              t_emag_disk, t_emag_diskonly, t_mflx4, t_mflx_disk, t_massflxz, t_massflxz_tot, t_massflxz_lh, t_massflxz_lh_tot
#ifdef MASS_COMPENS
   real   :: all_mcmpadd
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
   real   :: all_dmass_stars
#endif /* SN_DISTRIBUTION */
#ifdef SNE_DISTR
   real   :: sfrh_ldt = -1.0
   real   :: sfrp_ldt = -1.0
#endif /* SNE_DISTR */
   real, dimension(2) :: tot_massflxz = 0.0
   real               :: massflxz_ldt = 0.0

contains

!>
!! \brief Routine to set parameter values from namelist IO_PARAMS
!!
!! \n \n
!! @b IO_PARAMS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>emir</td><td> 0.0</td><td>real value</td><td>\copydoc iosupport::emir</td></tr>
!! <tr><td>emor</td><td>20.0</td><td>real value</td><td>\copydoc iosupport::emor</td></tr>
!! <tr><td>emoh</td><td> 0.5</td><td>real value</td><td>\copydoc iosupport::emoh</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_iosupport

      use dataio_pub, only: nh ! QA_WARN required for diff_nml
      use mpisetup,   only: ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast
      use units,      only: kpc

      implicit none

      namelist /IO_PARAMS/ emir, emor, emoh, zflx_layer, t_mcmp, t_mcmp_tot, t_emag, t_emag_tot, t_encr, t_encr_tot, t_dmass_stars, t_dmass_stars_tot, &
                           t_emag_disk, t_emag_diskonly, t_mflx4, t_mflx_disk, t_massflxz, t_massflxz_tot, t_massflxz_lh, t_massflxz_lh_tot

      emir =  0.0
      emor = 20.0
      emoh =  0.5
      zflx_layer = 3

      t_mcmp            = .true.
      t_mcmp_tot        = .true.
      t_emag            = .true.
      t_emag_tot        = .true.
      t_encr            = .false.
      t_encr_tot        = .false.
      t_dmass_stars     = .true.
      t_dmass_stars_tot = .false.
      t_emag_disk       = .true.
      t_emag_diskonly   = .true.
      t_mflx4           = .false.
      t_mflx_disk       = .true.
      t_massflxz        = .true.
      t_massflxz_tot    = .true.
      t_massflxz_lh     = .false.
      t_massflxz_lh_tot = .false.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=IO_PARAMS)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=IO_PARAMS, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "IO_PARAMS")
         read(nh%cmdl_nml,nml=IO_PARAMS, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "IO_PARAMS", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=IO_PARAMS)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = emir
         rbuff(2) = emor
         rbuff(3) = emoh

         ibuff(1) = zflx_layer

         lbuff(1)  = t_mcmp
         lbuff(2)  = t_mcmp_tot
         lbuff(5)  = t_emag
         lbuff(6)  = t_emag_tot
         lbuff(7)  = t_encr
         lbuff(8)  = t_encr_tot
         lbuff(13) = t_emag_disk
         lbuff(14) = t_emag_diskonly
         lbuff(15) = t_mflx4
         lbuff(16) = t_mflx_disk
         lbuff(17) = t_massflxz
         lbuff(18) = t_massflxz_tot
         lbuff(19) = t_massflxz_lh
         lbuff(20) = t_massflxz_lh_tot
         lbuff(21) = t_dmass_stars
         lbuff(22) = t_dmass_stars_tot
      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         emir = rbuff(1)
         emor = rbuff(2)
         emoh = rbuff(3)

         zflx_layer        = ibuff(1)

         t_mcmp            = lbuff(1)
         t_mcmp_tot        = lbuff(2)
         t_emag            = lbuff(5)
         t_emag_tot        = lbuff(6)
         t_encr            = lbuff(7)
         t_encr_tot        = lbuff(8)
         t_emag_disk       = lbuff(13)
         t_emag_diskonly   = lbuff(14)
         t_mflx4           = lbuff(15)
         t_mflx_disk       = lbuff(16)
         t_massflxz        = lbuff(17)
         t_massflxz_tot    = lbuff(18)
         t_massflxz_lh     = lbuff(19)
         t_massflxz_lh_tot = lbuff(20)
         t_dmass_stars     = lbuff(21)
         t_dmass_stars_tot = lbuff(22)
      endif

      emir = emir * kpc
      emor = emor * kpc
      emoh = emoh * kpc

   end subroutine init_iosupport

   subroutine galdisk_redostep
#ifdef MASS_COMPENS
      use dodges,   only: sum_mcmpadd, rec_mcmpadd, tot_mcmpadd, ald_mcd
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
      use sndistr,  only: sum_dmass_stars, rec_dmass_stars, tot_dmass_stars, ald_dms
#endif /* SN_DISTRIBUTION */
      implicit none
#ifdef MASS_COMPENS
      if (ald_mcd) tot_mcmpadd = tot_mcmpadd - all_mcmpadd
      sum_mcmpadd = rec_mcmpadd
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
      if (ald_dms) tot_dmass_stars = tot_dmass_stars - all_dmass_stars
      sum_dmass_stars = rec_dmass_stars
#endif /* SN_DISTRIBUTION */
   end subroutine galdisk_redostep

!-----------------------------------------------------------------------------

   subroutine galdisk_vars_hdf5(var, tab, ierrh, cg)

      use grid_cont,        only: grid_container
#if defined ANY_LIMITS || defined SNE_DISTR
      use named_array,      only: p3
      use named_array_list, only: qna
#endif /* ANY_LIMITS || SNE_DISTR */
#ifdef SNE_DISTR
      use constants,        only: xdim, ydim, zdim
      use global,           only: t
      use named_array_list, only: wna
      use sndistr,          only: SFR_n, SFRh_n, DIP_n, DIPh_n
#endif /* SNE_DISTR */
#ifdef MASS_COMPENS
      use dodges,           only: MCD_n
#endif /* MASS_COMPENS */
#ifdef ANY_LIMITS
      use dodges,           only: MSD_n, MMD_n
#endif /* ANY_LIMITS */
#ifdef THERM
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: emag, ekin
      use thermal,          only: thermal_active
      use units,            only: kboltz, mH
#endif /* THERM */

      implicit none

      character(len=*),               intent(in)    :: var
      real, dimension(:,:,:),         intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg
#ifdef THERM
      real, dimension(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) :: eint, kin_ener, mag_ener
      class(component_fluid), pointer                        :: pfl
#endif /* THERM */

      ierrh = 0
      select case (trim(var))
#ifdef SNE_DISTR
         case ("sfrl")
            if (qna%exists(SFR_n)) tab(:,:,:) = real(cg%q(qna%ind(SFR_n) )%span(cg%ijkse), kind=4)
         case ("sfrh")
            if (qna%exists(SFRh_n)) then
               tab(:,:,:) = real(cg%q(qna%ind(SFRh_n))%span(cg%ijkse)/(t-sfrh_ldt), kind=4)
               cg%q(qna%ind(SFRh_n))%arr = 0.0
               sfrh_ldt = t
            endif
         case ("dplx")
            if (wna%exists(DIP_n)) then
               p3 => cg%w(wna%ind(DIP_n) )%span(xdim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dply")
            if (wna%exists(DIP_n)) then
               p3 => cg%w(wna%ind(DIP_n) )%span(ydim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dplz")
            if (wna%exists(DIP_n)) then
               p3 => cg%w(wna%ind(DIP_n) )%span(zdim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dphx")
            if (wna%exists(DIPh_n)) then
               p3 => cg%w(wna%ind(DIPh_n))%span(xdim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dphy")
            if (wna%exists(DIPh_n)) then
               p3 => cg%w(wna%ind(DIPh_n))%span(ydim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
         case ("dphz")
            if (wna%exists(DIPh_n)) then
               p3 => cg%w(wna%ind(DIPh_n))%span(zdim, cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
            endif
#endif /* SNE_DISTR */
#ifdef MASS_COMPENS
         case ("mcmp")
            if (qna%exists(MCD_n)) then
               p3 => cg%q(qna%ind(MCD_n))%span(cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
               cg%q(qna%ind(MCD_n))%arr(:,:,:) = 0.0
            endif
#endif /* MASS_COMPENS */
#ifdef ANY_LIMITS
         case ("mmss")
            if (qna%exists(MMD_n)) then
               p3 => cg%q(qna%ind(MMD_n))%span(cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
               cg%q(qna%ind(MMD_n))%arr(:,:,:) = 0.0
            endif
         case ("dmss")
            if (qna%exists(MSD_n)) then
               p3 => cg%q(qna%ind(MSD_n))%span(cg%ijkse) ; tab(:,:,:) = real(p3, kind=4)
               cg%q(qna%ind(MSD_n))%arr(:,:,:) = 0.0
            endif
#endif /* ANY_LIMITS */
#ifdef THERM
         case ("temp")
            if (thermal_active) then
               !> \warning ONLY ONE FLUID IS USED!!!
               pfl => flind%all_fluids(1)%fl
               if (pfl%has_energy) then
                  kin_ener = ekin(cg%w(wna%fi)%span(pfl%imx,cg%ijkse), cg%w(wna%fi)%span(pfl%imy,cg%ijkse), cg%w(wna%fi)%span(pfl%imz,cg%ijkse), cg%w(wna%fi)%span(pfl%idn,cg%ijkse))
                  if (pfl%is_magnetized) then
                     mag_ener = emag(cg%w(wna%bi)%span(xdim,cg%ijkse), cg%w(wna%bi)%span(ydim,cg%ijkse), cg%w(wna%bi)%span(zdim,cg%ijkse))
                     eint = cg%w(wna%fi)%span(pfl%ien,cg%ijkse) - kin_ener - mag_ener
                  else
                     eint = cg%w(wna%fi)%span(pfl%ien,cg%ijkse) - kin_ener
                  endif
                  tab(:,:,:) = (pfl%gam-1)*mH/kboltz*eint/cg%w(wna%fi)%span(pfl%idn,cg%ijkse)
               endif
            endif
#endif /* THERM */
         case default
            ierrh = -1
            if (.false.) tab(:,:,:) = real(cg%u(-ierrh,:,:,:),kind=4) ! suppress compiler warnings
      end select

   end subroutine galdisk_vars_hdf5

!-----------------------------------------------------------------------------

   subroutine galdisk_tsl(user_vars, tsl_names)

      use constants,   only: LO, HI
      use diagnostics, only: pop_vector
#ifdef MASS_COMPENS
      use dodges,      only: tot_mcmpadd, sum_mcmpadd, ald_mcd
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
      use sndistr,     only: tot_dmass_stars, sum_dmass_stars, ald_dms
#endif /* SN_DISTRIBUTION */

      implicit none

      real,             dimension(:), intent(inout), allocatable           :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names
      real                                                                 :: emag_disk
      real, dimension(2)                                                   :: massflxz
      real, dimension(5)                                                   :: mflx5

      if (present(tsl_names)) then
#ifdef MASS_COMPENS
         if (t_mcmp)            call pop_vector(tsl_names, len(tsl_names(1)), ["mcmpadd    "])                                              !   add to header
         if (t_mcmp_tot)        call pop_vector(tsl_names, len(tsl_names(1)), ["tot_mcmpadd"])                                              !   add to header
#endif /* MASS_COMPENS */
#ifdef SNE_DISTR
#ifdef SN_DISTRIBUTION
         if (t_dmass_stars)     call pop_vector(tsl_names, len(tsl_names(1)), ["dmass_stars"])                                              !   add to header
         if (t_dmass_stars_tot) call pop_vector(tsl_names, len(tsl_names(1)), ["dmass_stars_tot"])                                          !   add to header
#endif /* SN_DISTRIBUTION */
#endif /* SNE_DISTR */
         if (t_emag_disk)       call pop_vector(tsl_names, len(tsl_names(1)), ["emag_disk  "])                                              !   add to header
         if (t_emag_diskonly)   call pop_vector(tsl_names, len(tsl_names(1)), ["em_diskonly"])                                              !   add to header
         if (t_mflx4)           call pop_vector(tsl_names, len(tsl_names(4)), ["mflx1      ", "mflx2      ", "mflx3      ", "mflx4      "]) !   add to header
         if (t_mflx_disk)       call pop_vector(tsl_names, len(tsl_names(1)), ["mflx_disk  "])                                              !   add to header
         if (t_massflxz)        call pop_vector(tsl_names, len(tsl_names(1)), ["mass_flxz"])                                                !   add to header
         if (t_massflxz_tot)    call pop_vector(tsl_names, len(tsl_names(1)), ["mass_flxz_tot"])                                            !   add to header
         if (t_massflxz_lh)     call pop_vector(tsl_names, len(tsl_names(2)), ["mass_flxz_l", "mass_flxz_h"])                               !   add to header
         if (t_massflxz_lh_tot) call pop_vector(tsl_names, len(tsl_names(2)), ["mass_flxz_l_tot", "mass_flxz_h_tot"])                       !   add to header
      else
#ifdef MASS_COMPENS
         call sum_all_tot(sum_mcmpadd, all_mcmpadd, tot_mcmpadd, user_vars, t_mcmp, t_mcmp_tot) ; ald_mcd = .true.
#endif /* MASS_COMPENS */
#ifdef SNE_DISTR
#ifdef SN_DISTRIBUTION
         call sum_all_tot(sum_dmass_stars, all_dmass_stars, tot_dmass_stars, user_vars, t_dmass_stars, t_dmass_stars_tot) ; ald_dms = .true.
#endif /* SN_DISTRIBUTION */
#endif /* SNE_DISTR */
         if (t_emag_disk) then
            call count_emag_disk(emag_disk,.false.)
            call pop_vector(user_vars,[emag_disk])                                                                      !   pop value
         endif
         if (t_emag_diskonly) then
            call count_emag_disk(emag_disk,.true.)
            call pop_vector(user_vars,[emag_disk])                                                                      !   pop value
         endif
         if (t_mflx4 .or. t_mflx_disk) then
            call count_mflx_disk(mflx5)
            if (t_mflx4)     call pop_vector(user_vars,mflx5(1:4))                                                      !   pop value
            if (t_mflx_disk) call pop_vector(user_vars,[mflx5(5)])                                                      !   pop value
         endif
         if (t_massflxz .or. t_massflxz_tot .or. t_massflxz_lh .or. t_massflxz_lh_tot) then
            call count_massflxz(massflxz)
            if (t_massflxz)        call pop_vector(user_vars,[massflxz(HI)-massflxz(LO)])                               !   pop value
            if (t_massflxz_tot)    call pop_vector(user_vars,[tot_massflxz(HI)-tot_massflxz(LO)])                       !   pop value
            if (t_massflxz_lh)     call pop_vector(user_vars,massflxz(:))                                               !   pop value
            if (t_massflxz_lh_tot) call pop_vector(user_vars,tot_massflxz(:))                                           !   pop value
         endif
      endif

   end subroutine galdisk_tsl

   subroutine sum_all_tot(su, al, tot, user_vars, t_i, t_i_tot)
      use constants,   only: pSUM
      use diagnostics, only: pop_vector
      use mpisetup,    only: piernik_MPI_Allreduce
      implicit none
      real,                            intent(inout) :: su, al ,tot
      real, allocatable, dimension(:), intent(inout) :: user_vars
      logical,                         intent(in)    :: t_i, t_i_tot
      if (.not.(t_i .or. t_i_tot)) return
      call piernik_MPI_Allreduce(su, pSUM)
      al = su ; tot = tot + al ; su = 0.0                      !   init for the next one or a couple of steps
      if (t_i)     call pop_vector(user_vars, [al] )           !   pop value
      if (t_i_tot) call pop_vector(user_vars, [tot])           !   pop value
   end subroutine sum_all_tot

   subroutine count_massflxz(flx_on_extbnd)

      use cg_leaves,        only: leaves
      use constants,        only: pSUM, I_ONE, zdim, LO, HI, ndims
      use dataio_pub,       only: die
      use domain,           only: is_multicg
      use fluidindex,       only: iarr_all_mz
      use global,           only: dt, t
      use grid_cont,        only: grid_container
      use mpisetup,         only: piernik_MPI_Allreduce
      use named_array_list, only: wna
      implicit none

      real, dimension(LO:HI),     intent(out) :: flx_on_extbnd
      integer(kind=4), dimension(ndims,LO:HI) :: lh
      integer(kind=4), parameter              :: ifl = 1
      type(grid_container), pointer           :: cg

      if (is_multicg) call die("[iosupport:count_massflxz] multicg not supported")
      flx_on_extbnd = 0.0
      cg => leaves%first%cg
      lh = cg%ijkse

      if (cg%ext_bnd(zdim,LO)) then
      lh(zdim,:) = cg%ijkse(zdim,LO) + zflx_layer - I_ONE
      flx_on_extbnd(LO) = sum(cg%w(wna%fi)%span(iarr_all_mz(ifl),lh)) * cg%dx * cg%dy * dt
      endif

      if (cg%ext_bnd(zdim,HI)) then
      lh(zdim,:) = cg%ijkse(zdim,HI) - zflx_layer + I_ONE
      flx_on_extbnd(HI) = sum(cg%w(wna%fi)%span(iarr_all_mz(ifl),lh)) * cg%dx * cg%dy * dt
      endif

      call piernik_MPI_Allreduce(flx_on_extbnd,pSUM)
      tot_massflxz = tot_massflxz + flx_on_extbnd / dt * (t - massflxz_ldt)
      massflxz_ldt = t

   end subroutine count_massflxz

   subroutine count_emag_disk(ed, diskonly)

      use cg_leaves,        only: leaves
      use constants,        only: LO, HI, xdim, ydim, zdim, pSUM
      use domain,           only: dom
      use func,             only: emag
      use funcgaldisk,      only: aktugalpos
      use grid_cont,        only: grid_container
      use mpisetup,         only: piernik_MPI_Allreduce
      use named_array,      only: p3, p4
      use named_array_list, only: qna, wna

      implicit none

      real,    intent(out)          :: ed
      logical, intent(in)           :: diskonly
      real                          :: y2
      real,    dimension(2)         :: bll, blr
      integer, dimension(2)         :: ind, il, ih, iemir
      integer                       :: i, j, k
      type(grid_container), pointer :: cg

      cg => leaves%first%cg
      p3 => cg%q(qna%wai)%span(cg%ijkse)
      p4 => cg%w(wna%bi)%span(cg%ijkse)
      p3 = emag(p4(xdim,:,:,:),p4(ydim,:,:,:),p4(zdim,:,:,:))
      iemir    = ceiling(emir*cg%idl(xdim:ydim))
      bll(:) = cg%fbnd(xdim:ydim, LO) - dom%nb*cg%dl(xdim:ydim)
      blr(:) = cg%fbnd(xdim:ydim, HI) + dom%nb*cg%dl(xdim:ydim)
      if (all((aktugalpos(xdim:ydim)+emir > bll) .and. (aktugalpos(xdim:ydim)-emir <= blr))) then
         ind =dom%nb+ceiling((aktugalpos(xdim:ydim)-cg%fbnd(xdim:ydim, LO))*cg%idl(xdim:ydim))
         il = max(ind-iemir,int(cg%lhn(xdim:ydim,LO)))
         ih = min(ind+iemir,int(cg%lhn(xdim:ydim,HI)))
         cg%wa(il(xdim):ih(xdim),il(ydim):ih(ydim),cg%ks:cg%ke) = 0.0
      endif
      if (diskonly) then
         do k = cg%ks, cg%ke
            if ((cg%z(k)-aktugalpos(zdim) > emoh) .or. (cg%z(k)-aktugalpos(zdim) < -emoh)) then
               cg%wa(:,:,k) = 0.0
            else
               do j = cg%js, cg%je
                  y2 = (cg%y(j)-aktugalpos(ydim))**2
                  do i = cg%is, cg%ie
                     if (sqrt((cg%x(i)-aktugalpos(xdim))**2+y2) > emor) cg%wa(i,j,k) = 0.0
                  enddo
               enddo
            endif
         enddo
      endif
      ed = sum(p3) * cg%dvol
      call piernik_MPI_Allreduce(ed, pSUM)

   end subroutine count_emag_disk

   subroutine count_mflx_disk(mflx)

      use cg_leaves,        only: leaves
      use constants,        only: LO, HI, xdim, ydim, zdim, ndims, pSUM
      use funcgaldisk,      only: aktugalpos
      use grid_cont,        only: grid_container
      use mpisetup,         only: piernik_MPI_Allreduce
      use named_array_list, only: wna

      implicit none

      real, dimension(5), intent(out)         :: mflx
      integer(kind=4), dimension(ndims,LO:HI) :: ind
      integer(kind=4)                         :: dir, dir2
      type(grid_container), pointer           :: cg

      cg => leaves%first%cg
      mflx(:) = 0.0
      if (.not.((cg%fbnd(zdim,LO) > aktugalpos(zdim)+emoh) .or. (cg%fbnd(zdim,HI) < aktugalpos(zdim)-emoh))) then
         ind(zdim,LO) = max(  floor((aktugalpos(zdim) - emoh - cg%fbnd(zdim,LO))*cg%idl(zdim), kind=4) + cg%lh1(zdim,LO), cg%ks)
         ind(zdim,HI) = min(ceiling((aktugalpos(zdim) + emoh - cg%fbnd(zdim,LO))*cg%idl(zdim), kind=4) + cg%lh1(zdim,LO), cg%ke)
         do dir = xdim, ydim
            dir2 = ndims - dir
            if ((cg%fbnd(dir2,LO) <= aktugalpos(dir2)) .and. (cg%fbnd(dir2,HI) >= aktugalpos(dir2))) then
               ind(dir2,:) = floor((aktugalpos(dir2) - cg%fbnd(dir2,LO))/cg%dl(dir2), kind=4) + cg%lh1(dir2,LO)
               if (.not.((cg%fbnd(dir,LO) > aktugalpos(dir)+emor) .or. (cg%fbnd(dir,HI) < aktugalpos(dir)+emir))) then
                  ind(dir,LO) = max(  floor((aktugalpos(dir) + emir - cg%fbnd(dir,LO))*cg%idl(dir), kind=4) + cg%lh1(dir,LO), cg%ijkse(dir,LO))
                  ind(dir,HI) = min(ceiling((aktugalpos(dir) + emor - cg%fbnd(dir,LO))*cg%idl(dir), kind=4) + cg%lh1(dir,LO), cg%ijkse(dir,HI))
                  mflx(dir) = sum(cg%w(wna%bi)%span(dir2,ind)) * cg%dvol * cg%idl(dir2) * real(3-2*dir)
               endif
               if (.not.((cg%fbnd(dir,LO) > aktugalpos(dir)-emir) .or. (cg%fbnd(dir,HI) < aktugalpos(dir)-emor))) then
                  ind(dir,LO) = max(  floor((aktugalpos(dir) - emor - cg%fbnd(dir,LO))*cg%idl(dir), kind=4) + cg%lh1(dir,LO), cg%ijkse(dir,LO))
                  ind(dir,HI) = min(ceiling((aktugalpos(dir) - emir - cg%fbnd(dir,LO))*cg%idl(dir), kind=4) + cg%lh1(dir,LO), cg%ijkse(dir,HI))
                  mflx(dir+2) = sum(cg%w(wna%bi)%span(dir2,ind)) * cg%dvol * cg%idl(dir2) * real(2*dir-3)
               endif
            endif
         enddo
      endif
      call piernik_MPI_Allreduce(mflx, pSUM)
      mflx(5) = sum(mflx(1:4))

   end subroutine count_mflx_disk

!-----------------------------------------------------------------------------

   subroutine galdisk_attrs_pre

#if defined MC_OR_OVLP || defined SNE_DISTR
      use constants, only: pSUM
      use mpisetup,  only: piernik_MPI_Allreduce
#endif /* MC_OR_OVLP || SNE_DISTR */
#ifdef MASS_COMPENS
      use dodges,    only: sum_mcmpadd
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
      use sndistr,   only: sum_dmass_stars
#endif /* SN_DISTRIBUTION */

      implicit none

#ifdef MASS_COMPENS
      all_mcmpadd = sum_mcmpadd ; call piernik_MPI_Allreduce(all_mcmpadd, pSUM)
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
      all_dmass_stars = sum_dmass_stars ; call piernik_MPI_Allreduce(all_dmass_stars, pSUM)
#endif /* SN_DISTRIBUTION */

   end subroutine galdisk_attrs_pre

!-----------------------------------------------------------------------------

   subroutine write_global_to_restart(file_id)

      use hdf5,        only: HID_T, SIZE_T
      use h5lt,        only: h5ltset_attribute_double_f
#ifdef MASS_COMPENS
      use dodges,      only: init_mass, tot_mcmpadd
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
      use sndistr,     only: tot_dmass_stars
#endif /* SN_DISTRIBUTION */

      implicit none

      integer(HID_T), intent(in)      :: file_id
      integer(SIZE_T)                 :: bufsize
      integer(kind=4)                 :: error

      bufsize = 1
#ifdef MASS_COMPENS
      call h5ltset_attribute_double_f(file_id, "/", "init_mass",     [init_mass],      bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "tot_mcmpadd",   [tot_mcmpadd],    bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "mcmpadd",       [all_mcmpadd],    bufsize, error)
#endif /* MASS_COMPENS */
#ifdef SNE_DISTR
      call h5ltset_attribute_double_f(file_id, "/", "sfrhldt",       [sfrh_ldt],       bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "sfrpldt",       [sfrp_ldt],       bufsize, error)
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
      call h5ltset_attribute_double_f(file_id, "/", "dmass_stars",     [all_dmass_stars], bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "tot_dmass_stars", [tot_dmass_stars], bufsize, error)
#endif /* SN_DISTRIBUTION */
      call h5ltset_attribute_double_f(file_id, "/", "massflxz_ldt",  [massflxz_ldt],   bufsize, error)
      bufsize = 2
      call h5ltset_attribute_double_f(file_id, "/", "tot_massflxz",  tot_massflxz,     bufsize, error)

   end subroutine write_global_to_restart

!-----------------------------------------------------------------------------

   subroutine read_initial_fld_from_restart(file_id)

      use hdf5,     only: HID_T
      use h5lt,     only: h5ltget_attribute_double_f
      use mpisetup, only: master
#ifdef MASS_COMPENS
      use dodges,   only: init_mass, sum_mcmpadd, tot_mcmpadd
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
      use sndistr,  only: tot_dmass_stars, sum_dmass_stars
#endif /* SN_DISTRIBUTION */

      implicit none

      integer(HID_T), intent(in)      :: file_id
      integer(kind=4)                 :: error
      real, dimension(1)              :: buff
      real, dimension(2)              :: buff2

#ifdef MASS_COMPENS
         call h5ltget_attribute_double_f(file_id, "/", "init_mass",     buff,  error)
         init_mass      = buff(1)
         call h5ltget_attribute_double_f(file_id, "/", "tot_mcmpadd",   buff,  error)
         tot_mcmpadd    = buff(1)
#endif /* MASS_COMPENS */
#ifdef SN_DISTRIBUTION
         call h5ltget_attribute_double_f(file_id, "/", "tot_dmass_stars", buff, error)
         tot_dmass_stars = buff(1)
#endif /* SN_DISTRIBUTION */
         call h5ltget_attribute_double_f(file_id, "/", "massflxz_ldt",  buff,  error)
         massflxz_ldt = buff(1)
         call h5ltget_attribute_double_f(file_id, "/", "tot_massflxz",  buff2, error)
         tot_massflxz = buff2(1:2)

         if (master) then
#ifdef MASS_COMPENS
            call h5ltget_attribute_double_f(file_id, "/", "mcmpadd",    buff,  error)
            sum_mcmpadd = buff(1)
#endif /* MASS_COMPENS */
#ifdef SNE_DISTR
            call h5ltget_attribute_double_f(file_id, "/", "sfrhldt",    buff,  error)
            sfrh_ldt = buff(1)
            call h5ltget_attribute_double_f(file_id, "/", "sfrpldt",    buff,  error)
            sfrp_ldt = buff(1)
#endif /* SNE_DISTR */
#ifdef SN_DISTRIBUTION
            call h5ltget_attribute_double_f(file_id, "/", "dmass_stars",  buff,  error)
            sum_dmass_stars = buff(1)
#endif /* SN_DISTRIBUTION */
         endif

   end subroutine read_initial_fld_from_restart

end module iosupport
