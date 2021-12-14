! $Id: dodges.F90 7959 2014-08-19 09:38:08Z wolt $
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

module dodges
! pulled by ANY_LIMITS
   use constants, only: dsetnamelen
   implicit none
   private
   public :: init_dodges, MSD_n, msd, MMD_n, mmd, register_mass_vars, dump_mass_update, interpose_after_restart
#ifdef VZ_LIMITS
   public :: vzlimiter
#endif /* VZ_LIMITS */
#ifdef OVLP_ON_BNDS
   public :: overlap_on_bnds
#endif /* OVLP_ON_BNDS */
#ifdef INIT_STATE_REFERENCE
   public :: ndn0, nom0, nox0, noy0, fill_user_extant_arrays
#endif /* INIT_STATE_REFERENCE */
#ifdef MASS_COMPENS
   public :: init_mass, save_init_mass, mass_loss_compensate, rec_mcmpadd, sum_mcmpadd, tot_mcmpadd, MCD_n, mcd, ald_mcd
   real                                        :: mass_loss, init_mass, rec_mcmpadd, sum_mcmpadd, tot_mcmpadd
   character(len=dsetnamelen), parameter       :: MCD_n = "MC_D"
   logical                                     :: mcd, ald_mcd
#endif /* MASS_COMPENS */
#ifdef INIT_STATE_REFERENCE
   character(len=dsetnamelen), parameter       :: ndn0 = "den0", nom0 = "omxy", nox0 = "omx0", noy0 = "omy0"
#endif /* INIT_STATE_REFERENCE */
   real                                        :: rd_max, zd_max       !< radius and altitude of MASS_COMPENS action
   real                                        :: mc_reduce            !< density reduction for MASS_COMPENS action (meant to have zero, smalld or rhoa values)
   real                                        :: M_dot_in, zmax, yin, zin, dvx, dvz
   real                                        :: floor_vz
   logical                                     :: overlap_vz, overlap_cr, overlap_dn
   logical                                     :: msd, mmd
   logical                                     :: sngl_mass_incl
   character(len=dsetnamelen), parameter       :: MSD_n = "MS_D", MMD_n = "MM_D"
   contains

!>
!! \brief Routine to set parameter values from namelist FLOORCEILS
!!
!! \n \n
!! @b FLOORCEILS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>floor_vz</td><td>-1.e99  </td><td>real value</td><td>\copydoc dodges::floor_vz</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_dodges

      use dataio_pub, only: nh ! QA_WARN required for diff_nml
      use mpisetup,   only: lbuff, rbuff, master, slave, piernik_MPI_Bcast
      use units,      only: kpc
#ifdef VERBOSE
      use dataio_pub, only: printinfo
#endif /* VERBOSE */

      implicit none

      namelist /FLOORCEILS/ floor_vz, overlap_vz, overlap_cr, overlap_dn, rd_max, zd_max, mc_reduce, M_dot_in, zmax, yin, zin, dvx, dvz, sngl_mass_incl

#ifdef VERBOSE
      call printinfo("[dodges:init_dodges] Commencing dodges module initialization")
#endif /* VERBOSE */


      floor_vz     = -1.e99
      rd_max       =  1.0e5
      zd_max       =  1.0e5
      mc_reduce    =  0.0
      M_dot_in     =  0.0
      zmax         =  0.0
      yin          =  0.0
      zin          =  0.0
      dvx          =  0.0
      dvz          =  0.0

      overlap_vz   = .false.
      overlap_cr   = .false.
      overlap_dn   = .false.
      sngl_mass_incl = .false.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=FLOORCEILS)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=FLOORCEILS, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "FLOORCEILS")
         read(nh%cmdl_nml,nml=FLOORCEILS, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "FLOORCEILS", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=FLOORCEILS)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)   = floor_vz
         rbuff(4)   = rd_max
         rbuff(5)   = zd_max
         rbuff(6)   = mc_reduce
         rbuff(7)   = M_dot_in
         rbuff(8)   = zmax
         rbuff(9)   = yin
         rbuff(10)  = zin
         rbuff(11)  = dvx
         rbuff(12)  = dvz

         lbuff(1)   = overlap_vz
         lbuff(2)   = overlap_cr
         lbuff(3)   = overlap_dn
         lbuff(4)   = sngl_mass_incl

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         floor_vz     = rbuff(1)
         rd_max       = rbuff(4)
         zd_max       = rbuff(5)
         mc_reduce    = rbuff(6)
         M_dot_in     = rbuff(7)
         zmax         = rbuff(8)
         yin          = rbuff(9)
         zin          = rbuff(10)
         dvx          = rbuff(11)
         dvz          = rbuff(12)

         overlap_vz   = lbuff(1)
         overlap_cr   = lbuff(2)
         overlap_dn   = lbuff(3)
         sngl_mass_incl = lbuff(4)

      endif

      rd_max = rd_max * kpc
      zd_max = zd_max * kpc

#ifdef INIT_STATE_REFERENCE
      call register_user_var
#endif /* INIT_STATE_REFERENCE */

#ifdef MASS_COMPENS
      tot_mcmpadd = 0.0 ; sum_mcmpadd = 0.0
#endif /* MASS_COMPENS */

#ifdef VERBOSE
      call printinfo("[dodges:init_dodges] finished. \o/")
#endif /* VERBOSE */

   end subroutine init_dodges

!-----------------------------------------------------------------------------

   subroutine register_mass_vars

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use common_hdf5,      only: hdf_vars
      use constants,        only: AT_IGNORE
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      mmd = any(hdf_vars == 'mmss')
      if (mmd) call all_cg%reg_var(MMD_n, restart_mode = AT_IGNORE)
      msd = any(hdf_vars == 'dmss')
      if (msd) call all_cg%reg_var(MSD_n, restart_mode = AT_IGNORE)
#ifdef MASS_COMPENS
      mcd = any(hdf_vars == 'mcmp')
      if (mcd) call all_cg%reg_var(MCD_n, restart_mode = AT_IGNORE)
#endif /* MASS_COMPENS */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (mmd) cg%q(qna%ind(MMD_n))%arr = 0.0
         if (msd) cg%q(qna%ind(MSD_n))%arr = 0.0
#ifdef MASS_COMPENS
         if (mcd) cg%q(qna%ind(MCD_n))%arr = 0.0
#endif /* MASS_COMPENS */
         cgl => cgl%nxt
      enddo

   end subroutine register_mass_vars

   subroutine dump_mass_update

#ifdef INIT_STATE_REFERENCE
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use fluidindex,       only: iarr_all_dn
      use grid_cont,        only: grid_container
      use named_array_list, only: qna, wna
#endif /* INIT_STATE_REFERENCE */

      implicit none

#ifdef INIT_STATE_REFERENCE
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      if (mmd) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            cg%q(qna%ind(MMD_n))%arr(:,:,:) = sum(cg%u(iarr_all_dn,:,:,:) - cg%w(wna%ind(ndn0))%arr(iarr_all_dn,:,:,:),dim=1)
            cgl => cgl%nxt
         enddo
      endif
#endif /* INIT_STATE_REFERENCE */

   end subroutine dump_mass_update

#ifdef INIT_STATE_REFERENCE
   subroutine register_user_var

      use cg_list_global, only: all_cg
      use constants,      only: AT_OUT_B
      use fluidindex,     only: flind

      implicit none

      call all_cg%reg_var(ndn0, restart_mode = AT_OUT_B, dim4 = flind%fluids)
      call all_cg%reg_var(nom0, restart_mode = AT_OUT_B, dim4 = flind%fluids)
      call all_cg%reg_var(nox0,                          dim4 = flind%fluids)
      call all_cg%reg_var(noy0,                          dim4 = flind%fluids)

   end subroutine register_user_var
#endif /* INIT_STATE_REFERENCE */

!-----------------------------------------------------------------------------

#ifdef INIT_STATE_REFERENCE
   subroutine fill_user_extant_arrays

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, LO, HI
      use grid_cont,        only: grid_container
      use named_array_list, only: wna

      implicit none

      integer                           :: i
      type(cg_list_element),    pointer :: cgl
      type(grid_container),     pointer :: cg
      real, dimension(:,:,:,:), pointer :: pomxy, pomx0, pomy0, pden0

      call leaves%leaf_arr4d_boundaries(wna%ind(ndn0))
      call leaves%leaf_arr4d_boundaries(wna%ind(nom0))

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         pden0 => cg%w(wna%ind(ndn0))%arr
         pomxy => cg%w(wna%ind(nom0))%arr
         pomx0 => cg%w(wna%ind(nox0))%arr
         pomy0 => cg%w(wna%ind(noy0))%arr
         do i = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            pomx0(:,:,i,:) = -cg%y(i)*pomxy(:,:,i,:)
         enddo
         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            pomy0(:,i,:,:) =  cg%x(i)*pomxy(:,i,:,:)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine fill_user_extant_arrays
#endif /* INIT_STATE_REFERENCE */

!-----------------------------------------------------------------------------

   subroutine interpose_after_restart

      implicit none
#ifdef INIT_STATE_REFERENCE
      call fill_user_extant_arrays
#endif /* INIT_STATE_REFERENCE */
#ifdef DISTR_INFLOW
      !call disposable_mass_inclusion
#endif /* DISTR_INFLOW */
   end subroutine interpose_after_restart

!-----------------------------------------------------------------------------

#ifdef VZ_LIMITS
!>
!! \brief Routine to limit infall velocity (that is velocity z-component) in a case of xy-plane disk simulation
!<
   subroutine vzlimiter

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: zdim, LO, HI
      use fluidindex,  only: iarr_all_dn, iarr_all_mz
      use funcgaldisk, only: aktugalpos
      use grid_cont,   only: grid_container

      implicit none

      integer                        :: k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            if (cg%z(k) > aktugalpos(zdim)) then
               cg%u(iarr_all_mz(:),:,:,k) = max(cg%u(iarr_all_mz(:),:,:,k),  floor_vz*cg%u(iarr_all_dn(:),:,:,k))
            elseif (cg%z(k) < aktugalpos(zdim)) then
               cg%u(iarr_all_mz(:),:,:,k) = min(cg%u(iarr_all_mz(:),:,:,k), -floor_vz*cg%u(iarr_all_dn(:),:,:,k))
            endif
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine vzlimiter
#endif /* VZ_LIMITS */

#ifdef OVLP_ON_BNDS
   subroutine overlap_on_bnds(dir, side, cg, wn, qn, emfdir)

      use constants,             only: ndims, xdim, ydim, zdim, LO, HI, I_ONE, I_TWO, I_THREE, wcu_n
      use dataio_pub,            only: warn
      use domain,                only: dom
      use fluidindex,            only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid_cont,             only: grid_container
      use named_array_list,      only: wna, qna
#ifdef COSM_RAYS
      use initcosmicrays,        only: smallecr
      use fluidindex,            only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      integer(kind=4),               intent(in)    :: dir, side
      type(grid_container), pointer, intent(inout) :: cg
      integer(kind=4),     optional, intent(in)    :: wn, qn, emfdir
      integer(kind=4), dimension(ndims,LO:HI)      :: l, r
      real, dimension(:,:,:),   pointer            :: p3
      real, dimension(:,:,:,:), pointer            :: pden0, pomx0, pomy0
      integer(kind=4)                              :: ib, ssign, nbc, nf, sn
      logical                                      :: fluid_bnd, mag_bnd, emf_bnd
#ifndef ZERO_BND_EMF
      real, dimension(:,:,:), allocatable          :: dvar
      real, dimension(:,:,:),   pointer            :: p3a
#endif /* !ZERO_BND_EMF */

      fluid_bnd = .false. ; mag_bnd = .false. ; emf_bnd = .false.

      if (present(wn)) then
         if (wn == wna%fi) fluid_bnd = .true.
         if (wn == wna%bi) mag_bnd   = .true.
         if (.not.fluid_bnd .and. .not.mag_bnd) call warn("[dodges:overlap_on_bnds]: unrecognized case for cg%w boundaries")
      endif

      if (present(qn)) then
         if ((qn == qna%wai) .or. (qn == qna%ind(wcu_n))) emf_bnd = .true.
         if (.not.emf_bnd) call warn("[dodges:overlap_on_bnds]: unrecognized case for cg%q boundaries")
      endif

      if (fluid_bnd) then
         l = cg%lhn ; r = l
         r(dir,:) = cg%ijkse(dir,side)
         ssign = I_TWO*side-I_THREE
         do ib = I_ONE, dom%nb
            l(dir,:) = cg%ijkse(dir,side)+ssign*ib
            cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
         enddo
         l(dir,side) = cg%lhn(dir,side)
         l(dir,ndims-side) = cg%ijkse(dir,side)+ssign
#ifdef COSM_RAYS
         if (overlap_cr .and. (dir /= zdim)) then
            cg%u(iarr_all_crs,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = 0.0
         else
            cg%u(iarr_all_crs,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallecr
         endif
#endif /* COSM_RAYS */
         if (dir == zdim) then
            l(dir,:) = cg%ijkse(dir,side) - [dom%nb, I_ONE] +(dom%nb+I_ONE)*(side-LO)
            if (side == LO) then
               cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = min(cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)),0.0)
            else
               cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = max(cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)),0.0)
            endif
         else
            pden0 => cg%w(wna%ind(ndn0))%span(l)
            pomx0 => cg%w(wna%ind(nox0))%span(l)
            pomy0 => cg%w(wna%ind(noy0))%span(l)
            if (overlap_dn) cg%u(iarr_all_dn,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = pden0
                            cg%u(iarr_all_mx,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = pden0*pomx0
                            cg%u(iarr_all_my,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = pden0*pomy0
            if (overlap_vz) cg%u(iarr_all_mz,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = 0.0
         endif
      endif

      if (mag_bnd) then
         l = cg%lhn ; r = l
         l(dir,:) = cg%lhn(dir,side)
         r(dir,:) = cg%lhn(dir,side) + I_THREE - I_TWO*side
         cg%b(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%b(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
      endif

      if (emf_bnd) then
               l = cg%lhn ; r = l
               sn = I_TWO*side - I_THREE
               nf = cg%ijkse(dir,side) -sn*10
               if (.not.(emfdir == xdim .and. side == HI) .and. .not.(emfdir == ydim .and. side == LO)) nf = nf + sn
               nbc = abs(cg%lhn(dir,side) - nf)
#ifdef ZERO_BND_EMF
               l(dir,I_THREE-side) = nf + sn
               p3 => cg%q(qn)%span(l) ; p3 = 0.0
#else /* !ZERO_BND_EMF */
               l(dir,:) = 1 ; allocate(dvar(l(xdim,LO):l(xdim,HI) ,l(ydim,LO):l(ydim,HI), l(zdim,LO):l(zdim,HI)))
               l = cg%lhn ; r = l
               r(dir,:) = nf ; p3a => cg%q(qn)%span(r)
               l(dir,:) = nf - sn
               dvar(:,:,:) = p3a - cg%q(qn)%span(l)
               do ib = 1, nbc
                  l(dir,:) = nf + sn*ib ; p3 => cg%q(qn)%span(l)
                  p3 = p3a + real(ib)*dvar
               enddo
               deallocate(dvar)
#endif /* ZERO_BND_EMF */
      endif

   end subroutine overlap_on_bnds
#endif /* OVLP_ON_BNDS */

#ifdef DISTR_INFLOW
!>
!! \brief assume mass inflow onto the galaxy
!<
   subroutine distribute_inflow_mass(forward)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, LO, HI
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_my
      use global,           only: smalld, dt
      use grid_cont,        only: grid_container
      use named_array_list, only: wna
#ifndef ISO
      use constants,        only: half
      use fluidindex,       only: flind, iarr_all_en
#endif /* !ISO */

      implicit none

      logical, intent(in)               :: forward
      integer                           :: i,j
      real                              :: dmass, mcmpadd
      real, dimension(:,:,:,:), pointer :: pden0, pomx0, pomy0
      type(cg_list_element),    pointer :: cgl
      type(grid_container),     pointer :: cg


      dmass = M_dot_in/init_mass * dt

      mcmpadd = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         pden0 => cg%w(wna%ind(ndn0))%arr
         pomx0 => cg%w(wna%ind(nox0))%arr
         pomy0 => cg%w(wna%ind(noy0))%arr

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               cg%u(iarr_all_dn,i,j,:) = cg%u(iarr_all_dn,i,j,:) + dmass * pden0(:,i,j,:)
               cg%u(iarr_all_mx,i,j,:) = cg%u(iarr_all_mx,i,j,:) + dmass * pomx0(:,i,j,:)*(pden0(:,i,j,:)-smalld)
               cg%u(iarr_all_my,i,j,:) = cg%u(iarr_all_my,i,j,:) + dmass * pomy0(:,i,j,:)*(pden0(:,i,j,:)-smalld)
#ifndef ISO
               cg%u(iarr_all_en,i,j,:) = cg%u(iarr_all_en,i,j,:) + dmass * pden0(:,i,j,:)*flind%ion%cs2/flind%ion%gam_1
               cg%u(iarr_all_en,i,j,:) = cg%u(iarr_all_en,i,j,:) + dmass * pden0(:,i,j,:)*half*(pomx0(:,i,j,:)**2+pomy0(:,i,j,:)**2)
#endif /* !ISO */
            enddo
         enddo
         mcmpadd = mcmpadd +  dmass * sum(pden0(:,cg%is:cg%ie,cg%js:cg%je,cg%ks:cg%ke)) * cg%dvol
         nullify(pden0,pomx0,pomy0)

         cgl => cgl%nxt
      enddo
      if (forward) rec_mcmpadd = sum_mcmpadd
      sum_mcmpadd = sum_mcmpadd + mcmpadd
      ald_mcd = .false.

   end subroutine  distribute_inflow_mass

!---------------------------------------------------------------------------------------------------------------------------------

!>
!! \brief single mass inclusion after reading restart file
!<
   subroutine disposable_mass_inclusion

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, LO, HI
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_my
      use global,           only: smalld, dt
      use grid_cont,        only: grid_container
      use named_array_list, only: wna
#ifndef ISO
      use constants,        only: half
      use fluidindex,       only: flind, iarr_all_en
#endif /* !ISO */

      implicit none

      integer                           :: i,j
      real                              :: dmass, mcmpadd
      real, dimension(:,:,:,:), pointer :: pden0, pomx0, pomy0
      type(cg_list_element),    pointer :: cgl
      type(grid_container),     pointer :: cg

      if (.not.sngl_mass_incl) return

      ! single mass inclusion should be specified in following lines:
      dmass = M_dot_in/init_mass * dt

      mcmpadd = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         pden0 => cg%w(wna%ind(ndn0))%arr
         pomx0 => cg%w(wna%ind(nox0))%arr
         pomy0 => cg%w(wna%ind(noy0))%arr

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               cg%u(iarr_all_dn,i,j,:) = cg%u(iarr_all_dn,i,j,:) + dmass * pden0(:,i,j,:)
               cg%u(iarr_all_mx,i,j,:) = cg%u(iarr_all_mx,i,j,:) + dmass * pomx0(:,i,j,:)*(pden0(:,i,j,:)-smalld)
               cg%u(iarr_all_my,i,j,:) = cg%u(iarr_all_my,i,j,:) + dmass * pomy0(:,i,j,:)*(pden0(:,i,j,:)-smalld)
#ifndef ISO
               cg%u(iarr_all_en,i,j,:) = cg%u(iarr_all_en,i,j,:) + dmass * pden0(:,i,j,:)*flind%ion%cs2/flind%ion%gam_1
               cg%u(iarr_all_en,i,j,:) = cg%u(iarr_all_en,i,j,:) + dmass * pden0(:,i,j,:)*half*(pomx0(:,i,j,:)**2+pomy0(:,i,j,:)**2)
#endif /* !ISO */
            enddo
         enddo
         mcmpadd = mcmpadd +   sum(pden0(:,cg%is:cg%ie,cg%js:cg%je,cg%ks:cg%ke)) * cg%dvol
         nullify(pden0,pomx0,pomy0)

         cgl => cgl%nxt
      enddo
      rec_mcmpadd = sum_mcmpadd
      sum_mcmpadd = sum_mcmpadd + mcmpadd
      ald_mcd = .false.

   end subroutine  disposable_mass_inclusion

!---------------------------------------------------------------------------------------------------------------------------------

!>
!! \brief Input mass cloud
!<
   subroutine cloud_inflow_mass(forward)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, pSUM
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_mz
      use global,           only: smalld, dt
      use grid_cont,        only: grid_container
      use mpisetup,         only: piernik_MPI_Allreduce
      use named_array_list, only: wna
#ifndef ISO
      use constants,        only: half
      use fluidindex,       only: flind, iarr_all_en
#endif /* !ISO */

      implicit none

      logical, intent(in)               :: forward
      integer                           :: j,k
      real                              :: dmass, mcmpadd, ncells
      type(cg_list_element),    pointer :: cgl
      type(grid_container),     pointer :: cg

      
      
      mcmpadd = 0.0
      ncells  = 0.0
      dmass   = 0.0

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if ( cg%x(cg%ijkse(xdim,HI)) < 48000.) exit ! TO IMPROVE!
         do j = cg%ijkse(ydim,LO), cg%ijkse(ydim,HI)
            if ( cg%y(j) > yin + zmax .or. cg%y(j) < yin - zmax ) cycle 
            do k = cg%ijkse(zdim,LO), cg%ijkse(zdim,HI)
               if ( cg%z(k) > zin + zmax .or. cg%z(k) < zin - zmax ) cycle
               ncells = ncells+1
            end do
         end do
         cgl => cgl%nxt
       end do

       call piernik_MPI_Allreduce(ncells, pSUM)

       cgl => leaves%first
       do while (associated(cgl))
         cg => cgl%cg
         if ( ncells > 0.0 ) dmass = M_dot_in / (cg%dy * cg%dz * abs(dvx) * ncells)
         if ( cg%x(cg%ijkse(xdim,HI)) < 48000.) exit  ! TO IMPROVE!
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            if ( cg%y(j) > yin+zmax .or. cg%y(j) < yin-zmax ) cycle 
            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               if ( cg%z(k) > zin + zmax .or. cg%z(k) < zin - zmax ) cycle

               cg%u(iarr_all_dn,cg%lhn(xdim,HI)-6,j,k) =  dmass
               cg%u(iarr_all_mx,cg%lhn(xdim,HI)-6,j,k) =  dmass * dvx
               cg%u(iarr_all_mz,cg%lhn(xdim,HI)-6,j,k) =  dmass * dvz
#ifndef ISO
               cg%u(iarr_all_en,cg%lhn(xdim,HI)-6,j,k) =  dmass * flind%ion%cs2/flind%ion%gam_1
               cg%u(iarr_all_en,cg%lhn(xdim,HI)-6,j,k) =  dmass * half * (dvx**2 + dvz**2)
#endif /* !ISO */
               mcmpadd = mcmpadd + dmass * cg%dvol
               
            enddo
         enddo

         cgl => cgl%nxt
      enddo
      if (forward) rec_mcmpadd = sum_mcmpadd
      sum_mcmpadd = sum_mcmpadd + mcmpadd
      ald_mcd = .false.

   end subroutine  cloud_inflow_mass
#endif /* DISTR_INFLOW */
!----------------------------------------------------------------------------------------------------------------------------------------
   
#ifdef MASS_COMPENS
   subroutine mass_loss_compensate(forward)

#ifdef DISTR_INFLOW
      use dataio_pub,   only: last_hdf_time
      implicit none
      logical, intent(in) :: forward
      if (last_hdf_time<1000.) call cloud_inflow_mass(forward)
      !call distribute_inflow_mass(forward)
#endif /* DISTR_INFLOW */

   end subroutine mass_loss_compensate

!-----------------------------------------------------------------------------

   subroutine save_init_mass

      implicit none

      mass_loss = 0.0
      init_mass = 0.0

      call total_mass(init_mass)

   end subroutine save_init_mass

!-----------------------------------------------------------------------------

   subroutine total_mass(tmass)

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, pSUM
      use fluidindex,  only: iarr_all_dn
      use funcgaldisk, only: aktugalpos
      use grid_cont,   only: grid_container
      use mpisetup,    only: piernik_MPI_Allreduce

      implicit none

      real, intent(out)              :: tmass
      real                           :: xx, x2, rr
      integer                        :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      tmass = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do i = cg%is, cg%ie
            xx = cg%x(i) - aktugalpos(xdim)
            if (abs(xx) <= rd_max) then
               x2 = xx**2
               do j = cg%js, cg%je
                  rr = sqrt(x2 + (cg%y(j) - aktugalpos(ydim))**2)
                  if (rr <= rd_max) then
                     do k = cg%ks, cg%ke
                        if (abs(cg%z(k) - aktugalpos(zdim)) <= zd_max) tmass = tmass + sum(cg%u(iarr_all_dn, i, j, k)) * cg%dvol
                     enddo
                  endif
               enddo
            endif
         enddo

         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(tmass, pSUM)

   end subroutine total_mass
#endif /* MASS_COMPENS */

end module dodges
