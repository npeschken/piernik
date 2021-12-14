! $Id: gravity_user.F90 7040 2012-11-13 13:28:03Z wolt $
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

module gravity_user

   implicit none
   private
   public :: init_galdisk_grav, grav_hdf5

   real, allocatable      :: gpotdisk(:,:,:), gpothalo(:,:,:), gpotbulge(:,:,:)
   real                   :: Mhalo, Mbulge, Mdisk, ahalo, bbulge, adisk, bdisk

   contains

   subroutine galdisk_grav_settings

      use units, only: gmu, kpc

      implicit none

      Mhalo = 4615.0*gmu               ! 24.3e10*Msun    !Milky Way
      Mbulge= 0.0          !606.0*gmu  !  0.8e10*Msun    !Milky Way (estimation)
      Mdisk = 3690.0*gmu               !  3.7e10*Msun    !Milky Way
      ahalo =   12.0*kpc               !    35. *kpc
      bbulge=    3.0*kpc   !0.3873*kpc !  2100. * pc
      adisk = 5.3178*kpc               !     4.9*kpc
      bdisk = 0.2500*kpc               !   150. * pc

   end subroutine galdisk_grav_settings

   subroutine init_galdisk_grav

      use dataio_pub,  only: die, warn
      use gravity,     only: get_gprofs, gprofs_target, grav_pot_3d, grav_type, user_grav, variable_gp

      implicit none

      !grav_type => grav_allen_santillan
      grav_type => grav_hdf5_interp
      !call galdisk_grav_settings
      grav_pot_3d => grav_pot_user

      !if (.not.user_grav) call die("[gravity_user:init_galdisk_grav] user_grav has been found set .false.! Gravity is not present at all!")
      if (.not.associated(get_gprofs) .and. (gprofs_target /= 'extgp')) call die("[gravity_user:init_galdisk_grav] get_gprofs and gprofs_target configuration is not proper for galdisk!")
      if (variable_gp) then
         call warn("[vtp_convert:switch_to_vgpot]: variable_gp has been found set .true. and time could be wasted on counting it every step! Switch to .false. value.")
         variable_gp = .false.
      endif

   end subroutine init_galdisk_grav

   subroutine grav_pot_user

      use axes_M,    only: axes
      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use grid_cont, only: grid_container
      use constants, only: cwdlen

      implicit none

      type(axes)                     :: ax
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      character(len=cwdlen)          :: filename

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         !call ax%allocate_axes(cg%lhn)
         !ax%x = cg%x
         !ax%y = cg%y
         !ax%z = cg%z
         !call grav_allen_santillan(cg%gp, ax, cg%lhn)
         filename = 'previous.h5'
         call grav_hdf5(cg, filename)
         !call ax%deallocate_axes

         cgl => cgl%nxt
      enddo

   end subroutine grav_pot_user

   subroutine grav_allen_santillan(gp, ax, lhn, flatten)

      use axes_M,      only: axes
      use constants,   only: ndims, xdim, ydim, LO, HI
      use dataio_pub,  only: die
      use diagnostics, only: ma1d, my_allocate, my_deallocate
      use func,        only: operator(.notequals.)
      use units,       only: kpc, newtong

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten

      logical                                             :: grav_log
      integer                                             :: i,j
      real                                                :: rgcs2, mahalo, khalo
      real, allocatable, dimension(:)                     :: gpdisk, gphalo, gpblg, gpot, rgsv, rgsv2, irgsv, frgsv

      if (present(flatten)) grav_log = flatten ! stop compilator warnings
      grav_log = .false.

! galactic case as in vollmer'01 (gravitational potential of Allen & Santillan '01)
      ma1d = [int(size(ax%z), kind=4)]
      call my_allocate(gpdisk, ma1d)
      call my_allocate(gphalo, ma1d)
      call my_allocate(gpblg,  ma1d)
      call my_allocate(gpot,   ma1d)
      call my_allocate(rgsv,   ma1d)
      call my_allocate(rgsv2,  ma1d)
      call my_allocate(irgsv,  ma1d)
      call my_allocate(frgsv,  ma1d)
      mahalo = Mhalo/1.02/ahalo
      khalo  = 1.0 + (100.0*kpc/ahalo)**1.02

      if (.not. associated(gp)) call die("[gravity_user:grav_allen_santillan] gp is not associated")
      do i = lhn(xdim,LO), lhn(xdim,HI)
         do j = lhn(ydim,LO), lhn(ydim,HI)
            rgcs2 = (ax%x(i)**2+ax%y(j)**2)
            rgsv2 = rgcs2 + ax%z(:)**2
            rgsv  = sqrt(rgsv2)

            where (rgsv .notequals. 0.0)
               irgsv = 1./rgsv
            elsewhere
               irgsv = 0.0       ! first term of gphalo approaches 0 for rgsv close to 0
            endwhere
            rgsv  = rgsv / ahalo
            frgsv = 1. + rgsv**1.02

            gpdisk = -Mdisk / sqrt(rgcs2 + (adisk + sqrt(rgsv2 - rgcs2 + bdisk**2))**2)
            gphalo = -irgsv * Mhalo * rgsv**2.02 / frgsv
            gphalo = gphalo - mahalo * (-1.02 / khalo + log(khalo))
            gphalo = gphalo + mahalo * (-1.02 / frgsv + log(frgsv))
            gpblg  = -Mbulge / sqrt(rgsv2 + bbulge**2)
            gpot = (gpdisk + gpblg + gphalo) * newtong
            if (grav_log) then
               gpotdisk(i,j,:) = gpdisk * newtong
               gpothalo(i,j,:) = gphalo * newtong
               gpotbulge(i,j,:)= gpblg  * newtong
            endif
            gp(i,j,:) = gpot
         enddo
      enddo
      call my_deallocate(gpdisk)
      call my_deallocate(gphalo)
      call my_deallocate(gpblg)
      call my_deallocate(gpot)
      call my_deallocate(rgsv)
      call my_deallocate(rgsv2)
      call my_deallocate(irgsv)
      call my_deallocate(frgsv)
   end subroutine grav_allen_santillan

   subroutine grav_hdf5(cg, filename)

      use axes_M,             only: axes
      use common_hdf5,        only: data_gname, cg_cnt_aname, n_cg_name
      use constants,          only: ndims, xdim, ydim, zdim, LO, HI, cwdlen
      use dataio_pub,         only: die, msg, restarted_sim
      use diagnostics,        only: ma1d, my_allocate, my_deallocate
      use domain,             only: dom
      use func,               only: operator(.notequals.)
      use grid_cont,          only: grid_container
      use hdf5,               only: HID_T, HSIZE_T, H5F_ACC_RDONLY_F, H5T_NATIVE_DOUBLE, h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5gopen_f, h5gclose_f, h5dopen_f, h5dclose_f, h5dread_f
      use named_array_list,   only: qna, wna
      use read_attr,          only: read_attribute
      use restart_hdf5_v2,    only: cg_essentials
      use set_get_attributes, only: get_attr
      use units,              only: kpc, newtong

      implicit none

      type(grid_container),  pointer                      :: cg
      logical                                             :: grav_log
      integer                                             :: i,j, ia
      real                                                :: rgcs2, mahalo, khalo
      real, allocatable, dimension(:)                     :: gpdisk, gphalo, gpblg, gpot, rgsv, rgsv2, irgsv, frgsv
      character(len=cwdlen)                               :: filename
      logical                                             :: file_exist
      integer(HID_T)                                      :: file_id              !< File identifier
      integer(kind=4)                                     :: error
      real,                     dimension(:), allocatable :: rbuf
      integer(kind=4),          dimension(:), allocatable :: ibuf
      integer(HID_T)                                      :: cgl_g_id,  cg_g_id   !< cg list and cg group identifiers
      type(cg_essentials), dimension(:), allocatable      :: cg_res
      integer, allocatable, dimension(:)                  :: qr_lst
      integer(HSIZE_T), dimension(:), allocatable         :: dims
      integer(kind=8), dimension(xdim:zdim)               :: own_off, restart_off, o_size
      real, dimension(:,:,:),   allocatable               :: a3d
      integer(HID_T)                                      :: dset_id

      if (restarted_sim) then
         cg%gp = 0.0
         return
      endif

      inquire(file = trim(filename), exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[grav_hdf5] No previous run to start from'
         call die(msg)
      endif

      call h5open_f(error)

      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

      !call get_attr(file_id, "piernik", rbuf)

      !if (size(rbuf) /= 1) call die("[grav_hdf5] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")

      call h5gopen_f(file_id, data_gname, cgl_g_id, error)       ! open "/data"

      allocate(ibuf(1))
      call read_attribute(cgl_g_id, cg_cnt_aname, ibuf)          ! open "/data/cg_count"
      if (ibuf(1) <= 0) call die("[grav_hdf5] Empty cg list")

      allocate(cg_res(ibuf(1)))
      deallocate(ibuf)

      do ia = lbound(cg_res, dim=1), ubound(cg_res, dim=1)
         if (ia .gt. 1) then
            write(msg, '(3a)') '[grav_hdf5] Fow now only one cg is accepted in input'
            call die(msg)
         endif
         call h5gopen_f(cgl_g_id, n_cg_name(ia), cg_g_id, error) ! open "/data/grid_%08d, ia"
         allocate(dims(ndims))
         dims= (/ dom%n_d + 2*dom%nb, dom%n_d + 2*dom%nb, dom%n_d + 2*dom%nb /)
         call h5dopen_f(cg_g_id, "gpot", dset_id, error)
         allocate(a3d(dims(xdim), dims(ydim), dims(zdim)))
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, a3d, dims(:), error)
         cg%gp = a3d(cg%lhn(xdim,LO)+dom%nb+1:cg%lhn(xdim,HI)+dom%nb+1, cg%lhn(ydim,LO)+dom%nb+1:cg%lhn(ydim,HI)+dom%nb+1, cg%lhn(zdim,LO)+dom%nb+1:cg%lhn(zdim,HI)+dom%nb+1)
         deallocate(a3d)
         deallocate(dims)
         call h5dclose_f(dset_id, error)
         call h5gclose_f(cg_g_id, error)
      end do
      deallocate(cg_res)
      call h5gclose_f(cgl_g_id, error)
      call h5fclose_f(file_id, error)

    end subroutine grav_hdf5

   subroutine grav_hdf5_interp(gp, ax, lhn, flatten)

      use axes_M,             only: axes
      use constants,          only: ndims, xdim, ydim, zdim, LO, HI, cwdlen
      use dataio_pub,         only: die
      use diagnostics,        only: ma1d, my_allocate, my_deallocate
      use func,               only: operator(.notequals.)
      use units,              only: kpc, newtong
      use dataio_pub,         only: die, msg
      use hdf5,               only: HID_T, HSIZE_T, H5F_ACC_RDONLY_F, H5T_NATIVE_DOUBLE, h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5gopen_f, h5gclose_f, h5dopen_f, h5dclose_f, h5dread_f
      use set_get_attributes, only: get_attr
      use common_hdf5,        only: data_gname, cg_cnt_aname, n_cg_name
      use restart_hdf5_v2,    only: cg_essentials
      use named_array_list,   only: qna, wna
      use read_attr,          only: read_attribute
      use gravity,            only: nsub
      use mpisetup,           only: proc
      use domain,             only: dom

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten

      logical                                             :: grav_log
      integer                                             :: i,j, k, l, m, ia
      real                                                :: rgcs2, mahalo, khalo
      real, allocatable, dimension(:)                     :: gpdisk, gphalo, gpblg, gpot, rgsv, rgsv2, irgsv, frgsv
      character(len=cwdlen)                             :: filename
      logical                                           :: file_exist
      integer(HID_T)                                    :: file_id              !< File identifier
      integer(kind=4)                                   :: error
      real,                     dimension(:), allocatable :: rbuf
      integer(kind=4),          dimension(:), allocatable :: ibuf
      integer(HID_T)                                    :: cgl_g_id,  cg_g_id   !< cg list and cg group identifiers
      type(cg_essentials), dimension(:), allocatable    :: cg_res
      integer, allocatable, dimension(:)           :: qr_lst
      integer(HSIZE_T), dimension(:), allocatable  :: dims
      integer(kind=8), dimension(xdim:zdim)        :: own_off, restart_off, o_size
      real, dimension(:,:,:),   allocatable        :: a3d
      integer(HID_T)                               :: dset_id
      real                                         :: z, dz

      if (.not. associated(gp)) call die("[gravity_user:grav_hdf5_interp] gp is not associated")

      filename = 'previous.h5'
      inquire(file = trim(filename), exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[grav_hdf5] No previous run to start from'
         call die(msg)
      endif
      

      call h5open_f(error)


      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

      !call get_attr(file_id, "piernik", rbuf)

      !if (size(rbuf) /= 1) call die("[grav_hdf5] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")

      call h5gopen_f(file_id, data_gname, cgl_g_id, error)       ! open "/data"

      allocate(ibuf(1))
      call read_attribute(cgl_g_id, cg_cnt_aname, ibuf)          ! open "/data/cg_count"
      if (ibuf(1) <= 0) call die("[grav_hdf5] Empty cg list")

      allocate(cg_res(ibuf(1)))
      deallocate(ibuf)

      do ia = lbound(cg_res, dim=1), ubound(cg_res, dim=1)
         if (ia .gt. 1) then
            write(msg, '(3a)') '[grav_hdf5] Fow now only one cg is accepted in input'
            call die(msg)
         endif
         call h5gopen_f(cgl_g_id, n_cg_name(ia), cg_g_id, error) ! open "/data/grid_%08d, ia"
         allocate(dims(ndims))
         dims= (/ dom%n_d + 2*dom%nb, dom%n_d + 2*dom%nb, dom%n_d + 2*dom%nb /)
         call h5dopen_f(cg_g_id, "gpot", dset_id, error)
         allocate(a3d(dims(xdim), dims(ydim), dims(zdim)))
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, a3d, dims(:), error)

         i = ( ( dom%L_(xdim)/2. + ax%x(1) )  / dom%L_(xdim) * dom%n_d(xdim) ) + dom%nb + 0.5
         j = ( ( dom%L_(ydim)/2. + ax%y(1) )  / dom%L_(ydim) * dom%n_d(ydim) ) + dom%nb + 0.5
         dz = dom%L_(zdim) / dom%n_d(zdim)
         l=0
         m = 1
         do k = lhn(zdim, LO), lhn(zdim, HI)
            z = (m - dom%n_d(xdim)/2 - dom%nb - 0.5)  * dom%L_(zdim) / dom%n_d(xdim)
            if ( (k .lt. lhn(zdim, LO) + nsub/2) .or. (k .ge. lhn(zdim, HI) - nsub/2)) then
               gp(1,1,k) = a3d(i, j, m)
               !print *, 'nope', i, j, m, gp(1,1,k)
            else
               gp(1,1,k) = a3d(i,j,m) * (1 - (ax%z(k) - z) / dz ) + a3d(i,j,m+1) * (ax%z(k) - z) / dz
               !print *, proc, i, j, m, ax%z(k), z, a3d(i,j,m), a3d(i,j,m+1), gp(1,1,k)
               l=l+1
               if (l .eq. nsub) then
                  m = m+1
                  l=0
               endif
            endif
         end do
            
         deallocate(a3d)
         call h5dclose_f(dset_id, error)
         call h5gclose_f(cg_g_id, error)
      end do

      deallocate(cg_res)
      call h5gclose_f(cgl_g_id, error)
      call h5fclose_f(file_id, error)

    end subroutine grav_hdf5_interp

   subroutine grav_nbody(cg)

      use named_array_list, only: qna
      use constants, only: sgp_n
      use grid_cont, only: grid_container
      use multigrid_gravity, only: multigrid_solve_grav

      implicit none

      type(grid_container),  pointer :: cg
      integer(kind=4), dimension(0)   :: i_sg_dens

      !print *, cg%q(qna%ind(sgp_n))%arr(0:3, 0:3, 0:3)
      call multigrid_solve_grav(i_sg_dens)
      !print *, cg%q(qna%ind(sgp_n))%arr(0:3, 0:3, 0:3)
      !print *, '               NBODY'
      cg%gp = cg%q(qna%ind(sgp_n))%arr

    end subroutine grav_nbody

    subroutine grav_nbody_interp(gp, ax, lhn, flatten)

      use axes_M,      only: axes
      use named_array_list, only: qna
      use constants, only: sgp_n, ndims, xdim, ydim, zdim, LO, HI
      use grid_cont, only: grid_container
      use multigrid_gravity, only: multigrid_solve_grav
      use gravity,         only: nsub
      use cg_list,   only: cg_list_element
      use func,     only: operator(.notequals.)
      use cg_leaves, only: leaves
      use dataio_pub,  only: die
      use mpisetup, only: proc

      implicit none

      real, dimension(:,:,:), pointer   :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(0)     :: i_sg_dens
      integer                           :: i, j, k, l, m
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real                           :: xi, yj
      logical, optional, intent(in)  :: flatten

      if (.not. associated(gp)) call die("[gravity_user:grav_nbody_interp] gp is not associated")
      !call multigrid_solve_grav(i_sg_dens)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         l=0
         m=-4
         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            if (cg%y(j) .notequals. ax%y(1)) cycle
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               if (cg%x(i) .notequals. ax%x(1)) cycle
               do k = lhn(zdim, LO), lhn(zdim, HI)
                  gp(1,1,k) = cg%q(qna%ind(sgp_n))%arr(i,j,m) * (1 - (ax%z(k)+ cg%dl(zdim)/2 -cg%z(m)) * cg%idl(zdim)) + cg%q(qna%ind(sgp_n))%arr(i,j,m+1) * (ax%z(k)+ cg%dl(zdim)/2 -cg%z(m)) * cg%idl(zdim)
                  print *, proc, ax%z(k) + cg%dl(zdim)/2, cg%z(m), cg%z(m+1), cg%q(qna%ind(sgp_n))%arr(i,j,m), cg%q(qna%ind(sgp_n))%arr(i,j,m+1), gp(1,1,k)
                  l=l+1
                  if (l .eq. nsub) then
                     m = m+1
                     l=0
                  endif
               end do
            end do
         end do

         cgl => cgl%nxt
      enddo

    end subroutine grav_nbody_interp

    
end module gravity_user
