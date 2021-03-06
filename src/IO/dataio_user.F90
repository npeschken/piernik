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
#include "piernik.h"
!>
!! \brief module providing hooks for user I/O addons
!<
module dataio_user

   implicit none

   public  ! QA_WARN there are no secrets here

#ifdef HDF5
   interface
      subroutine vars_hdf5(var, tab, ierrh, cg)

         use grid_cont,  only: grid_container

         implicit none

         character(len=*),               intent(in)    :: var
         real, dimension(:,:,:),         intent(inout) :: tab
         integer,                        intent(inout) :: ierrh
         type(grid_container), pointer,  intent(in)    :: cg

      end subroutine vars_hdf5
   end interface

   interface
      subroutine add_attr(file_id)

         use hdf5, only: HID_T

         implicit none

         integer(HID_T), intent(in) :: file_id

      end subroutine add_attr
   end interface

   interface
      subroutine plt_attr(gr2_id,vdname)

         use hdf5, only: HID_T

         implicit none

         integer(HID_T),   intent(in) :: gr2_id
         character(len=*), intent(in) :: vdname

      end subroutine plt_attr
   end interface
#endif /* HDF5 */

   interface
      subroutine postout(output,dump)

         use constants, only: RES, TSL

         implicit none

         integer(kind=4),             intent(in) :: output
         logical, dimension(RES:TSL), intent(in) :: dump

      end subroutine postout
   end interface

   interface
      subroutine tsl_out(user_vars, tsl_names)

         implicit none

         real,             dimension(:), intent(inout), allocatable           :: user_vars
         character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names

      end subroutine tsl_out
   end interface

   interface
      subroutine add_data
         implicit none
      end subroutine add_data
   end interface

#ifdef HDF5
   procedure(add_attr),  pointer :: user_attrs_rd         => Null() !< A hook to read user-defined attributes from the restart file
   procedure(add_attr),  pointer :: user_attrs_wr         => Null() !< A hook to write user-defined attributes to the restart file
   procedure(plt_attr),  pointer :: user_plt_attrs        => Null()
   procedure(vars_hdf5), pointer :: user_vars_hdf5        => Null()
#endif /* HDF5 */
   procedure(add_data),  pointer :: user_attrs_pre        => Null()
   procedure(tsl_out),   pointer :: user_tsl              => Null()
   procedure(add_data),  pointer :: user_reg_var_restart  => Null() !< A hook to have user-defined fields in the restart files
   procedure(postout),   pointer :: user_post_write_data  => Null()

end module dataio_user
