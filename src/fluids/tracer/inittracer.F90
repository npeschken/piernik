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
!! \brief Initialization of tracer fluids
!!
!! In this module, the following variables are defined:
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>variable</b></td><td width="135pt"><b>type</b></td><td width="200pt"><b>description</b></td></tr>
!! <tr><td>ntracers</td><td>integer</td><td>Number of tracer fluids enabled in the simulation. Allowed values: 0-10. Default value: 0 (no tracer fluids enabled).</td></tr>
!! <tr><td>iarr_trc</td><td>integer array</td><td>Array of indices of tracer fluids in the main data array.</td></tr>
!! </table>
!! \n \n
!<

module inittracer

   use constants, only: INT4
   implicit none

   private
   public :: init_tracer, cleanup_tracer, tracer_index, iarr_trc, trace_fluid, tracers_max, ntracers

   integer(kind=INT4), dimension(:), allocatable :: iarr_trc     !< Array of indices to tracers
   integer(kind=INT4), dimension(:), allocatable :: trace_fluid  !< Which fluid to trace? Currently possible: ionized, neutral, dust
   integer(kind=INT4), protected :: ntracers                     !< Number of tracers
   integer(kind=INT4), parameter :: tracers_max = 10_INT4        !< Maximum allowed number of tracer fluids (arbitrary, some hardcoded bits in data_hdf5 depend on it)

   ! Note: trace_fluid is allocatable and listed in FLUID_TRACER. It MUST be allocated
   ! before any read(nml=FLUID_TRACER) statement, otherwise the read will target an
   ! unallocated array and crash. See init_tracer.
   namelist /FLUID_TRACER/ ntracers, trace_fluid

contains

!>
!! \brief Read tracer fluid parameters from namelist and broadcast via MPI
!!
!! Reads the FLUID_TRACER namelist from the configuration file and broadcasts
!! the parameters to all MPI processes.
!!
!! \n \n
!! @b FLUID_TRACER
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>ntracers</td><td>0</td><td>1-10</td><td>\copydoc inittracer::ntracers</td></tr>
!! <tr><td>trace_fluid</td><td>1</td><td>1-3</td><td>\copydoc inittracer::trace_fluid</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_tracer(num_fluids)

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: V_INFO
      use dataio_pub, only: warn, nh, die, printinfo, msg
      use mpisetup,   only: master, slave, ibuff

      implicit none

      integer, intent(in) :: num_fluids

      integer :: i

      if (tracers_max > ubound(ibuff, 1) - 1) call die("[inittracer:init_tracer] Too big tracers_max value - the trace fluid array would not fit ibuff array.")
      allocate(trace_fluid(tracers_max))  ! must precede namelist read; see module-level note

      ntracers = 0_INT4  ! no tracer fluid by default
      trace_fluid = 1_INT4  ! trace first fluid by default

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=FLUID_TRACER)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=FLUID_TRACER, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "FLUID_TRACER")
         read(nh%cmdl_nml,nml=FLUID_TRACER, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "FLUID_TRACER", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=FLUID_TRACER)
         close(nh%lun)
         call nh%compare_namelist()

         ibuff(1)   = ntracers
         ibuff(2:1+tracers_max) = trace_fluid
      endif

      call piernik_MPI_Bcast(ibuff)

      if (slave) then

         ntracers    = int(ibuff(1), kind=INT4)
         trace_fluid = int(ibuff(2:1+tracers_max), kind=INT4)

      endif

      if (ntracers < 0) then
         call warn("[inittracer:init_tracer] Negative number of tracers provided. Reset to 0")
         ntracers = 0
      elseif (ntracers > tracers_max) then
         write(msg, "(A,I0)") "[inittracer:init_tracer] Too many tracers provided. Reset to maximum allowed value: ", tracers_max
         call warn(msg)
         ntracers = tracers_max
      endif

      if (any(trace_fluid(1:ntracers) < 1) .or. any(trace_fluid(1:ntracers) > num_fluids)) then
         call warn("[inittracer:init_tracer] Some values of trace_fluid(:) are invalid. Reset them to the default value (1).")
         where (trace_fluid(1:ntracers) < 1 .or. trace_fluid(1:ntracers) > num_fluids) trace_fluid(1:ntracers) = 1_INT4
      endif

      if (ntracers == 0) then
         msg = "[inittracer:init_tracer] No tracer fluids are enabled"
      else
         write(msg, "(A,I0, A)") "[inittracer:init_tracer] Number of tracer fluids enabled: ", ntracers, ". Traced fluids: "
         do i = 1, ntracers
            write(msg(len_trim(msg)+1:), "(' ',I0)") trace_fluid(i)
         enddo
      endif
      if (master) call printinfo(msg, V_INFO)

      allocate(iarr_trc(ntracers))

   end subroutine init_tracer

!> \brief Deallocate tracer fluid arrays
!!
!! Deallocates arrays allocated during initialization.
   subroutine cleanup_tracer

      implicit none

      if (allocated(iarr_trc)) deallocate(iarr_trc)
      if (allocated(trace_fluid)) deallocate(trace_fluid)

   end subroutine cleanup_tracer

!> \brief Set indices of tracer fluids in the main data array
!!
!! Assigns index ranges for tracer fluid variables and computes the iarr_trc array.
   subroutine tracer_index(flind)

      use constants,  only: I_ONE
      use fluidtypes, only: var_numbers

      implicit none

      type(var_numbers), intent(inout) :: flind
      integer :: i

      flind%trc%beg  = flind%all + I_ONE
      flind%trc%all  = ntracers
      flind%all      = flind%all + flind%trc%all
      flind%trc%end  = flind%all

      if (ntracers > 0) iarr_trc = int([(i, i = 0, ntracers-1)], kind=INT4) + flind%trc%beg

      ! Tracer index position is set to -1 since tracers are initialized last
      ! and no other components depend on their index position
!      flind%components = flind%components + 1
!      flind%trc%pos    = flind%components
      flind%trc%pos    = -1

   end subroutine tracer_index

end module inittracer
