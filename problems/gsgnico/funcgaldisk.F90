! $Id: funcgaldisk.F90 7063 2012-11-15 10:55:24Z wolt $
#include "user_macro.h"

module funcgaldisk

   use constants, only: ndims

   implicit none

   private
   public :: take_care_about_bnd_indx, update_b_with_A, v0, g0, aktugalpos, aktugalvel

   real, dimension(ndims) :: v0, g0, aktugalpos, aktugalvel !< initial and current velocity and position of galaxy

   contains

   subroutine take_care_about_bnd_indx(ii,nn,iiu,iid,sfq)

      use constants, only: LO, HI

      implicit none

      integer                           :: ii, iiu, iid
      integer(kind=4), dimension(LO:HI) :: nn
      real                              :: sfq

      iiu = ii+1
      iid = ii-1
      sfq = 0.5
      if (ii == nn(LO)) iid = ii
      if (ii == nn(HI)) iiu = ii
      if (any(ii == nn)) sfq = 1.0

   end subroutine take_care_about_bnd_indx

   subroutine update_b_with_A(cg,an,bn,initpb)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use domain,     only: dom
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: an, bn
      logical,             optional, intent(in) :: initpb
      integer, dimension(ndims)                 :: nl, nu
      real, dimension(:,:,:,:), pointer         :: pA, pB

      nl(:) = cg%lhn(:,LO) + dom%D_(:)
      nu(:) = cg%lhn(:,HI) - dom%D_(:)
      pA => cg%w(an)%arr
      pB => cg%w(bn)%arr
      if (present(initpb)) pB = 0.0
      pB(xdim,:,:nu(ydim),:nu(zdim)) = pB(xdim,:,:nu(ydim),:nu(zdim)) + (pA(zdim,:,nl(ydim):cg%lhn(ydim,HI),:nu(zdim)) - pA(zdim,:,:nu(ydim),:nu(zdim)))*cg%idl(ydim) &
         &                                                            - (pA(ydim,:,:nu(ydim),nl(zdim):cg%lhn(zdim,HI)) - pA(ydim,:,:nu(ydim),:nu(zdim)))*cg%idl(zdim)
      pB(ydim,:nu(xdim),:,:nu(zdim)) = pB(ydim,:nu(xdim),:,:nu(zdim)) + (pA(xdim,:nu(xdim),:,nl(zdim):cg%lhn(zdim,HI)) - pA(xdim,:nu(xdim),:,:nu(zdim)))*cg%idl(zdim) &
         &                                                            - (pA(zdim,nl(xdim):cg%lhn(xdim,HI),:,:nu(zdim)) - pA(zdim,:nu(xdim),:,:nu(zdim)))*cg%idl(xdim)
      pB(zdim,:nu(xdim),:nu(ydim),:) = pB(zdim,:nu(xdim),:nu(ydim),:) + (pA(ydim,nl(xdim):cg%lhn(xdim,HI),:nu(ydim),:) - pA(ydim,:nu(xdim),:nu(ydim),:))*cg%idl(xdim) &
         &                                                            - (pA(xdim,:nu(xdim),nl(ydim):cg%lhn(ydim,HI),:) - pA(xdim,:nu(xdim),:nu(ydim),:))*cg%idl(ydim)

   end subroutine update_b_with_A

end module funcgaldisk
