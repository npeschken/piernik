! Code Copyright (C) 2006 Michal Hanasz
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
!  References:
!
!  A multi-state HLLD approximate Riemann solver for ideal magnetohydrodynamics.
!  Takahiro Miyoshi, Kanya Kusano
!  Journal of Computational Physics 208 (2005) 315-344
!
!  ->Solve one dimensional Riemann problem using adiabatic HLLD scheme
!
!  Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!  Dr. Artur Gawryszczak, CAMK, Warszawa.
!---------------------------------------------------------------------------------------------------------------------------


#include "piernik.h"
!>
!!  \brief This module implements HLLD Riemann solver following the work of Miyoshi & Kusano (2005)
!<
module hlld
! pulled by RIEMANN

   implicit none

   private
   public :: riemann_wrap

contains
!>
!! \brief Wrapper for the Riemann solver that takes care of fluid differences
!!
!! OPT: check if passing pointers here will improve performance
!<

   subroutine riemann_wrap(ql, qr, b_cc_l, b_cc_r, psi_l, psi_r, flx, mag_cc, psi_cc)

      use constants,  only: xdim, zdim, DIVB_HDC
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use global,     only: divB_0_method

      implicit none

      real, dimension(:,:), target, intent(in)  :: ql, qr          ! left and right fluid states
      real, dimension(:,:), target, intent(in)  :: b_cc_l, b_cc_r  ! left and right magnetic field states (relevant only for IONIZED fluid)
      real, dimension(:,:), target, intent(in)  :: psi_l, psi_r    ! left and right psi field states (relevant only for GLM method)
      real, dimension(:,:), target, intent(out) :: flx, mag_cc, psi_cc ! output fluxes: fluid, magnetic field and psi

      integer :: i
      class(component_fluid), pointer :: fl

      real, dimension(size(b_cc_l,1), size(b_cc_l,2)), target :: b0, bf0
      real, dimension(size(psi_l, 1), size(psi_l, 2)), target :: p0, pf0
      real, dimension(:,:), pointer :: p_flx, p_ql, p_qr
      real, dimension(:,:), pointer :: p_bcc, p_bccl, p_bccr
      real, dimension(:,:), pointer :: p_psif, p_psi_l, p_psi_r

      do i = 1, flind%fluids
         fl    => flind%all_fluids(i)%fl
         p_flx => flx(:, fl%beg:fl%end)
         p_ql  => ql(:, fl%beg:fl%end)
         p_qr  => qr(:, fl%beg:fl%end)
         if (fl%is_magnetized) then
            p_bccl => b_cc_l(:, xdim:zdim)
            p_bccr => b_cc_r(:, xdim:zdim)
            p_bcc  => mag_cc(:, xdim:zdim)
            if (divB_0_method == DIVB_HDC) then
               p_psi_l => psi_l(:,:)
               p_psi_r => psi_r(:,:)
               p_psif  => psi_cc(:,:)
            else  ! CT
               p0 = 0.
               p_psi_l => p0
               p_psi_r => p0
               p_psif  => pf0
            endif
         else ! ignore all magnetic field
            b0 = 0.
            p_bccl => b0
            p_bccr => b0
            p_bcc  => bf0
            p0 = 0.
            p_psi_l => p0
            p_psi_r => p0
            p_psif  => pf0
         endif

         call riemann_hlld(p_flx, p_ql, p_qr, p_bcc, p_bccl, p_bccr, p_psi_l, p_psi_r, p_psif, fl%gam)
      enddo

   end subroutine riemann_wrap

   subroutine riemann_hlld(f,ul,ur,b_cc,b_ccl,b_ccr,psil,psir,psi,gamma)

    ! external procedures

    use constants,  only: half, zero, one, xdim, ydim, zdim, idn, imx, imy, imz, ien, DIVB_HDC
    use func,       only: operator(.notequals.), operator(.equals.)
    use global,     only: divB_0_method
    use hdc,        only: chspeed

    ! arguments

    implicit none

    real, dimension(:,:), pointer, intent(out)   :: f
    real, dimension(:,:), pointer, intent(in)    :: ul, ur
    real, dimension(:,:), pointer, intent(out)   :: b_cc
    real, dimension(:,:), pointer, intent(in)    :: b_ccl, b_ccr
    real, dimension(:,:), pointer, intent(in)    :: psil, psir
    real, dimension(:,:), pointer, intent(out)   :: psi
    real,                          intent(in)    :: gamma

    ! Local variables

    integer                                      :: i
    real, parameter                              :: four = 4.0
    real                                         :: sm, sl, sr
    real                                         :: alfven_l, alfven_r, c_fastm, gampr_l, gampr_r
    real                                         :: slsm, srsm, slvxl, srvxr, smvxl, smvxr, srmsl, dn_l, dn_r
    real                                         :: b_lr, b_lrgam, magprl, magprr, prt_star, b_sig, enl, enr
    real                                         :: coeff_1, dn_lsqt, dn_rsqt, add_dnsq, mul_dnsq
    real                                         :: vb_l, vb_starl, vb_r, vb_starr, vb_2star
    real                                         :: prl, prr


    ! Local arrays

    real, dimension(size(f, 2))                  :: fl, fr
    real, dimension(size(f, 2))                  :: u_starl, u_starr, u_2starl, u_2starr
    real, dimension(xdim:zdim)                   :: v_2star, v_starl, v_starr
    real, dimension(xdim:zdim)                   :: b_cclf, b_ccrf
    real, dimension(xdim:zdim)                   :: b_starl, b_starr, b_2star
    logical                                      :: has_energy
    real                                         :: ue

#ifdef ISO
#  error Isothermal EOS is not implemented yet in this Riemann solver.
#endif /* ISO */

    ! SOLVER

    ! suppress complains caused by -Wmaybe-uninitialized
    b_cclf = 0.
    b_ccrf = 0.
    !if (divB_0_method /= DIVB_HDC) b_cc(xdim,:) = 0.
    has_energy = (ubound(ul, dim=2) >= ien)
    ue = 0.

    ! Eq. 42, Dedner et al.
    if (divB_0_method .eq. DIVB_HDC) then
       psi(:, 1)     = chspeed * chspeed * half*((b_ccr(:, xdim)+b_ccl(:, xdim)) - (psir(:, 1)-psil(:, 1))/chspeed)
       b_cc(:, xdim) = half*((psir(:, 1)+psil(:, 1)) - chspeed*(b_ccr(:, xdim)-b_ccl(:, xdim)))
    else
       b_cc(:, xdim) = 0.
    endif

    do i = 1, size(f, 1)

       ! Left and right states of magnetic pressure

       magprl  =  half*sum(b_ccl(i, xdim:zdim)*b_ccl(i, xdim:zdim))
       magprr  =  half*sum(b_ccr(i, xdim:zdim)*b_ccr(i, xdim:zdim))

       ! Left and right states of total pressure
       ! From fluidupdate.F90, utoq() (1) is used in hydro regime and (2) in MHD regime. In case of vanishing magnetic fields the magnetic components do not contribute and hydro results are obtained trivially.

       if (has_energy) then

          prl = ul(i, ien) + magprl ! ul(i, ien) is the left state of gas pressure
          prr = ur(i, ien) + magprr ! ur(i, ien) is the right state of gas pressure

          ! Left and right states of energy Eq. 2.

          enl = (ul(i, ien)/(gamma -one)) + half*ul(i, idn)*sum(ul(i, imx:imz)**2) + half*sum(b_ccl(i, xdim:zdim)**2)
          enr = (ur(i, ien)/(gamma -one)) + half*ur(i, idn)*sum(ur(i, imx:imz)**2) + half*sum(b_ccr(i, xdim:zdim)**2)

          ! Left and right states of gamma*p_gas

          gampr_l = gamma*ul(i, ien)
          gampr_r = gamma*ur(i, ien)

       else ! this is for DUST (presureless fluid)

          ! check if it is consistent
          prl = magprl
          prr = magprr

          enl = half*ul(i, idn)*sum(ul(i, imx:imz)**2) + half*sum(b_ccl(i, xdim:zdim)**2)
          enr = half*ur(i, idn)*sum(ur(i, imx:imz)**2) + half*sum(b_ccr(i, xdim:zdim)**2)

          gampr_l = 0.
          gampr_r = 0.

       endif

       ! Left and right states of fast magnetosonic waves Eq. 3

       c_fastm = sqrt(half*max( &
             ((gampr_l+sum(b_ccl(i, xdim:zdim)**2)) + sqrt((gampr_l+sum(b_ccl(i, xdim:zdim)**2))**2-(four*gampr_l*b_ccl(i, xdim)**2)))/ul(i, idn), &
             ((gampr_r+sum(b_ccr(i, xdim:zdim)**2)) + sqrt((gampr_r+sum(b_ccr(i, xdim:zdim)**2))**2-(four*gampr_r*b_ccr(i, xdim)**2)))/ur(i, idn)) )

       ! Estimates of speed for left and right going waves Eq. 67

       sl  =  min(ul(i, imx) ,ur(i, imx)) - c_fastm
       sr  =  max(ul(i, imx), ur(i, imx)) + c_fastm

       ! Left flux

       fl(idn) = ul(i, idn)*ul(i, imx)
       fl(imx) = ul(i, idn)*ul(i, imx)**2 + prl - b_ccl(i, xdim)**2  ! Total left state of pressure, so prl
       fl(imy:imz) = ul(i, idn)*ul(i, imy:imz)*ul(i, imx) - b_ccl(i, xdim)*b_ccl(i, ydim:zdim)
       if (has_energy) fl(ien) = (enl + prl)*ul(i, imx) - b_ccl(i, xdim)*(sum(ul(i, imx:imz)*b_ccl(i, xdim:zdim))) ! Total left state of pressure, so prl
       b_cclf(ydim:zdim) = b_ccl(i, ydim:zdim)*ul(i, imx) - b_ccl(i, xdim)*ul(i, imy:imz)

       ! Right flux

       fr(idn) = ur(i, idn)*ur(i, imx)
       fr(imx) = ur(i, idn)*ur(i, imx)**2 + prr - b_ccr(i, xdim)**2  ! Total right state of pressure, so prl
       fr(imy:imz) = ur(i, idn)*ur(i, imy:imz)*ur(i, imx) - b_ccr(i, xdim)*b_ccr(i, ydim:zdim)
       if (has_energy) fr(ien) = (enr + prr)*ur(i, imx) - b_ccr(i, xdim)*(sum(ur(i, imx:imz)*b_ccr(i, xdim:zdim)))  ! Total right state of pressure, so prl
       b_ccrf(ydim:zdim) = b_ccr(i, ydim:zdim)*ur(i, imx) - b_ccr(i, xdim)*ur(i, imy:imz)

       ! HLLD fluxes

       if (sl .ge.  zero) then
          f(i,:)  =  fl
          b_cc(i, ydim:zdim) = b_cclf(ydim:zdim)
       else if (sr .le.  zero) then
          f(i,:)  =  fr
          b_cc(i, ydim:zdim) = b_ccrf(ydim:zdim)
       else

          ! Speed of contact discontinuity Eq. 38
          ! Total left and right states of pressure, so prr and prl sm_nr/sm_dr

          if ((sr - ur(i, imx))*ur(i, idn) .equals. (sl - ul(i, imx))*ul(i, idn)) then
             sm = (sl + sr) / 2.
          else
             sm =   ( ((sr - ur(i, imx))*ur(i, idn)*ur(i, imx) - prr) - &
                  &   ((sl - ul(i, imx))*ul(i, idn)*ul(i, imx) - prl) ) / &
                  &   ((sr - ur(i, imx))*ur(i, idn) - &
                  &    (sl - ul(i, imx))*ul(i, idn))
          endif

          ! Speed differences

          slsm  =  sl - sm
          srsm  =  sr - sm

          slvxl  =  sl - ul(i, imx)
          srvxr  =  sr - ur(i, imx)

          smvxl  =  sm - ul(i, imx)
          smvxr  =  sm - ur(i, imx)

          srmsl  =  sr - sl

          ! Co-efficients

          dn_l     =  ul(i, idn)*slvxl
          dn_r     =  ur(i, idn)*srvxr
          b_lr     =  b_ccl(i, xdim)*b_ccr(i, xdim)
          b_lrgam  =  b_lr/gamma

          ! Pressure of intermediate state Eq. (23)

          prt_star  =  half*((prl+dn_l*smvxl) + (prr+dn_r*smvxr))  !< Check for 0.5. Total left and right states of pressure, so prr and prl

          ! Normal components of velocity and magnetic field

          v_starl(xdim)  =  sm
          v_starr(xdim)  =  sm

          b_starl(xdim)  =  b_ccl(i, xdim)
          b_starr(xdim)  =  b_ccr(i, xdim)

          ! Transversal components of magnetic field for left states (Eq. 45 & 47), taking degeneracy into account
          coeff_1  =  dn_l*slsm - b_lr
          if (has_energy) ue = ul(i, ien)
          if ((coeff_1 .notequals. zero) .and. b_lrgam .le. ue) then  ! Left state of gas pressure, so ul(i, ien)
             b_starl(ydim:zdim) = b_ccl(i, ydim:zdim) * (dn_l*slvxl - b_lr)/coeff_1
          else
             ! Calculate HLL left states
             b_starl(ydim:zdim) = ((sr*b_ccr(i, ydim:zdim) - sl*b_ccl(i, ydim:zdim)) - (b_ccrf(ydim:zdim) - b_cclf(ydim:zdim)))/srmsl
          endif

          coeff_1  =  dn_r*srsm - b_lr
          if (has_energy) ue = ur(i, ien)
          if ((coeff_1 .notequals. zero) .and. b_lrgam .le. ue) then  ! Right state of gas pressure, so ur(i, ien)
             b_starr(ydim:zdim) = b_ccr(i, ydim:zdim) * (dn_r*srvxr - b_lr)/coeff_1
          else
             ! Calculate HLL right states
             b_starr(ydim:zdim) = ((sr*b_ccr(i, ydim:zdim) - sl*b_ccl(i, ydim:zdim)) - (b_ccrf(ydim:zdim) - b_cclf(ydim:zdim)))/srmsl
          endif

          ! Transversal components of velocity Eq. 42
          v_starl(ydim:zdim) = ul(i, imy:imz)
          if (b_ccl(i, xdim) .notequals. 0.) v_starl(ydim:zdim) = v_starl(ydim:zdim) + b_ccl(i, xdim)/dn_l * (b_ccl(i, ydim:zdim) - b_starl(ydim:zdim))
          v_starr(ydim:zdim) = ur(i, imy:imz)
          if (b_ccr(i, xdim) .notequals. 0.) v_starr(ydim:zdim) = v_starr(ydim:zdim) + b_ccr(i, xdim)/dn_r * (b_ccr(i, ydim:zdim) - b_starr(ydim:zdim))

          ! Dot product of velocity and magnetic field

          vb_l  =  sum(ul(i, imx:imz)*b_ccl(i, :))
          vb_r  =  sum(ur(i, imx:imz)*b_ccr(i, :))
          vb_starl  =  sum(v_starl*b_starl)
          vb_starr  =  sum(v_starr*b_starr)


          ! Left intermediate state conservative form

          u_starl(idn)  =  dn_l/slsm
          u_starl(imx:imz)  =  u_starl(idn)*v_starl

          ! Right intermediate state conservative form

          u_starr(idn)  =  dn_r/srsm
          u_starr(imx:imz)  =  u_starr(idn)*v_starr

          ! Total energy of left and right intermediate states Eq. (48)

          if (has_energy) then
             u_starl(ien) = (slvxl*enl - prl*ul(i, imx) + prt_star*sm + b_ccl(i, xdim)*(vb_l - vb_starl))/slsm  ! Total left state of pressure
             u_starr(ien) = (srvxr*enr - prr*ur(i, imx) + prt_star*sm + b_ccr(i, xdim)*(vb_r - vb_starr))/srsm  ! Total right state of pressure
          endif

          ! Cases for B_x .ne. and .eq. zero

          if (abs(b_ccl(i, xdim)) > zero) then

             ! Left and right Alfven waves velocity Eq. 51

             dn_lsqt  =  sqrt(u_starl(idn))
             dn_rsqt  =  sqrt(u_starr(idn))

             alfven_l  =  sm - abs(b_ccl(i, xdim))/dn_lsqt
             alfven_r  =  sm + abs(b_ccr(i, xdim))/dn_rsqt

             ! Intermediate discontinuities

             if (alfven_l > zero) then

                ! Left intermediate flux Eq. 64

                f(i, :) = fl + sl*(u_starl - [ ul(i, idn), ul(i, idn)*ul(i, imx:imz), enl ] )
                b_cc(i, ydim:zdim) = b_cclf(ydim:zdim) + sl*(b_starl(ydim:zdim) - b_ccl(i, ydim:zdim))

             else if (alfven_r < zero) then

                ! Right intermediate flux Eq. 64

                f(i, :) = fr + sr*(u_starr - [ ur(i, idn), ur(i, idn)*ur(i, imx:imz), enr ] )
                b_cc(i, ydim:zdim) = b_ccrf(ydim:zdim) + sr*(b_starr(ydim:zdim) - b_ccr(i, ydim:zdim))

             else ! alfven_l .le. zero .le. alfven_r

                ! Arrange for sign of normal component of magnetic field

                if (b_ccl(i, xdim) .ge. zero) then

                   b_sig = one

                else

                   b_sig = -one

                endif

                ! Sum and product of density square-root

                add_dnsq  =  dn_lsqt + dn_rsqt
                mul_dnsq  =  dn_lsqt*dn_rsqt

                ! Components of velocity Eq. 39, 59, 60 and magnetic field Eq. 61, 62

                v_2star(xdim)      = sm
                v_2star(ydim:zdim) = ((dn_lsqt*v_starl(ydim:zdim) + dn_rsqt*v_starr(ydim:zdim)) + b_sig*(b_starr(ydim:zdim) - b_starl(ydim:zdim)))/add_dnsq

                b_2star(xdim)      = half * (b_ccl(i, xdim) + b_ccr(i, xdim))

                b_2star(ydim:zdim) = ((dn_lsqt*b_starr(ydim:zdim) + dn_rsqt*b_starl(ydim:zdim)) + b_sig*mul_dnsq*(v_starr(ydim:zdim) - v_starl(ydim:zdim)))/add_dnsq

                ! Dot product of velocity and magnetic field

                vb_2star  =  sum(v_2star*b_2star)

                ! Choose right Alfven wave according to speed of contact discontinuity

                if (sm >= zero) then
                   ! Conservative variables for left Alfven intermediate state
                   u_2starl(idn)  =  u_starl(idn)
                   u_2starl(imx:imz)  =  u_starl(idn)*v_2star

                   ! Energy of Alfven intermediate state Eq. 63
                   if (has_energy) u_2starl(ien)  =  u_starl(ien) - b_sig*dn_lsqt*(vb_starl - vb_2star)

                endif

                if (sm <= zero) then
                   ! Conservative variables for right Alfven intermediate state
                   u_2starr(idn)  =  u_starr(idn)
                   u_2starr(imx:imz)  =  u_starr(idn)*v_2star

                   ! Energy of Alfven intermediate state Eq. 63
                   if (has_energy) u_2starr(ien)  =  u_starr(ien) + b_sig*dn_rsqt*(vb_starr - vb_2star)

                endif

                if (sm > zero) then
                   ! Left Alfven intermediate flux Eq. 65
                   f(i, :) = fl + alfven_l*u_2starl - (alfven_l - sl)*u_starl - sl* [ ul(i, idn), ul(i, idn)*ul(i, imx:imz), enl ]
                   b_cc(i, ydim:zdim) = b_cclf(ydim:zdim) + alfven_l*b_2star(ydim:zdim) - (alfven_l - sl)*b_starl(ydim:zdim) - sl*b_ccl(i, ydim:zdim)

                else if (sm < zero) then
                   ! Right Alfven intermediate flux Eq. 65
                   f(i, :) = fr + alfven_r*u_2starr - (alfven_r - sr)*u_starr - sr* [ ur(i, idn), ur(i, idn)*ur(i, imx:imz), enr ]
                   b_cc(i, ydim:zdim) = b_ccrf(ydim:zdim) + alfven_r*b_2star(ydim:zdim) - (alfven_r - sr)*b_starr(ydim:zdim) - sr*b_ccr(i, ydim:zdim)

                else ! sm = 0

                   ! Left and right Alfven intermediate flux Eq. 65
                   f(i, :) = half*( &
                        (fl + alfven_l*u_2starl - (alfven_l - sl)*u_starl - sl* [ ul(i, idn), ul(i, idn)*ul(i, imx:imz), enl ]) + &
                        (fr + alfven_r*u_2starr - (alfven_r - sr)*u_starr - sr* [ ur(i, idn), ur(i, idn)*ur(i, imx:imz), enr ]))

                   b_cc(i, ydim:zdim) = half*( &
                        (b_cclf(ydim:zdim) + alfven_l*b_2star(ydim:zdim) - (alfven_l - sl)*b_starl(ydim:zdim) - sl*b_ccl(i, ydim:zdim)) + &
                        (b_ccrf(ydim:zdim) + alfven_r*b_2star(ydim:zdim) - (alfven_r - sr)*b_starr(ydim:zdim) - sr*b_ccr(i, ydim:zdim)))

                endif  ! sm = 0

             endif  ! alfven_l .le. 0 and alfven_r .ge. 0

          else ! B_x = 0

             ! Intermediate state for B_x = 0

             if (sm > zero) then

                ! Left intermediate flux Eq. 64

                f(i, :)  =  fl + sl*(u_starl - [ ul(i, idn), ul(i, idn)*ul(i, imx:imz), enl ])
                b_cc(i, ydim:zdim) = b_cclf(ydim:zdim) + sl*(b_starl(ydim:zdim) - b_ccl(i, ydim:zdim))

             else if (sm < zero) then

                f(i, :)  =  fr + sr*(u_starr - [ ur(i, idn), ur(i, idn)*ur(i, imx:imz), enr ])
                b_cc(i, ydim:zdim) = b_ccrf(ydim:zdim) + sr*(b_starr(ydim:zdim) - b_ccr(i, ydim:zdim))

             else ! sm = 0

                ! Average left and right flux if both sm = 0 = B_x

                f(i, :) = half * (fl + sl*(u_starl - [ ul(i, idn), ul(i, idn)*ul(i, imx:imz), enl ]) + &
                     &           fr + sr*(u_starr - [ ur(i, idn), ur(i, idn)*ur(i, imx:imz), enr ]))
                b_cc(i, ydim:zdim) = half*(b_cclf(ydim:zdim) + sl*(b_starl(ydim:zdim) - b_ccl(i, ydim:zdim)) + &
                     &                    b_ccrf(ydim:zdim) + sr*(b_starr(ydim:zdim) - b_ccr(i, ydim:zdim)))

             endif  ! sm = 0

          endif     ! B_x = 0

       endif

    enddo

  end subroutine riemann_hlld

end module hlld
