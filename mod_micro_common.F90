!
! Copyright (c) 2024 Graziano Giuliani from UNESCO ICTP
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
module mod_micro_common

  use mod_constants
  use mod_dimensions

  implicit none

  private

  logical , public :: budget_compute = .false.
  integer(ik4) , public :: nssopt = 1
  integer(ik4) , public :: iautoconv = 1
  real(rkx) , public :: vfqr = 4.0_rkx
  real(rkx) , public :: vfqi = 0.15_rkx
  real(rkx) , public :: vfqs = 1.0_rkx
  real(rkx) , public :: auto_rate_khair = 0.355_rkx
  real(rkx) , public :: auto_rate_kessl = 1.e-3_rkx
  real(rkx) , public :: auto_rate_klepi = 0.5e-3_rkx
  real(rkx) , public :: rkconv = 1.666e-4_rkx ! 1.0/6000.0
  real(rkx) , public :: skconv = 1.0e-3_rkx
  real(rkx) , public :: rcldiff = 1.0e-6_rkx
  real(rkx) , public :: rcovpmin = 0.1_rkx
  real(rkx) , public :: rpecons = 5.547e-5_rkx

  type mod_2_micro
    real(rkx) , pointer , dimension(:) :: xlat
    real(rkx) , pointer , dimension(:) :: xlon
    real(rkx) , pointer , dimension(:) :: ps
    real(rkx) , pointer , dimension(:) :: ht
    real(rkx) , pointer , dimension(:,:) :: phs
    real(rkx) , pointer , dimension(:,:) :: pfs
    real(rkx) , pointer , dimension(:,:) :: delz
    real(rkx) , pointer , dimension(:,:) :: t
    real(rkx) , pointer , dimension(:,:) :: rho
    real(rkx) , pointer , dimension(:,:) :: pverv
    real(rkx) , pointer , dimension(:,:,:) :: qxx
    real(rkx) , pointer , dimension(:,:) :: qs
    real(rkx) , pointer , dimension(:,:) :: heatrt
    real(rkx) , pointer , dimension(:,:) :: qdetr
    real(rkx) , pointer , dimension(:,:) :: cldf
    integer(ik4) , pointer , dimension(:) :: ldmsk
  end type mod_2_micro

  type micro_2_mod
    real(rkx) , pointer , dimension(:) :: trrate
    real(rkx) , pointer , dimension(:) :: rainnc
    real(rkx) , pointer , dimension(:) :: lsmrnc
    real(rkx) , pointer , dimension(:) :: snownc
    real(rkx) , pointer , dimension(:,:) :: tten
    real(rkx) , pointer , dimension(:,:,:) :: qxten
  end type micro_2_mod

  public :: mod_2_micro , micro_2_mod , allocate_intf , deallocate_intf

  contains

  subroutine allocate_intf(mo2mi,mi2mo,ici1,ici2)
    implicit none
    integer(ik4) , intent(in) :: ici1 , ici2
    type(mod_2_micro) , intent(inout) :: mo2mi
    type(micro_2_mod) , intent(inout) :: mi2mo

    allocate(mo2mi%xlat(ici1:ici2))
    allocate(mo2mi%xlon(ici1:ici2))
    allocate(mo2mi%ps(ici1:ici2))
    allocate(mo2mi%ht(ici1:ici2))
    allocate(mo2mi%ldmsk(ici1:ici2))
    allocate(mo2mi%phs(ici1:ici2,kz))
    allocate(mo2mi%pfs(ici1:ici2,kzp1))
    allocate(mo2mi%delz(ici1:ici2,kz))
    allocate(mo2mi%t(ici1:ici2,kz))
    allocate(mo2mi%rho(ici1:ici2,kz))
    allocate(mo2mi%pverv(ici1:ici2,kz))
    allocate(mo2mi%qxx(ici1:ici2,kz,nqx))
    allocate(mo2mi%qs(ici1:ici2,kz))
    allocate(mo2mi%heatrt(ici1:ici2,kz))
    allocate(mo2mi%qdetr(ici1:ici2,kz))
    allocate(mo2mi%cldf(ici1:ici2,kz))

    allocate(mi2mo%trrate(ici1:ici2))
    allocate(mi2mo%rainnc(ici1:ici2))
    allocate(mi2mo%lsmrnc(ici1:ici2))
    allocate(mi2mo%snownc(ici1:ici2))
    allocate(mi2mo%tten(ici1:ici2,kz))
    allocate(mi2mo%qxten(ici1:ici2,kz,nqx))

  end subroutine allocate_intf

  subroutine deallocate_intf(mo2mi,mi2mo)
    implicit none
    type(mod_2_micro) , intent(inout) :: mo2mi
    type(micro_2_mod) , intent(inout) :: mi2mo
    deallocate(mo2mi%xlat)
    deallocate(mo2mi%xlon)
    deallocate(mo2mi%ps)
    deallocate(mo2mi%ht)
    deallocate(mo2mi%ldmsk)
    deallocate(mo2mi%phs)
    deallocate(mo2mi%pfs)
    deallocate(mo2mi%delz)
    deallocate(mo2mi%t)
    deallocate(mo2mi%rho)
    deallocate(mo2mi%pverv)
    deallocate(mo2mi%qxx)
    deallocate(mo2mi%qs)
    deallocate(mo2mi%heatrt)
    deallocate(mo2mi%qdetr)
    deallocate(mo2mi%cldf)

    deallocate(mi2mo%trrate)
    deallocate(mi2mo%rainnc)
    deallocate(mi2mo%lsmrnc)
    deallocate(mi2mo%snownc)
    deallocate(mi2mo%tten)
    deallocate(mi2mo%qxten)
  end subroutine deallocate_intf

end module mod_micro_common

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
