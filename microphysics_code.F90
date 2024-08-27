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
program microphysics_code
  use mod_constants
  use mod_dimensions
  use mod_micro_common
  use mod_micro_nogtom
  use iso_fortran_env

  implicit none

  integer :: year_start = 1950
  integer :: year_end   = 2024

  character(len=*) , parameter :: e5_base = 'data'
  character(len=*) , parameter :: cmip6_base = 'data'
  character(len=*) , parameter :: macv2_base = 'data/MACV2'
  character(len=*) , parameter :: e5_mask = e5_base//'/fixed/landmask.nc'
  character(len=*) , parameter :: e5_geop = e5_base//'/fixed/geopot.nc'
  character(len=maxpath) :: e5_ps , e5_t , e5_qv , e5_qc , e5_qi
  character(len=maxpath) :: e5_cld , e5_pvv
  character(len=16) :: timepart
  type(mod_2_micro) :: mo2mi
  type(micro_2_mod) :: mi2mo
  integer(ik4) :: rgrid , reduced_points , nlev
  integer(ik4) :: n , k , is , it , ith , ntimes
  integer(ik4) :: n1 , n2
  integer(ik4) :: iyear , imonth , mday
  integer(ik8) :: istep
  real(rkx) , dimension(:) , pointer :: am , bm , ai , bi
  real(rkx) , dimension(:,:) , pointer :: zq
  real(rkx) :: lwc , iwc , tcels , tempc , fsr , desr , aiwc , biwc
  real(rkx) :: kabsi , kabsl , kabs , arg , cldemis
  real(rk8) :: eccen , obliq , mvelp , obliqr , lambm0 , mvelpp , calday
  real(rk8) :: declin , eccf
  real(rk8) :: dt
  integer(ik8) :: tstart , tstop , trate

  call system_clock(count_rate=trate)
  call read_dimensions(e5_geop,rgrid,nlev,reduced_points)

  write(output_unit,*) 'Number of reduced grid points: ', rgrid
  write(output_unit,*) 'Number of latitudes          : ', reduced_points
  write(output_unit,*) 'Number of vertical levels    : ', nlev

  call init_dimensions(rgrid,nlev)

  write(output_unit,*) 'Allocating vertical grids....'
  allocate(am(nlev),bm(nlev))
  allocate(ai(nlev+1),bi(nlev+1))

  write(output_unit,*) 'Allocating types....'

  ! This must be for MPI par....
  n1 = 271000
  n2 = 271080

  dt = 120.0_rkx ! This should work for a 25km....

  call allocate_intf(mo2mi,mi2mo,n1,n2)
  call allocate_nogtom(n1,n2)
  allocate(zq(n1:n2,nlev+1))

  call read_geolocation(e5_mask,mo2mi%xlat,mo2mi%xlon)

  print *, 'Latitude extrema are  :', mo2mi%xlat(n1), mo2mi%xlat(n2)
  print *, 'Longitude extrema are :', mo2mi%xlon(n1), mo2mi%xlon(n2)

  call read_static_2di(e5_mask,'lsm',mo2mi%ldmsk)
  call read_vertical_grid(e5_geop,am,bm,ai,bi)
  call read_static_2d(e5_geop,'z',mo2mi%ht)

  do concurrent ( n = n1:n2 )
    mo2mi%ht(n) = mo2mi%ht(n)/egrav
  end do

  ! Call microphysics code

  call init_nogtom(n1,n2,mo2mi%ldmsk)

  iyear = year_start
  istep = 0

  do iyear = year_start , year_end

    calday = 0.0_rk8

    do imonth = 1 , 12

      write(timepart,'(i0.4,a,i0.2)') iyear, '_', imonth
      e5_ps = e5_base//'/mlevel/lnps_'//trim(timepart)//'.nc'
      e5_t = e5_base//'/mlevel/t_'//trim(timepart)//'.nc'
      e5_qv = e5_base//'/mlevel/q_'//trim(timepart)//'.nc'
      e5_qc = e5_base//'/mlevel/clwc_'//trim(timepart)//'.nc'
      e5_qi = e5_base//'/mlevel/ciwc_'//trim(timepart)//'.nc'
      e5_cld = e5_base//'/mlevel/cc_'//trim(timepart)//'.nc'
      e5_pvv = e5_base//'/mlevel/w_'//trim(timepart)//'.nc'

      ntimes = ntimes_in_file(e5_ps)

      do it = 1 , ntimes

        mday = (it-1)/4

        print *, '########################################################'
        print *, iyear , imonth, mday+1, calday

        ! Time in hourly file consistent with 6h in mlev
        ith = (it-1) * 6 + 1

        call read_field_2d(e5_ps,'lnsp',mo2mi%ps,1,it)
        call read_field_3d(e5_t,'t',mo2mi%t,it)
        call read_field_3d(e5_qv,'q',mo2mi%qxx(:,:,iqqv),it)
        call read_field_3d(e5_qc,'clwc',mo2mi%qxx(:,:,iqql),it)
        call read_field_3d(e5_qi,'ciwc',mo2mi%qxx(:,:,iqqi),it)
        call read_field_3d(e5_cld,'cc',mo2mi%cldf,it)
        call read_field_3d(e5_pvv,'w',mo2mi%pverv,it)

        do concurrent ( n = n1:n2 )
          mo2mi%ps(n) = exp(mo2mi%ps(n))
        end do
        do is = 1 , nqx
          do k = 1 , nlev
            do n = n1 , n2
              mo2mi%qxx(n,k,is) = mo2mi%qxx(n,k,is)/(1.0_rkx-mo2mi%qxx(n,k,is))
            end do
          end do
        end do
        do k = 1 , nlev
          do n = n1 , n2
            mo2mi%phs(n,k) = am(k)+ bm(k)*mo2mi%ps(n)
            mo2mi%rho(n,k) = (mo2mi%phs(n,k))/(mo2mi%t(n,k)*rgas)
            mo2mi%qs(n,k) = pfwsat(mo2mi%t(n,k),mo2mi%phs(n,k))
          end do
        end do
        do n = n1 , n2
          mo2mi%pfs(n,1) = 1.0e-10_rkx
        end do
        do k = 2 , nlev+1
          do n = n1 , n2
            mo2mi%pfs(n,k) = ai(k)+ bi(k)*mo2mi%ps(n)
          end do
        end do

        do n = n1 , n2
          zq(n,nlev+1) = mo2mi%ht(n)
        end do
        do k = nlev, 1 , -1
          do n = n1 , n2
            zq(n,k) = zq(n,k) + rovg*mo2mi%t(n,k)* &
              log(mo2mi%pfs(n,k+1)/mo2mi%pfs(n,k))
          end do
        end do
        do k = 1 , nlev
          do n = n1 , n2
            mo2mi%delz(n,k) =  zq(n,k)- zq(n,k+1)
            ! Set for now to zero
            mo2mi%heatrt(n,k) = 0.0_rkx
            mo2mi%qdetr(n,k) = 0.0_rkx
            mo2mi%pverv(n,k) = 0.0_rkx
          end do
        end do

        mi2mo%rainnc = 0.0_rkx
        mi2mo%snownc = 0.0_rkx
        mi2mo%trrate = 0.0_rkx
        mi2mo%lsmrnc = 0.0_rkx ! This cumulates, but we don't need here it to
        mi2mo%tten = 0.0_rkx
        mi2mo%qxten = 0.0_rkx

        call system_clock(tstart)
        call nogtom(n1,n2,dt,mo2mi,mi2mo)
        call system_clock(tstop)
        print *, '########################################################'
        print *, 'Elapsed time : ',real(tstop-tstart,rk8)/real(trate,rk8), &
          ' s'
        print *, '########################################################'
        print *, 'Large scale rain (max)       : ', maxval(mi2mo%rainnc)
        print *, 'Large scale snow (max)       : ', maxval(mi2mo%snownc)
        print *, 'Large scale T tendency (max) : ', maxval(mi2mo%tten)
        print *, 'Large scale T tendency (min) : ', minval(mi2mo%tten)
        print *, 'Large scale QV tendency (max): ', &
          maxval(mi2mo%qxten(:,:,iqqv))
        print *, 'Large scale QV tendency (min): ', &
          minval(mi2mo%qxten(:,:,iqqv))
        print *, 'Large scale QC tendency (max): ', &
          maxval(mi2mo%qxten(:,:,iqql))
        print *, 'Large scale QC tendency (min): ', &
          minval(mi2mo%qxten(:,:,iqql))
        print *, '########################################################'

        calday = calday + 1.0_rkx/4.0_rkx
        istep = istep + 1
      end do
    end do
  end do

  call deallocate_intf(mo2mi,mi2mo)
  call deallocate_nogtom( )
  deallocate(am,bm)
  deallocate(ai,bi)
  deallocate(zq)

  contains

  integer function ntimes_in_file(fname) result(nt)
    use netcdf
    character(len=*) , intent(in) :: fname
    integer :: ncid , ncstat , idimid
    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_dimid(ncid, 'time', idimid)
    if ( ncstat /= nf90_noerr ) stop 'No time dimension in '//trim(fname)
    ncstat = nf90_inquire_dimension(ncid, idimid, len=nt)
    if ( ncstat /= nf90_noerr ) stop 'rgrid dim time error in '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end function ntimes_in_file

  subroutine read_geolocation(fname,lat,lon)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname
    real(rkx) , dimension(n1:n2) , intent(out) :: lat , lon
    real(rk8) , dimension(reduced_points) :: rlat
    real(rk8) , dimension(nr) :: flat , flon
    integer , dimension(reduced_points) :: nlons
    real(rkx) , dimension(reduced_points) :: dlons
    integer :: ncid , ncvar , ncstat
    integer :: i , j , n , istart

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'reduced_points',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable reduced_points in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,nlons)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot read variable reduced_points in file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'lat',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable lat in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,rlat)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot read variable lat in file '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
    do i = 1 , reduced_points
      dlons(i) = 360_rkx / real(nlons(i),rkx)
    end do
    istart = 1
    do i = 1 , reduced_points
      do j = istart , istart+nlons(i)-1
        flat(j) = rlat(i)
        !flon(j) = 0.5_rkx*dlons(i) + (j-istart) * dlons(i)
        ! According to EMOS, 0 is always present (WHY, THOUGH?)
        flon(j) = (j-istart) * dlons(i)
      end do
      istart = istart + nlons(i)
    end do
    do n = n1 , n2
      lat(n) = flat(n)
      lon(n) = flon(n)
    end do
  end subroutine read_geolocation

  subroutine read_vertical_grid(fname,am,bm,ai,bi)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname
    real(rkx) , dimension(nlev) :: am , bm
    real(rkx) , dimension(nlev+1) :: ai , bi
    integer :: ncid , ncvar , ncstat
    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'hyam',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable hyam in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,am)
    if ( ncstat /= nf90_noerr ) stop 'error read var hyam in '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'hybm',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable hybm in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,bm)
    if ( ncstat /= nf90_noerr ) stop 'error read var hybm in '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'hyai',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable hyai in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,ai)
    if ( ncstat /= nf90_noerr ) stop 'error read var hyai in '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'hybi',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable hybi in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,bi)
    if ( ncstat /= nf90_noerr ) stop 'error read var hybi in '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_vertical_grid

  subroutine read_static_2di(fname,vname,f)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname , vname
    integer(ik4) , dimension(n1:n2) , intent(out) :: f
    integer :: ncid , ncvar , ncstat
    integer , dimension(2) :: istart , icount

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,vname,ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable '//trim(vname)//' in file '//trim(fname)
    istart(1) = n1
    istart(2) = 1
    icount(1) = n2-n1+1
    icount(2) = 1
    ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_static_2di

  subroutine read_static_2d(fname,vname,f)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname , vname
    real(rkx) , dimension(n1:n2) , intent(out) :: f
    integer :: ncid , ncvar , ncstat
    integer , dimension(2) :: istart , icount

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,vname,ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable '//trim(vname)//' in file '//trim(fname)
    istart(1) = n1
    istart(2) = 1
    icount(1) = n2-n1+1
    icount(2) = 1
    ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_static_2d

  subroutine read_field_2d(fname,vname,f,il,it)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname , vname
    real(rkx) , dimension(n1:n2) , intent(out) :: f
    integer(ik4) , intent(in) :: il, it
    integer :: ncid , ncvar , ncstat

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,vname,ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable '//trim(vname)//' in file '//trim(fname)
    if ( il == 0 ) then
      notime : block
        integer , dimension(2) :: istart , icount
        istart(1) = n1
        istart(2) = it
        icount(1) = n2-n1+1
        icount(2) = 1
        ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
        if ( ncstat /= nf90_noerr ) &
          stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
      end block notime
    else
      hastime : block
        integer , dimension(3) :: istart , icount
        istart(1) = n1
        istart(2) = il
        istart(3) = it
        icount(1) = n2-n1+1
        icount(2) = 1
        icount(3) = 1
        ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
        if ( ncstat /= nf90_noerr ) &
          stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
      end block hastime
    end if
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_field_2d

  subroutine read_field_3d(fname,vname,f,it)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname , vname
    real(rkx) , dimension(n1:n2,nlev) , intent(out) :: f
    integer(ik4) , intent(in) :: it
    integer :: ncid , ncvar , ncstat
    integer :: i , k

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,vname,ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable '//trim(vname)//' in file '//trim(fname)
    if ( it == 0 ) then
      notime : block
        integer , dimension(3) :: istart , icount
        istart(1) = n1
        istart(2) = 1
        istart(3) = 1
        icount(1) = n2-n1+1
        icount(2) = nlev
        icount(3) = 1
        ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
        if ( ncstat /= nf90_noerr ) &
          stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
      end block notime
    else
      hastime : block
        integer , dimension(3) :: istart , icount
        istart(1) = n1
        istart(2) = 1
        istart(3) = it
        icount(1) = n2-n1+1
        icount(2) = nlev
        icount(3) = 1
        ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
        if ( ncstat /= nf90_noerr ) &
          stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
      end block hastime
    end if
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_field_3d

  subroutine read_dimensions(fname,rgrid,nlev,reduced_points)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname
    integer(ik4) , intent(out) :: rgrid , nlev , reduced_points
    integer :: ncid , ncstat , idimid

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_dimid(ncid, 'rgrid', idimid)
    if ( ncstat /= nf90_noerr ) stop 'No rgrid dimension in '//trim(fname)
    ncstat = nf90_inquire_dimension(ncid, idimid, len=rgrid)
    if ( ncstat /= nf90_noerr ) stop 'rgrid dim read error in '//trim(fname)
    ncstat = nf90_inq_dimid(ncid, 'nhym', idimid)
    if ( ncstat /= nf90_noerr ) stop 'No nhym dimension in '//trim(fname)
    ncstat = nf90_inquire_dimension(ncid, idimid, len=nlev)
    if ( ncstat /= nf90_noerr ) stop 'nhym dim read error in '//trim(fname)
    ncstat = nf90_inq_dimid(ncid, 'reduced_points', idimid)
    if ( ncstat /= nf90_noerr ) stop 'No lat dimension in '//trim(fname)
    ncstat = nf90_inquire_dimension(ncid, idimid, len=reduced_points)
    if ( ncstat /= nf90_noerr ) stop 'lat dim read error in '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_dimensions

  ! Computes saturation pressurre
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
!DIR$ ATTRIBUTES FORCEINLINE :: pfesat
  pure elemental real(rkx) function pfesat(t,p) result(es)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t , p ! Temperature (K) , Pressure (Pa)

    real(rk8) :: td , t_limit , esat
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8) , parameter :: a0 =  0.611213476e+03_rk8
    real(rk8) , parameter :: a1 =  0.444007856e+02_rk8
    real(rk8) , parameter :: a2 =  0.143064234e+01_rk8
    real(rk8) , parameter :: a3 =  0.264461437e-01_rk8
    real(rk8) , parameter :: a4 =  0.305903558e-03_rk8
    real(rk8) , parameter :: a5 =  0.196237241e-05_rk8
    real(rk8) , parameter :: a6 =  0.892344772e-08_rk8
    real(rk8) , parameter :: a7 = -0.373208410e-10_rk8
    real(rk8) , parameter :: a8 =  0.209339997e-13_rk8
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8) , parameter :: c0 =  0.611123516e+03_rk8
    real(rk8) , parameter :: c1 =  0.503109514e+02_rk8
    real(rk8) , parameter :: c2 =  0.188369801e+01_rk8
    real(rk8) , parameter :: c3 =  0.420547422e-01_rk8
    real(rk8) , parameter :: c4 =  0.614396778e-03_rk8
    real(rk8) , parameter :: c5 =  0.602780717e-05_rk8
    real(rk8) , parameter :: c6 =  0.387940929e-07_rk8
    real(rk8) , parameter :: c7 =  0.149436277e-09_rk8
    real(rk8) , parameter :: c8 =  0.262655803e-12_rk8

    t_limit = t - tzero
    if ( t_limit > 100.0_rk8 ) t_limit = 100.0_rk8
    if ( t_limit < -75.0_rk8 ) t_limit = -75.0_rk8
    td = t_limit
    if ( td >= 0.0_rk8 ) then
      esat = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    else
      esat = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    end if
    es = real(min(esat,0.15_rk8*p),rkx) ! pa
  end function pfesat

!DIR$ ATTRIBUTES FORCEINLINE :: pfwsat
  pure elemental real(rkx) function pfwsat(t,p,e) result(ws) ! In kg/kg
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) :: es ! , qs
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t,p)
    end if
    ws = ep2 * (es / (p - es))
    ! Bolton 1980
    ! qs = ep2 * es / (p - (0.378_rkx * es))
    ! ws = qs * ( d_one - qs)
  end function pfwsat

end program microphysics_code

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
