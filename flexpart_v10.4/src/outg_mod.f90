! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module outg_mod

  use par_mod, only: dep_prec, sp

  implicit none

  real,allocatable, dimension (:) :: outheight
  real,allocatable, dimension (:) :: outheighthalf
  real,allocatable, dimension (:,:) :: oroout
  real,allocatable, dimension (:,:) :: orooutn
  real,allocatable, dimension (:,:) :: area
  real,allocatable, dimension (:,:) :: arean
  real,allocatable, dimension (:,:,:) :: volume
  real,allocatable, dimension (:,:,:) :: volumen
  real,allocatable, dimension (:,:,:) :: areaeast
  real,allocatable, dimension (:,:,:) :: areanorth
  real,allocatable, dimension (:,:,:) :: densityoutgrid
  real,allocatable, dimension (:,:,:) :: densitydrygrid ! added RLT 
  real,allocatable, dimension (:,:,:) :: factor_drygrid ! added RLT 
  real,allocatable, dimension (:,:,:) :: factor3d
  real,allocatable, dimension (:,:,:) :: grid
  real(sp),allocatable, dimension (:,:) :: wetgrid
  real(sp),allocatable, dimension (:,:) :: drygrid
  real,allocatable, dimension (:,:,:) :: gridsigma
  real(dep_prec),allocatable, dimension (:,:) :: drygridsigma
  real(dep_prec),allocatable, dimension (:,:) :: wetgridsigma
  real,allocatable, dimension (:) :: sparse_dump_r
  real,allocatable, dimension (:) :: sparse_dump_u
  integer,allocatable, dimension (:) :: sparse_dump_i

  ! To print concentration within mixing height, added by ZW
  ! -------------------- start --------------------
  real,allocatable, dimension (:,:,:) :: conc_3d_byz, conc_3d_byzn
  real,allocatable, dimension (:,:) :: conc_2d_lasttime, conc_2d_lasttimen, conc_2d_lasttime_allrelease
  real,allocatable, dimension (:,:) :: conc_2d, conc_2dn, hmix_2d, hmix_2dn

  ! read in external hmix
  real,allocatable, dimension (:,:,:) :: data_in
  real,allocatable, dimension (:,:) :: hmix_readin_2d, data_in_onehour
  real,allocatable, dimension (:)  :: time_in
  integer :: nlon_in, nlat_in, ntime_in
  integer ::  read_month = 0 ! identify which month of hmix data has been read
  integer ::  index_within_month = 0 ! identify hour index of this month for hmix data  
  ! -------------------- end --------------------
end module outg_mod
