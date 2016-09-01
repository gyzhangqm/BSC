program expt
  implicit none

  integer :: nodu, nodl, io,ioo, i, j, k, ierr ,il, counterstart, counterend, counter, l, ndp
  integer :: pair(2),lul(3)
  real, dimension(:), allocatable :: xu, yu, xl, yl, bump_pos, xunew, yunew, xlnew, ylnew, dispu, displ
  integer, dimension(:), allocatable :: unodes, lnodes
  real, dimension(:), allocatable :: dpa, dpb
  real (kind = 8) :: three(3), pairreal(2), h, scale, sum, m
  character(256) :: input
  integer, parameter :: t_b = 4
  real(kind = 8), parameter :: pi = acos(-1.)

!-------------------------------------------------------------------------------
!---- Finds out the nodes belonging to the upper and lowers surfaces -----------
  nodu = 0
  nodl = 0

  open(1,file='cylinder.fix.nod')
  do
    read(1,*,iostat=io) pair
    if ( io == 0 ) then
      if (pair(2) == 4) then
        nodu = nodu + 1
      endif
      if (pair(2) == 3) then
        nodl = nodl + 1
      endif
    endif
    if(io==5010) then
    else if (io/=0) then
      exit
    end if
  enddo
  allocate(unodes(nodu))
  allocate(lnodes(nodl))
  rewind(1)
  i = 1
  j = 1
  do
    read(1,*,iostat=io) pair
    if ( io == 0 ) then
      if (pair(2) == 4) then
        unodes(i) = pair(1)
         i = i + 1
       endif
      if (pair(2) == 3) then
        lnodes(j) = pair(1)
        j = j +1
      endif
    endif

    if(io==5010) then

    else if (io/=0) then
      exit
    end if
  enddo

!-------------------------------------------------------------------------------
!Extracts the coordinates of the grid points which lie on the upper and lower surface

allocate(xu(nodu))
allocate(yu(nodu))
allocate(xl(nodl))
allocate(yl(nodl))

counter = 0
i = 1
j = 1
open(1,file='cylinder.geo.dat',status='old',iostat=ierr)
do while (ierr.eq.0)
  read(1,'(A)',iostat=ierr) input
  counter = counter + 1
  il=len_trim(input)
  if ( input(1:il) == 'COORDINATES' ) then
    counterstart = counter + 1
  end if
  if ( input(1:il) == 'END_COORDINATES' ) then
    counterend = counter - 1
  end if
enddo

rewind(1)

do l = 1,counterstart-1
  read(1,*)
end do
do l = counterstart,counterend
  read(1,*)three
  do k = 1,nodu
      if ( unodes(k) == three(1) ) then
        xu(i) = three(2)
        yu(i) = three(3)
        i = i + 1
      end if
    enddo
    do k = 1,nodl
      if ( lnodes(k) == three(1) ) then
        xl(j) = three(2)
        yl(j) = three(3)
        j = j + 1
      end if
    enddo
end do
close(1)
!-------------------------------------------------------------------------------
!--- Finding Number of design parameters for des_var file and reading them -----
ndp = 0
open(1,file='des_var',status='old',iostat=ierr)

  do while (ierr.eq.0)
    read(1,*,iostat=ierr)
    ndp = ndp + 1
  end do
  ndp = (ndp - 1)/2
  allocate(dpa(ndp))
  allocate(dpb(ndp))

rewind(1)
  do i = 1,ndp
    read(1,*,iostat=ierr) pairreal
    dpa(i) = pairreal(2)
  enddo
  do i = 1,ndp
    read(1,*,iostat=ierr) pairreal
    dpb(i) = pairreal(2)
  enddo
close(1)
!-------------------------------------------------------------------------------
!----------- Hicks-Henne Update ------------------------------------------------
scale = 1.
allocate(bump_pos(ndp))
allocate(xunew(nodu))
allocate(yunew(nodu))
allocate(xlnew(nodl))
allocate(ylnew(nodl))
allocate(dispu(nodl))
allocate(displ(nodl))

h = 1./(ndp+1)
do i = 1,ndp
  bump_pos(i) = i*h
enddo
xunew = xu/scale
yunew = yu/scale
xlnew = xl/scale
ylnew = yl/scale

do i = 1,nodu
  sum = 0
  do j = 1,ndp
      m = log(0.5)/log(bump_pos(j))
      sum = sum + dpa(j)*(sin(pi*xu(i)**m)**t_b)
      !gradu(i,j) = sin(pi*xu(i)**m)**t_b;
  enddo
  yunew(i) = yunew(i) + sum
  yunew(i) = yunew(i) * scale
  dispu(i) = yunew(i) - yu(i)
enddo

do i = 1,nodl
  sum = 0
  do j = 1,ndp
      m = log(0.5)/log(bump_pos(j))
      sum = sum + dpb(j)*(sin(pi*xl(i)**m)**t_b)
      !gradu(i,j) = sin(pi*xu(i)**m)**t_b;
  enddo
  ylnew(i) = ylnew(i) + sum
  ylnew(i) = ylnew(i) * scale
  displ(i) = ylnew(i) - yl(i)
enddo
!-------------------------------------------------------------------------------
open(1,file='mesh_disp.txt',status='unknown',iostat=ierr)

do while (ierr==0)
  read(1,iostat=ierr) three
  print *,ierr
  print *,lul
  !if ( three(1) == 20 ) then
    !write(1,*) three(1), h, h
  !end if
enddo


!-------------------------------------------------------------------------------
open(unit = 1,file = "dump/dumptest.txt")
do i=1,ndp
  write(1,*) i, dpa(i)
enddo
do i=1,ndp
  write(1,*) i+ndp, dpb(i)
enddo
close(1)

end program expt
