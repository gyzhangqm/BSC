program expt
  use iso_fortran_env, only : IOSTAT_EOR
  implicit none

  integer :: nodu, nodl, io,ioo, i, j, k, flag1, flag2, ierr ,il, counterstart, counterend, counter, l
  integer ::pair(2)
  real, dimension(:), allocatable :: xu, yu, xl, yl
  integer, dimension(:), allocatable :: unodes, lnodes
  real (kind = 8) :: three(3)
!  character(len = :), allocatable :: input, output
!  character(len = 1) :: buffer
!  character(len = 1) :: input
  character(256) :: input


!-------------------------------------------------------------------------------
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

allocate(xu(nodu))
allocate(yu(nodu))
allocate(xl(nodl))
allocate(yl(nodl))

counter = 0
i = 1
j = 1
flag1 = 1
flag2 = 0
open(1,file='cylinder.geo.dat',status='old',iostat=ierr)

do while (ierr.eq.0)
  read(1,'(A)',iostat=ierr) input
  counter = counter + 1
  il=len_trim(input)
  !print *,input(1:il)
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

print *,counter
print *,counterstart
print *,counterend
open(unit = 1,file = "dump/dumptest.txt")
do i=1,nodu
      write(1,*) unodes(i), xu(i), yu(i)
enddo

end program expt
