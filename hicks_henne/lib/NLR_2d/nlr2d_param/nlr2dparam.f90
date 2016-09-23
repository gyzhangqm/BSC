program param
  implicit none

  type intreal
    integer,dimension(1)   :: ints
    real,dimension(2)      :: floats
  endtype intreal

  type(intreal) :: mesh_disp, geoinput, geoinput1, geoinput2, temp
  type(intreal), dimension(:), allocatable :: panelgeoinput


  integer :: nodu, nodl, nodf, i, j, k, p, ierr ,il, counterstart, counterend, counter, l, ndp, totnodes
  integer :: pair(2), flag, pos
  real, dimension(:), allocatable :: xu, yu, xl, yl, xf, yf, bump_pos, xunew, yunew, xlnew, ylnew, xfnew, yfnew
  integer, dimension(:), allocatable :: unodes, lnodes, fnodes
  real, dimension(:), allocatable :: dpa, dpb
  real (kind = 8) :: three(3), pairreal(2), h, scaleu, scalel, sum, m, xutrans, xltrans, xflap, yflap, flapangle
  real (kind = 8) :: xref, yref
  integer, parameter :: t_b = 4
  real(kind = 8), parameter :: pi = acos(-1.)
  character(256) :: input, cmd, geo, fixnod, code
  integer :: foo, flagupper, flaglower, flagflap

  call get_command_argument(1, cmd)
  geo = trim(cmd)//".geo.dat"
  fixnod = trim(cmd)//".fix.nod"
  code = trim(cmd)//".codes"
!-------------------------------------------------------------------------------
!-------- Reading the flags for upper and lower surfaces from .codes file ------
  open(1,file = code)
  read(1,*) input
  read(1,'(a)') input
  input = trim(input)
  pos = index(input,":")
  read(input(pos+1:), *) flagupper
  read(1,'(a)') input
  input = trim(input)
  pos = index(input,":")
  read(input(pos+1:), *) flaglower
  read(1,'(a)') input
  input = trim(input)
  pos = index(input,":")
  read(input(pos+1:), *) flagflap
  close(1)

!-------------------------------------------------------------------------------
!---- Finds out the nodes belonging to the upper and lowers surfaces -----------
  nodu = 0
  nodl = 0
  nodf = 0

  open(1,file=fixnod,status='old',iostat=ierr)
  do
    read(1,*,iostat=ierr) pair
    if ( ierr == 0 ) then
      if (pair(2) == flagupper) then
        nodu = nodu + 1
      endif
      if (pair(2) == flaglower) then
        nodl = nodl + 1
      endif
      if (pair(2) == flagflap) then
        nodf = nodf + 1
      endif
    endif
    if(ierr==5010) then
    else if (ierr/=0) then
      exit
    end if
  enddo
  allocate(unodes(nodu))
  allocate(lnodes(nodl))
  allocate(fnodes(nodf))
  rewind(1)
  i = 1
  j = 1
  k = 1
  do
    read(1,*,iostat=ierr) pair
    if ( ierr == 0 ) then
      if (pair(2) == flagupper) then
        unodes(i) = pair(1)
         i = i + 1
       endif
      if (pair(2) == flaglower) then
        lnodes(j) = pair(1)
        j = j + 1
      endif
      if (pair(2) == flagflap) then
        fnodes(k) = pair(1)
         k = k + 1
       endif
    endif

    if(ierr==5010) then

    else if (ierr/=0) then
      exit
    end if
  enddo

!-------------------------------------------------------------------------------
!Extracts the coordinates of the grid points which lie on the upper and lower surface

counter = 0
open(2,file='dumpallu.txt')
open(3,file='dumpalll.txt')
open(4,file='dumpflap.txt')

open(1,file=geo,status='old',iostat=ierr)
do while (ierr.eq.0)
  read(1,'(A)',iostat=ierr) input
  counter = counter + 1
  il=len_trim(input)
  if ( input(1:il) == '  COORDINATES' ) then
    counterstart = counter + 1
  end if
  if ( input(1:il) == '  END_COORDINATES' ) then
    counterend = counter - 1
  end if
enddo

rewind(1)

do l = 1,counterstart-1
  read(1,*)
end do
do l = counterstart,counterend
  read(1,*) geoinput
  do k = 1,nodu
      if ( unodes(k) == geoinput%ints(1) ) then
        write(2,*) geoinput
      end if
    enddo
    do k = 1,nodl
      if ( lnodes(k) == geoinput%ints(1) ) then
        write(3,*) geoinput
      end if
    enddo
    do k = 1,nodf
      if ( fnodes(k) == geoinput%ints(1) ) then
        write(4,*) geoinput
      end if
    enddo
end do
close(1)
close(2)
close(3)
close(4)
totnodes = counterend - counterstart +1


!-------------------------------------------------------------------------------
!---------- Sorting data wrt x -------------------------------------------------

allocate(panelgeoinput(nodu))
open(1,file='dumpallu.txt')
do k = 1,nodu
  read(1,*) panelgeoinput(k)
enddo
close(1)
do i = 1,nodu
  do j = 1,nodu
    if ( i == j ) then
    else
      if ( panelgeoinput(i)%floats(1) .lt. panelgeoinput(j)%floats(1) ) then
        temp = panelgeoinput(i)
        panelgeoinput(i) = panelgeoinput(j)
        panelgeoinput(j) = temp
      end if
    end if
  enddo
enddo
open(1,file='dumpallu.txt')
do k = 1,nodu
  write(1,*) panelgeoinput(k)
enddo
deallocate(panelgeoinput)

allocate(panelgeoinput(nodl))
open(1,file='dumpalll.txt')
do k = 1,nodl
  read(1,*) panelgeoinput(k)
enddo
close(1)
do i = 1,nodl
  do j = 1,nodl
    if ( i == j ) then
    else
      if ( panelgeoinput(i)%floats(1) .gt. panelgeoinput(j)%floats(1) ) then
        temp = panelgeoinput(i)
        panelgeoinput(i) = panelgeoinput(j)
        panelgeoinput(j) = temp
      end if
    end if
  enddo
enddo
open(1,file='dumpalll.txt')
do k = 1,nodl
  write(1,*) panelgeoinput(k)
enddo
deallocate(panelgeoinput)


!-------------------------------------------------------------------------------
!--- Finding Number of design parameters for des_var file and reading them -----
ndp = 0
open(1,file='des_vars.dat',status='old',iostat=ierr)

  do
    read(1,*,iostat=ierr)
    if ( ierr .ne. 0 ) then
      exit
    end if
    ndp = ndp + 1
  end do
  ndp = (ndp - 3)/2
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
  read(1,*) pairreal
  xflap = pairreal(2)
  read(1,*) pairreal
  yflap = pairreal(2)
  read(1,*) pairreal
  flapangle = pairreal(2)
close(1)

!-------------------------------------------------------------------------------
!----------- Hicks-Henne Update ------------------------------------------------

allocate(xu(nodu))
allocate(yu(nodu))
allocate(xl(nodl))
allocate(yl(nodl))

open(1,file='dumpallu.txt')
do k = 1,nodu
  read(1,*) geoinput
  xu(k) = geoinput%floats(1)
  yu(k) = geoinput%floats(2)
enddo
close(1)

open(1,file='dumpalll.txt')
do k = 1,nodl
  read(1,*) geoinput
  xl(k) = geoinput%floats(1)
  yl(k) = geoinput%floats(2)
enddo
close(1)


scaleu = abs(maxval(xu) - minval(xu))
scalel = abs(maxval(xl) - minval(xl))

allocate(bump_pos(ndp))
allocate(xunew(nodu))
allocate(yunew(nodu))
allocate(xlnew(nodl))
allocate(ylnew(nodl))


h = 1./(ndp+1)
do i = 1,ndp
  bump_pos(i) = i*h
enddo
xunew = xu/scaleu
yunew = yu

xlnew = xl/scalel
ylnew = yl
xutrans = minval(xunew)
xltrans = minval(xlnew)
xunew = xunew - xutrans
xlnew = xlnew - xltrans


do i = 1,nodu
  sum = 0.
  do j = 1,ndp
      m = log(0.5)/log(bump_pos(j))
      sum = sum + dpa(j)*(sin(pi*xunew(i)**m)**t_b)
  enddo
  yunew(i) = yunew(i) + sum
  xunew(i) = xunew(i) + xutrans
  xunew(i) = xunew(i) * scaleu
enddo

do i = 1,nodl
  sum = 0.
  do j = 1,ndp
      m = log(0.5)/log(bump_pos(j))
      sum = sum + dpb(j)*(sin(pi*xlnew(i)**m)**t_b)
  enddo
  ylnew(i) = ylnew(i) + sum
  xlnew(i) = xlnew(i) + xltrans
  xlnew(i) = xlnew(i) * scalel
enddo

open(1,file='dumpupdatedu.txt')
open(2,file='dumpallu.txt')
do k = 1,nodu
  read(2,*) geoinput
  write(1,*) geoinput%ints(1), xunew(k), yunew(k)
enddo
close(1)
close(2)

open(1,file='dumpupdatedl.txt')
open(2,file='dumpalll.txt')
do k = 1,nodl
  read(2,*) geoinput
  write(1,*) geoinput%ints(1), xlnew(k), ylnew(k)
enddo
close(1)
close(2)

!-------------------------------------------------------------------------------
!---------------------- Flap update --------------------------------------------
allocate(xf(nodf))
allocate(yf(nodf))
allocate(xfnew(nodf))
allocate(yfnew(nodf))

open(1,file='dumpflap.txt')
do k = 1,nodf
  read(1,*) geoinput
  xf(k) = geoinput%floats(1)
  yf(k) = geoinput%floats(2)
enddo
close(1)
xref = maxval(xf)
open(1,file='dumpflap.txt')
do k = 1,nodf
  read(1,*) geoinput
  if ( xref == geoinput%floats(1) ) then
    yref = geoinput%floats(2)
  end if
enddo
close(1)

do k = 1,nodf
  xfnew(k) = ((xf(k) - xref)*cos(flapangle)) - ((yf(k) - yref)*sin(flapangle)) + xref + xflap
  yfnew(k) = ((xf(k) - xref)*sin(flapangle)) + ((yf(k) - yref)*cos(flapangle)) + yref + yflap
enddo

open(1,file='dumpupdatedflap.txt')
open(2,file='dumpflap.txt')
do k = 1,nodf
  read(2,*) geoinput
  write(1,*) geoinput%ints(1), xfnew(k), yfnew(k)
enddo

!-------------------------------------------------------------------------------
!-------------------- Concatenation --------------------------------------------

open(1,file='baseline.txt')
open(2,file='dumpallu.txt')
open(3,file='dumpalll.txt')
open(4,file='dumpflap.txt')
do k = 1,nodu
  read(2,*) geoinput
  write(1,*) geoinput
enddo
do k = 1,nodl
  read(3,*) geoinput
  write(1,*) geoinput
enddo
do k = 1,nodf
  read(4,*) geoinput
  write(1,*) geoinput
enddo
close(1)
close(2)
close(3)
close(4)

open(1,file='new.txt')
open(2,file='dumpupdatedu.txt')
open(3,file='dumpupdatedl.txt')
open(4,file='dumpupdatedflap.txt')
do k = 1,nodu
  read(2,*) geoinput
  write(1,*) geoinput
enddo
do k = 1,nodl
  read(3,*) geoinput
  write(1,*) geoinput
enddo
do k = 1,nodf
  read(4,*) geoinput
  write(1,*) geoinput
enddo
close(1)
close(2)
close(3)
close(4)

!-------------------------------------------------------------------------------
!------- Writing into mesh_disp.txt file ---------------------------------------

open(unit = 1, file = 'mesh_disp.dat', status = 'unknown', iostat = ierr)
open(2,file = 'baseline.txt')
open(3,file = 'new.txt')
do k = 1,nodu+nodl+nodf
  read(2,*) geoinput1
  read(3,*) geoinput2
  write(1,*) geoinput1%ints(1), geoinput2%floats(1)-geoinput1%floats(1), geoinput2%floats(2)-geoinput1%floats(2)
enddo
close(1)
close(2)
close(3)

call system('rm dumpallu.txt')
call system('rm dumpalll.txt')
call system('rm dumpflap.txt')
call system('rm dumpupdatedu.txt')
call system('rm dumpupdatedl.txt')
call system('rm dumpupdatedflap.txt')
!call system('rm baseline.txt')
!call system('rm new.txt')


end program param
