program param
  implicit none

  type intreal
    integer,dimension(1)   :: ints
    real,dimension(3)      :: floats
  endtype intreal
  type(intreal) :: geoinput, temp, geoinput1, geoinput2
  type(intreal), dimension(:), allocatable :: panelgeoinput, panelgeoinputu, panelgeoinputl, uleft, uright, lleft, lright

  integer, parameter :: t_b = 4
  real(kind = 8), parameter :: pi = acos(-1.), gamma = 10.
  character(256) :: input, cmd, geo, fixnod, code, adjoint, gradtitle
  integer :: flagupper, flaglower, pos, nodu, nodl,  i, j, k, l, ierr, counterstart, counterend, counter, totnodes
  integer :: nx, nl, totalpanels, pair(2), panels, ndp, leftbound, rightbound, flag, dummy
  integer, dimension(:), allocatable :: unodes, lnodes, panelflag, masterpanelno, derivativenodenou, derivativenodenol
  real, dimension(:), allocatable :: zpos, zpanelpos, dpa(:,:), dpb(:,:), areaold, areanew, x, y, derivative, lambda(:,:)
  real, dimension(:), allocatable :: xu, yu, xl, yl, bump_pos, xunew, yunew, xlnew, ylnew, areagrad(:,:), graddata(:,:)
  real, dimension(:), allocatable :: volumegrad, graddataleft(:,:), graddataright(:,:)
  real(kind = 8) :: buffer, pairreal(2), h, scaleu, scalel, sum, m, xutrans, xltrans, firstterm, secondterm, area
  real(kind = 8) :: volumeold, volumenew

  type data
    integer, dimension(:), allocatable :: nodenos
    real, dimension(:), allocatable :: x,y,z
  endtype data
  type(data) :: datau, datal


  call get_command_argument(1, cmd)
  geo = trim(cmd)//".geo.dat"
  fixnod = trim(cmd)//".fix.nod"
  code = trim(cmd)//".codes"
  adjoint = trim(cmd)//".ensi.DISPM-000001"
  gradtitle = trim(cmd)//"_grad.dat"

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
  close(1)

!-------------------------------------------------------------------------------
!---- Finds out the nodes belonging to the upper and lowers surfaces -----------
  nodu = 0
  nodl = 0

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
    endif
    if(ierr==5010) then
    else if (ierr/=0) then
      exit
    end if
  enddo
  allocate(unodes(nodu))
  allocate(lnodes(nodl))
  rewind(1)
  i = 1
  j = 1
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
    endif

    if(ierr==5010) then

    else if (ierr/=0) then
      exit
    end if
  enddo

!-------------------------------------------------------------------------------
!Extracts the coordinates of the grid points which lie on the upper and lower surface

counter = 0
i = 1
j = 1

open(2,file = 'dumpallu.txt')
open(3,file = 'dumpalll.txt')
open(1,file=geo,status='old',iostat=ierr)
do while (ierr.eq.0)
  read(1,'(A)',iostat=ierr) input
  counter = counter + 1
  if ( trim(input) == ' COORDINATES' ) then
    counterstart = counter + 1
  end if
  if ( trim(input) == ' END_COORDINATES' ) then
    counterend = counter - 1
  end if
enddo
totnodes = counterend - counterstart +1

rewind(1)

do l = 1,counterstart-1
  read(1,*)
end do
do l = counterstart,counterend
  read(1,*) geoinput
  do k = 1,nodu
      if ( unodes(k) == geoinput%ints(1) ) then
        write(2,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), geoinput%floats(3)
      end if
    enddo
    do k = 1,nodl
      if ( lnodes(k) == geoinput%ints(1) ) then
        write(3,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), geoinput%floats(3)
      end if
    enddo
end do

close(1)
close (2)
close(3)
!-------------------------------------------------------------------------------
!---------------- Reading panel_data.dat ---------------------------------------

open(1,file = 'panel_data.dat')
read(1,'(a)') input
input = trim(input)
pos = index(input,":")
read(input(pos+1:), *) panels
read(1,'(a)') input
allocate(zpos(panels))
do i = 1,panels
  read(1,*) zpos(i)
enddo
close(1)

open(1,file = 'dumpallu.txt', iostat=ierr)
counter = 0
do while (ierr .eq. 0)
  read(1,*,iostat=ierr) geoinput
  if ( geoinput%floats(3) == zpos(1) ) then
    counter = counter + 1
  end if
enddo
close(1)
nx = counter
totalpanels = nodu/nx
open(1,file = 'dumpalll.txt', iostat=ierr)
counter = 0
do while (ierr .eq. 0)
  read(1,*,iostat=ierr) geoinput
  if ( geoinput%floats(3) == zpos(1) ) then
    counter = counter + 1
  end if
enddo
close(1)
nl = counter

do i = 1,panels
  do j = 1,panels
    if ( i == j ) then
    else
      if ( zpos(j) .gt. zpos(i) ) then
        buffer = zpos(i)
        zpos(i) = zpos(j)
        zpos(j) = buffer
      end if
    end if
  enddo
enddo

!-------------------------------------------------------------------------------
!---------------- Extracting data ----------------------------------------------

allocate(zpanelpos(totalpanels))
allocate(panelflag(totalpanels))


open(1,file = 'dumpallu.txt', iostat=ierr)
read(1,*,iostat=ierr) geoinput
zpanelpos(1) = geoinput%floats(3)
i = 2
do while (ierr .eq. 0)
  pos  = 0
  read(1,*,iostat=ierr) geoinput
  do j = 1,i-1
    if ( geoinput%floats(3) == zpanelpos(j) ) then
      pos = 1
      exit
    end if
  enddo
  if ( pos .eq. 0 ) then
    zpanelpos(i) = geoinput%floats(3)
    i = i + 1
  end if
enddo
close(1)

do i = 1,totalpanels
  do j = 1,totalpanels
    if ( i == j ) then
    else
      if ( zpanelpos(j) .gt. zpanelpos(i) ) then
        buffer = zpanelpos(i)
        zpanelpos(i) = zpanelpos(j)
        zpanelpos(j) = buffer
      end if
    end if
  enddo
enddo


open(1,file = 'dumpallu.txt', iostat=ierr)
open(2,file = 'dumpall_u.txt')


do i = 1,totalpanels
  do
    j = 1
    read(1,*,iostat=ierr) geoinput
    if ( ierr/=0 ) then
      exit
    end if
    if ( geoinput%floats(3) == zpanelpos(i) ) then
      write(2,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), zpanelpos(i)
      j = j + 1
    end if
  enddo
  rewind(1)
  ierr = 0
enddo
close(1)
close(2)

open(1,file = 'dumpalll.txt', iostat=ierr)
open(2,file = 'dumpall_l.txt')

do i = 1,totalpanels
  do
    j = 1
    read(1,*,iostat=ierr) geoinput
    if (ierr/=0) then
      exit
    endif
    if ( geoinput%floats(3) == zpanelpos(i) ) then
      write(2,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), zpanelpos(i)
      j = j + 1
    end if
  enddo
  rewind(1)
  ierr = 0
enddo
close(1)
close(2)


call system('rm dumpallu.txt')
call system('rm dumpalll.txt')



do i = 1,totalpanels
  panelflag(i) = 0
  do j = 1,panels
    if ( zpanelpos(i) == zpos(j) ) then
      panelflag(i) = 1
    end if
  enddo
enddo

allocate(masterpanelno(panels))
j = 1
do i = 1,totalpanels
  if ( panelflag(i) == 1 ) then
    masterpanelno(j) = i
    j = j + 1
  end if
enddo


!-------------------------------------------------------------------------------
!--------------- Sorting all the data with respect to x ------------------------

allocate(panelgeoinput(nx))

open(1,file = 'dumpall_u.txt', iostat=ierr)
open(2,file = 'dumpallu.txt')

do i = 1,totalpanels
  do j = 1,nx
    read(1,*,iostat=ierr) panelgeoinput(j)
  enddo

  do k = 1,nx
    do l = 1,nx
      if ( k == l ) then
      else
        if ( panelgeoinput(l)%floats(1) .lt. panelgeoinput(k)%floats(1) ) then
          temp = panelgeoinput(l)
          panelgeoinput(l) = panelgeoinput(k)
          panelgeoinput(k) = temp
        end if
      end if
    enddo
  enddo

  do k = 1,nx
    write(2,*) panelgeoinput(k)
  enddo
enddo
close(1)
close(2)

deallocate(panelgeoinput)
allocate(panelgeoinput(nl))

open(1,file = 'dumpall_l.txt', iostat=ierr)
open(2,file = 'dumpalll.txt')

do i = 1,totalpanels
  do j = 1,nl
    read(1,*,iostat=ierr) panelgeoinput(j)
  enddo

  do k = 1,nl
    do l = 1,nl
      if ( k == l ) then
      else
        if ( panelgeoinput(l)%floats(1) .gt. panelgeoinput(k)%floats(1) ) then
          temp = panelgeoinput(l)
          panelgeoinput(l) = panelgeoinput(k)
          panelgeoinput(k) = temp
        end if
      end if
    enddo
  enddo

  do k = 1,nl
    write(2,*) panelgeoinput(k)
  enddo
enddo
close(1)
close(2)

call system('rm dumpall_u.txt')
call system('rm dumpall_l.txt')

!-------------------------------------------------------------------------------
!--------------- Reading des_vars.dat ------------------------------------------

ndp = 0
open(1,file='des_vars.dat',status='old',iostat=ierr)

  do while (ierr.eq.0)
    read(1,*,iostat=ierr)
    ndp = ndp + 1
  end do

  ndp = (ndp - 1)/(2*panels)
  allocate(dpa(panels,ndp))
  allocate(dpb(panels,ndp))

rewind(1)
  do j = 1,panels
    do i = 1,ndp
      read(1,*,iostat=ierr) pairreal
      dpa(j,i) = pairreal(2)
    enddo
    do i = 1,ndp
      read(1,*,iostat=ierr) pairreal
      dpb(j,i) = pairreal(2)
    enddo
  enddo
close(1)


!-------------------------------------------------------------------------------
!----------------------- Hicks henne update ------------------------------------
allocate(panelgeoinputu(nx))
allocate(panelgeoinputl(nl))
allocate(bump_pos(ndp))

allocate(xu(nx))
allocate(yu(nx))
allocate(xl(nl))
allocate(yl(nl))
allocate(xunew(nx))
allocate(yunew(nx))
allocate(xlnew(nl))
allocate(ylnew(nl))


allocate(derivativenodenou(nx))
allocate(derivativenodenol(nl))
allocate(derivative(2*ndp*panels))



h = 1./(ndp+1)
do i = 1,ndp
  bump_pos(i) = i*h
enddo

open(1,file='dumpallu.txt',status='old',iostat=ierr)
open(2,file='tobeupdatedpanelsu.txt')

do while (ierr.eq.0)
  read(1,*,iostat=ierr) geoinput
  do i = 1,panels
    if ( geoinput%floats(3) == zpos(i) ) then
      write(2,*) geoinput
    end if
  enddo
enddo

close(1)
close(2)

open(1,file='dumpalll.txt',status='old',iostat=ierr)
open(2,file='tobeupdatedpanelsl.txt')

do while (ierr.eq.0)
  read(1,*,iostat=ierr) geoinput
  do i = 1,panels
    if ( geoinput%floats(3) == zpos(i) ) then
      write(2,*) geoinput
    end if
  enddo
enddo

close(1)
close(2)

open(1,file='tobeupdatedpanelsu.txt')
open(2,file='tobeupdatedpanelsl.txt')
open(3,file='updatedpanelu.txt')
open(4,file='updatedpanell.txt')

open(5,file='gradu.txt')
open(6,file='gradl.txt')

do k = 1,panels
  do j = 1,nx
    read(1,*,iostat=ierr) panelgeoinputu(j)
    xu(j) = panelgeoinputu(j)%floats(1)
    yu(j) = panelgeoinputu(j)%floats(2)
    derivativenodenou(j) = panelgeoinputu(j)%ints(1)

  enddo
  do j = 1,nl
    read(2,*) panelgeoinputl(j)
    xl(j) = panelgeoinputl(j)%floats(1)
    yl(j) = panelgeoinputl(j)%floats(2)
    derivativenodenol(j) = panelgeoinputl(j)%ints(1)
  enddo

  scaleu = abs(maxval(xu) - minval(xu))
  scalel = abs(maxval(xl) - minval(xl))

  xunew = xu/scaleu
  yunew = yu
  xlnew = xl/scalel
  ylnew = yl
  xutrans = minval(xunew)
  xltrans = minval(xlnew)
  xunew = xunew - xutrans
  xlnew = xlnew - xltrans

  do i = 1,nx
    sum = 0.
    derivative = 0.
    do j = 1,ndp
        m = log(0.5)/log(bump_pos(j))
        sum = sum + dpa(k,j)*(sin(pi*xunew(i)**m)**t_b)
        derivative((2*(k-1)*ndp)+j) = (sin(pi*xunew(i)**m)**t_b)
    enddo
    write(5,*) derivativenodenou(i), derivative
    yunew(i) = yunew(i) + sum
    xunew(i) = xunew(i) + xutrans
    xunew(i) = xunew(i) * scaleu
  enddo

  do i = 1,nl
    sum = 0.
    derivative = 0.
    do j = 1,ndp
        m = log(0.5)/log(bump_pos(j))
        sum = sum + dpb(k,j)*(sin(pi*xlnew(i)**m)**t_b)
        derivative((2*(k-1)*ndp)+ ndp + j) = (sin(pi*xlnew(i)**m)**t_b)
    enddo
    write(6,*) derivativenodenol(i), derivative
    ylnew(i) = ylnew(i) + sum
    xlnew(i) = xlnew(i) + xltrans
    xlnew(i) = xlnew(i) * scalel
  enddo

  do j = 1,nx
    panelgeoinputu(j)%floats(1) = xunew(j)
    panelgeoinputu(j)%floats(2) = yunew(j)
    write(3,*) panelgeoinputu(j)
  enddo
  do j = 1,nl
    panelgeoinputl(j)%floats(1) = xlnew(j)
    panelgeoinputl(j)%floats(2) = ylnew(j)
    write(4,*) panelgeoinputl(j)
  enddo
enddo
close(1)
close(2)
close(3)
close(4)
close(5)
close(6)

call system('rm tobeupdatedpanelsu.txt')
call system('rm tobeupdatedpanelsl.txt')


!-------------------------------------------------------------------------------
!-------------- Interpolation --------------------------------------------------

allocate(uleft(nx))
allocate(uright(nx))
allocate(lleft(nl))
allocate(lright(nl))

open(1,file='dumpallu.txt',status='old',iostat=ierr)
open(2,file='updatedpanelu.txt')
open(3,file='dumpupdatedu.txt')

do k = 1,nx*masterpanelno(1)
  read(1,*)
enddo
do k = 1,nx
  read(2,*) uleft(k)
  write(3,*) uleft(k)
enddo
do k = 1,nx
  read(2,*) uright(k)
enddo

do j = 1,panels-1
  do i = masterpanelno(j)+1,masterpanelno(j+1)-1
    do k = 1,nx
      read(1,*) geoinput
      buffer = uleft(k)%floats(2) + (((uright(k)%floats(2)-uleft(k)%floats(2))/(uright(k)%floats(3)-uleft(k)%floats(3)))*(geoinput%floats(3) - uleft(k)%floats(3)))
      write(3,*) geoinput%ints(1), geoinput%floats(1), buffer, geoinput%floats(3)
    enddo
  enddo

  do k = 1,nx
    read(1,*)
    write(3,*) uright(k)
  enddo
  uleft = uright
  do k = 1,nx
    read(2,*, iostat = ierr) uright(k)
  enddo
enddo


close(1)
close(2)
close(3)

open(1,file='dumpalll.txt',status='old',iostat=ierr)
open(2,file='updatedpanell.txt')
open(3,file='dumpupdatedl.txt')

do k = 1,nl*masterpanelno(1)
  read(1,*)
enddo
do k = 1,nl
  read(2,*) lleft(k)
  write(3,*) lleft(k)
enddo
do k = 1,nl
  read(2,*) lright(k)
enddo

do j = 1,panels-1
  do i = masterpanelno(j)+1,masterpanelno(j+1)-1
    do k = 1,nl
      read(1,*) geoinput
      buffer = lleft(k)%floats(2) + (((lright(k)%floats(2)-lleft(k)%floats(2))/(lright(k)%floats(3)-lleft(k)%floats(3)))*(geoinput%floats(3) - lleft(k)%floats(3)))
      write(3,*) geoinput%ints(1), geoinput%floats(1), buffer, geoinput%floats(3)
    enddo
  enddo

  do k = 1,nl
    read(1,*)
    write(3,*) lright(k)
  enddo
  lleft = lright
  do k = 1,nl
    read(2,*, iostat = ierr) lright(k)
  enddo
enddo


close(1)
close(2)
close(3)

call system('rm updatedpanelu.txt')
call system('rm updatedpanell.txt')
!-------------------------------------------------------------------------------
!------------------ Extrapolation ----------------------------------------------
open(1,file='dumpallu.txt',status='old',iostat=ierr)
open(2,file='dumpupdatedu.txt')
open(3,file='extrapolatedu.txt')

do k = 1,nx
  read(2,*) uleft(k)
enddo
do k = 1,nx
  read(2,*) uright(k)
enddo
do i = 1,masterpanelno(1)-1
  do k = 1,nx
    read(1,*) geoinput
    buffer = uleft(k)%floats(2) + (((uright(k)%floats(2)-uleft(k)%floats(2))/(uright(k)%floats(3)-uleft(k)%floats(3)))*(geoinput%floats(3) - uleft(k)%floats(3)))
    write(3,*) geoinput%ints(1), geoinput%floats(1), buffer, geoinput%floats(3)
  enddo
enddo
rewind(2)
rewind(1)
do
  read(2,*,iostat = ierr) geoinput
  if ( ierr .ne. 0 ) then
    exit
  end if
  write(3,*) geoinput
enddo
do k = 1,nx*masterpanelno(panels)
  read(1,*)
enddo
rewind(2)
do k = 1,masterpanelno(panels)-masterpanelno(1)-1
  do i = 1,nx
    read(2,*)
  enddo
enddo
do k = 1,nx
  read(2,*) uleft(k)
enddo
do k = 1,nx
  read(2,*) uright(k)
enddo
do i = masterpanelno(panels)+1,totalpanels
  do k = 1,nx
    read(1,*) geoinput
    buffer = uleft(k)%floats(2) + (((uright(k)%floats(2)-uleft(k)%floats(2))/(uright(k)%floats(3)-uleft(k)%floats(3)))*(geoinput%floats(3) - uleft(k)%floats(3)))
    write(3,*) geoinput%ints(1), geoinput%floats(1), buffer, geoinput%floats(3)
  enddo
enddo


close(1)
close(2)
close(3)



open(1,file='dumpalll.txt',status='old',iostat=ierr)
open(2,file='dumpupdatedl.txt')
open(3,file='extrapolatedl.txt')

do k = 1,nl
  read(2,*) lleft(k)
enddo
do k = 1,nl
  read(2,*) lright(k)
enddo
do i = 1,masterpanelno(1)-1
  do k = 1,nl
    read(1,*) geoinput
    buffer = lleft(k)%floats(2) + (((lright(k)%floats(2)-lleft(k)%floats(2))/(lright(k)%floats(3)-lleft(k)%floats(3)))*(geoinput%floats(3) - lleft(k)%floats(3)))
    write(3,*) geoinput%ints(1), geoinput%floats(1), buffer, geoinput%floats(3)
  enddo
enddo
rewind(2)
rewind(1)
do
  read(2,*,iostat = ierr) geoinput
  if ( ierr .ne. 0 ) then
    exit
  end if
  write(3,*) geoinput
enddo
do k = 1,nl*masterpanelno(panels)
  read(1,*)
enddo
rewind(2)
do k = 1,masterpanelno(panels)-masterpanelno(1)-1
  do i = 1,nl
    read(2,*)
  enddo
enddo
do k = 1,nl
  read(2,*) lleft(k)
enddo
do k = 1,nl
  read(2,*) lright(k)
enddo
do i = masterpanelno(panels)+1,totalpanels
  do k = 1,nl
    read(1,*) geoinput
    buffer = lleft(k)%floats(2) + (((lright(k)%floats(2)-lleft(k)%floats(2))/(lright(k)%floats(3)-lleft(k)%floats(3)))*(geoinput%floats(3) - lleft(k)%floats(3)))
    write(3,*) geoinput%ints(1), geoinput%floats(1), buffer, geoinput%floats(3)
  enddo
enddo


close(1)
close(2)
close(3)


!-------------------------------------------------------------------------------
!---------------- Concatenation ------------------------------------------------

open(1,file='baseline.txt')
open(2,file='dumpallu.txt')
open(3,file='dumpalll.txt')

do i = 1,totalpanels
  do k = 1,nx
    read(2,*) geoinput
    write(1,*) geoinput
  enddo
  do k = 1,nl
    read(3,*) geoinput
    write(1,*) geoinput
  enddo
enddo

close(1)
close(2)
close(3)

open(1,file='new.txt')
open(2,file='extrapolatedu.txt')
open(3,file='extrapolatedl.txt')

do i = 1,totalpanels
  do k = 1,nx
    read(2,*) geoinput
    write(1,*) geoinput
  enddo
  do k = 1,nl
    read(3,*) geoinput
    write(1,*) geoinput
  enddo
enddo

close(1)
close(2)
close(3)

call system('rm dumpupdatedu.txt')
call system('rm dumpupdatedl.txt')



!-------------------------------------------------------------------------------
!------------------ Area Calculation of each panel -----------------------------


allocate(areaold(totalpanels))
allocate(areanew(totalpanels))

allocate(x(nx+nl))
allocate(y(nx+nl))

open(1,file='baseline.txt')
open(2,file='new.txt')

do j = 1,totalpanels
  !---- old areas
  do k = 1,nx+nl
    read(1,*)geoinput
    x(k) = geoinput%floats(1)
    y(k) = geoinput%floats(2)
  enddo

  firstterm = 0.
  secondterm = 0.
  do i = 1,nx+nl-1
    firstterm = firstterm + x(i)*y(i+1)
    secondterm = secondterm + x(i+1)*y(i)
  enddo
  area = 0.5*abs(firstterm - secondterm + (x(nx+nl)*y(1)) - (x(1)*y(nx+nl)))
  areaold(j) = area

  !----- new areas
  do k = 1,nx+nl
    read(2,*)geoinput
    x(k) = geoinput%floats(1)
    y(k) = geoinput%floats(2)
  enddo

  firstterm = 0.
  secondterm = 0.
  do i = 1,nx+nl-1
    firstterm = firstterm + x(i)*y(i+1)
    secondterm = secondterm + x(i+1)*y(i)
  enddo
  area = 0.5*abs(firstterm - secondterm + (x(nx+nl)*y(1)) - (x(1)*y(nx+nl)))
  areanew(j) = area

enddo

close(1)
close(2)


!-------------------------------------------------------------------------------
!---------------- Calculating volumes ------------------------------------------

volumeold = 0.
volumenew = 0.
do i = 1,totalpanels-1
  volumeold = volumeold + (((areaold(i+1)+areaold(i))/2) * (zpanelpos(i+1) - zpanelpos(i)))
  volumenew = volumenew + (((areanew(i+1)+areanew(i))/2) * (zpanelpos(i+1) - zpanelpos(i)))
enddo


!-------------------------------------------------------------------------------
!------------- Gradient Interpolation ------------------------------------------

allocate(graddataleft(nx,panels*2*ndp))
allocate(graddataright(nx,panels*2*ndp))


open(1,file='dumpallu.txt',status='old',iostat=ierr)
open(2,file='gradu.txt')
open(3,file='gradallu.txt')

do k = 1,nx*masterpanelno(1)
  read(1,*)
enddo
do k = 1,nx
  read(2,*) dummy, graddataleft(k,:)
  write(3,*) dummy, graddataleft(k,:)
enddo
do k = 1,nx
  read(2,*) dummy, graddataright(k,:)
enddo

do j = 1,panels-1
  do i = masterpanelno(j)+1,masterpanelno(j+1)-1
    do k = 1,nx
      read(1,*) geoinput
      derivative = 0.
      derivative = graddataleft(k,:) + (graddataright(k,:) - graddataleft(k,:))*((geoinput%floats(3) - zpanelpos(masterpanelno(j)))/(zpanelpos(masterpanelno(j+1)) - zpanelpos(masterpanelno(j))))
      write(3,*) geoinput%ints(1), derivative
    enddo
  enddo

  do k = 1,nx
    read(1,*) geoinput
    write(3,*) geoinput%ints(1), graddataright(k,:)
    graddataleft(k,:) = graddataright(k,:)
  enddo
  do k = 1,nx
    read(2,*, iostat = ierr) dummy, graddataright(k,:)
  enddo
enddo

close(1)
close(2)
close(3)

deallocate(graddataleft)
deallocate(graddataright)
allocate(graddataleft(nl,panels*2*ndp))
allocate(graddataright(nl,panels*2*ndp))

open(1,file='dumpalll.txt',status='old',iostat=ierr)
open(2,file='gradl.txt')
open(3,file='gradalll.txt')

do k = 1,nl*masterpanelno(1)
  read(1,*)
enddo
do k = 1,nl
  read(2,*) dummy, graddataleft(k,:)
  write(3,*) dummy, graddataleft(k,:)
enddo
do k = 1,nl
  read(2,*) dummy, graddataright(k,:)
enddo

do j = 1,panels-1
  do i = masterpanelno(j)+1,masterpanelno(j+1)-1
    do k = 1,nl
      read(1,*) geoinput
      derivative = 0.
      derivative = graddataleft(k,:) + (graddataright(k,:) - graddataleft(k,:))*((geoinput%floats(3) - zpanelpos(masterpanelno(j)))/(zpanelpos(masterpanelno(j+1)) - zpanelpos(masterpanelno(j))))
      write(3,*) geoinput%ints(1), derivative
    enddo
  enddo

  do k = 1,nl
    read(1,*) geoinput
    write(3,*) geoinput%ints(1), graddataright(k,:)
    graddataleft(k,:) = graddataright(k,:)
  enddo
  do k = 1,nl
    read(2,*, iostat = ierr) dummy, graddataright(k,:)
  enddo
enddo

close(1)
close(2)
close(3)

deallocate(graddataleft)
deallocate(graddataright)


!-------------------------------------------------------------------------------
!---------- Gradient Extrapolation ---------------------------------------------

allocate(graddataright(nx,panels*2*ndp))
allocate(graddataleft(nx,panels*2*ndp))

open(1,file='dumpallu.txt',status='old',iostat=ierr)
open(2,file='gradallu.txt')
open(3,file='extrapolatedgradu.txt')
open(4,file='gradu.txt')

do k = 1,nx
  read(4,*) dummy, graddataleft(k,:)
enddo
do k = 1,nx
  read(4,*) dummy, graddataright(k,:)
enddo
do i = 1,masterpanelno(1)-1
  do k = 1,nx
    read(1,*) geoinput
    derivative = 0.
    derivative = graddataleft(k,:) + (graddataright(k,:) - graddataleft(k,:))*((geoinput%floats(3) - zpanelpos(masterpanelno(1)))/(zpanelpos(masterpanelno(2)) - zpanelpos(masterpanelno(1))))
    write(3,*) geoinput%ints(1), derivative
  enddo
enddo

do
  read(2,*,iostat = ierr) dummy, derivative
  if ( ierr .ne. 0 ) then
    exit
  end if
  read(1,*)
  write(3,*) dummy, derivative
enddo
rewind(4)
do i = 1,panels-2
  do k = 1,nx
    read(4,*)
  enddo
enddo
do k = 1,nx
  read(4,*) dummy, graddataleft(k,:)
enddo
do k = 1,nx
  read(4,*) dummy, graddataright(k,:)
enddo
do i = masterpanelno(panels)+1,totalpanels
  do k = 1,nx
    read(1,*) geoinput
    derivative = 0.
    derivative = graddataleft(k,:) + (graddataright(k,:) - graddataleft(k,:))*((geoinput%floats(3) - zpanelpos(masterpanelno(panels)-1))/(zpanelpos(masterpanelno(panels)) - zpanelpos(masterpanelno(panels)-1)))
    write(3,*) geoinput%ints(1), derivative
  enddo
enddo


close(1)
close(2)
close(3)
close(4)
deallocate(graddataleft)
deallocate(graddataright)

allocate(graddataright(nl,panels*2*ndp))
allocate(graddataleft(nl,panels*2*ndp))

open(1,file='dumpalll.txt',status='old',iostat=ierr)
open(2,file='gradalll.txt')
open(3,file='extrapolatedgradl.txt')
open(4,file='gradl.txt')

do k = 1,nl
  read(4,*) dummy, graddataleft(k,:)
enddo
do k = 1,nl
  read(4,*) dummy, graddataright(k,:)
enddo
do i = 1,masterpanelno(1)-1
  do k = 1,nl
    read(1,*) geoinput
    derivative = 0.
    derivative = graddataleft(k,:) + (graddataright(k,:) - graddataleft(k,:))*((geoinput%floats(3) - zpanelpos(masterpanelno(1)))/(zpanelpos(masterpanelno(2)) - zpanelpos(masterpanelno(1))))
    write(3,*) geoinput%ints(1), derivative
  enddo
enddo
do
  read(2,*,iostat = ierr) dummy, derivative
  if ( ierr .ne. 0 ) then
    exit
  end if
  read(1,*)
  write(3,*) dummy, derivative
enddo
rewind(4)
do i = 1,panels-2
  do k = 1,nl
    read(4,*)
  enddo
enddo
do k = 1,nl
  read(4,*) dummy, graddataleft(k,:)
enddo
do k = 1,nl
  read(4,*) dummy, graddataright(k,:)
enddo
do i = masterpanelno(panels)+1,totalpanels
  do k = 1,nl
    read(1,*) geoinput
    derivative = 0.
    derivative = graddataleft(k,:) + (graddataright(k,:) - graddataleft(k,:))*((geoinput%floats(3) - zpanelpos(masterpanelno(panels)-1))/(zpanelpos(masterpanelno(panels)) - zpanelpos(masterpanelno(panels)-1)))
    write(3,*) geoinput%ints(1), derivative
  enddo
enddo


close(1)
close(2)
close(3)
close(4)
deallocate(graddataleft)
deallocate(graddataright)

!-------------------------------------------------------------------------------
!---------------- Gradient Concatenation ---------------------------------------

open(1,file='allgrad.txt')
open(2,file='extrapolatedgradu.txt')
open(3,file='extrapolatedgradl.txt')

do i = 1,totalpanels
  do k = 1,nx
    read(2,*) dummy,derivative
    write(1,*) dummy,derivative
  enddo
  do k = 1,nl
    read(3,*) dummy,derivative
    write(1,*) dummy,derivative
  enddo
enddo

close(1)
close(2)
close(3)

!-------------------------------------------------------------------------------
!-------------- Areagrad Calculation of each panel -----------------------------
!-------------- and Volume gradient Calculation --------------------------------

allocate(areagrad(totalpanels,panels*2*ndp))
allocate(graddata(nx+nl,panels*2*ndp))

open(1,file='allgrad.txt')
open(2,file='baseline.txt')

do i = 1,totalpanels
  do j = 1,nx+nl
    read(1,*)dummy, derivative
    graddata(j,:) = derivative
    read(2,*) geoinput
    x(j) = geoinput%floats(1)
  enddo
  do k = 1,panels*2*ndp
    firstterm = 0.
    secondterm = 0.
    do j = 1,nx+nl-1
      firstterm = firstterm + x(j)*graddata(j+1,k)
      secondterm = secondterm + x(j+1)*graddata(j,k)
    enddo
    area = 0.5*(firstterm - secondterm + (x(nx+nl)*graddata(1,k)) - (x(1)*graddata(nx+nl,k)))
    areagrad(i,k) = area
  enddo
enddo

close(1)
close(2)

allocate(volumegrad(panels*2*ndp))


do k = 1,panels*2*ndp
  buffer = 0.
  do i = 1,totalpanels-1
    buffer = buffer + (((areagrad(i+1,k)+areagrad(i,k))/2) * (zpanelpos(i+1) - zpanelpos(i)))
  enddo
  volumegrad(k) = buffer
enddo




!-------------------------------------------------------------------------------
!------------- Reading Adjoint vector and writing {casename}_grad.dat ----------

allocate(lambda(totnodes,3))
open(1,file=adjoint,status='old',iostat=ierr)
do i = 1,4
  read(1,*) input
enddo
do i = 1,totnodes
  read(1,*) buffer
  lambda(i,1) = buffer
enddo
do i = 1,totnodes
  read(1,*) buffer
  lambda(i,2) = buffer
enddo
do i = 1,totnodes
  read(1,*) buffer
  lambda(i,3) = buffer
enddo
close(1)

open(1,file='allgrad.txt',iostat = ierr)
open(2,file='afterdottingwithlambda.txt')
do
  read(1,*,iostat = ierr) dummy, derivative
  if ( ierr/=0 ) then
    exit
  end if
  do i = 1,panels*2*ndp
    derivative(i) = 0.*lambda(dummy,1) + derivative(i)*lambda(dummy,2) + 0.*lambda(dummy,3)
  enddo
  write(2,*) dummy,derivative
enddo
close(1)
close(2)



open(1,file='afterdottingwithlambda.txt',iostat = ierr)
open(2,file=gradtitle)
write(2,*) 'START'
do i = 1,panels*2*ndp
  sum = 0.
  do
    read(1,*,iostat = ierr) dummy, derivative
    if ( ierr/=0 ) then
      exit
    end if
    sum = sum + derivative(i)
  enddo
  write(2,*)i,sum + gamma*(volumeold-0.5)*volumegrad(i)
  rewind(1)
enddo
write(2,*) 'END'
close(1)
close(2)

open(unit = 1, file = 'functional.dat')
open(unit = 2, file = 'functional_all.dat')
read(1,*)buffer
write(2,*)buffer
write(2,*)volumeold
rewind(1)
write(1,*)buffer + (gamma*0.5*(volumeold-0.5)**2)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
print *,(volumenew-volumeold)/1e-8
print *,volumegrad(9)
do k = 1,totalpanels
  print *,areagrad(k,:)
enddo

!call system('rm dumpallu.txt')
!call system('rm dumpalll.txt')
!call system('rm gradu.txt')
!call system('rm gradl.txt')
!call system('rm gradallu.txt')
!call system('rm gradalll.txt')
!call system('rm afterdottingwithlambda.txt')
!call system('rm allgrad.txt')
!call system('rm baseline.txt')
!call system('rm new.txt')
!call system('rm extrapolatedu.txt')
!call system('rm extrapolatedl.txt')
!call system('rm extrapolatedgradu.txt')
!call system('rm extrapolatedgradl.txt')




end program param
