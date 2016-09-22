program param
  implicit none

  type intreal
    integer,dimension(1)   :: ints
    real,dimension(2)      :: floats
  endtype intreal

  type(intreal) :: mesh_disp, geoinput, geoinput1, geoinput2, temp
  type(intreal), dimension(:), allocatable :: panelgeoinput



  integer :: nodu, nodl, nodf, i, j, k, p, ierr ,il, counterstart, counterend, counter, l, ndp, totnodes
  integer :: pair(2), flag, pos, dummy
  real, dimension(:), allocatable :: xu, yu, xl, yl, xf, yf, bump_pos, derivativeu(:,:)
  real, dimension(:), allocatable :: derivativel(:,:), derivativefx(:,:), derivativefy(:,:), derivative, zeroderivative
  integer, dimension(:), allocatable :: unodes, lnodes, fnodes
  real, dimension(:), allocatable :: dpa, dpb, x, y, graddata(:,:), areagrad, lambda(:,:)
  real (kind = 8) :: three(3), pairreal(2), h, scaleu, scalel, sum, m, xutrans, xltrans, xflap, yflap, flapangle
  real (kind = 8) :: xref, yref,firstterm, secondterm, area, buffer
  integer, parameter :: t_b = 4
  real(kind = 8), parameter :: pi = acos(-1.), gamma = 0.
  character(256) :: input, cmd, geo, fixnod, code, adjoint, gradtitle
  integer :: foo, flagupper, flaglower, flagflap

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
      if ( panelgeoinput(l)%floats(1) .lt. panelgeoinput(k)%floats(1) ) then
        temp = panelgeoinput(l)
        panelgeoinput(l) = panelgeoinput(k)
        panelgeoinput(k) = temp
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
      if ( panelgeoinput(l)%floats(1) .gt. panelgeoinput(k)%floats(1) ) then
        temp = panelgeoinput(l)
        panelgeoinput(l) = panelgeoinput(k)
        panelgeoinput(k) = temp
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

allocate(derivativeu(nodu,(2*ndp)+3))
allocate(derivativel(nodl,(2*ndp)+3))

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

h = 1./(ndp+1)
do i = 1,ndp
  bump_pos(i) = i*h
enddo
xu = xu/scaleu
yu = yu

xl = xl/scalel
yl = yl
xutrans = minval(xu)
xltrans = minval(xl)
xu = xu - xutrans
xl = xl - xltrans


do i = 1,nodu
  derivativeu(i,:) = 0.
  do j = 1,ndp
      m = log(0.5)/log(bump_pos(j))
      derivativeu(i,j) = (sin(pi*xu(i)**m)**t_b)
  enddo
  xu(i) = xu(i) + xutrans
  xu(i) = xu(i) * scaleu
enddo

do i = 1,nodl
  derivativel(i,:) = 0.
  do j = 1,ndp
      m = log(0.5)/log(bump_pos(j))
      derivativel(i,ndp+j) = (sin(pi*xl(i)**m)**t_b)
  enddo
  xl(i) = xl(i) + xltrans
  xl(i) = xl(i) * scalel
enddo

open(2,file='dumpallu.txt')
open(3,file='gradu.txt')
do k = 1,nodu
  read(2,*) geoinput
  write(3,*) geoinput%ints(1), derivativeu(k,:)
enddo
close(2)
close(3)

open(2,file='dumpalll.txt')
open(3,file='gradl.txt')
do k = 1,nodl
  read(2,*) geoinput
  write(3,*) geoinput%ints(1), derivativel(k,:)
enddo
close(2)
close(3)

!-------------------------------------------------------------------------------
!---------------------- Flap update --------------------------------------------
allocate(xf(nodf))
allocate(yf(nodf))

allocate(derivativefy(nodf,(2*ndp)+3))
allocate(derivativefx(nodf,(2*ndp)+3))

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
  derivativefx(k,:) = 0.
  derivativefy(k,:) = 0.
  derivativefx(k,(2*ndp)+1) = 1.
  derivativefx(k,(2*ndp)+3) = -((xf(k) - xref)*sin(flapangle)) - ((yf(k) - yref)*cos(flapangle))
  derivativefy(k,(2*ndp)+2) = 1.
  derivativefy(k,(2*ndp)+3) = ((xf(k) - xref)*cos(flapangle)) - ((yf(k) - yref)*sin(flapangle))
enddo

open(2,file='dumpflap.txt')
open(3,file='gradfx.txt')
open(4,file='gradfy.txt')
do k = 1,nodf
  read(2,*) geoinput
  write(3,*) geoinput%ints(1), derivativefx(k,:)
  write(4,*) geoinput%ints(1), derivativefy(k,:)
enddo
close(2)
close(3)
close(4)
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


!-------------------------------------------------------------------------------
!------------ Area Calculation -------------------------------------------------

allocate(x(nodu+nodl))
allocate(y(nodu+nodl))

open(1,file='baseline.txt')
do k = 1,nodu+nodl
  read(1,*) geoinput
  x(k) = geoinput%floats(1)
  y(k) = geoinput%floats(2)
enddo
close(1)
firstterm = 0.
secondterm = 0.
do i = 1,nodu+nodl-1
  firstterm = firstterm + x(i)*y(i+1)
  secondterm = secondterm + x(i+1)*y(i)
enddo
area = 0.5*abs(firstterm - secondterm + (x(nodu+nodl)*y(1)) - (x(1)*y(nodu+nodl)))


!-------------------------------------------------------------------------------
call system('rm dumpallu.txt')
call system('rm dumpalll.txt')
call system('rm dumpflap.txt')


!-------------------------------------------------------------------------------
!------------ Gradient Concatenation -------------------------------------------

allocate(derivative((2*ndp)+3))
allocate(zeroderivative((2*ndp)+3))
open(1,file='gradu.txt')
open(2,file='gradl.txt')
open(3,file='gradfy.txt')
open(4,file='gradally.txt')
open(5,file='gradallx.txt')
open(6,file='gradfx.txt')
zeroderivative = 0.
do k = 1,nodu
  read(1,*) dummy, derivative
  write(4,*) dummy,derivative
  write(5,*) dummy, zeroderivative
enddo
do k = 1,nodl
  read(2,*) dummy, derivative
  write(4,*) dummy,derivative
  write(5,*) dummy, zeroderivative
enddo
do k = 1,nodf
  read(3,*) dummy, derivative
  write(4,*) dummy,derivative
  read(6,*) dummy, derivative
  write(5,*) dummy,derivative
enddo

close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
!-------------------------------------------------------------------------------

call system('rm gradu.txt')
call system('rm gradl.txt')
call system('rm gradfx.txt')
call system('rm gradfy.txt')

!-------------------------------------------------------------------------------
!------------ areagrad calculation ---------------------------------------------
allocate(graddata(nodu+nodl, (2*ndp)+3))
allocate(areagrad((2*ndp)+3))

open(1,file='gradally.txt')
open(2,file='baseline.txt')

do j = 1,nodu+nodl
  read(1,*)dummy, graddata(j,:)
  read(2,*) geoinput
  x(j) = geoinput%floats(1)
enddo

do k = 1,(2*ndp)+3
  firstterm = 0.
  secondterm = 0.
  do j = 1,nodu+nodl-1
    firstterm = firstterm + x(j)*graddata(j+1,k)
    secondterm = secondterm + x(j+1)*graddata(j,k)
  enddo
  areagrad(k) = 0.5*(firstterm - secondterm + (x(nodu+nodl)*graddata(1,k)) - (x(1)*graddata(nodu+nodl,k)))
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

open(1,file='gradally.txt',iostat = ierr)
open(2,file='afterdottingwithlambda.txt')
open(3,file='gradallx.txt')
do
  read(1,*,iostat = ierr) dummy, derivative
  read(3,*,iostat = ierr) dummy, zeroderivative
  if ( ierr/=0 ) then
    exit
  end if
  do i = 1,(2*ndp)+3
    derivative(i) = zeroderivative(i)*lambda(dummy,1) + derivative(i)*lambda(dummy,2) + 0.*lambda(dummy,3)
  enddo
  write(2,*) dummy,derivative
enddo
close(1)
close(2)
close(3)

open(1,file='afterdottingwithlambda.txt',iostat = ierr)
open(2,file=gradtitle)
write(2,*) 'START'
do i = 1,(2*ndp)+3
  sum = 0.
  do
    read(1,*,iostat = ierr) dummy, derivative
    if ( ierr/=0 ) then
      exit
    end if
    sum = sum + derivative(i)
  enddo
  write(2,*)sum + gamma*(area-0.5)*areagrad(i)
  rewind(1)
enddo
write(2,*) 'END'
close(1)
close(2)

open(unit = 1, file = 'functional.dat')
open(unit = 2, file = 'functional_all.dat')
read(1,*)buffer
write(2,*)buffer
write(2,*)area
rewind(1)
write(1,*)buffer + (gamma*0.5*(area-0.5)**2)

!-------------------------------------------------------------------------------

call system('rm gradallx.txt')
call system('rm gradally.txt')
call system('rm afterdottingwithlambda.txt')
call system('rm baseline.txt')

end program param
