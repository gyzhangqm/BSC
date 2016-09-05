program param
  implicit none

  integer :: nodu, nodl, i, j, k, ierr ,il, counterstart, counterend, counter, l, ndp, totnodes
  integer :: pair(2), flag, pos
  real, dimension(:), allocatable :: xu, yu, xl, yl, bump_pos, xunew, yunew, xlnew, ylnew, dispu, displ, gradu(:,:), gradl(:,:), grad(:,:), lambda(:,:), gradout(:,:), adjout
  integer, dimension(:), allocatable :: unodes, lnodes
  real, dimension(:), allocatable :: dpa, dpb
  real (kind = 8) :: three(3), pairreal(2), h, scaleu, scalel, sum, m, xutrans, xltrans, dummy
  integer, parameter :: t_b = 4
  real(kind = 8), parameter :: pi = acos(-1.)
  character(256) :: input, cmd, geo, fixnod, adjoint, code
  integer :: foo, flagupper, flaglower

  call get_command_argument(1, cmd)
  geo = trim(cmd)//".geo.dat"
  fixnod = trim(cmd)//".fix.nod"
  adjoint = trim(cmd)//".ensi.DISPM-000001"
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
        j = j +1
      endif
    endif

    if(ierr==5010) then

    else if (ierr/=0) then
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
open(1,file=geo,status='old',iostat=ierr)
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
totnodes = counterend - counterstart +1
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
scaleu = abs(maxval(xu) - minval(xu))
scalel = abs(maxval(xl) - minval(xl))

allocate(bump_pos(ndp))
allocate(xunew(nodu))
allocate(yunew(nodu))
allocate(xlnew(nodl))
allocate(ylnew(nodl))
allocate(dispu(nodl))
allocate(displ(nodl))
allocate(gradu(nodu,ndp))
allocate(gradl(nodl,ndp))


h = 1./(ndp+1)
do i = 1,ndp
  bump_pos(i) = i*h
enddo
xunew = xu/scaleu
yunew = yu/scaleu
xlnew = xl/scalel
ylnew = yl/scalel
xutrans = minval(xunew)
xltrans = minval(xlnew)
xunew = xunew - xutrans
xlnew = xlnew - xltrans


do i = 1,nodu
  sum = 0.
  do j = 1,ndp
      m = log(0.5)/log(bump_pos(j))
      sum = sum + dpa(j)*(sin(pi*xunew(i)**m)**t_b)
      gradu(i,j) = sin(pi*xunew(i)**m)**t_b;
  enddo
  yunew(i) = yunew(i) + sum
  yunew(i) = yunew(i) * scaleu
  dispu(i) = yunew(i) - yu(i)
  xunew(i) = xunew(i) + xutrans
  xunew(i) = xunew(i) * scaleu
enddo

do i = 1,nodl
  sum = 0.
  do j = 1,ndp
      m = log(0.5)/log(bump_pos(j))
      sum = sum + dpb(j)*(sin(pi*xlnew(i)**m)**t_b)
      gradl(i,j) = sin(pi*xlnew(i)**m)**t_b;
  enddo
  ylnew(i) = ylnew(i) + sum
  ylnew(i) = ylnew(i) * scalel
  displ(i) = ylnew(i) - yl(i)
  xlnew(i) = xlnew(i) + xltrans
  xlnew(i) = xlnew(i) * scalel

enddo
!-------------------------------------------------------------------------------
!---- Making Grad matrix -------------------------------------------------------
allocate(grad(totnodes,2*ndp))
do k=1,totnodes
  flag = 0
  do i = 1, nodu
    if ( unodes(i) == k ) then
      !grad(k,:) = gradu(i,:)
      do j = 1,ndp
        grad(k,j) = gradu(i,j)
        grad(k,ndp+j) = 0.
      enddo
      flag = 1
    end if
  enddo
  do i = 1, nodl
    if ( lnodes(i) == k ) then
      grad(k,:) = gradl(i,:)
      do j = 1,ndp
        grad(k,ndp+j) = gradl(i,j)
        grad(k,j) = 0.
      enddo
      flag = 1
    end if
  enddo
  if ( flag == 0 ) then
    grad(k,:) = 0.
  end if
enddo


!-------------------------------------------------------------------------------
!------ Reading Lambda from Adjoint Alya  and writing gradout.txt---------------
allocate(lambda(totnodes,2))
allocate(gradout(totnodes,2*ndp))
allocate(adjout(2*ndp))
open(1,file=adjoint,status='old',iostat=ierr)
do i = 1,4
  read(1,*) input
enddo
do i = 1,totnodes
  read(1,*) dummy
  lambda(i,1) = dummy
enddo
do i = 1,totnodes
  read(1,*) dummy
  lambda(i,2) = dummy
enddo
close(1)

do k = 1,totnodes
  do i = 1,2*ndp
    gradout(k,i) = 0.*lambda(k,1) + grad(k,i)*lambda(k,2)
  enddo
enddo

do i = 1,2*ndp
  adjout(i) = 0.
  do k = 1,totnodes
    adjout(i) = adjout(i) + gradout(k,i)
  enddo
enddo

open(unit = 1, file = 'grad.txt', status = 'unknown', iostat = ierr)
write(1,*)"START"
do k = 1,2*ndp
  !write(1,*) (grad(k,j), j=1,2*ndp)
  write(1,*) adjout(k)
enddo
write(1,*)"END"
close(1)


end program param
