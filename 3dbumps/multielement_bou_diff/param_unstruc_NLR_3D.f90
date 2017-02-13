program param
  implicit none

  type intreal
    integer,dimension(1)   :: ints
    real,dimension(3)      :: floats
  endtype intreal
  type(intreal), dimension(:), allocatable :: datau, datal, dataunew, datalnew, dataf, datafnew
  type(intreal) :: geoinput

  real(kind = 8), parameter :: pi = acos(-1.)
  integer, parameter :: bumps = 5, panels = 6, t_b = 4
  integer :: pos, flagupper, flaglower, flagflap, nodu, nodl,nodf, i, j, k, l,p , ierr, pair(2), counter, counterstart, counterend, totnodes, ndp
  integer, dimension(:), allocatable :: unodes, lnodes,fnodes, flagpanel
  character(256) :: input, cmd, geo, fixnod, code
  real, dimension(:), allocatable :: temp, temp2, des_varsu, des_varsl, y_pos, bump_pos, func(:,:), spanwise, chordwise
  real(kind = 8) :: rootchord, tipchord, rooty, tipy, rootlex, roottex, tiplex, tiptex, rootflapx, rootflapz, tipflapx, tipflapz, pairreal(2), m, xtrans, scale, sum, x ,y, xflap, zflap, flapangle, xref,zref

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
  l = 1
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
        fnodes(l) = pair(1)
        l = l + 1
      endif
    endif

    if(ierr==5010) then

    else if (ierr/=0) then
      exit
    end if
  enddo
  close(1)
!-------------------------------------------------------------------------------
!Extracts the coordinates of the grid points which lie on the upper and lower surface

allocate(datau(nodu))
allocate(datal(nodl))
allocate(dataunew(nodu))
allocate(datalnew(nodl))
allocate(dataf(nodf))
allocate(datafnew(nodf))

  counter = 0
  i = 1
  j = 1
  p = 1

  open(1,file=geo,status='old',iostat=ierr)
  do while (ierr.eq.0)
    read(1,'(A)',iostat=ierr) input
    counter = counter + 1
    if ( trim(input) == '  COORDINATES' ) then
      counterstart = counter + 1
    end if
    if ( trim(input) == '  END_COORDINATES' ) then
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
          datau(i) = geoinput
          i = i + 1
        end if
      enddo
      do k = 1,nodl
        if ( lnodes(k) == geoinput%ints(1) ) then
          datal(j) = geoinput
          j = j  + 1
        end if
      enddo
      do k = 1,nodf
        if ( fnodes(k) == geoinput%ints(1) ) then
          dataf(p) = geoinput
          p = p  + 1
        end if
      enddo
  end do

  close(1)

  open(1,file = 'baseline.txt')
  do i = 1,nodu
    write(1,*) datau(i)
  enddo
  do i = 1,nodl
    write(1,*) datal(i)
  enddo
  do i = 1,nodf
    write(1,*) dataf(i)
  enddo
  close(1)

!-------------------------------------------------------------------------------
!----- Finding root chord length and tip chord length --------------------------
allocate(temp(nodu))
do i = 1,nodu
  temp(i) = datau(i)%floats(2)
enddo
rooty = minval(temp)
tipy = maxval(temp)
deallocate(temp)

counter = 0
do i = 1,nodu
  if ( datau(i)%floats(2) == rooty ) then
    counter = counter + 1
  end if
enddo
allocate(temp(counter))
j = 1
do i = 1,nodu
  if ( datau(i)%floats(2) == rooty ) then
    temp(j) = datau(i)%floats(1)
    j = j + 1
  end if
enddo
rootchord =abs(maxval(temp) - minval(temp))
rootlex = minval(temp)
roottex = maxval(temp)
deallocate(temp)

counter = 0
do i = 1,nodu
  if ( datau(i)%floats(2) == tipy ) then
    counter = counter + 1
  end if
enddo
allocate(temp(counter))
j = 1
do i = 1,nodu
  if ( datau(i)%floats(2) == tipy ) then
    temp(j) = datau(i)%floats(1)
    j = j + 1
  end if
enddo
tipchord = abs(maxval(temp) - minval(temp))
tiplex = minval(temp)
tiptex = maxval(temp)
deallocate(temp)

counter = 0
do i = 1,nodf
  if ( dataf(i)%floats(2) == rooty ) then
    counter = counter + 1
  end if
enddo
allocate(temp(counter))
allocate(temp2(counter))
j = 1
do i = 1,nodf
  if ( dataf(i)%floats(2) == rooty ) then
    temp(j) = dataf(i)%floats(1)
    temp2(j) = dataf(i)%floats(3)
    j = j + 1
  end if
enddo
rootflapx = minval(temp)
pos = minloc(temp,1)
rootflapz = temp2(pos)
deallocate(temp)
deallocate(temp2)

counter = 0
do i = 1,nodf
  if ( dataf(i)%floats(2) == tipy ) then
    counter = counter + 1
  end if
enddo
allocate(temp(counter))
allocate(temp2(counter))
j = 1
do i = 1,nodf
  if ( dataf(i)%floats(2) == tipy ) then
    temp(j) = dataf(i)%floats(1)
    temp2(j) = dataf(i)%floats(3)
    j = j + 1
  end if
enddo
tipflapx = minval(temp)
pos = minloc(temp,1)
tipflapz = temp2(pos)
deallocate(temp)
deallocate(temp2)

!print *, rooty, 'rooty'
!print *, tipy, 'tipy'
!print *, rootchord, 'rootchord'
!print *, tipchord, 'tipchord'
!print *, rootlex, 'rootlex'
!print *, tiplex, 'tiplex'
!print *, roottex, 'roottex'
!print *, tiptex, 'tiptex'

!-------------------------------------------------------------------------------
!------- Reading des_vars ------------------------------------------------------
ndp = (bumps)*panels
allocate(des_varsu(ndp))
allocate(des_varsl(ndp))
open(1,file = 'des_vars.dat')
do i = 1,ndp
  read(1,*) pairreal
  des_varsu(i) = pairreal(2)
enddo
do i = 1,ndp
  read(1,*) pairreal
  des_varsl(i) = pairreal(2)
enddo
read(1,*) pairreal
xflap = pairreal(2)
read(1,*) pairreal
zflap = pairreal(2)
read(1,*) pairreal
flapangle = pairreal(2)
close(1)

!-------------------------------------------------------------------------------
!------------- Update ----------------------------------------------------------
allocate(y_pos(panels))
allocate(bump_pos(bumps))
allocate(func(panels,bumps))
allocate(spanwise(panels))
allocate(chordwise(bumps))
allocate(flagpanel(panels-1))

do i = 1,bumps
  bump_pos(i) = i*1./(bumps+1)
enddo

do i = 1,panels
  y_pos(i) = rooty + ((i-1)*(tipy-rooty)/(panels-1))
enddo


do i = 1,nodu
  x = datau(i)%floats(1)
  y = datau(i)%floats(2)
  xtrans = rootlex + ((tiplex-rootlex)*(y-rooty)/(tipy-rooty))
  scale = rootchord + ((tiptex-roottex)*(y-rooty)/(tipy-rooty)) - ((tiplex-rootlex)*(y-rooty)/(tipy-rooty))
  x = x - xtrans
  x = x/scale
  if ( x .le. 0.0 ) then
    x = 0.
  end if

  do j = 1,bumps
    m = log(0.5)/log(bump_pos(j))
    chordwise(j) = (sin(pi*x**m)**t_b)
  enddo


  flagpanel = 0
  do j = 1,panels-1
    if ( y .ge. y_pos(j) .and. y .le. y_pos(j+1) ) then
      flagpanel(j) = 1
    end if
  enddo

  spanwise = 0.
  do j = 1,panels-1
    if ( flagpanel(1) == 1 ) then
      spanwise(1) = (y_pos(2)-y)/(y_pos(2)-y_pos(1))
      spanwise(2) = (y-y_pos(1))/(y_pos(2)-y_pos(1))
      exit
    elseif ( flagpanel(panels-1) == 1 ) then
      spanwise(panels) = (y-y_pos(panels-1))/(y_pos(panels)-y_pos(panels-1))
      spanwise(panels-1) = 1 + (y_pos(panels-1)-y)/(y_pos(panels)-y_pos(panels-1))
      exit
    elseif ( flagpanel(j) == 1 ) then
      spanwise(j) = 1 + (y_pos(j)-y)/(y_pos(j+1)-y_pos(j))
      spanwise(j+1) = (y-y_pos(j))/(y_pos(j+1)-y_pos(j))
    end if
  enddo

  sum = 0.
  p = 1
  do j = 1,panels
    do k = 1,bumps
      func(j,k) = spanwise(j)*chordwise(k)
      sum = sum + des_varsu(p)*func(j,k)
      p = p + 1
    enddo
  enddo

  dataunew(i) = datau(i)
  dataunew(i)%floats(3) = datau(i)%floats(3) + sum
  !if ( ISNAN(sum) ) then
  !  print *, 'upper', datau(i)%ints(1),x, datau(i)%floats(1), xtrans, scale
  !end if
enddo

do i = 1,nodl
  x = datal(i)%floats(1)
  y = datal(i)%floats(2)
  xtrans = rootlex + ((tiplex-rootlex)*(y-rooty)/(tipy-rooty))
  scale = rootchord + ((tiptex-roottex)*(y-rooty)/(tipy-rooty)) - ((tiplex-rootlex)*(y-rooty)/(tipy-rooty))
  x = x - xtrans
  x = x/scale
  if ( x .le. 0.0 ) then
    x = 0.
  end if

  do j = 1,bumps
    m = log(0.5)/log(bump_pos(j))
    chordwise(j) = (sin(pi*x**m)**t_b)
  enddo



  flagpanel = 0
  do j = 1,panels-1
    if ( y .ge. y_pos(j) .and. y .le. y_pos(j+1) ) then
      flagpanel(j) = 1
    end if
  enddo

  spanwise = 0.
  do j = 1,panels-1
    if ( flagpanel(1) == 1 ) then
      spanwise(1) = (y_pos(2)-y)/(y_pos(2)-y_pos(1))
      spanwise(2) = (y-y_pos(1))/(y_pos(2)-y_pos(1))
      exit
    elseif ( flagpanel(panels-1) == 1 ) then
      spanwise(panels) = (y-y_pos(panels-1))/(y_pos(panels)-y_pos(panels-1))
      spanwise(panels-1) = 1 + (y_pos(panels-1)-y)/(y_pos(panels)-y_pos(panels-1))
      exit
    elseif ( flagpanel(j) == 1 ) then
      spanwise(j) = 1 + (y_pos(j)-y)/(y_pos(j+1)-y_pos(j))
      spanwise(j+1) = (y-y_pos(j))/(y_pos(j+1)-y_pos(j))
    end if
  enddo

  sum = 0.
  p = 1
  do j = 1,panels
    do k = 1,bumps
      func(j,k) = spanwise(j)*chordwise(k)
      sum = sum + des_varsl(p)*func(j,k)
      p = p + 1
    enddo
  enddo


  datalnew(i) = datal(i)
  datalnew(i)%floats(3) = datal(i)%floats(3) + sum
  !if ( ISNAN(sum) ) then
  !  print *, 'lower', datal(i)%ints(1), x, datal(i)%floats(1), xtrans, scale
  !end if

enddo

do i = 1,nodf
  xref = rootflapx + (((dataf(i)%floats(2)-rooty)/(tipy-rooty))*(tipflapx-rootflapx))
  zref = rootflapz + (((dataf(i)%floats(2)-rooty)/(tipy-rooty))*(tipflapz-rootflapz))
  datafnew(i)%ints(1) = dataf(i)%ints(1)
  datafnew(i)%floats(1) = ((dataf(i)%floats(1) - xref)*cos(flapangle)) - ((dataf(i)%floats(3) - zref)*sin(flapangle)) + xref + xflap
  datafnew(i)%floats(2) = dataf(i)%floats(2)
  datafnew(i)%floats(3) = ((dataf(i)%floats(1) - xref)*sin(flapangle)) + ((dataf(i)%floats(3) - zref)*cos(flapangle)) + zref + zflap
enddo

print *,rootflapx,'rootflapx'
print *,rootflapz,'rootflapz'
print *,tipflapx,'tipflapx'
print *,tipflapz,'tipflapz'

open(1,file = 'new.txt')
do i = 1,nodu
  write(1,*) dataunew(i)
enddo
do i = 1,nodl
  write(1,*) datalnew(i)
enddo
do i = 1,nodf
  write(1,*) datafnew(i)
enddo
close(1)

open(1,file = 'mesh_disp.dat')
do i = 1,nodu
  write(1,*) datau(i)%ints(1), dataunew(i)%floats(1)-datau(i)%floats(1), dataunew(i)%floats(2)-datau(i)%floats(2), dataunew(i)%floats(3)-datau(i)%floats(3)
enddo
do i = 1,nodl
  write(1,*) datal(i)%ints(1), datalnew(i)%floats(1)-datal(i)%floats(1), datalnew(i)%floats(2)-datal(i)%floats(2), datalnew(i)%floats(3)-datal(i)%floats(3)
enddo
do i = 1,nodf
  write(1,*) dataf(i)%ints(1), datafnew(i)%floats(1)-dataf(i)%floats(1), datafnew(i)%floats(2)-dataf(i)%floats(2), datafnew(i)%floats(3)-dataf(i)%floats(3)
enddo
close(1)


end program param
