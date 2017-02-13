program gsd
  implicit none
  character(256) :: fixnod, input, cmd, bou ,geo,code, set
  integer :: flagupper, flaglower, flagflap,flagte, flagtip, pos, ierr, pair(2), nbe, i, j ,k, l, m, n, counter, flag, nodu, nodl, nodf, nodte, nodtip, bouall(6)
  integer, dimension(:), allocatable :: unodes, lnodes, fnodes, geoboudata(:,:), fixboudata(:,:)
  call get_command_argument(1, cmd)
  geo = trim(cmd)//".geo.dat"
  fixnod = trim(cmd)//".fix.nod"
  code = trim(cmd)//".codes"
  bou = trim(cmd)//".fix.bou"
  set = trim(cmd)//".set.dat"

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
    read(1,'(a)') input
    input = trim(input)
    pos = index(input,":")
    read(input(pos+1:), *) flagte
    read(1,'(a)') input
    input = trim(input)
    pos = index(input,":")
    read(input(pos+1:), *) flagtip
    close(1)

    !-------------------------------------------------------------------------------
    !------- Read .fix.bou file ----------------------------------------------------

    nodu = 0
    nodl = 0
    nodf = 0
    nodte = 0
    nodtip = 0
    open(1, file = bou, iostat = ierr)
    read(1,*)
    do
      read(1,*,iostat = ierr) bouall
      if ( ierr .ne. 0 ) then
        exit
      end if
      if ( bouall(6) == flagupper ) then
        nodu = nodu + 1
      elseif ( bouall(6) == flaglower) then
        nodl = nodl + 1
      elseif ( bouall(6) == flagflap) then
        nodf = nodf + 1
      elseif ( bouall(6) == flagte) then
        nodte = nodte + 1
      elseif ( bouall(6) == flagtip) then
        nodtip = nodtip + 1
      end if
    enddo
    close(1)

    !print *,nodu,'nodu'
    !print *,nodl,'nodl'
    !print *,nodf,'nodf'
    !print *,nodte,'nodte'
    !print *,nodtip,'nodtip'

    nbe = 0
    open(1, file = bou, iostat = ierr)
    read(1,*)
    do
      read(1,*,iostat = ierr) bouall
      if( ierr .ne. 0) exit
      nbe = nbe + 1
    enddo
    allocate(fixboudata(nodu+nodl+nodf+nodte+nodtip,6))

    rewind(1)
    i = 1
    read(1,*)
    do k = 1,nbe
      read(1,*) bouall
      if ( bouall(6) == flagupper .or. bouall(6) == flaglower .or. bouall(6) == flagflap .or. bouall(6) == flagte .or. bouall(6) == flagtip ) then
        fixboudata(i,:) = bouall
        i = i + 1
      end if
    enddo
    close(1)

    counter = 0
    allocate(geoboudata(nbe,5))
    open(1, file = geo, iostat = ierr)
    do
      read(1,'(a)',iostat = ierr) input
      if( ierr .ne. 0) exit
      counter = counter + 1
      if ( trim(input) == '  BOUNDARIES, ELEMENT' ) then
        exit
      end if
    enddo
    rewind(1)
    do k = 1,counter
      read(1,*)
    enddo
    do k = 1,nbe
      read(1,*) geoboudata(k,:)
    enddo
    close(1)

!-------------------------------------------------------------------------------
open(1,file=set)
write(1,*) 'START_BOUNDARIES'
do i = 1,nodu+nodl+nodf+nodte+nodtip
  do j = 1,nbe
    if ( fixboudata(i,3) == geoboudata(j,2) .and. fixboudata(i,4) == geoboudata(j,3) .and. fixboudata(i,5) == geoboudata(j,4)) then
      write(1,*) geoboudata(j,1), fixboudata(i,6)
      exit
    end if
  enddo
enddo
write(1,*) 'END_BOUNDARIES'
close(1)

end program gsd
