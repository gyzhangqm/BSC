program gfn
  implicit none
  character(256) :: fixnod, input, cmd, bou ,geo,code
  integer :: flagupper, flaglower, flagflap, pos, ierr, pair(2), nbe, i, j ,k, counter, flag
  integer, dimension(:), allocatable :: boudata(:,:), geoboudata(:,:)


  call get_command_argument(1, cmd)
  geo = trim(cmd)//".geo.dat"
  fixnod = trim(cmd)//".fix.nod"
  code = trim(cmd)//".codes"
  bou = trim(cmd)//".fix.bou"
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
!------- Read .fix.bou file ----------------------------------------------------

nbe = 0
open(1, file = bou, iostat = ierr)
read(1,*)
do
  read(1,*,iostat = ierr) pair
  if( ierr .ne. 0) exit
  nbe = nbe + 1
enddo
close(1)
allocate(boudata(nbe,2))
open(1, file = bou, iostat = ierr)
read(1,*)
do k = 1,nbe
  read(1,*) boudata(k,:)
enddo


!-------------------------------------------------------------------------------
!------------------------- Read .geo.dat ---------------------------------------
counter = 0
allocate(geoboudata(nbe,4))
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

counter = 0
open(1,file=fixnod)
write(1,*) 'ON_NODES'
do k = 1,nbe
  if ( boudata(k,2) == flagupper ) then
    write(1,*) geoboudata(k,2), flagupper
    counter = counter + 1
  end if
  if ( boudata(k,2) == flaglower ) then
    write(1,*) geoboudata(k,2), flaglower
    counter = counter + 1
  end if
  if ( boudata(k,2) == flagflap ) then
    write(1,*) geoboudata(k,2), flagflap
    counter = counter + 1
  end if
enddo
write(1,*) 'END_ON_NODES'
close(1)


end program gfn
