program mc
  implicit none

  type intreal
    integer,dimension(1)   :: ints
    real,dimension(3)      :: floats
  endtype intreal
  type(intreal) :: geoinput


  integer, parameter :: zprecision = 6
  character(256) :: input, cmd, geo, fixnod, code
  integer :: flagupper, flaglower, pos, nodu, nodl,  i, j, k, l, ierr ,il, counterstart, counterend, counter, totnodes, pair(2), buffer, flag
  integer, dimension(:), allocatable :: unodes, lnodes
  real (kind = 8) :: bufferreal, epsilon


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
!print *,nodu
!print *,nodl

!-------------------------------------------------------------------------------
!Extracts the coordinates of the grid points which lie on the upper and lower surface
!---- and multiplyig z coordinate with 10^precision ----------------------------

counter = 0
i = 1
j = 1

open(2,file = 'dump/dumpall.txt')
open(1,file=geo,status='old',iostat=ierr)
open(3,file = 'corrected.geo.dat')
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
  read(1,'(a)') input
  write(3,*) trim(input)
end do
do l = counterstart,counterend
  flag = 0
  read(1,*) geoinput
  do k = 1,nodu
      if ( unodes(k) == geoinput%ints(1) ) then
        buffer = geoinput%floats(3)*10**zprecision
        bufferreal = buffer
        bufferreal = bufferreal/10**zprecision
        write(2,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), bufferreal
        write(3,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), bufferreal
        flag = 1
      end if
    enddo
    do k = 1,nodl
      if ( lnodes(k) == geoinput%ints(1) ) then
        buffer = geoinput%floats(3)*10**zprecision
        bufferreal = buffer
        bufferreal = bufferreal/10**zprecision
        write(2,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), bufferreal
        write(3,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), bufferreal
        flag = 1
      end if
    enddo
    if ( flag == 0 ) then
      write(3,*) geoinput%ints(1), geoinput%floats(1), geoinput%floats(2), geoinput%floats(3)
    end if
end do

ierr = 0
do while (ierr.eq.0)
  read(1,'(A)',iostat=ierr) input
  write(3,*) trim(input)
enddo


close(1)
close(2)
close(3)
totnodes = counterend - counterstart +1

cmd = 'rm '//geo
call execute_command_line(cmd)
cmd = 'mv corrected.geo.dat '//geo
call execute_command_line(cmd)





end program mc
