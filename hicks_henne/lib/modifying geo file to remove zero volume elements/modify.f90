program param
  implicit none

  type intreal
    integer,dimension(1)   :: ints
    integer,dimension(4)      :: floats
  endtype intreal

  TYPE(intreal) :: geoinput

  integer :: nodu, nodl, i, j, k, ierr ,il, counterstart, counterend, counter, l, ndp, totnodes
  integer :: pair(2), flag, pos
  real, dimension(:), allocatable :: xu, yu, xl, yl, bump_pos, xunew, yunew, xlnew, ylnew, dispu, displ
  integer, dimension(:), allocatable :: unodes, lnodes
  real, dimension(:), allocatable :: dpa, dpb
  real (kind = 8) :: three(3), pairreal(2), h, scaleu, scalel, sum, m, xutrans, xltrans
  integer, parameter :: t_b = 4
  real(kind = 8), parameter :: pi = acos(-1.)
  character(256) :: input, cmd, geo, fixnod, code, geonew
  integer :: foo, flagupper, flaglower

  call get_command_argument(1, cmd)
  geo = trim(cmd)//".geo.dat"
  fixnod = trim(cmd)//".fix.nod"
  code = trim(cmd)//".codes"
  geonew = trim(cmd)//"_new.geo.dat"

  counter = 0
  open(1,file=geo,status='old',iostat=ierr)
  open(2,file=geonew)
  do while (ierr.eq.0)
    read(1,'(A)',iostat=ierr) input
    counter = counter + 1
    il=len_trim(input)
    if ( input(1:il) == '  ELEMENTS' ) then
      counterstart = counter + 1
    end if
    if ( input(1:il) == '  END_ELEMENTS' ) then
      counterend = counter - 1
    end if
  enddo
  totnodes = counterend - counterstart +1
  rewind(1)

  j = 1
  do l = 1,counterstart-1
    read(1,'(a)') input
    write(2,'(a)') trim(input)
  end do
  do l = counterstart,counterend
    read(1,*,iostat = ierr) geoinput
    if ( geoinput%ints(1) .eq. 117806) then
    elseif (geoinput%ints(1) .eq. 160101) then
    elseif (geoinput%ints(1) .eq. 256413) then
    elseif (geoinput%ints(1) .eq. 283919) then
    elseif (geoinput%ints(1) .eq. 320652) then
    elseif (geoinput%ints(1) .eq. 361450) then
    elseif (geoinput%ints(1) .eq. 402873) then
    elseif (geoinput%ints(1) .eq. 431422) then
    elseif (geoinput%ints(1) .eq. 461622) then
    else
      write(2,*) j, geoinput%floats
      j = j + 1
    end if
  end do
  do
    if(ierr .ne. 0) exit
    read(1,'(a)',iostat = ierr) input
    write(2,'(a)') trim(input)
  enddo
  close(1)
  close(2)



end program param
