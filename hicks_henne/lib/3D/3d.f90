program param
  implicit none

  type intreal
    integer,dimension(1)   :: ints
    real,dimension(3)      :: floats
  endtype intreal
  type(intreal) :: geoinput, temp
  type(intreal), dimension(:), allocatable :: panelgeoinput


  integer, parameter :: t_b = 4
  real(kind = 8), parameter :: pi = acos(-1.)
  character(256) :: input, cmd, geo, fixnod, code
  integer :: flagupper, flaglower, pos, nodu, nodl,  i, j, k, l, ierr, counterstart, counterend, counter, totnodes, pair(2), panels
  integer :: nx, totalpanels
  integer, dimension(:), allocatable :: unodes, lnodes, panelflag
  real, dimension(:), allocatable :: zpos, zpanelpos
  real(kind = 8) :: buffer

  type data
    integer, dimension(:), allocatable :: nodenos
    real, dimension(:), allocatable :: x,y,z
  endtype data
  type(data) :: datau, datal


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
      if ( zpanelpos(j) .lt. zpanelpos(i) ) then
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
  do while (ierr .eq. 0)
    j = 1
    read(1,*,iostat=ierr) geoinput
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
  do while (ierr .eq. 0)
    j = 1

    read(1,*,iostat=ierr) geoinput
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


call execute_command_line('rm dumpallu.txt')
call execute_command_line('rm dumpalll.txt')

do i = 1,totalpanels
  panelflag(i) = 0
  do j = 1,panels
    if ( zpanelpos(i) == zpos(j) ) then
      panelflag(i) = 1
    end if
  enddo
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
        if ( panelgeoinput(l)%floats(1) .le. panelgeoinput(k)%floats(1) ) then
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

open(1,file = 'dumpall_l.txt', iostat=ierr)
open(2,file = 'dumpalll.txt')

do i = 1,totalpanels
  do j = 1,nx
    read(1,*,iostat=ierr) panelgeoinput(j)
  enddo

  do k = 1,nx
    do l = 1,nx
      if ( k == l ) then
      else
        if ( panelgeoinput(l)%floats(1) .le. panelgeoinput(k)%floats(1) ) then
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

call execute_command_line('rm dumpall_u.txt')
call execute_command_line('rm dumpall_l.txt')





end program param
