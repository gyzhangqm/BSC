program gfn
  implicit none
  character(256) :: fixnod, input, cmd, bou ,geo,code
  integer :: flagupper, flaglower, pos, ierr, pair(2), nbe, i, j ,k, counter, flag, nodu, nodl
  integer, dimension(:), allocatable :: boudata(:,:), geoboudata(:,:), boudataupper(:,:), boudatalower(:,:), unodes, lnodes


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

i = 0
j = 0

do k = 1,nbe
  if ( boudata(k,2) == flagupper ) then
    i = i + 1
  end if
  if ( boudata(k,2) == flaglower ) then
    j = j + 1
  end if
enddo
nodu = i
nodl = j
allocate(boudataupper(i,5))
allocate(boudatalower(j,5))

i = 1
j = 1
do k = 1,nbe
  if ( boudata(k,2) == flagupper ) then
    boudataupper(i,:) = geoboudata(k,:)
    i = i + 1
  end if
  if ( boudata(k,2) == flaglower ) then
    boudatalower(j,:) = geoboudata(k,:)
    j = j + 1
  end if
enddo
allocate(unodes(nodu*3))
allocate(lnodes(nodl*3))
k = 1
do i = 1, nodu
  unodes(k) = boudataupper(i,2)
  unodes(k+1) = boudataupper(i,3)
  unodes(k+2) = boudataupper(i,4)
  k = k + 3
enddo
k = 1
do i = 1, nodl
  lnodes(k) = boudatalower(i,2)
  lnodes(k+1) = boudatalower(i,3)
  lnodes(k+2) = boudatalower(i,4)
  k = k + 3
enddo

call Unique1DArray_D(unodes)
call Unique1DArray_D(lnodes)
nodu = size(unodes)
nodl = size(lnodes)

open(1,file=fixnod)
write(1,*) 'ON_NODES'
do i = 1,nodu
  write(1,*) unodes(i), flagupper
enddo
do i = 2,nodl
  write(1,*) lnodes(i), flaglower
enddo
write(1,*) 'END_NODES'
close(1)



contains
  subroutine Unique1DArray_D(Arr_a)
    implicit none
    integer,dimension(:),allocatable::Arr_a,Arr_b
    logical,dimension(:), allocatable::mask
    integer,dimension(:),allocatable::index_vector,indexSos
    integer::i,j,num
    num=size(Arr_a);  allocate(mask(num)); mask = .true.
    do i=num,2,-1
      mask(i)=.not.(any(Arr_a(:i-1)==Arr_a(i)))
    end do
    allocate(indexSos(size(PACK([(i,i=1,num)],mask))))
    allocate(index_vector(size(indexSos))); index_vector=PACK([(i,i=1,num)],mask)
    allocate(Arr_b(size(index_vector)))
    Arr_b=Arr_a(index_vector)
    call move_alloc(Arr_b,Arr_a)
  end subroutine Unique1DArray_D





end program gfn
