program gfn
  implicit none
  character(256) :: fixnod, input, cmd, bou ,geo,code
  integer :: flagupper, flaglower, pos, ierr, pair(2), nbe, i, j ,k
  integer :: counter, flag, foilcode, totnodes
  integer, dimension(:), allocatable :: boudata(:,:), geoboudata(:,:),nodes
  real, dimension(:), allocatable :: geocoorddata(:,:), sample(:,:)


  call get_command_argument(1, cmd)
  geo = trim(cmd)//".geo.dat"
  fixnod = trim(cmd)//".fix.nod"
  code = trim(cmd)//".codes"
  bou = trim(cmd)//".fix.bou"
!-------------------------------------------------------------------------------
!-------- Reading the flags for upper and lower surfaces from .codes file ------
  open(1,file = code)
  read(1,'(a)') input
  read(1,'(a)') input
  pos = index(trim(input),":")
  read(input(pos+1:), *) foilcode
  read(1,'(a)') input
  input = trim(input)
  pos = index(input,":")
  read(input(pos+1:), *) flagupper
  read(1,'(a)') input
  input = trim(input)
  pos = index(input,":")
  read(input(pos+1:), *) flaglower

!-------------------------------------------------------------------------------
!------- Read .fix.bou file ----------------------------------------------------

nbe = 0
open(1, file = bou, iostat = ierr)

do
  read(1,*,iostat = ierr) pair
  if( ierr .ne. 0) then
    exit
  else
    nbe = nbe + 1
  endif
enddo
close(1)
allocate(boudata(nbe,2))
open(1, file = bou, iostat = ierr)
do k = 1,nbe
  read(1,*) boudata(k,:)
enddo


!-------------------------------------------------------------------------------
!------------------------- Read .geo.dat ---------------------------------------
counter = 0
allocate(geoboudata(nbe,3))
open(1, file = geo, iostat = ierr)
do
  read(1,'(a)',iostat = ierr) input
  if( ierr .ne. 0) exit
  counter = counter + 1
  if ( trim(input) == ' boundaries' ) then
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
do k = 1,nbe
  if ( boudata(k,2) == foilcode ) then
    counter = counter + 2
  end if
enddo
allocate(nodes(counter))
counter = 1
do k = 1,nbe
  if ( boudata(k,2) == foilcode ) then
    nodes(counter) = geoboudata(k,2)
    nodes(counter+1) = geoboudata(k,3)
    counter = counter + 2
  end if
enddo

call Unique1DArray_D(nodes)

!------------------------------------------------------------------------------
!----- Reading coordinates to find if node lies on upper surface or lower -----
allocate(sample(1,3))
counter = 0
open(1, file = geo, iostat = ierr)
do
  read(1,'(a)',iostat = ierr) input

  counter = counter + 1
  if ( trim(input) == ' COORDINATES' ) then
    exit
  end if
enddo
rewind(1)
totnodes = 0
do k = 1,counter
  read(1,*)
enddo
do
  read(1,*,iostat = ierr) sample
  if( ierr .ne. 0) then
    exit
  else
    totnodes = totnodes + 1
  endif
enddo
rewind(1)
allocate(geocoorddata(totnodes,3))
do k = 1,counter
  read(1,*)
enddo
do k = 1,totnodes
  read(1,*) geocoorddata(k,:)
enddo
close(1)

open(1,file=fixnod)
write(1,*) 'ON_NODES'
do k = 1,size(nodes)
  if ( geocoorddata(nodes(k),3) .gt. 0.0 ) then
    write(1,*) nodes(k), flagupper
  end if
  if ( geocoorddata(nodes(k),3) .lt. 0.0  ) then
    write(1,*) nodes(k), flaglower
  end if
enddo
write(1,*) 'END_ON_NODES'
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
    allocate(index_vector(size(indexSos)));
    index_vector=PACK([(i,i=1,num)],mask)
    allocate(Arr_b(size(index_vector)))
    Arr_b=Arr_a(index_vector)
    call move_alloc(Arr_b,Arr_a)
  end subroutine Unique1DArray_D



end program gfn
