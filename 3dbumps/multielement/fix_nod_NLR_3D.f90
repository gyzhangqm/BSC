program gfn
  implicit none
  character(256) :: fixnod, input, cmd, bou ,geo,code
  integer :: flagupper, flaglower, flagflap, pos, ierr, pair(2), nbe, i, j ,k, l, counter, flag, nodu, nodl, nodf, bouall(6)
  integer, dimension(:), allocatable :: unodes, lnodes, fnodes


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

nodu = 0
nodl = 0
nodf = 0
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
  end if
enddo

allocate(unodes(nodu*3))
allocate(lnodes(nodl*3))
allocate(fnodes(nodf*3))
rewind(1)
read(1,*)
i = 1
j = 1
k = 1
do
  read(1,*,iostat = ierr) bouall
  if ( ierr .ne. 0 ) then
    exit
  end if
  if ( bouall(6) == flagupper ) then
    unodes(i) = bouall(3)
    unodes(i+1) = bouall(4)
    unodes(i+2) = bouall(5)
    i = i + 3
  elseif ( bouall(6) == flaglower ) then
    lnodes(j) = bouall(3)
    lnodes(j+1) = bouall(4)
    lnodes(j+2) = bouall(5)
    j = j + 3
  elseif ( bouall(6) == flagflap ) then
    fnodes(k) = bouall(3)
    fnodes(k+1) = bouall(4)
    fnodes(k+2) = bouall(5)
    k = k + 3
  end if
enddo

close(1)

call Unique1DArray_D(unodes)
call Unique1DArray_D(lnodes)
call Unique1DArray_D(fnodes)
nodu = size(unodes)
nodl = size(lnodes)
nodf = size(fnodes)

open(1,file=fixnod)
write(1,*) 'ON_NODES'
do i = 1,nodu
  write(1,*) unodes(i), flagupper
enddo
do i = 2,nodl
  write(1,*) lnodes(i), flaglower
enddo
do i = 2,nodf
  write(1,*) fnodes(i), flagflap
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
