! Program to implement Hicks-Henne Bump parameterization

program hickshenne
implicit none

integer, parameter :: t_b = 4
real(kind = 8), parameter :: pi = acos(-1.)
real, dimension(:), allocatable :: xu, ubase,unew, xl, lbase,lnew, file(:,:), gradu(:,:), gradl(:,:), dpa, dpb, bump_pos, dp(:,:)
integer :: i, j, Nx, Ny, io, n
real(kind = 8) :: suma, sumb, m, h, pair(2)

n = 0
open(1,file="design_parameters/dp.txt")

  do
    read(1,*,iostat=io) pair
    if (io/=0) exit
    n = n + 1
  end do
  !write(*,*) 'The number of design parameters: '
  !print *,n
  rewind(1)

  allocate(dp(n,2))
  do i=1,n
    read(1,*) dp(i,:)
  end do
  close(1)
  allocate(dpa(n))
  allocate(dpb(n))
  allocate(bump_pos(n))
  dpa = dp(:,1)
  dpb= dp(:,2)
  deallocate(dp)


Nx = 0
open(1,file="airfoils/naca0012upper.txt")

  do
    read(1,*,iostat=io) pair
    if (io/=0) exit
    Nx = Nx + 1
  end do
  !write(*,*) 'The number of coordinates on upper curve: '
  !print *,Nx
  rewind(1)

  allocate(file(Nx,2))

  do i=1,Nx
    read(1,*) file(i,:)
  end do

  close(1)
  allocate(xu(Nx))
  allocate(ubase(Nx))
  allocate(unew(Nx))
  allocate(gradu(Nx,n))
  xu = file(:,1)
  ubase = file(:,2)
  deallocate(file)

  Ny = 0
  open(1,file="airfoils/naca0012lower.txt")

    do
      read(1,*,iostat=io) pair
      if (io/=0) exit
      Ny = Ny + 1
    end do
    !write(*,*) 'The number of coordinates on lower curve: '
    !print *,Ny
    rewind(1)

    allocate(file(Ny,2))

    do i=1,Ny
      read(1,*) file(i,:)
    end do

    close(1)
    allocate(xl(Ny))
    allocate(lbase(Ny))
    allocate(lnew(Ny))
    allocate(gradl(Ny,n))
    xl = file(:,1)
    lbase = file(:,2)
h = 1./(n+1)
do i = 1,n
  bump_pos(i) = i*h
enddo


do i = 1,Nx
    suma = 0;
    do j = 1,n
        m = log(0.5)/log(bump_pos(j));
        suma = suma + dpa(j)*(sin(pi*xu(i)**m)**t_b);
        gradu(i,j) = sin(pi*xu(i)**m)**t_b;
    enddo
    unew(i) = ubase(i) + suma;
enddo
do i = 1,Ny
    sumb = 0;
    do j = 1,n
        m = log(0.5)/log(bump_pos(j));
        sumb = sumb + dpb(j)*(sin(pi*xl(i)**m)**t_b);
        gradl(i,j) = sin(pi*xl(i)**m)**t_b;
    enddo
    lnew(i) = lbase(i) + sumb;
enddo


open(unit = 1,file = "updatedu.txt")
do i=1,Nx
      write(1,*) xu(i), unew(i)
enddo
open(unit = 1,file = "updatedl.txt")
do i=1,Ny
      write(1,*) xl(i), lnew(i)
enddo
open(unit = 1,file = "grad.txt")
do i=1,Nx
    write(1,*) (gradu(i,j), j=1,n)
enddo
do i=1,Ny
    write(1,*) (gradl(i,j), j=1,n)
enddo


end program hickshenne
