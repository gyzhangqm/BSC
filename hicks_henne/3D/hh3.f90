! Program to implement Hicks-Henne Bump parameterization in 3D

program hickshenne
implicit none

integer, parameter :: t_b = 4, panels = 5, hh = 4
real(kind = 8), parameter :: pi = acos(-1.), span = 1196.3, tspan = 2734.5, xspan = 1578.7, sweep = 26.7*pi/180, leangle = 60*pi/180, teangle = 74.219120*pi/180
real, dimension(:), allocatable :: xu, ubase, xl, lbase, ubasepanel(:,:), lbasepanel(:,:), unewpanel(:,:), lnewpanel(:,:), xubasepanel(:,:), xlbasepanel(:,:), xunewpanel(:,:), xlnewpanel(:,:), file(:,:), dpas, dpbs
integer :: i, j, k, Nu, Nl, io, foo
real(kind = 8), dimension(panels) :: xpanelpos, ypanelpos, chordlen, scales, dispx, dispy, twist, twistnew
real(kind = 8), dimension(hh) :: bump_pos
real(kind = 8), dimension(panels,hh) :: dpa, dpb
real(kind = 8) :: sum, m, h, pair(2)

Nu = 0
open(1,file="airfoils/oneraup.txt")
do
  read(1,*,iostat=io) pair
  if(io/=0) exit
  Nu = Nu+1
enddo
rewind(1)
allocate(file(Nu,2))
do i=1,Nu
  read(1,*) file(i,:)
enddo
close(1)
allocate(xu(Nu))
allocate(ubase(Nu))
allocate(ubasepanel(panels,Nu))
allocate(unewpanel(panels,Nu))
allocate(xubasepanel(panels,Nu))
allocate(xunewpanel(panels,Nu))
xu = file(:,1)
ubase = file(:,2)
deallocate(file)

Nl = 0
open(1,file="airfoils/oneradown.txt")
do
  read(1,*,iostat=io) pair
  if(io/=0) exit
  Nl = Nl+1
enddo
rewind(1)
allocate(file(Nu,2))
do i=1,Nl
  read(1,*) file(i,:)
enddo
close(1)
allocate(xl(Nl))
allocate(lbase(Nl))
allocate(lbasepanel(panels,Nl))
allocate(lnewpanel(panels,Nl))
allocate(xlbasepanel(panels,Nl))
allocate(xlnewpanel(panels,Nl))
xl = file(:,1)
lbase = file(:,2)
deallocate(file)

foo = hh*panels
allocate(file(foo,2))
allocate(dpas(foo))
allocate(dpbs(foo))
open(1,file="dp/dp.txt")
do i=1,foo
  read(1,*) file(i,:)
  dpas(i) = file(i,1)
  dpbs(i) = file(i,2)
enddo
close(1)
deallocate(file)
do i = 1,panels
    do j = 1,hh
        dpa(i,j) = dpas((i-1)*hh + j)
        dpb(i,j) = dpbs((i-1)*hh + j)
    enddo
enddo

h = 1./(hh+1)
open(1,file="dp/dispx.txt")
open(2,file="dp/dispy.txt")
open(3,file="dp/scale.txt")
open(4,file="dp/twist.txt")
do i = 1,panels
  read(1,*) dispx(i)
  read(2,*) dispy(i)
  read(3,*) scales(i)
  read(4,*) twistnew(i)
  twist(i) = 0
  ypanelpos(i) = (span/(panels-1))*(i-1)
enddo
close(1)
close(2)
close(3)
close(4)
do i = 1,hh
  bump_pos(i) = 1*h
enddo


do i = 1,panels
    ! Baseline Config Generation start
    xpanelpos(i) = (tspan-ypanelpos(i))/tan(0.5*pi-sweep)
    chordlen(i) = ((tspan-ypanelpos(i))/tan(leangle)) - ((tspan-ypanelpos(i))/tan(teangle))
    !Scaling
    ubasepanel(i,:) = chordlen(i)*ubase
    lbasepanel(i,:) = chordlen(i)*lbase
    xubasepanel(i,:) = chordlen(i)*xu
    xlbasepanel(i,:) = chordlen(i)*xl
    !Translating
    xubasepanel(i,:) = xubasepanel(i,:) - (tspan-ypanelpos(i))/tan(leangle)
    xlbasepanel(i,:) = xlbasepanel(i,:) - (tspan-ypanelpos(i))/tan(leangle)
    !Twisting
    xubasepanel(i,:) = (xubasepanel(i,:) - xubasepanel(i,Nu))*cos(twist(i)*pi/180) + (ubasepanel(i,:) - ubasepanel(i,Nu))*sin(twist(i)*pi/180) + xubasepanel(i,Nu);
    ubasepanel(i,:) = -(xubasepanel(i,:) - xubasepanel(i,Nu))*sin(twist(i)*pi/180) + (ubasepanel(i,:) - ubasepanel(i,Nu))*cos(twist(i)*pi/180) + ubasepanel(i,Nu);
    xlbasepanel(i,:) = (xlbasepanel(i,:) - xlbasepanel(i,Nl))*cos(twist(i)*pi/180) + (lbasepanel(i,:) - lbasepanel(i,Nl))*sin(twist(i)*pi/180) + xlbasepanel(i,Nl);
    lbasepanel(i,:) = -(xlbasepanel(i,:) - xlbasepanel(i,Nl))*sin(twist(i)*pi/180) + (lbasepanel(i,:) - lbasepanel(i,Nl))*cos(twist(i)*pi/180) + lbasepanel(i,Nl);

    !Baseline Config Generation end

    !Updating

    !Untwisting
    xunewpanel(i,:) = (xubasepanel(i,:) - xubasepanel(i,Nu))*cos(-twist(i)*pi/180) + (ubasepanel(i,:) - ubasepanel(i,Nu))*sin(-twist(i)*pi/180) + xubasepanel(i,Nu);
    unewpanel(i,:) = -(xubasepanel(i,:) - xubasepanel(i,Nu))*sin(-twist(i)*pi/180) + (ubasepanel(i,:) - ubasepanel(i,Nu))*cos(-twist(i)*pi/180) + ubasepanel(i,Nu);
    xlnewpanel(i,:) = (xlbasepanel(i,:) - xlbasepanel(i,Nl))*cos(-twist(i)*pi/180) + (lbasepanel(i,:) - lbasepanel(i,Nl))*sin(-twist(i)*pi/180) + xlbasepanel(i,Nl);
    lnewpanel(i,:) = -(xlbasepanel(i,:) - xlbasepanel(i,Nl))*sin(-twist(i)*pi/180) + (lbasepanel(i,:) - lbasepanel(i,Nl))*cos(-twist(i)*pi/180) + lbasepanel(i,Nl);

    !Translating back
    xunewpanel(i,:) = xunewpanel(i,:) + (tspan-ypanelpos(i))/tan(leangle)
    xlnewpanel(i,:) = xlnewpanel(i,:) + (tspan-ypanelpos(i))/tan(leangle)

    !Scaling down
    unewpanel(i,:) = unewpanel(i,:)/chordlen(i)
    lnewpanel(i,:) = lnewpanel(i,:)/chordlen(i)
    xunewpanel(i,:) = xunewpanel(i,:)/chordlen(i)
    xlnewpanel(i,:) = xlnewpanel(i,:)/chordlen(i)

    !Hicks-Henne Update
    do k = 1,Nu
        sum = 0
        do j = 1,hh
            m = log(0.5)/log(bump_pos(j))
            sum = sum + dpa(i,j)*(sin(pi*xunewpanel(i,k)**m)**t_b)
            !gradu(i,j) = sin(pi*xunewpanel(i,k)^m)^t_b;
        enddo
        unewpanel(i,k) = unewpanel(i,k) + sum;
    enddo
    do k = 1,Nl
        sum = 0
        do j = 1,hh
            m = log(0.5)/log(bump_pos(j))
            sum = sum + dpb(i,j)*(sin(pi*xlnewpanel(i,k)**m)**t_b)
            !gradu(i,j) = sin(pi*xunewpanel(i,k)^m)^t_b
        enddo
        lnewpanel(i,k) = lnewpanel(i,k) + sum
    enddo

    !Scaling
    unewpanel(i,:) = scales(i)*unewpanel(i,:)
    lnewpanel(i,:) = scales(i)*lnewpanel(i,:)
    xunewpanel(i,:) = scales(i)*xunewpanel(i,:)
    xlnewpanel(i,:) = scales(i)*xlnewpanel(i,:)

    !Translating
    xunewpanel(i,:) = xunewpanel(i,:) - dispx(i)
    xlnewpanel(i,:) = xlnewpanel(i,:) - dispx(i)
    unewpanel(i,:) = unewpanel(i,:) + dispy(i)
    lnewpanel(i,:) = lnewpanel(i,:) + dispy(i)

    !twisting
    xunewpanel(i,:) = (xunewpanel(i,:) - xunewpanel(i,Nu))*cos(twistnew(i)*pi/180) + (unewpanel(i,:) - unewpanel(i,Nu))*sin(twistnew(i)*pi/180) + xunewpanel(i,Nu);
    unewpanel(i,:) = -(xunewpanel(i,:) - xunewpanel(i,Nu))*sin(twistnew(i)*pi/180) + (unewpanel(i,:) - unewpanel(i,Nu))*cos(twistnew(i)*pi/180) + unewpanel(i,Nu);
    xlnewpanel(i,:) = (xlnewpanel(i,:) - xlnewpanel(i,Nl))*cos(twistnew(i)*pi/180) + (lnewpanel(i,:) - lnewpanel(i,Nl))*sin(twistnew(i)*pi/180) + xlnewpanel(i,Nl);
    lnewpanel(i,:) = -(xlnewpanel(i,:) - xlnewpanel(i,Nl))*sin(twistnew(i)*pi/180) + (lnewpanel(i,:) - lnewpanel(i,Nl))*cos(twistnew(i)*pi/180) + lnewpanel(i,Nl);



enddo

open(unit=1,file="updatedu.txt")
open(unit=2,file="updatedl.txt")
do i = 1,panels
  do j = 1,Nu
    write(1,*) xunewpanel(i,j), unewpanel(i,j), ypanelpos(i)
    write(2,*) xlnewpanel(i,j), lnewpanel(i,j), ypanelpos(i)
  enddo
enddo
close(1)
close(2)




end program hickshenne
