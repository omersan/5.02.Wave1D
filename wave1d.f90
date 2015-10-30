!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Methods for wave equation (1D advection equation)
!     du/dt + a*du/dx = 0, where a=1.0d0
!     Periodic bc
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) 
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 13, 2015
!-----------------------------------------------------------------------------!

program wave1d
implicit none
integer::i,k,nx,nt,ns,nf,np
real*8 ::dx,dt,x0,xL,pi,t,Tmax,a,c,xx,uu
real*8,allocatable ::u(:),x(:)

!Domain
x0 = 0.0d0 !left
xL = 1.0d0 !right

!number of points
nx = 40

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

!maximum time desired
Tmax = 1.0d0

!convective constant:
a = 1.0d0

!CFL number
c = 0.5d0

!time step
dt = c*dx/a

!number of points in time
nt = nint(Tmax/dt)

!number of snapshots to plot
ns = nt

!frequency for plotting
nf = nint(dfloat(nt)/dfloat(ns))

!u: convected variable 
allocate(u(0:nx))

!initial condition
pi = 4.0d0*datan(1.0d0)
t = 0.0d0
do i=0,nx
u(i) = dsin(2.0d0*pi*x(i))
end do


!Plot initial condition
open(18,file='u.plt')
write(18,*) 'variables ="x","u"'
write(18,100)'zone f=point i=',nx+1,',t="t =',t,'"'
do i=0,nx
write(18,*) x(i),u(i)
end do


!time integration
do k=1,nt
  
    !FTCS scheme (unconditionally unstable)
    !call FTCS(nx,dx,dt,a,u)
       
    !upwind scheme
    !call upwind(nx,dx,dt,a,u)
    
    !Lax-Wendroff scheme
    !call LW(nx,dx,dt,a,u)

    !Lax scheme
    !call Lax(nx,dx,dt,a,u)
   
   	!MacCormack scheme
    !call MC(nx,dx,dt,a,u)
     
    !RK4 scheme
    call RK4(nx,dx,dt,a,u)
         
    !update t
    t = t+dt 
 


    !plot field
    if (mod(k,nf).eq.0) then
	write(18,100)'zone f=point i=',nx+1,',t="t =',t,'"'
	do i=0,nx
	write(18,*) x(i),u(i)
	end do
    end if

    
	print*,k,t,maxval(u)


end do
  
close(18)


!Plot at t=Tmax
open(12, file="numerical.plt")
write(12,*)'variables ="x","u"'
do i=0,nx
	write(12,*) dfloat(i)*dx,u(i)
end do
close(12)
   
! Writing exact solution using np points
np = 2000 !use 2000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="x","u"'
	do i=0,np
		xx = dfloat(i)/dfloat(np)
		uu = dsin(2.0d0*pi*(xx-a*t))
		write(12,*) xx,uu
	end do
close(12)




100 format(a16,i8,a10,f8.4,a3)
end


!-----------------------------------------------------------------------------!
!FTCS scheme
!-----------------------------------------------------------------------------!
subroutine FTCS(nx,dx,dt,a,u)
implicit none
integer::nx,i
real*8 ::dx,dt,a,c
real*8 ::u(0:nx),v(-1:nx+1)

c = a*dt/dx

do i=0,nx
v(i) = u(i)
end do
v(-1) = v(nx-1)   !periodic
v(nx+1) = v(1)    !periodic

do i=0,nx
	u(i) = v(i) - c*(v(i+1)-v(i-1))  
end do

end 

!-----------------------------------------------------------------------------!
!upwind scheme
!-----------------------------------------------------------------------------!
subroutine upwind(nx,dx,dt,a,u)
implicit none
integer::nx,i
real*8 ::dx,dt,a,c
real*8 ::u(0:nx),v(-1:nx+1)

c = a*dt/dx

do i=0,nx
v(i) = u(i)
end do
v(-1) = v(nx-1)   !periodic
v(nx+1) = v(1)    !periodic

do i=0,nx
	u(i) = v(i) - c*(v(i)-v(i-1))  
end do

end 


!-----------------------------------------------------------------------------!
!Lax-Wendroff scheme
!-----------------------------------------------------------------------------!
subroutine LW(nx,dx,dt,a,u)
implicit none
integer::nx,i
real*8 ::dx,dt,a,c
real*8 ::u(0:nx),v(-1:nx+1)

c = a*dt/dx

do i=0,nx
v(i) = u(i)
end do
v(-1) = v(nx-1)   !periodic
v(nx+1) = v(1)    !periodic

do i=0,nx
	u(i) = v(i) - 0.5d0*c*(v(i+1)-v(i-1)) + 0.5d0*c*c*(v(i+1)-2.0d0*v(i)+v(i-1))  
end do

end 

!-----------------------------------------------------------------------------!
!Lax scheme
!-----------------------------------------------------------------------------!
subroutine Lax(nx,dx,dt,a,u)
implicit none
integer::nx,i
real*8 ::dx,dt,a,c
real*8 ::u(0:nx),v(-1:nx+1)

c = a*dt/dx

do i=0,nx
v(i) = u(i)
end do
v(-1) = v(nx-1)   !periodic
v(nx+1) = v(1)    !periodic

do i=0,nx
	u(i) = 0.5d0*(v(i+1)+v(i-1)) - 0.5d0*c*(v(i+1)-v(i-1))  
end do

end 

!-----------------------------------------------------------------------------!
!MacCormack scheme
!-----------------------------------------------------------------------------!
subroutine MC(nx,dx,dt,a,u)
implicit none
integer::nx,i
real*8 ::dx,dt,a,c
real*8 ::u(0:nx),v(-1:nx+1),y(-1:nx+1)

c = a*dt/dx

do i=0,nx
v(i) = u(i)
end do
v(-1) = v(nx-1)   !periodic
v(nx+1) = v(1)    !periodic

!Predictor
do i=0,nx
	y(i) = v(i) - c*(v(i+1)-v(i))  
end do
	y(-1) = y(nx-1) !periodic

!corrector
do i=0,nx
	u(i) = 0.5d0*(u(i)+y(i)) - 0.5d0*c*(y(i)-y(i-1))  
end do

end 

!-----------------------------------------------------------------------------!
!Runge-Kutta 4 scheme + 2nd-order central scheme
!-----------------------------------------------------------------------------!
subroutine RK4(nx,dx,dt,a,u)
implicit none
integer::nx,i
real*8 ::dx,dt,a
real*8 ::u(0:nx)
real*8 ::k1(0:nx),k2(0:nx),k3(0:nx),k4(0:nx),r(0:nx)

    call RHS(nx,dx,a,u,r)
		do i=0,nx
		k1(i) = dt*r(i)      
    	end do
	
    call RHS(nx,dx,a,u+k1/2.0d0,r)
		do i=0,nx
		k2(i) = dt*r(i)      
    	end do

    call RHS(nx,dx,a,u+k2/2.0d0,r)
		do i=0,nx
		k3(i) = dt*r(i)      
    	end do

    call RHS(nx,dx,a,u+k3,r)
		do i=0,nx
		k4(i) = dt*r(i)      
    	end do
   
do i=0,nx
u(i) = u(i) + (k1(i)+2.0d0*(k2(i)+k3(i))+k4(i))/6.0d0     
end do

end 


!-----------------------------------------------------------------------------!
!Right Hand Side (RHS) for RK scheme
!2nd order central scheme
!4th order central scheme
!-----------------------------------------------------------------------------!
subroutine RHS(nx,dx,a,u,r)
implicit none
integer::nx,i
real*8 ::dx,a
real*8 ::u(0:nx),r(0:nx),v(-2:nx+2)

do i=0,nx
v(i) = u(i)
end do
v(-1) = v(nx-1)   !periodic
v(-2) = v(nx-2)   !periodic
v(nx+1) = v(1)    !periodic
v(nx+2) = v(2)    !periodic

do i=0,nx
	!r(i) = - a*(v(i+1)-v(i-1))/(2.0d0*dx)
    r(i) = - a*(-v(i+2)+8.0d0*v(i+1)-8.0d0*v(i-1)+v(i-2))/(12.0d0*dx)
end do

end 










