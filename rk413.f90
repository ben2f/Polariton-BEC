!-------------------------------------------------------------------------------------------------------------------------------------------
!										  RK4 METHOD FOR 1ST ORDER DIFFERENTIAL EQUATION
!-------------------------------------------------------------------------------------------------------------------------------------------

function f1(x,y,z,t,gme,g,g1,omega)
	implicit none
	double precision::x,y,z,f1,t,gme,omega,g1,g
	f1=2*(gme-z)*x-2.5D0/sqrt(2D0)*g1*(x*x)
end function f1

function f2(x,y,z,t,gme,g,g1,omega)
	implicit none
	double precision::x,y,z,f2,t,g1,g,omega,gme
	f2=4D0*z*y+(g1*x*y)/sqrt(2D0)
end function f2

function f3(x,y,z,t,gme,g,g1,omega)
	implicit none
	double precision::x,y,z,f3,t,g,omega,g1,gme,dt
	dt=1D0
	f3=-0.5D0*(omega*omega)*dt+0.5D0/(y*y)+(g*x)/(2D0*sqrt(2D0)*y)-2D0*z*z
end function f3

!function f4(x,y,z,t,gme,g,g1,omega)
!	implicit none
!	double precision::x,y,z,f4,t,g1,g,omega,gme,w
!	w=(5D0/(4D0*sq!!rt(2D0)))
!	f4=-1D0/(2D0*(y**2.0))-w*g*(x**2.0)
!end function f4

program rk4
	implicit none
	
	double precision::f1,f2,f3,t,t2,tf,x,y,z,h,y2,x2,z2
	double precision::k,k1,k2,k3,k4,l,l1,l2,l3,l4,m,m1,m2,m3,m4
	double precision::gme,g1,g,omega,D,c,d1
	
	g=5.0D0
	g1=0.01D0
	omega=0.01D0
	
	gme=0.001D0

	open(16,file="LS002.txt")
	open(17,file="initialb.txt")
	
	t=0.0D0
	tf=2500.0D0
	h=0.001D0
	
	x=sqrt(2.0)*(gme/g1)
	z=-gme/4.0
	c=(gme*g1)**2.0+(g*gme)**2.0+4*g1*g1*omega*omega
    d1=g1*(gme**2)+4*g1*omega*omega
    y=(2*gme*g + 2*sqrt(c))/d1
		
	write(17,*) "INITIAL PARAMETERS"
	write(17,*)"t	","Amplitude(A)	","	Width(R)	","Chirp(beta)","	Final t	","	dt	"
	write(17,*) t,"					",sqrt(x),"						",sqrt(y),"					",z,"				",tf,"		",h,"	"
	
	do while (abs(t)<abs(tf))
		
		!print*,t,x,y,z
		
		k1=h*f1(x,y,z,t,gme,g,g1,omega)
		l1=h*f2(x,y,z,t,gme,g,g1,omega)
		m1=h*f3(x,y,z,t,gme,g,g1,omega)
		
		t2=t+0.5D0*h
		
		x2=x+0.5D0*k1
		y2=y+0.5D0*l1
		z2=z+0.5D0*m1
		
		
		k2=h*f1(x2,y2,z2,t2,gme,g,g1,omega)
		l2=h*f2(x2,y2,z2,t2,gme,g,g1,omega)
		m2=h*f3(x2,y2,z2,t2,gme,g,g1,omega)
		
		x2=x+k2*0.5D0
		y2=y+l2*0.5D0
		z2=z+m2*0.5D0
		
		k3=h*f1(x2,y2,z2,t2,gme,g,g1,omega)
		l3=h*f2(x2,y2,z2,t2,gme,g,g1,omega)
		m3=h*f3(x2,y2,z2,t2,gme,g,g1,omega)
		
		t2=t+h
		x2=x+k3
		y2=y+l3
		z2=z+m3
		
		k4=h*f1(x2,y2,z2,t2,gme,g,g1,omega)
		l4=h*f2(x2,y2,z2,t2,gme,g,g1,omega)
		m4=h*f3(x2,y2,z2,t2,gme,g,g1,omega)
		
		D=1D0/6D0
		k=D*(k1+2*k2+2*k3+k4)
		l=D*(l1+2*l2+2*l3+l4)
		m=D*(m1+2*m2+2*m3+m4)
		
		x=x+k
		y=y+l
		z=z+m
		t=t2
		
		
		write(16,*)t,sqrt(x),sqrt(y)
		
	end do
	!write(17,*) t,"	",x,"			",y,"			",z,"	",tf,"		",h,"	"
	
	close(17)
	close(16)
	print*,"COMPLETED"
end program
	