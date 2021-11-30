from math import*
from numpy import*
from matplotlib.pyplot import*
from mpl_toolkits import mplot3d


xmin = -float(input("Donner les dim de l'écran l= ")); l=-xmin
ymin = -float(input("Donner les dim de l'écran L= ")); L=-ymin
lambd = eval(input("Donner la longueur d'onde lambda= "))
Nb = 2500
Rep = int(input('Combien de fente(s)? '))
f = eval(input('Donner la distance f= '))

Incx0 = float(input("Angle d'incidence en x (si pas d'incidence x = 0): ")); Incx0 = Incx0*180/pi
Incy0 = float(input("Angle d'incidence en y (si pas d'incidence y = 0): ")); Incy0 = Incy0*180/pi

x = linspace(xmin,abs(xmin),Nb)
y = linspace(xmin,abs(ymin),Nb)
t = linspace(0,2*pi,Nb)

def TFD2c(F,X,Y):
	def f(X,Y):
		return eval(F)
	Fonc = abs(fft.fft2(f(X,Y)))
	Ny = int(len(Fonc[:,0])); Nx = int(len(Fonc[0,:]))
	P = zeros((Nx,Ny))
	Carre1 = Fonc[0:int(Nx/2),0:int(Ny/2)]
	Carre2 = Fonc[0:int(Nx/2),int(Ny/2):Ny]
	Carre3 = Fonc[int(Nx/2):Nx,0:int(Ny/2)]
	Carre4 = Fonc[int(Nx/2):Nx,int(Ny/2):Ny]
	P[0:int(Nx/2),0:int(Ny/2)] = Carre4
	P[0:int(Nx/2),int(Ny/2):Ny] = Carre3
	P[int(Nx/2):Nx,0:int(Ny/2)] = Carre2
	P[int(Nx/2):Nx,int(Ny/2):Ny] = Carre1
	return P

def Porte2D(x0,y0,a,b,X,Y,xmin,ymin,n):
	Porte=zeros((n,n))
	dx = 2*abs(xmin)/n; dy = 2*abs(ymin)/n
	for i in range(n):
		for j in range(n):
			if (x0-a/2<=(xmin+i*dx)<=x0+a/2):
				if (y0-b/2<=(ymin+j*dy)<=y0+b/2):
					Porte[i,j]=1
			else:
				Porte[i,j]=0
	return Porte

def Elipse2D(x0,y0,a,b,X,Y,xmin,ymin,n):
	Elipse=zeros((n,n))
	for i in range(n):
		dx = xmin + i*2*abs(xmin)/n
		for j in range(n):
			dy = ymin + j*2*abs(ymin)/n
			if (((dx-x0)/a)**2+((dy-y0)/b)**2<=1):
				Elipse[i,j]=1
			else:
				Elipse[i,j]=0
	return Elipse

def Utilisateur2Dc(x0,y0,F,R,X,Y,xmin,ymin,n):
	Util=zeros((n,n))
	def Alpha(x,y):
		return eval(F)
	def Omega(x,y):
		return eval(R)
	for i in range(n):
		dx = xmin + i*2*abs(xmin)/n
		for j in range(n):
			dy = ymin + j*2*abs(ymin)/n
			if (Alpha(dx-x0,dy-y0)<=Omega(dx-x0,dy-y0)):
				Util[i,j]=1
			else:
				Util[i,j]=0
	return Util

def Utilisateur2Dp(x0,y0,F,R,x,t,xmin,ymin,n):
	Util=zeros((n,n))
	def Xx(t):
		return eval(F)
	def Yy(t):
		return eval(R)
	for i in range(n):
		dx = (xmin + i*2*abs(xmin)/n)
		for j in range(n):
			dy = (ymin + j*2*abs(ymin)/n)
			if (Alpha(dx-x0,dy-y0)<=Omega(dx-x0,dy-y0)):
				Util[i,j]=1
			else:
				Util[i,j]=0
	return Util

X, Y = meshgrid(x,y)

for i in range(0,Nb):
	dx = xmin + i*2*abs(xmin)/Nb
	for j in range(0,Nb):
		dy = ymin + j*2*abs(ymin)/Nb
		X[i,j] = (1/lambd)*(dx/f-sin(Incx0)*pi/180)
		Y[i,j] = (1/lambd)*(dy/f-sin(Incy0)*pi/180)


Fonc = zeros((Nb,Nb))
Unit = ones((Nb,Nb))
for i in range(Rep):
	Cas1 = input('fente rectangulaire, ellipsoïdale ou un équation (r/e/u): ')
	if (Cas1 == 'e'):
		x0 = float(input('Décalage en x= '))
		y0 = float(input('Décalage en y= '))
		a = eval(input('Valeur du grand axe, a= '))
		b = eval(input('Valeur du petit axe, b= '))
		Fonc += Elipse2D(x0,y0,a,b,X,Y,xmin,ymin,Nb)
	if (Cas1 == 'r'):
		x0 = float(input('Décalage en x= '))
		y0 = float(input('Décalage en y= '))
		a = eval(input('Largeur, a= '))
		b = eval(input('Hauteur, b= '))
		Fonc += Porte2D(x0,y0,a,b,X,Y,xmin,ymin,Nb)
	if (Cas1 == 'u'):
		Cas2 = input('Equation polaire ou cartésienne ? (p/c): ')
		if (Cas2 == 'p'):
			x0 = float(input('Décalage en x= '))
			y0 = float(input('Décalage en y= '))
			Fc = str(input('donner votre équation x(t)= '))
			R = str(input('donner votre équation, y(t)= '))
			Fonc += Utilisateur2Dp(x0,y0,Fc,R,X,Y,xmin,ymin,Nb)
		if (Cas2 == 'c'):
			x0 = float(input('Décalage en x= '))
			y0 = float(input('Décalage en y= '))
			Fc = str(input('donner votre équation f(x,y)= '))
			R =  str(input('donner votre condition, R(x,y)= '))
			Fonc += Utilisateur2Dc(x0,y0,Fc,R,X,Y,xmin,ymin,Nb)

nx = len(Fonc[0,:]); ny = len(Fonc[:,0])
for i in range(nx):
	for j in range(ny):
		if Fonc[i,j]>1:
			Fonc[i,j]=1
		else:0

def Masque(F,X,Y,Rep,Nb):
	return Fonc

def filtre(X,Y):
	fil = input("Donner l'équation du filtre F(X,Y) = ")
	return eval(fil)

Mas = Fonc


choix = input("Mettre un filtre ? (o/n) ")
if (choix == 'o'):
	filtt = TFD2c('filtre(X,Y)',X,Y)
	Fonc = ((cos(Incx0)+cos(Incy0))/2)*TFD2c('Masque(F,X,Y,Rep,Nb)',X,Y)*filtt
if (choix == 'n'):
	Fonc = ((cos(Incx0)+cos(Incy0))/2)*TFD2c('Masque(F,X,Y,Rep,Nb)',X,Y)
	filtt = zeros((Nb,Nb))

Foncx = Fonc[:,int(len(Fonc[:,0])/2)]
Foncy = Fonc[int(len(Fonc[0,:])/2),:]

def Max(F):
	nx = len(Fonc[0,:]); ny = len(Fonc[:,0]); m=0
	for i in range(nx):
		for j in range(ny):
			if F[i,j]>m:
				m = F[i,j]
	return m

fig = subplots(figsize=(10, 10))
imshow(filtt.T/Max(filtt),cmap='hot', origin="lower", extent=[-l/2, l/2, -L/2, L/2])
colorbar(label='Luminance')
xlabel("Distance x en (m)")
ylabel("Distance y en (m)")
title('Figure d interférence')

fig = subplots(figsize=(10, 10))
imshow(Fonc.T/Max(Fonc),cmap='hot', origin="lower", extent=[-l/2, l/2, -L/2, L/2])
colorbar(label='Luminance')
xlabel("Distance x en (m)")
ylabel("Distance y en (m)")
title('filtre')

fig = subplots(figsize=(10, 10))
imshow(transpose(Mas),cmap='hot', origin="lower", extent=[-l/2, l/2, -L/2, L/2])
xlabel("Distance x en (m)")
ylabel("Distance y en (m)")
title('Masque')

fig = subplots(figsize=(10, 10))
subplot(211)
plot(x,Foncx/Foncx.max())
xlabel("Distance x en (m)")
ylabel("Luminance")
title('I(x,0)')
subplot(212)
plot(y,Foncy/Foncy.max())
xlabel("Distance y en (m)")
ylabel("Luminance")
title('I(0,y)')
show()
