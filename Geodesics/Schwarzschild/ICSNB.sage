#Initial conditions
b_c=3^(3/2).n()
#Deflection in degrees
d=10
#Impact parameter
b=4/d*180/pi.n()
#b=b_c
print('b= ',b)
#Angular Momentum per mass
L_0=5
print('L_0= ',L_0)
#Energy per mass
E_0=L_0/b
r_0=100
dotr_0=-sqrt( E_0^2 + L_0^2/r_0^2*(2/r_0-1) ).n()
print('dotr_0= ',dotr_0)
phi_0=arcsin(b/r_0).n()
print('phi_0=',phi_0)
dotphi_0=(L_0/r_0^2).n()
print('dotphi_0=',dotphi_0)
t_0=0
dott_0=abs(E_0)*r_0/(r_0-2)

#Initial conditions vector
x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

#Speed of light in vacuum
c=1
#Schwarzschild radius
r_s=2

h=0.0001 #Step size
tau_0=0
tau_f=1_000
count=10^4