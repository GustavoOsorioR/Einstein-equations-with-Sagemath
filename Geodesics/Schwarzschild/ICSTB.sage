#Initial conditions
#Angular Momentum per mass
L_0=sqrt(24).n()
print('L_0= ',L_0)
r_0=15
#r maximum of effective potential
#r_0=(L_0^2-(L_0^4-12*L_0^2)^(1/2))/2
print('r_0= ',r_0)
dotr_0=0
print('dotr_0= ',dotr_0)
phi_0=0
print('phi_0_0= ',phi_0)
dotphi_0=L_0/r_0^2
print('dotphi_0= ',dotphi_0)
t_0=0
print('t_0= ',t_0)
#Energy per mass
E_0=sqrt( dotr_0^2 + (L_0^2/r_0^2+1)*(1-2/r_0) ).n()
print('E_0=',E_0)
dott_0=abs(E_0)*r_0/(r_0-2)
print('dott_0=',dott_0)

#Initial conditions vector
x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

#Speed of light in vacuum
c=1
#Schwarzschild radius
r_s=2

h=0.0001 #Step size
tau_0=0
tau_f=3_820
count=10^4