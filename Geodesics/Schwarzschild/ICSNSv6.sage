#Impact parameter
b=R_sun/r_c
print('b= ',b,'r_c')
#Initial conditions in characteristic units
dott_0=1
print('dott_0= ',dott_0)
r_0=AU/r_c
#Initial Energy per mass
E_0=dott_0*c^2*(1-r_s/r_0)
print('E_0= ',E_0)
#Angular Momentum per mass
L_0=b*E_0/c
print('L_0= ',L_0)
dotr_0=-sqrt( E_0^2/c^2 + L_0^2/r_0^2*(r_s/r_0-1) ).n()
print('dotr_0= ',dotr_0)
phi_0=arcsin(b/r_0).n()
print('phi_0= ',phi_0)
dotphi_0=(L_0/r_0^2).n()
print('dotphi_0= ',dotphi_0)
t_0=0
print('t_0= ',t_0)

#Initial conditions vector
x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

h=0.0001 #Step size
tau_0=0
tau_f=600
count=10^3