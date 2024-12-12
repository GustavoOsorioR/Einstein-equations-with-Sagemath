#RN null
L_0=L[0]
print('L_0= ',L_0)
E_0=E[0]
print('E_0= ',E_0)
r_0=r[0]
print('r_0= ',r_0)
#Charge
q=0.2
print('q=',q)
delta = r_0^2 - 2*r_0 + q^2
dotr_0=-(E[0]^2 - L[0]^2/r[0]^4*delta)^(1/2)
print('dotr_0= ',dotr_0)
phi_0=phi[0]
print('phi_0=',phi_0)
dotphi_0=dotphi[0]
print('dotphi_0=',dotphi_0)
t_0=0
dott_0=E_0*r_0^2/delta
print('dott_0=',dott_0)
#Initial conditions vector
x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

h=0.0001
tau_0=0
tau_f=1_000
count=10^4