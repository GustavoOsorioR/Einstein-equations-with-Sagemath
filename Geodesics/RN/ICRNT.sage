#RN Timelike
L_0=L[0]
print('L_0= ',L_0)
r_0=r[0]
print('r_0= ',r_0)
dotr_0=dotr[0]
print('dotr_0= ',dotr_0)
phi_0=phi[0]
print('phi_0=',phi_0)
dotphi_0=dotphi[0]
print('dotphi_0=',dotphi_0)
t_0=0
#Charge
q=0.2
delta = r_0^2 - 2*r_0 + q^2
E_0 = (dotr_0^2 + delta/r_0^2*(L_0^2/r_0^2 + 1))^(1/2)
dott_0=E_0*r_0^2/delta
print('dott_0=',dott_0)
print('q=',q)
#Initial conditions vector
x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

h=0.0001
tau_0=0
tau_f=1_540
count=10^4