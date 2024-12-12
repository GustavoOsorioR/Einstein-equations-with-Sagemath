k=eval(input('Select a value for k.\nk=1 for timelike geodesics,\nk=0 for null geodesics.\nk: '))
while k!=0 and k!=1:
    k=eval(input('No valid value for k.\nSelect a valid value for k: '))
print('k=',k)

def solution(r,rdot,L,a):
    A = 1/2*(1 + a^2/r^2*(1 + 2/r))
    B = -2*a*L/r^3
    Delta=r^2-2*r+a^2
    C= L^2/(2*r^2)*(2/r-1)-k*Delta/(2*r^2) - rdot^2/2
    X=solve(A*x^2 + B*x + C,x)
    sol=[X[i].right().n() for i in range(len(X))]
    return sol

#Angular momentum
a=-0.2
print('a=',a)
L_0=L[0]
print('L_0= ',L_0)
r_0=r[0]
print('r_0= ',r_0)
dotr_0=dotr[0]
print('dotr_0= ',dotr_0)
phi_0=phi[0]
print('phi_0=',phi_0)
E_0=solution(r_0,dotr_0,L_0,a)[1]
print('E_0= ',E_0)
delta = r_0^2 - 2*r_0 + a^2
dotphi_0 = 1/delta*((1-2/r_0)*L_0 + 2*a*E_0/r_0)
print('dotphi_0=',dotphi_0)
t_0=0
dott_0 = 1/delta*( (r_0^2+(1+2/r_0)*a^2)*E_0 - 2*a*L_0/r_0)
print('dott_0=',dott_0)
#Initial conditions vector
x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

h=0.0001
tau_0=0
tau_f=1_468
count=10^4