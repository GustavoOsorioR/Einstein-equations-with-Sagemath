#Solution
print('a=',a)
[tiempos,sol]=RK4v3(Kerr,x_0,tau_0,tau_f,h,count)

N=len(tiempos)

r=[sol[i][0] for i in range(N)]
dotr=[sol[i][1] for i in range(N)]
phi=[sol[i][2] for i in range(N)]
dotphi=[sol[i][3] for i in range(N)]
t=[sol[i][4] for i in range(N)]
dott=[sol[i][5] for i in range(N)]
    
#Kerr Energy
E=[(1-2/r[i])*dott[i] + 2*a*dotphi[i]/r[i] for i in range(N)]
#initial energy
E_0=E[0]
if E_0!=0:
    #Relative error of energy
    RError=[abs((E[i]-E_0)/E_0 ) for i in range(1,N)]

#Kerr Angular momentum
L=[-2*a*dott[i]/r[i] + (r[i]^2+a^2*(1+2/r[i]))*dotphi[i] for i in range(N)]
#initial angular momentum
L_0=L[0]
if L_0!=0:
    #Relative error of angular momentum
    RLError=[abs((L[i]-L_0)/L_0) for i in range(1,N)]