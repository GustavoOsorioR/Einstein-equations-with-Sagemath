x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

#solution
#[tiempos,sol]=RK4v6(SO,x_0,tau_0,tau_f,h,count,T)
[tiempos,sol]=RK4v3(SO,x_0,tau_0,tau_f,h,count)
#[tiempos,sol]=RK4v8(SO,x_0,tau_0,tau_f,h,T,ap)

N=len(tiempos)

r=[sol[i][0] for i in range(N)]
dotr=[sol[i][1] for i in range(N)]
phi=[sol[i][2] for i in range(N)]
dotphi=[sol[i][3] for i in range(N)]
t=[sol[i][4] for i in range(N)]
dott=[sol[i][5] for i in range(N)]

#Energy
E=[dott[i]*(1-2/r[i]) for i in range(N)]
#initial energy
E_0=E[0]
if E_0!=0:
    #relative error of energy
    RError=[abs((E[i]-E_0)/E_0 ) for i in range(1,N)]

#angular momentum
L=[r[i]^2*dotphi[i] for i in range(N)]
#initial angular momentum
L_0=L[0]
if L_0!=0:
    #relative error of angular momentum
    RLError=[abs((L[i]-L_0)/L_0) for i in range(1,N)]