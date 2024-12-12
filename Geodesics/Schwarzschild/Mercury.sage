#Mercury 
e=0.205623
#perihelion
ph_si=0.30750*UA
#aphelion
ap_si=0.46670*UA
ap=ap_si/r_c
#phi's derivative in SI
dotphi_si=1.28207200482977e-6
#Sidereal orbit period in years
T=87.969/365
#perihelion precesion of Mercury in arcs per century
pr=43.1

#initial conditions in r_c and t_c units
r_0=ph_si/r_c
dotr_0=0
phi_0=0
dotphi_0=dotphi_si*t_c
t_0=0
#Newtonian Energy E_N per unit of mass
E_N=dotr_0^2/2 + r_0^2*dotphi_0^2/2 - 1/r_0
dott_0=abs(E_N)*r_0/((c*t_c/r_c)^2*(r_0-r_sc))

#initial conditions vector
x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

h=0.0001
tau_0=0
tau_f=20_000
count=4150

#solution
[tiempos,sol]=RK4v6(SO,x_0,tau_0,tau_f,h,count,T)
#[tiempos,sol]=RK4v8(SO,x_0,tau_0,tau_f,h,T,ap)

N=len(tiempos)

r=[sol[i][0] for i in range(N)]
dotr=[sol[i][1] for i in range(N)]
phi=[sol[i][2] for i in range(N)]
dotphi=[sol[i][3] for i in range(N)]
t=[sol[i][4] for i in range(N)]
dott=[sol[i][5] for i in range(N)]

#Energy in r_c and t_c units
E=[(c*t_c/r_c)^2*dott[i]*(1-r_sc/r[i]) for i in range(N)]
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