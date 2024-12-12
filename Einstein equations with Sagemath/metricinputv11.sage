print('Version 11')
#variables
var('t,r,theta,phi,M,a,k,q')
x=[t,r,theta,phi]

m=eval(input('Select an integer value for m: 1 to 4.\nm=1 for Schwarzschild,\nm=2 for FRW,\nm=3 for RN,\nor m=4 for Kerr.\nm: '))
while m!=1 and m!=2 and m!=3 and m!=4:
    m=eval(input('No valid value for m.\nSelect a valid value for m: '))
if m==1:
    print('You chose Schwarzschild metric.')
elif m==2:
    print('You chose FRW metric.')
elif m==3:
    print('You chose RN metric.')
elif m==4:
    print('You chose Kerr metric.')

R_1=0
R_2=0
K=0

if m==1:
    #Schwarzschild
    g=matrix(4,4,[factor(-(1-(2*M)/r)),0,0,0,0,factor(1/(1-(2*M)/r)),0,0,0,0,r^2,0,0,0,0,r^2*sin(theta)^2])
    #Stress-energy tensor
    #Vacuum
    #T_{ij} = 0
    #T_{ij} = T[i,j]
    T=matrix(4,4,0)    
elif m==2:
    #FRW
    t=var('t')
    a = function("a")(t)
    g=matrix(4,4,[-1,0,0,0,0,factor(a^2/(1-k*r^2)),0,0,0,0,a^2*r^2,0,0,0,0,a^2*r^2*sin(theta)^2])
    #Perfect fluid
    #4-velocity
    U=[-1,0,0,0]
    var('rho,P')
    #Stress-energy tensor
    #T_{ij} = (rho+P)U_iU_j + Pg_{ij}
    T=matrix(4,4,[[(rho+P)*U[i]*U[j]+P*g[i,j] for j in range(4)] for i in range(4)])
elif m==3:
    #Reissner-Nordstrom
    g=matrix(4,4,[-(1-2*M/r+q^2/r^2),0,0,0,0,factor((1-2*M/r+q^2/r^2)^(-1)),0,0,0,0,r^2,0,0,0,0,r^2*sin(theta)^2])
    g_inv=g.inverse()
    #Electromagnetic stress tensor
    #Electromagnetic four-potencial
    A=[q/r,0,0,0]
    #Faraday tensor
    F=matrix(4,4,[[diff(A[j],x[i])-diff(A[i],x[j]) for j in range(4)] for i in range(4)])
    #F_1_{ij} = F_{ia}F^{ja}
    F_1=matrix(4,4,[[sum(F[i,k]*F[j,l]*g_inv[l,k] for k in range(4) for l in range(4) ) for j in range(4)] for i in range(4)] )
    #F_2 = F_{ab}F^{ab}
    F_2=sum(F[i,j]*F[k,l]*g_inv[k,i]*g_inv[l,j] for k in range(4) for l in range(4) for j in range(4) for i in range(4))
    #T_{ij} = 1/4pi(F_{ia}F_{j}^{a} - 1/4g_{ij}F_{ab}F^{ab})
    T_triang=[[factor( 1/(4*pi)*(F_1[i,j]-1/4*g[i,j]*F_2) ) for j in range(i,4)] for i in range(4)]
    T=matrix(4,4,[[T_triang[i][j-i] if j>=i else  T_triang[j][i-j]  for j in range(4)] for i in range(4)])
elif m==4:
    #Kerr
    g=matrix(4,4,[factor( -(1-(2*M*r)/(r^2+a^2*cos(theta)^2)) ),0,0,-(2*M*a*r*sin(theta)^2)/(r^2+a^2*cos(theta)^2),0,(r^2+a^2*cos(theta)^2)/(r^2-2*M*r+a^2),0,0, 0,0,r^2+a^2*cos(theta)^2,0,-(2*M*a*r*sin(theta)^2)/(r^2+a^2*cos(theta)^2),0,0,factor((r^2+a^2+(2*M*a^2*r*sin(theta)^2)/(r^2+a^2*cos(theta)^2))*sin(theta)^2)])
    #Stress-energy tensor
    #Vacuum
    #T_{ij} = 0
    T=matrix(4,4,0)