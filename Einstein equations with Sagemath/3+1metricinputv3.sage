print('Version 3')
#variables
var('t,r,theta,phi,M,a,k,q')
#3-vector
x=[r,theta,phi]

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

if m==1:
    #Schwarszchild
    #lapse function
    N=sqrt((1-(2*M)/r))
    #Shift vector
    # N_{i} = N_v[i]
    N_v=[0,0,0]
    #3-metric
    h=matrix(3,3,[factor(1/(1-(2*M)/r)),0,0,0,r^2,0,0,0,r^2*sin(theta)^2])
    #3-momentum density vector
    #J_{i} = J[i]
    J = [0 for i in range(4)]
    #3-energy-momentum tensor
    #S_{ij} = S_m[i,j]
    S_m = matrix(4,4,0)
    #trace of 3-energy-momentum tensor
    S=0
    #total energy density
    E=0
elif m==2:
    #FRW
    #lapse function
    N=1
    #Shift vector
    # N_{i} = N_v[i]
    N_v=[0,0,0]
    # N^{i} = N_u[i]
    N_u=[0,0,0]
    t=var('t')
    a = function("a")(t)
    #metric
    g=matrix(4,4,[-1,0,0,0,0,factor(a^2/(1-k*r^2)),0,0,0,0,a^2*r^2,0,0,0,0,a^2*r^2*sin(theta)^2])
    g_inv=g.inverse()
    #3-metric
    h=matrix(3,3,[[g[i+1,j+1] for j in range(3)] for i in range(3)])
    #Perfect fluid
    #4-velocity
    U=[-1,0,0,0]
    var('rho,P')
    #Stress-energy tensor
    #T_{ij} = (rho+P)U_iU_j + Pg_{ij}
    T=matrix(4,4,[[(rho+P)*U[i]*U[j]+P*g[i,j] for j in range(4)] for i in range(4)])
    #T^{ij}=T_u[i,j]
    T_u=matrix(4,4,[[sum(T[l,m]*g_inv[i,l]*g_inv[j,m] for l in range(4) for m in range(4)) for j in range(4)] for i in range(4)])
    #4-velocity
    #n_{i} = n[i]
    n=[-N,0,0,0]
    #h_{ij} = g_{ij}+n_in_j
    #^4h_{ij}=H[i,j]
    H=matrix(4,4,[[g[i,j]+n[i]*n[j] for j in range(4)] for i in range(4)])
    #3-energy-momentum tensor
    #S_{ij}=h_{il}h_{jm}T^{lm}
    #S_{ij}=S_m[i,j]
    S_m=matrix(4,4,[[sum(H[i,l]*H[j,m]*T_u[l,m] for l in range(4) for m in range(4)) for j in range(4)] for i in range(4)])
    #3-momentum density vector
    #J_{i}=-h_{il}T^{lm}n_m
    #J_{i}=J[i]
    J=[-sum(H[i,l]*T_u[l,m]*n[m] for l in range(4) for m in range(4)) for i in range(4)]
    #total energy density
    #E=T^{lm}n_ln_m
    E=sum(T_u[l,m]*n[l]*n[m] for l in range(4) for m in range(4))
    #trace of 3-energy-momentum tensor
    #S=g^{ij}S_{ij}
    S=sum(S_m[l,m]*g_inv[l,m] for l in range(4) for m in range(4))
elif m==3:
    #RN
    #lapse function
    N=sqrt(1-2*M/r+q^2/r^2).factor()
    #Shift vector
    # N_{i} = N_v[i]
    N_v=[0,0,0]
    #metric
    g=matrix(4,4,[-(1-2*M/r+q^2/r^2),0,0,0,0,factor((1-2*M/r+q^2/r^2)^(-1)),0,0,0,0,r^2,0,0,0,0,r^2*sin(theta)^2])
    g_inv=g.inverse()
    #3-metric
    h=matrix(3,3,[[g[i+1,j+1] for j in range(3)] for i in range(3)])
    #Electromagnetic stress tensor
    y=[t,r,theta,phi]
    #Electromagnetic four-potencial
    A=[q/r,0,0,0]
    #Faraday tensor
    F=matrix(4,4,[[diff(A[j],y[i])-diff(A[i],y[j]) for j in range(4)] for i in range(4)])
    #F_1_{ij} = F_{ia}F^{ja}
    F_1=matrix(4,4,[[sum(F[i,k]*F[j,l]*g_inv[l,k] for k in range(4) for l in range(4) ) for j in range(4)] for i in range(4)] )
    #F_2 = F_{ab}F^{ab}
    F_2=sum(F[i,j]*F[k,l]*g_inv[k,i]*g_inv[l,j] for k in range(4) for l in range(4) for j in range(4) for i in range(4))
    #T_{ij} = 1/4pi(F_{ia}F_{j}^{a} - 1/4g_{ij}F_{ab}F^{ab})
    T_triang=[[factor( 1/(4*pi)*(F_1[i,j]-1/4*g[i,j]*F_2) ) for j in range(i,4)] for i in range(4)]
    T=matrix(4,4,[[T_triang[i][j-i] if j>=i else  T_triang[j][i-j]  for j in range(4)] for i in range(4)])
    #T^{ij}=T_u[i,j]
    T_u=matrix(4,4,[[sum(T[l,m]*g_inv[i,l]*g_inv[j,m] for l in range(4) for m in range(4)) for j in range(4)] for i in range(4)])
    #4-velocity
    #n_{i} = n[i]
    n=[-N,0,0,0]
    #h_{ij} = g_{ij}+n_in_j
    #^4h_{ij}=H[i,j]
    H=matrix(4,4,[[g[i,j]+n[i]*n[j] for j in range(4)] for i in range(4)])
    #3-energy-momentum tensor
    #S_{ij}=h_{il}h_{jm}T^{lm}
    #S_{ij}=S_m[i,j]
    S_m=matrix(4,4,[[sum(H[i,l]*H[j,m]*T_u[l,m] for l in range(4) for m in range(4)) for j in range(4)] for i in range(4)])
    #3-momentum density vector
    #J_{i}=-h_{il}T^{lm}n_m
    #J_{ij}=J[i,j]
    J=[-sum(H[i,l]*T_u[l,m]*n[m] for l in range(4) for m in range(4)) for i in range(4)]
    #total energy density
    #E=T^{lm}n_ln_m
    E=sum(T_u[l,m]*n[l]*n[m] for l in range(4) for m in range(4)).factor()
    #trace of 3-energy-momentum tensor
    #S=h^{ij}S_{ij}
    S=(sum(S_m[l,m]*g_inv[l,m] for l in range(4) for m in range(4))).full_simplify()
elif m==4:
    #Kerr
    rho=r^2+a^2*cos(theta)^2
    delta=r^2-2*M*r+a^2
    #3-metric
    h=matrix(3,3,[(r^2+a^2*cos(theta)^2)/(r^2-2*M*r+a^2),0,0,0,r^2+a^2*cos(theta)^2,0,0,0,factor((r^2+a^2+(2*M*a^2*r*sin(theta)^2)/(r^2+a^2*cos(theta)^2))*sin(theta)^2)])
    #lapse function
    N=sqrt(delta*rho/((r^2+a^2)^2-delta*a^2*sin(theta)^2))
    #Shift vector
    # N_{i} = N_v[i]
    N_v=[0,0,(2*M*a*r*sin(theta)^2)/(r^2+a^2*cos(theta)^2)]
    #3-momentum density vector
    #J_{i} = J[i]
    J = [0 for i in range(4)]
    #S_{ij} = S_m[i,j]
    #3-energy-momentum tensor
    S_m = matrix(4,4,0)
    #trace of 3-energy-momentum tensor
    S=0
    #total energy density
    E=0