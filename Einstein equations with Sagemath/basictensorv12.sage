print('Version 12')
print('Metric tensor')
print(g)
g_inv=g.inverse()
print(g_inv)

#g_{ij,t} = dg[i][j][0]; g_{ij,r} = dg[i][j][1]; g_{ij,theta} = dg[i][j][2]; g_{ij,phi} = dg[i][j][3]
#upper triangular matrix
dg_triang=[[[diff(g[i,j],x[k]) for k in range(4)] for j in range(i,4)] for i in range(4)]
dg=[[[dg_triang[i][j-i][k] if j>=i else dg_triang[j][i-j][k] for k in range(4)] for j in range(4)] for i in range(4)]

#show(dg)

#Christoffel symbols
#gamma^i_{jk} = 1/2*g^{il}(g_{lj,k}+g_{kl,j}-g_{jk,l})
gamma_triang=[[[sum((1/2)*g_inv[i,l]*(dg[l][j][k] + dg[k][l][j] - dg[j][k][l]) for l in range(4)) for k in range(j,4)] for j in range(4)]  for i in range(4)]
#gamma^i_{jk} = gamma[i][j][k]
gamma=[[[gamma_triang[i][j][k-j] if k>=j else  gamma_triang[i][k][j-k] for k in range(4)]for j in range(4)] for i in range(4)]

#show(gamma)

#Riemann tensor
#R_{ijk}^l = gamma^l_{ik,j}-gamma^l_{jk,i}+gamma^a_{ik}gamma^l_{aj}-gamma^a_{jk}gamma^l_{ai}
R_triang=[[[[diff(gamma[l][i][k],x[j]) - diff(gamma[l][j][k],x[i]) + sum(gamma[a][i][k]*gamma[l][a][j] - gamma[a][j][k]*gamma[l][a][i] for a in range(4)) for l in range(4)] for k in range(4)] for j in range(i,4)] for i in range(4)]
#R_{ijk}^l = R[i][j][k][l]
R=[[[[ R_triang[i][j-i][k][l] if j>=i else   -R_triang[j][i-j][k][l] for l in range(4)] for k in range(4)] for j in range(4)] for i in range(4)]

#show(R)

#R_d Riemann tensor with the indices downwards
def Riemann_1(g,R):
    #R_{ijkl} = g_{lm}R_{ijk}^m
    R_1_triang=[[[[sum(g[l,m]*R[i][j][k][m] for m in range(4)) for l in range(k,4)] for k in range(4)] for j in range(i,4)] for i in range(4)]
    #R_{ijkl} = R_d[i][j][k][l]
    return [[[[R_1_triang[i][j-i][k][l-k] if j>=i and l>=k else  -R_1_triang[j][i-j][k][l-k] if j<i and l>=k else  -R_1_triang[i][j-i][l][k-l] if j>=i and l<k else   R_1_triang[j][i-j][l][k-l] for l in range(4)] for k in range(4)] for j in range(4)] for i in range(4)]
     
#R_u Riemann tensor with the indices upwards
def Riemann_2(g_inv,R):
    #R^{ijkl} = g^{ia}g^{jb}g^{kc}R_{abc}^l
    #R^{ijkl} = R_u[i][j][k][l]
    return [[[[sum( g_inv[i,a]*g_inv[j,b]*g_inv[k,c]*R[a][b][c][l] for a in range(4) for b in range(4) for c in range(4) ) for l in range(4)] for k in range(4)] for j in range(4)] for i in range(4)]

def Kr(R_1,R_2):
#Kretschmann
    #K = R_{ijkl}R^{ijkl}
    Kretsch=sum(R_1[i][j][k][l]*R_2[i][j][k][l] for i in range(4) for j in range(4) for k in range(4) for l in range(4) )
    Kretsch=Kretsch.factor().full_simplify()
    return Kretsch

if K==1:
    R_d=Riemann_1(g,R)
    R_u=Riemann_2(g_inv,R)
    Kr=Kr(R_d,R_u)
    if R_1==1 and R_2==0:
        print('Riemann tensor with the indices downwards')
        print(R_d)
    elif R_1==0 and R_2==1:
        print('Riemann tensor with the indices upwards')
        print(R_u)
    elif R_1==1 and R_2==1:
        print('Riemann tensor with the indices downwards')
        print(R_d)
        print('Riemann tensor with the indices upwards')
        print(R_u)
    print('Kretschmann scalar')
    print(Kr)
else:
    if R_1==1 and R_2==0:
        R_d=Riemann_1(g,R)
        print('Riemann tensor with the indices downwards')
        print(R_d)
    elif R_1==0 and R_2==1:
        R_u=Riemann_2(g_inv,R)
        print('Riemann tensor with the indices upwards')
        print(R_u)
    elif R_1==1 and R_2==1:
        R_d=Riemann_1(g,R)
        print('Riemann tensor with the indices downwards')
        print(R_d)
        R_u=Riemann_2(g_inv,R)
        print('Riemann tensor with the indices upwards')
        print(R_u)

#Ricci tensor
#R_{ij} = R_{imj}^m
TR_triang=[[sum( R[i][m][j][m] for m in range(4)) for j in range(i,4)] for i in range(4)]
#R_{ij} = TR[i,j]
TR=matrix(4,4,[[TR_triang[i][j-i] if j>=i else  TR_triang[j][i-j]  for j in range(4)] for i in range(4)])
     
#show(TR)
     
#Ricci scalar
#R = g^{ij}R_{ij}
#R = ER
ER=sum(g_inv[i,j]*TR[i,j] for i in range(4) for j in range(4))

#show(ER)

#Einstein Tensor
#G_{ij} = R_{ij} - g_{ij}R/2
G_triang=[[(TR[i,j]-g[i,j]*ER/2).factor().full_simplify() for j in range(i,4)] for i in range(4)]
#G_{ij} = G[i,j]
G=matrix(4,4,[[G_triang[i][j-i] if j>=i else  G_triang[j][i-j]  for j in range(4)] for i in range(4)])

#show(G)

#Einstein field equations
#E_{ij} = G_{ij} == 8piT_{ij}
E_triang=[[(G[i,j] == 8*pi*T[i,j]).full_simplify() for j in range(i,4)] for i in range(4)]
#E_{ij} = E[i,j]
E=matrix(4,4,[[E_triang[i][j-i] if j>=i else  E_triang[j][i-j]  for j in range(4)] for i in range(4)])