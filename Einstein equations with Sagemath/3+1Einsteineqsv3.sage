print('Version 3')
print('3-Metric tensor')
print(h)
h_inv=h.inverse()
print(h_inv)

#h_{ij,r} = dh[i][j][0]; h_{ij,theta} = dh[i][j][1]; h_{ij,phi} = dh[i][j][2]
#upper triangular matrix
dh_triang=[[[diff(h[i,j],x[k]) for k in range(3)] for j in range(i,3)] for i in range(3)]
dh=[[[dh_triang[i][j-i][k] if j>=i else dh_triang[j][i-j][k] for k in range(3)] for j in range(3)] for i in range(3)]

#show(dh)

#Christoffel symbols
#gamma^{i}_{jk} = 1/2*h^{il}(h_{lj,k}+h_{kl,j}-h_{jk,l})
gamma_triang=[[[sum((1/2)*h_inv[i,l]*(dh[l][j][k] + dh[k][l][j] - dh[j][k][l]) for l in range(3)) for k in range(j,3)] for j in range(3)]  for i in range(3)]
#gamma^i_{jk} = gamma[i][j][k]
gamma=[[[gamma_triang[i][j][k-j] if k>=j else  gamma_triang[i][k][j-k] for k in range(3)]for j in range(3)] for i in range(3)]

#show(gamma)

#Riemann tensor
#R_{ijk}^{l} = gamma^{l}_{ik,j}-gamma^{l}_{jk,i}+gamma^{a}_{ik}gamma^{l}_{aj}-gamma^{a}_{jk}gamma^{l}_{ai}
R_triang=[[[[diff(gamma[l][i][k],x[j]) - diff(gamma[l][j][k],x[i]) + sum(gamma[a][i][k]*gamma[l][a][j] - gamma[a][j][k]*gamma[l][a][i] for a in range(3)) for l in range(3)] for k in range(3)] for j in range(i,3)] for i in range(3)]
#R_{ijk}^{l} = R[i][j][k][l]
R=[[[[R_triang[i][j-i][k][l] if j>=i else -R_triang[j][i-j][k][l] for l in range(3)] for k in range(3)] for j in range(3)] for i in range(3)]

#show(R)

#Ricci tensor
#R_{ij} = R_{imj}^m
TR_triang=[[sum( R[i][m][j][m] for m in range(3)) for j in range(i,3)] for i in range(3)]
#R_{ij} = TR[i,j]
TR=matrix(3,3,[[TR_triang[i][j-i] if j>=i else  TR_triang[j][i-j]  for j in range(3)] for i in range(3)])
     
#show(TR)
     
#Ricci scalar
#R = h^{ij}R_{ij}
#R = ER
ER=sum(h_inv[i,j]*TR[i,j] for i in range(3) for j in range(3))

#show(ER)

#D_{i}N_{j}
DN=matrix(3,3,[[diff(N_v[j],x[i]) - sum(N_v[l]*gamma[l][i][j] for l in range(3)) for j in range(3)] for i in range(3)])
#Extrinsic curvature
#K_{ij}=-1/(2N)(h_{ij,t}+D_jN_i+D_iN_j)
#K_{ij} = K_m[i,j]
K_m=matrix(3,3,[[-1/(2*N)*(diff(h[i,j],t) + DN[j,i] + DN[i,j]) for j in range(3)] for i in range(3)])
#trace of the extrinsic curvature
#K = h^{ij}K_{ij}
K=sum(h_inv[i,j]*K_m[i,j] for j in range(3) for i in range(3))
#K_u[i,j] = K^{ij}
K_u=matrix(3,3,[[sum(K_m[l,m]*h_inv[l,i]*h_inv[m,j] for m in range(3) for l in range(3)) for j in range(3)] for i in range(3)])
#KK = K_{ij}K^{ij}
KK=sum(K_m[i,j]*K_u[i,j] for j in range(3) for i in range(3))
#Hamiltonian constraint
#^{3}R+K^2-K_{ij}K^{ij}=16pi*G_0E
Hc= (ER+K^2-KK).full_simplify() == 16*pi*E

#K^{i}_{j} = Kij[i,j]
Kij=matrix(3,3,[[sum(K_m[l,j]*h_inv[l,i] for l in range(3)) for j in range(3)] for i in range(3)])
#D_{l}K^{l}_{i} = DK[i]
DK=[sum(diff(Kij[l,i],x[l]) for l in range(3))+ sum(- Kij[l,m]*gamma[m][l][i] + Kij[m,i]*gamma[l][m][l] for l in range(3) for m in range(3)) for i in range(3)]
#Momentum constraints
#D_lK^l_i-D_iK=8\pi G_0 J_i
Mc=[(DK[i]-diff(K,x[i])).full_simplify() == 8*pi*J[i+1] for i in range(3)]

#N^i = N_u[i]
N_u=[sum(N_v[l]*h_inv[l,i] for l in range(3)) for i in range(3)]
#D_{i}N =dN[i]
dN=[diff(N,x[i]) for i in range(3)]
# D_{i}D_{j}N = DDN[i,j]
DDN=matrix(3,3,[[diff(dN[j],x[i])-sum(dN[l]*gamma[l][i][j] for l in range(3)) for j in range(3)] for i in range(3)])
#dynamic Einstein equations
#K_{ij,t}+N^lK_{ij,l}+K_{li}N^l_{,j}+K_{lj}N^l_{,i}+ D_iD_jN-N[^3R_{ij} + KK_{ij} - 2K_{ij}K^l_{j}]
#=4\pi N[(S-E)h_{ij}-2S_{ij}]
DEEL=matrix(3,3,[[( diff(K_m[i,j],t) + sum(N_u[l]*diff(K_m[i,j],x[l]) + K_m[l,i]*diff(N_u[l],x[j]) + K_m[l,j]*diff(N_u[l],x[i]) for l in range(3)) + DDN[i,j] - N*( TR[i,j] + K*K_m[i,j] - sum( 2*K_m[i,l]*Kij[l,j] for l in range(3) ) ) ).full_simplify() for j in range(3)] for i in range(3)])
DEE=matrix(3,3,[[DEEL[i,j] == 4*pi*N*((S - E)*h[i,j] - 2*S_m[i+1,j+1]) for j in range(3)] for i in range(3)])
#The 3+1 formulation of Eintein equations was taken from the article: Marcelo Salgado. The cauchy problem of scalar-tensor theories of gravity.Classical and Quantum Gravity, 23(14):4719-4741, jul 2006