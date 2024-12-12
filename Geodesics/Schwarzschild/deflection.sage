xs=[r[i]*cos(phi[i]) for i in range(len(tiempos))]
ys=[r[i]*sin(phi[i]) for i in range(len(tiempos))]
phi0=[(phi[i],xs[i],ys[i]) for i in range(len(phi)) if abs(phi[i]-phi[0])/phi[0]<=10^-1]
#print(len(phi0))
phif=[(phi[i],xs[i],ys[i]) for i in range(len(phi)) if abs(phi[i]-phi[-1])/phi[-1]<=10^-1]
#print(len(phif))
#Initial angle ang0
#xi = phi0[0][1], xf = phi0[-1][1], #yi = phi0[0][2], yf = phi0[-1][2]
xi = phi0[0][1] 
xf = phi0[-1][1]
yi = phi0[0][2]
yf = phi0[-1][2]
ang0=arctan( (yf-yi)/(xf-xi) )*180/pi.n()
print(ang0)
#final angle angf
#xi = phif[0][1], xf = phif[-1][1], #yi = phif[0][2], yf = phif[-1][2]
xi = phif[0][1] 
xf = phif[-1][1]
yi = phif[0][2]
yf = phif[-1][2]
angf=arctan( (yf-yi)/(xf-xi) )*180/pi.n()
print(angf)
#deflection of light in degrees
deflection_N=angf-ang0
print('deflection_N= ',deflection_N)