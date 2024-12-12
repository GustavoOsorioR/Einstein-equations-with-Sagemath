#Mercury's orbital eccentricity 
e=0.205623
#Perihelion
ph_si=0.30750*AU
#Aphelion
ap_si=0.46670*AU
ap=ap_si/r_c
#phi's derivative in SI
dotphi_si=1.28207200482977e-6
#Sidereal orbit period in years
T=87.969/365
#Perihelion precesion of Mercury in arcs per century
pr=43.1

#Initial conditions in characteristic units
r_0=ph_si/r_c
dotr_0=0
phi_0=0
dotphi_0=dotphi_si*t_c
t_0=0
L_0=r_0^2*dotphi_0
#Schwarzschild Energy E per unit of mass in characteristic units
E=c*sqrt( dotr_0^2+(L_0^2/r_0^2+c^2)*(1-r_s/r_0) )
dott_0=E/((1-r_s/r_0)*c^2)
print('dott_0=',dott_0)

#Initial conditions vector
x_0=vector([r_0,dotr_0,phi_0,dotphi_0,t_0,dott_0])

h=0.0001 #Step size
tau_0=0
tau_f=20_000
count=4150