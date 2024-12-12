#Astronomical unit in SI
AU=149_597_870_700
print('AU=',AU,"m")
#Newtonian constant of gravitation in SI 
G=6.674_3*10^(-11)
print('G=',G,"Nm^2/kg^2")
#Speed of light in vacuum
c=299_792_458
print('c=',c,"m/s")
#Solar radius in SI
R_sun=695_700_000
print('R_sun=',R_sun,"m")
#Solar mass
M=1.988_47*10^30
print("M=",M,"kg")
#Characteristic time (t_p=t*t_c) 1 year
t_c=60*60*24*365
#Characteristic time (t_p=t*t_c)
#t_c=R_sun/c.n()
print('t_c=',t_c,"s")
#Characteristic r (r_p=r*r_c)
r_c=(t_c^2*G*M)^(1/3)
print("r_c=",r_c,"m")
#Schwarzschild radius in SI
r_s=2*G*M/c^2
print("r_s=",r_s,"m")
#Schwarzschild radius in characteristic units
r_s=r_s/r_c
print("r_s=",r_s,"r_c")
#Speed of light in vacuum in characteristic units
c=(c*t_c/r_c)
print('c=',c,"r_c/t_c")