def ph_pre(pr,tiempos,phi):
    delang=[arctan( tan(phi[i+1]) )-arctan( tan(phi[i]) ) for i in range(len(tiempos)-1)]
    deltiempos=[tiempos[i+1]-tiempos[i] for i in range(len(tiempos)-1)]
    precession=[delang[i]/deltiempos[i] for i in range(len(tiempos)-1)]
    prom_precession=mean(precession)
    print('Precession= ',prom_precession,'rad/y')
    pres=prom_precession*3_600*18*10^3/(pi.n())
    print('Precession= ',pres,'arcs/cy')
    #Porcent relative error
    pre=abs(pres-pr)/pr*100
    print('Error= ',pre,'%')