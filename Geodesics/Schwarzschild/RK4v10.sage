def RK4v10(f,x_0,t_0,t_f,h,T):
    print('RK4v10\n')
    #perihelion
    r=x_0[0]
    print(r)
    ts=t_0
    times=[ts]
    xs=x_0
    sol=[xs]
    csol=[]
    ctsol=[]
    absdotr=[]
    while ts<=t_f+T/2+T*5/100:
        k1 = f(xs,ts)
        k2 = f(xs + (h/2)*k1,ts + h/2)
        k3 = f(xs + (h/2)*k2,ts + h/2)
        k4 = f(xs + h*k3,ts + h)
        xs = xs + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        ts+=h
        if abs(xs[1])<=10^-3 and t_f-(T/2+T*5/100)<=ts and r-(r*5/100)<xs[0]<r+r*5/100:
            csol.append(xs)
            ctsol.append(ts)
            #absolute value of radial velocity
            absdotr.append(abs(xs[1]))
    indexper=absdotr.index(min(absdotr))
    sol.append(csol[indexper])
    times.append(ctsol[indexper])
    return (times,sol)