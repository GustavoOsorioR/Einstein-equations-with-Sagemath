def RK4v6(f,x_0,t_0,t_f,h,count,T):
    print('RK4v6\n')
    #Perihelion
    r=x_0[0]
    ts=t_0
    times=[ts]
    xs=x_0
    sol=[xs]
    j=0
    solper=[]
    per=[]
    tper=[]
    while ts<=t_f:
        k1 = f(xs,ts)
        k2 = f(xs + (h/2)*k1,ts + h/2)
        k3 = f(xs + (h/2)*k2,ts + h/2)
        k4 = f(xs + h*k3,ts + h)
        xs = xs + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        ts+=h
        if abs(r-xs[0])/r<=10^-5 and ts>t_0+T/2:
            solper.append(xs)
            per.append(xs[0])
            tper.append(ts)
            if ts-tper[0]>T/2:
                j+=1
                if Mod(j,count)==0:
                    per.pop()
                    p=min(per)
                    sol.append(solper[per.index(p)])
                    times.append(tper[per.index(p)])
                solper=[xs]
                per=[xs[0]]
                tper=[ts]
    return (times,sol)