def RK4v3(f,x_0,t_0,t_f,h,count):
    print('RK4v3\n')
    ts=t_0
    times=[ts]
    xs=x_0
    sol=[xs]
    i=0
    while ts<=t_f:
        i+=1
        k1 = f(xs,ts)
        k2 = f(xs + (h/2)*k1,ts + h/2)
        k3 = f(xs + (h/2)*k2,ts + h/2)
        k4 = f(xs + h*k3,ts + h)
        xs = xs + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        ts+=h
        if Mod(i,count)==0:
            sol.append(xs)
            times.append(ts)
    return (times,sol)