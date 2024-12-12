#Kerr Geodesic Differential Equations
def Kerr(X,t):
    #X=(w_1, w_2, w_3, w_4, w_5, w_6)
    Delta = X[0]^2 - 2*X[0] + a^2
    dotw2 = -Delta*X[5]^2/X[0]^4 + 2*a*Delta*X[5]*X[3]/X[0]^4 + (X[0]-a^2)*X[1]^2/(X[0]*Delta) + Delta*(X[0]^3-a^2)*X[3]^2/X[0]^4
    dotw4 = -2*a*X[5]*X[1]/(Delta*X[0]^2) + 2*(2*X[0]^2+a^2-X[0]^3)*X[1]*X[3]/(Delta*X[0]^2)
    dotw6 = -2*(X[0]^2+a^2)*X[1]*X[5]/(Delta*X[0]^2) + 2*a*(3*X[0]^2+a^2)*X[1]*X[3]/(X[0]^2*Delta)
    return vector([X[1],dotw2,X[3],dotw4,X[5],dotw6])