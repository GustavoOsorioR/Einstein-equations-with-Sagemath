k=eval(input('Select a value for k.\nk=1 for timelike geodesics,\nk=0 for null geodesics.\nk: '))
while k!=0 and k!=1:
    k=eval(input('No valid value for k.\nSelect a valid value for k: '))
print('k=',k)
#RN Geodesic Differential Equations
def RN(X,t):
    #X=(w_1, w_2, w_3, w_4, w_5, w_6)
    dotw2=X[0]*X[3]^2*( 1-3/X[0]+2*q^2/X[0]^2 )+k/X[0]^2*(q^2/X[0]-1 )
    dotw4=-2*X[1]*X[3]/X[0]
    dotw6=2*( q^2-X[0] )*X[1]*X[5]/( X[0]*( q^2+X[0]*( X[0]-2 ) ) )
    return vector([X[1],dotw2,X[3],dotw4,X[5],dotw6])