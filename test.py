
def convert(L):
    H=L.copy()
    for i  in range(len(L)):
      H[i]=tuple(L[i])
    return H

L=[[i,i] for i in range(4)]
U=convert(L)
print(U)
