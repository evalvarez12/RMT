

def H(self) :
  return self.conjugate().T

def GUE(dim) :
  A=randn(dim,dim) + 1j*randn(dim,dim)
  return sqrt(1/(2.*dim))*( A + H(A)) 
  
  
  