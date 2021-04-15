

function basic_linalg_problem(;Ns=1000,Nc=100)

  A = Diagonal(2*ones(Ns))
  C = Diagonal(ones(Nc))
  B2 = zeros(Nc,Ns)
  for j in 1:min(Nc,Ns)
      B2[j,j] = 1.0
  end
  B1t = B2'
  rhs1v = ones(Ns)
  rhs2v = 2*ones(Nc)

  Abig = [A B1t;B2 C]
  rhsbig = [rhs1v;rhs2v]
  solex = similar(rhsbig)
  solex .= Abig\rhsbig;

  return A,B2,B1t,C,rhs1v,rhs2v, solex

end
