

function right2left(rho_0,C,N,D,sites)


    print(" R to L ")

    fe = 0
    V_LL = 1.0
    VL_table = []

    orthogonalize!(C, N)
    for j = 1:N -1
      orthogonalize!(rho_0, j+1)
      A = copy(prime(rho_0))
      B = copy(rho_0)
      V_LL = V_LL*A[j]
      V_LL = V_LL*B[j]
      V_LL = V_LL*(C[j]*delta(dag(sites[j]''),sites[j]))
      V_LL = V_LL*delta(dag(sites[j]''),sites[j]''')
      push!(VL_table, V_LL)
    end

    V_R = 1.0
    for j_oc = N-1:-1:1

      orthogonalize!(rho_0, j_oc) 
      orthogonalize!(C, j_oc) 

      A = copy(prime(rho_0))
      B = copy(rho_0)

      if j_oc - 1 > 0
        V_L = copy(VL_table[j_oc-1])
      elseif j_oc -1 == 0
        V_L = 1.
      end

      if j_oc+2 <= N
        V_R =A[j_oc+2]*V_R
        V_R = B[j_oc+2]*V_R
        V_R = (C[j_oc+2]*delta(dag(sites[j_oc+2]''),sites[j_oc+2]))*V_R
        V_R = delta(dag(sites[j_oc+2]''),sites[j_oc+2]''')*V_R 
      end
      
      C_E = (((((V_L*A[j_oc])*B[j_oc])*A[j_oc+1])*B[j_oc+1])*V_R)

      if j_oc == 1
        U,S,V = svd(C_E,sites[j_oc],sites[j_oc]'';maxdim=D);
      else
        U,S,V = svd(C_E,sites[j_oc],sites[j_oc]'',commonind(C_E,C[j_oc-1]);maxdim=D);
      end

      
      fe = log(norm(S))
      S = S/norm(S) 

    
      C[j_oc] = U*S
      C[j_oc] = C[j_oc]*delta(dag(sites[j_oc]''),sites[j_oc]''')
      C[j_oc] = C[j_oc]*delta(sites[j_oc],dag(sites[j_oc]''))
      C[j_oc] = dag(C[j_oc])
      
      C[j_oc+1] = V
      C[j_oc+1] = C[j_oc+1]*delta(dag(sites[j_oc+1]''),sites[j_oc+1]''')
      C[j_oc+1] = C[j_oc+1]*delta(sites[j_oc+1],dag(sites[j_oc+1]''))
      C[j_oc+1] = dag(C[j_oc+1])

    end

    fe = fe/(N-1)

    return C, rho_0, fe

end