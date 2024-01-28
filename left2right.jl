

function left2right(rho_0,C,N,D,sites)

        ## ============
        ## == L to R ==
        ## ============


        print(" L to R ")
        V_RR = 1.0
        VR_table = []

        orthogonalize!(C,1)
        for j = N:-1:3
          # free_energy +=  log(norm(rho_0))
          orthogonalize!(rho_0,j-2)
          A = prime(rho_0)
          V_RR = A[j]*V_RR
          V_RR = rho_0[j]*V_RR
          V_RR = (C[j]*delta(dag(sites[j]''),sites[j]))*V_RR
          V_RR = delta(dag(sites[j]''),sites[j]''')*V_RR 
          push!(VR_table, V_RR)
        end
        

        
        V_L = 1. 
        for j_oc = 1:1:N-1
          
          orthogonalize!(rho_0, j_oc) 
          orthogonalize!(C, j_oc) 

          A = deepcopy(prime(rho_0))
          B = deepcopy(rho_0)

          if j_oc > 1
            V_L = V_L*A[j_oc-1]
            V_L = V_L*rho_0[j_oc-1]
            V_L = V_L*(C[j_oc-1]*delta(dag(sites[j_oc-1]''),sites[j_oc-1]))
            V_L = V_L*delta(dag(sites[j_oc-1]''),sites[j_oc-1]''')  
          end
        
          if N - j_oc - 1 > 0
            V_R = VR_table[N-j_oc-1]
          elseif  N - j_oc - 1 == 0
            V_R = 1.
          end

          C_E = (((((V_L*A[j_oc])*B[j_oc])*A[j_oc+1])*B[j_oc+1])*V_R)

          if j_oc == 1
            U,S,V = svd(C_E,sites[j_oc],sites[j_oc]''; maxdim = D);
          else
            U,S,V = svd(C_E,sites[j_oc],sites[j_oc]'',commonind(C_E,C[j_oc-1]); maxdim = D);
          end


          S = S/norm(S) 

          C[j_oc] = U
          C[j_oc] = C[j_oc]*delta(dag(sites[j_oc]''),sites[j_oc]''')
          C[j_oc] = C[j_oc]*delta(sites[j_oc],dag(sites[j_oc]''))
          C[j_oc] = dag(C[j_oc])

          C[j_oc+1] = S*V
          C[j_oc+1] = C[j_oc+1]*delta(dag(sites[j_oc+1]''),sites[j_oc+1]''')
          C[j_oc+1] = C[j_oc+1]*delta(sites[j_oc+1],dag(sites[j_oc+1]''))
          C[j_oc+1] = dag(C[j_oc+1])

        end

        return C, rho_0

end