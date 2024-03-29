using ITensors
using DelimitedFiles
using ITensors.HDF5
using LinearAlgebra
using CSV
using DataFrames
using MAT


include("left2right.jl")
include("right2left.jl")

let

  if isdir("Results")
    rm("Results", recursive=true)
  end
  mkdir("Results")

  
  N = 4

  sites = siteinds( "S=1/2", N, conserve_qns=true)

  
  ## ========= constructing the Hamiltonian ==================================
 
  osH = OpSum()
  J1 = 1.

  ## reading from file 
  # nn = readdlm("bondsN"*string(N)*".dat", Int)
  
  for j=1:N-1
    # osH += (J1,"Sz",j,"Sz",j+1)
    osH += (0.5*J1,"S+",j,"S-",j+1)
    osH += (0.5*J1,"S-",j,"S+",j+1)
  end
  
  osI = OpSum()
  osI += ("Id",1)

  nz = 6

  for iz = 0:nz-1

    ## ===== constructing the density matrix MPO, rho_0, at tau_0=beta_0/2 ====

    tau_0 = 1e-3
    beta_0 = (2.0^(iz/(nz-1)))*tau_0
    rho_0 = MPO(osI-beta_0*osH, sites)
    ham = MPO(osH,sites)
    free_energy = 0
    
    nmax = 15 # Maximum value to beta to be reached is tau_0*2^nmax

    for it = 1:nmax 

      beta_0 = 2*beta_0
      print("beta =", beta_0, " ")
      
      D = 50
      C = dag(prime(prime(rho_0))) 
      orthogonalize!(C, N)

      nsweep = 4 # Number recommended in PRX

      for isweep = 1:nsweep

        C, rho_0 = left2right(rho_0,C,N,D,sites)
        C, rho_0 = right2left(rho_0,C,N,D,sites)

      end 

      for j = 1:N
        rho_0[j]=(dag(C[j])*delta(dag(sites[j]'''),sites[j]'))*delta(sites[j]'',dag(sites[j]))
      end


      Z_n = 1.0
      for j = 1:N
        Z_n = Z_n*(rho_0[j]*delta(dag(sites[j]'),sites[j]))
      end
      
      u = 1.0
      for j = 1:N
        u = (u*rho_0[j])*dag(ham[j])
      end

      free_energy +=  log(Z_n[])

      u = u[]/Z_n[]/N

      println(u[])

      io = open("Results/beta_list.dat", "a") 
      writedlm(io,[beta_0 u free_energy])
      close(io)

    end

  end

  
path = "Results/"
df1 = DataFrame(CSV.File(path*"beta_list.dat"))

df1 = Matrix(df1)
beta = df1[:,1]
ener = df1[:,2]
free_energy = df1[:,3]

file = matopen(path*"XTRG.mat", "w")
write(file, "beta", beta)
write(file, "ener", ener)
write(file, "free_energy", free_energy)
write(file, "N", N)
close(file)



end



