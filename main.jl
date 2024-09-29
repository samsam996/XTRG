using DelimitedFiles
using LinearAlgebra
using CSV
using DataFrames
using MAT
using ITensors
using ITensors.HDF5

include("src/XTRG_code.jl")


let 
 

    # Maximum value to beta to be reached is tau_0*2^nmax
    
    number_site = 5
    max_bond_dimension = 20
    nz = 10 # parallelisation du problème
    nmax = 10
    tau_0 = 1e-4 # initial beta

    XTRG(number_site, max_bond_dimension, nz, nmax, tau_0)
    
end