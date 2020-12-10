include("/home/heinen/workcopies/HartreeFock_in_Julia/functions.jl")
include("/home/heinen/workcopies/HartreeFock_in_Julia/book_keeping.jl")
include("/home/heinen/workcopies/HartreeFock_in_Julia/computing_integrals.jl")

using LinearAlgebra


global pi = 3.14159

# xyz filen ame
filename = "HHe.xyz"

# Extract number of atoms, atom types, and atom coordinates
numAtoms, atom_type, atom_coordinates = read_xyz(filename)

##println("Number of Atoms: $numAtoms")
##println("Atom types: $atom_type")
##println("Coordinates: $atom_coordinates")


B = 0

for atom in atom_type
  global B += max_quantum_number[atom]
end

##println("New B = $B")

# Number of electrons
N = 2

# Dictionary of atomic charges
charge_dict = Dict("H" => 1, "He" => 2, "C" => 6)


########################################################################################################################
# MAIN
########################################################################################################################
# CALCULATE INTEGRALS S, T, V, and DIAGONALISATION
########################################################################################################################

# Initailize matrices
S = zeros(Float64, B, B)
T = zeros(Float64, B, B)
V = zeros(Float64, B, B)

multi_electron_tensor = zeros(Float64, B, B, B, B)

S, T, V, multi_electron_tensor = compute_integrals(S, T, V, multi_electron_tensor)

Hcore = T + V

# Symmetric Orthogonalisation of basis (p144)
evalS, U = eigen(S)
ev_idx = sortperm(evalS, rev=true)
evalS = evalS[ev_idx]
U = U[:,ev_idx]
diagS = transpose(U)* (S*U)

diagS_minhalf = add_diag(diagS)
X = U * (diagS_minhalf * transpose(U))


########################################################################################################################
# HARTREE FOCK ALGORITHM
#########################################################################################################################

# Initial guess at P
P          = zeros(Float64, B, B)
P_previous = zeros(Float64, B, B)
P_list = []

threshhold = 100

while threshhold > 10^-4
  # Calculate Fock Matrix
  G = zeros(Float64, B,B)
  for i in 1:B
    for j in 1:B
      for x in 1:B
        for y in 1:B
          G[i,j] += P[x,y]*(multi_electron_tensor[i,j,y,x]-0.5*multi_electron_tensor[i,x,y,j])
        end
      end
    end
  end
  Fock = Hcore + G

  # Calculate Fock matrix in orthogonalised base
  Fockprime = *(transpose(X), *(Fock, X))
  evalFockprime, Cprime = eigen(Fockprime)

  # Correct ordering of eigenvalues and eigenvectors (sorting from ground MO as first column of C, else we get the wrong P)
  idx = sortperm(evalFockprime)
  global evalFockprime = evalFockprime[idx]
  Cprime = Cprime[:,idx]

  global C = *(X, Cprime)

  # From new P (note, we only sum over electron pairs - we DON'T sum over the entire basis set)
  for i in 1:B
    for j in 1:B
      for a in 1:trunc(Int, N/2)
        P[i,j] = 2*C[i,a]*C[j,a]
      end
    end
  end
  push!(P_list, P)

  global threshhold = SD_successive_density_matrix_elements(P_previous, P, B)
  global P_previous = copy(P)
end

# Get some Results
numSteps = size(P_list)[1]
orbEnergies1 = evalFockprime[1]
orbEnergies2 = evalFockprime[2]

println("\n")
println("STO3G Restricted Closed Shell HF algorithm took $numSteps iterations to convergence")
println("\n")
println("The orbital energies are $orbEnergies1 and $orbEnergies2")
println("\n")
println("The orbital matrix is:\n\n $C")
println("\n")
println("The density/bond order matrix is:\n\n $P")
println("\n")

