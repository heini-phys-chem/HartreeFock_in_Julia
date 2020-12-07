include("/home/heinen/workcopies/HartreeFock_in_Julia/functions.jl")
include("/home/heinen/workcopies/HartreeFock_in_Julia/book_keeping.jl")
include("/home/heinen/workcopies/HartreeFock_in_Julia/computing_integrals.jl")

using LinearAlgebra


global pi = 3.14159

# xyz filen ame
filename = "H2.xyz"

# Extract number of atoms, atom types, and atom coordinates
numAtoms, atom_type, atom_coordinates = read_xyz(filename)

println("Number of Atoms: $numAtoms")
println("Atom types: $atom_type")
println("Coordinates: $atom_coordinates")


B = 0

for atom in atom_type
  global B += max_quantum_number[atom]
end

println("New B = $B")

# Number of electrons
N = 2

# Dictionary of atomic charges
charge_dict = Dict("H" => 1, "C" => 6)


########################################################################################################################
# MAIN
########################################################################################################################

# Initailize matrices
S = zeros(Float64, B, B)
T = zeros(Float64, B, B)
V = zeros(Float64, B, B)

multi_electron_tensor = zeros(Float64, B, B, B, B)

S, T, V, multi_electron_tensor = compute_integrals(S, T, V, multi_electron_tensor)

println("S: $S")
println("T: $T")
println("V: $V")
println("multi elec: $multi_electron_tensor")

Hcore = T + V
println("\nHcore: $Hcore")

# Symmetric Orthogonalisation of basis (p144)
evalS = eigvals(S)
U = eigvecs(S)
println("\nevalsS: $evalS\nU: $U")

#println("SU: ", dot(S, U))
#SU = zeros(Float64, 1, 1)
#SU[1] = dot(S[1], U[1])
#SU[2] = dot(S[2], U[2])
#
#println("SU: $SU")
#
#
#diagS = dot(U, SU)
#println("diagS: $diagS")





