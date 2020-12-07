include("/home/heinen/workcopies/HartreeFock_in_Julia/functions.jl")
include("/home/heinen/workcopies/HartreeFock_in_Julia/book_keeping.jl")


global pi = 3.14159

# xyz filen ame
filename = "H2.xyz"

# Extract number of atoms, atom types, and atom coordinates
numAtoms, atom_type, atom_coordinates = read_xyz(filename)

println("Number of Atoms: $numAtoms")
println("Atom types: $atom_type")
println("Coordinates: $atom_coordinates")


B = 0
println("B = $B")

for atom in atom_type
  global B += max_quantum_number[atom]
  println(atom)
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

# Iterate through atoms
for (idx_a, val_a) in enumerate(atom_type)
  Za = charge_dict[val_a]
  Ra = atom_coordinates[idx_a]

  # Iterate through quantum numbers (1s, 2s, etc.)
  for m in 0:max_quantum_number[val_a]
    println("Quantum number m: $m")
  end
end
