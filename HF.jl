
global pi = 3.14159

# Reading in xyz file
function read_xyz(f)

  lines = readlines(f)

  numAtoms         = parse(Int, lines[1])
  atom_type        = [] 
  atom_coordinates = [] 

  for i in 1:numAtoms+2
    if i < 3
      continue
    end
    push!(atom_type, split(lines[i], ' ')[1])
    push!(atom_coordinates, [parse(Float64, split(lines[i], ' ')[2]),
                             parse(Float64, split(lines[i], ' ')[3]),
                             parse(Float64, split(lines[i], ' ')[4])] )
  end

  return numAtoms, atom_type, atom_coordinates

end

# Integrals between gaussian orbitals (p410)
function gaussian_product(A, B)
  a, Ra = A
  b, Rb = B
  p = a + b

  diff = abs2(Ra-Rb)                # squared diff of the two centers
  N    = (4*a*b/(pi^2))*0.75    # normalisation
  K    = N*exp(-a*b/p*diff)         # new prefactor
  Rp   = (a*Ra + b*Rb)/p            # new center

  return p, diff, K, Rp
end


# Overlapp integral (p411)
function overlap(A, B)
  p, diff, K, Rp = gaussian_product(A, B)
  predactor = (pi/p)^1.5

  return prefactor*K
end

# Kinetic integra
function kinetic(A, B)
  p, diff, K, Rp = gaussian_product(A, B)
  prefactor = (pi/p)^1.5

  a, Ra = A
  b, Rb = B

  reduced_exponent = a*b/p

  return reduced_exponent*(3-2*reduced_exponent*diff)*prefactor*K
end

# F0 function (Boys function) for calculating potential and e-e repulsion terms (p414)
# TODO add higher orbitals (2p, 3d) -> only works for s orbitals
function F0(t)
  if t == 0
    return 1
  else
    return (0.5*(pi/t)^0.5)*erf(t^0.5)
  end
end

# Nuclear-electron integral (p412)
function potential(A, B, atom_idx)
  p, diff, K, Rp = gaussian_product(A, B)
  Rc = atom_coordinates[atom_idx]
  Zc = charge_dict[atoms[atom_idx]]

  return (-2*pi*Zc/p)*K*F0(p*abs2(Rp-Rc))
end

# (ab|cd) integral (p413)
function multi(A, B, C, D)
  p, diff_ab, K_ab, Rp = gaussian_product(A, B)
  q, diff_cd, K_cd, Rq = gaussian_product(C, D)

  multi_prefactor = 2*pi^2.5*(p*q*(p+q)^0.5)^-1

  return multi_prefactor*K_ab*K_cd*F0(p*q/(p+q)*abs2(Rp-Rq))
end

########################################################################################################################
# MAIN
########################################################################################################################


gaussian_product([1,1], [2,2])


# xyz filen ame
filename = "H2.xyz"

# Extract number of atoms, atom types, and atom coordinates
numAtoms, atom_type, atom_coordinates = read_xyz(filename)

println(numAtoms)
println(atom_type)
println(atom_coordinates)


# Basis set variables
# -------------------

# Number of gaussians used to form a contracted gaussian orbital (p153)
STOnG = 3

# dictionary of zeta values (p159-160, 170)
zeta_dict = Dict("H" => [1.24], "C" => [5.64, 1.72])

# Dictionary containing the maximal quantum numbers of each atom for a minimal STO3G basis set
max_quantum_number = Dict("H" => 1, "C" => 2)

# Gaussian contraction coefficients (p157), going up to 2s orbitals
D = [ [0.4444635, 0.533328, 0.154329],
      [0.7000115, 0.399513, -0.099672] ]

# Gaussian orbital exponents (p153)
alpha = [ [0.109818, 0.450771, 2.22766]
          [0.0751386, 0.231031, 0.994203] ]

# Basis set size TODO: add loop B+=max_quantum_numer[atom] for atom in atoms
B = 0

# Number of electrons
N = 2

# Dictionary of atomic charges
charge_dict = Dict("H" => 1, "C" => 6)


