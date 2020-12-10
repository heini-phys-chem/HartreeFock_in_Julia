using LinearAlgebra
using SpecialFunctions

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

  diff = norm(Ra-Rb)^2          # squared diff of the two centers
  N    = (4*a*b/(pi^2))^0.75    # normalisation
  K    = N*exp(-a*b/p*diff)     # new prefactor
  Rp   = (a*Ra + b*Rb)/p        # new center

  return p, diff, K, Rp
end


# Overlapp integral (p411)
function overlap(A, B)
  p, diff, K, Rp = gaussian_product(A, B)
  prefactor = (pi/p)^1.5

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
  Zc = charge_dict[atom_type[atom_idx]]

  return (-2*pi*Zc/p)*K*F0(p*norm(Rp-Rc)^2)
end

# (ab|cd) integral (p413)
function multi(A, B, C, D)
  p, diff_ab, K_ab, Rp = gaussian_product(A, B)
  q, diff_cd, K_cd, Rq = gaussian_product(C, D)

  multi_prefactor = 2*pi^2.5*(p*q*(p+q)^0.5)^-1

  return multi_prefactor*K_ab*K_cd*F0(p*q/(p+q)*norm(Rp-Rq)^2)
end

# check difference of two most recent guesses of the density matrix
function SD_successive_density_matrix_elements(Ptilde, P, B)
  x = 0

  for i in 1:B
    for j in 1:B
      x += B^-2*( Ptilde[i,j] - P[i,j])^2
    end
  end

  return x^0.5
end

# Multiply 0.5 to diagonal of a matrix M
function add_diag(M)
  diag_idx = diagind(M)

  for i in diag_idx
    M[i] = M[i]^-0.5
  end

  return M
end
