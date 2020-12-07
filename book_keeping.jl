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
alpha = [ [0.109818, 0.450771, 2.22766],
          [0.0751386, 0.231031, 0.994203] ]


