include("/home/heinen/workcopies/HartreeFock_in_Julia/functions.jl")
include("/home/heinen/workcopies/HartreeFock_in_Julia/book_keeping.jl")

using ProgressBars

function compute_integrals(S, T, V, multi_electron_tensor)
  # Iterate through atoms
  println("\n[-] Calculating integrals")
  #for (idx_a, val_a) in tqdm(enumerate(atom_type))
  for (idx_a, val_a) in enumerate(atom_type)
    Za = charge_dict[val_a]
    Ra = atom_coordinates[idx_a]
  
    # Iterate through quantum numbers (1s, 2s, etc.)
    for m in 1:max_quantum_number[val_a]
      d_vec_m = D[m]
      zeta = zeta_dict[val_a][m]
      alpha_vec_m = alpha[m]*zeta^2
  
      # Iterate over the contraction coefficients
      for p in 1:STOnG
  
        # Iterate through atoms once again
        for (idx_b, val_b) in enumerate(atom_type)
          Zb = charge_dict[val_b]
          Rb = atom_coordinates[idx_b]
  
          for n in 1:max_quantum_number[val_b]
            d_vec_n = D[n]
            zeta = zeta_dict[val_b][n]
            alpha_vec_n = alpha[n]*zeta^2
  
            for q in 1:STOnG
              a = (idx_a)*(m)
              b = (idx_b)*(n)
  
              # Generate the overlap, kinetic and potential matrices
              S[a,b] += d_vec_m[p]*d_vec_n[q]*overlap((alpha_vec_m[p], Ra), (alpha_vec_n[q], Rb))
              T[a,b] += d_vec_m[p]*d_vec_n[q]*kinetic((alpha_vec_m[p], Ra), (alpha_vec_n[q], Rb))
  
              for i in 1:numAtoms
                V[a,b] += d_vec_m[p]*d_vec_n[q]*potential((alpha_vec_m[p], Ra), (alpha_vec_n[q], Rb), i)
              end
  
              # Two more iterations to get the multi-electron-tensor
              for (idx_c, val_c) in enumerate(atom_type)
                Zc = charge_dict[val_c]
                Rc = atom_coordinates[idx_c]
  
                for k in 1:max_quantum_number[val_c]
                  d_vec_k = D[k]
                  zeta = zeta_dict[val_c][k]
                  alpha_vec_k = alpha[k]*zeta^2
  
                  for r in 1:STOnG
                    for (idx_d, val_d) in enumerate(atom_type)
                      Zd = charge_dict[val_d]
                      Rd = atom_coordinates[idx_d]
                      for l in 1:max_quantum_number[val_d]
                        d_vec_l = D[l]
                        zeta = zeta_dict[val_d][l]
                        alpha_vec_l = alpha[l]*zeta^2
  
                        for s in 1:STOnG
                          c = idx_c*k
                          d = idx_d*l
  
                          multi_electron_tensor[a,b,c,d] += d_vec_m[p]*d_vec_n[q]*d_vec_k[r]*d_vec_l[s]*(
                                                          multi((alpha_vec_m[p], Ra),
                                                                (alpha_vec_n[q], Rb),
                                                                (alpha_vec_k[r], Rc),
                                                                (alpha_vec_l[s], Rd))
                                                         )
                        end
                      end
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
  println("[+] Calculating integrals\n")

  return S, T, V, multi_electron_tensor

end
