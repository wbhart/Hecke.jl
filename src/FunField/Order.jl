export FunFieldOrd, basis_mat

type FunFieldOrd{S <: PolyElem}

  funfield::FunField{S}
  basis_ff::Array{FunFieldElem{S}, 1}
  basis_ord::Array{RingElem, 1}
  basis_mat::Mat{Fraction{S}}
  basis_mat_inv::Mat{Fraction{S}}
  disc::S

  function FunFieldOrd(F::FunField{S}, A::Mat{Fraction{S}})
    z = new{S}()
    z.funfield = F
    z.basis_mat = A
    return z
  end
end

function Order{S}(F::FunField{S}, A::Mat{Fraction{S}})
  z = FunFieldOrd{S}(F, A)
  return z
end

funfield(O::FunFieldOrd) = O.funfield

degree(O::FunFieldOrd) = degree(O.funfield)

function show(io::IO, O::FunFieldOrd)
  print(io, "Order in $(funfield(O)) with basis matrix\n$(basis_mat(O))")
end

function basis_mat(O::FunFieldOrd)
  if isdefined(O, :basis_mat)
    return O.basis_mat
  else
    @assert isdefined(O, :basis_ff)
    z = MatrixSpace(base_field(funfield(O)), degree(O), degree(O))
    for i in 1:degree(O)
      v = coeff_vector(basis_ff(O)[i])
      for j in 1:degree(O)
        z[i, j] = v[j]
      end
    end
    O.basis_mat = z
  end
end

## some helper functions

inv(A::Mat{Fraction}) = solve(A, one(parent(A)))
