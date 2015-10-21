export FunFieldOrd, basis_mat

type Test{T <: Real}
  x::Array{T, 1}
end

function +{T <: Real}(x::Test{T}, y::Test{T})
  return x.x[1] + y.x[1]
end

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

function hnf{T}(A::Mat{T})
  m = rows(A)
  n = cols(A)
  W = parent(A)()
# Step 1
  i = m
  j = n
  k = n
  if m <= n
    l = 1
  elseif
    l = m - n + 1
  end
# Step 2
  if j != 1
    j = j - 1
    while A[i,j] != 0
      j = j - 1
    end
# Step 3
    u, v, d = gcdx(A[i, k], A[i, j])
    B = u * A[k] + v * A[j]
    A[j] = A[i,k]/d * A[j] - a[i,j]*A[k]
    A[k] = B
  end
# Step 4
  b = a[i,k]
  if b == 0
    k = k + 1
  elseif j > k
    q = divrem(A[i,j], b)[1]
    A[j] = A[i] - q*A[k]
  end
# Step 5
  if i == l
    for j = 1:(n-k+1)
      W[j] = A[j+k-1]
    end
  else
    i = i - 1
    k = k - 1
    j = k
    goto 2
  end
  return W
end
