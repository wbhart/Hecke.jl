export FunField, FunFieldElem, FunctionField

import Nemo: show_minus_one

# type for the function field
#
# base_field = F_q(x)
# poly_ring = base_field[t]
# field = base_field[t]/(poly)

type FunField{S <: PolyElem}

  base_field::FractionField{S}
  poly_ring::PolynomialRing{Fraction{S}}
  field::ResidueRing{Poly{Fraction{S}}}
  poly::Poly{Fraction{S}}

  function FunField(f::Poly{Fraction{S}})
    z = new()
    z.base_field = base_ring(f)
    z.poly_ring = parent(f)
    z.field = ResidueRing(z.poly_ring, f)
    z.poly = f
    return z
  end
end

function FunctionField{S}(f::Poly{Fraction{S}})
  z = FunField{S}(f)
  return z
end

base_field(F::FunField) = F.base_field

poly(F::FunField) = F.poly

degree(F::FunField) = degree(poly(F))

function show(io::IO, F::FunField)
  print(io, "Function field with defining polynomial $(poly(F))")
end

### now the type for the elements

type FunFieldElem{S}
  parent::FunField
  data::Residue{Poly{Fraction{S}}}

  function FunFieldElem()
    z = new()
    return z
  end
end

parent(x::FunFieldElem) = x.parent

### overloading the parent object

function call{S}(F::FunField{S})
  z = FunFieldElem{S}()
  z.data = F.field(0)
  z.parent = F
  return z
end

function call{S}(F::FunField{S}, x)
  z = FunFieldElem{S}()
  z.data = F.field(x)
  z.parent = F
  return z
end

################################################################################
#
#  Coefficient vector
#
################################################################################

function coeff_vector{S}(x::FunFieldElem{S})
  z = Array(Fraction{S}, degree(parent(x)))
  for i in 1:degree(parent(x))
    z[i] = coeff(x.data.data, i - 1)
  end
  return z
end

################################################################################
#
#  I/O
#
################################################################################

function show(io::IO, x::FunFieldElem)
  print(io, x.data)
end

################################################################################
#
#  Binary operations
#
################################################################################

function +{S}(x::FunFieldElem{S}, y::FunFieldElem{S})
  x.parent != y.parent && error("Must have same parents")
  z = parent(x)()
  z.data = x.data + y.data
  z.parent = x.parent
  return z
end

function -{S}(x::FunFieldElem{S}, y::FunFieldElem{S})
  x.parent != y.parent && error("Must have same parents")
  z = parent(x)()
  z.data = x.data - y.data
  return z
end

function *{S}(x::FunFieldElem{S}, y::FunFieldElem{S})
  x.parent != y.parent && error("Must have same parents")
  z = parent(x)()
  z.data = x.data * y.data
  return z
end

function //{S}(x::FunFieldElem{S}, y::FunFieldElem{S})
  x.parent != y.parent && error("Must have same parents")
  z = parent(x)()
  z.data = divexact(x.data, y.data)
  return z
end

function ^(x::FunFieldElem, n::Int)
  z = parent(x)()
  z.data = (x.data)^n
  return z
end

### ignore this

mod{S <: RingElem}(f::Poly{S}, g::Poly{S}) = pseudorem(f, g)

#function ^{S <: RingElem}(f::ResidueElem{Poly{S}}, n::Int)
#  if n == 0
#    return one(parent(f))
#  elseif n == 1
#    return f
#  elseif mod(n, 2) == 0
#    return (f*f)^(div(n, 2))
#  elseif mod(n, 2) == 1
#    return f*f^(n-1)
#  end
#end
#
#function gcdinv{T}(a::PolyElem{T}, b::PolyElem{T})
#   check_parent(a, b)
#   if length(a) == 0
#      if length(b) == 0
#         return zero(base_ring(a)), zero(parent(a))
#      else
#         d = inv(lead(b))
#         return b*d, zero(parent(a))
#      end
#   end
#   if length(b) == 0
#      d = inv(lead(b))
#      return a*d, d
#   end
#   if length(a) < length(b)
#      a, b = b, a
#      u1, u2 = zero(parent(a)), one(parent(a))
#   else
#      u1, u2 = one(parent(a)), zero(parent(a))
#   end
#   lena = length(a)
#   lenb = length(b)
#   c1 = content(a)
#   c2 = content(b)
#   A = divexact(a, c1)
#   B = divexact(b, c2)
#   u1 *= inv(c1)
#   u2 *= inv(c2)
#   while lenb > 0
#      d = lena - lenb
#      (Q, B), A = pseudodivrem(A, B), B
#      lena = lenb
#      lenb = length(B)
#      u2, u1 = u1 - Q*u2, u2
#   end
#   d = gcd(c1, c2)
#   A, u1 = d*A, d*u1
#   d = inv(lead(A))
#   return d*A, d*u1
#end
#
#function inv{S <: Nemo.FiniteFieldElem}(x::PolyElem{S})
#  degree(x) != 0 && error("Element not invertible")
#  return parent(x)(inv(lead(x)))
#end

show_minus_one{T}(::Type{Fraction{T}}) = show_minus_one(T)
