
export rational_reconstruction, farey_lift, div, valence, leading_coefficient,
       trailing_coefficient, constant_coefficient

function PolynomialRing(R::Ring)
  return PolynomialRing(R, "_x")
end

function FiniteField(p::Integer)
  return ResidueRing(ZZ, p)
end

function FiniteField(p::fmpz)
  return ResidueRing(ZZ, p)
end

function fmpz(a::Residue{Nemo.fmpz})
  return a.data
end

function lift(R::FlintIntegerRing, a::Residue{Nemo.fmpz})
  return a.data
end

function Base.call(R::FlintIntegerRing, a::Residue{Nemo.fmpz})
  return a.data
end

## given some r/s = a mod b and deg(r) = n, deg(s) <= m find r,s
## a and b better be polynomials in the same poly ring.
## seems to work for Q (Qx) and Fp experimentally
#
# possibly should be rewritten to use an asymptotically (and practically)
# faster algorithm. For Q possibly using CRT and fast Fp techniques
# Algorithm copies from the bad-primes paper

function rational_reconstruction{S}(a::PolyElem{S}, b::PolyElem{S}, n::Int, m::Int)
  R = a.parent
  if degree(a) <= n return true, a, R(1); end

  M = MatrixSpace(R, 2, 2)()
  M[1,1] = b
  M[2,1] = a
  M[2,2] = R(1)

  T = MatrixSpace(R, 2, 2)()
  T[1,2] = R(1)
  T[2,1] = R(1)
  while degree(M[2,1]) > n
    q, r = divrem(M[1,1], M[2,1])
    T[2,2] = -q
    M = T*M
  end
  if degree(M[2,2]) <= m 
    return true, M[2,1], M[2,2]
  end

  return false, M[2,1], M[2,2]
end

function rational_reconstruction{T}(a::PolyElem{T}, b::PolyElem{T})
  return rational_reconstruction(a, b, div(degree(b), 2), div(degree(b), 2))
end

# to appease the Singular crowd...
farey_lift = rational_reconstruction


# in at least 2 examples produces the same result as Magma
# can do experiments to see if dedicated Berlekamp Massey would be
# faster as well as experiments if Berlekamp Massey yields faster 
# rational_reconstruction as well.
# Idea of using the same agorithm due to E. Thome
#

function berlekamp_massey{T}(a::Array{T, 1})
  Rx,x = PolynomialRing(parent(a[1]))
  f = Rx(a)
  xn= x^length(a)

  fl, n, d = rational_reconstruction(f, xn)
  if fl
    return true, d*(inv(trailing_coefficient(d)))
  else
    return false, Rx(0)
  end
end


function div(f::PolyElem, g::PolyElem)
  q,r = divrem(f,g)
  return q
end

# probably better off in c and faster
function valence(f::PolyElem)
  c = f(0)
  while c==0 && degree(f)>0
    f = div(f, gen(parent(f)))
    c = f(0)
  end
  return c
end

function leading_coefficient(f::PolyElem)
  return coeff(f, degree(f))
end

function trailing_coefficient(f::PolyElem)
  return coeff(f, 0)
end

constant_coefficient = trailing_coefficient

################################################################################
#
#  Functions for nmod_poly
#
################################################################################

#CF: div is neccessary for general Euc operations!
div(x::nmod_poly, y::nmod_poly) = divexact(x,y)

################################################################################
#
# Valuation
#
################################################################################
#CF TODO: use squaring for fast large valuation
#         use divrem to combine steps

function valuation(z::nmod_poly, p::nmod_poly)
  check_parent(z, p)
  z == 0 && error("Not yet implemented")
  v = 0
  while mod(z, p) == 0
    z = div(z, p)
    v += 1
  end
  return v, z
end 

function resultant(f::fmpz_poly, g::fmpz_poly, d::fmpz, nb::Int)
  z = fmpz()
  ccall((:fmpz_poly_resultant_modular_div, :libflint), Void, 
     (Ptr{fmpz}, Ptr{fmpz_poly}, Ptr{fmpz_poly}, Ptr{fmpz}, Int), 
     &z, &f, &g, &d, nb)
  return z
end


##############################################################
#
# Hensel
#
##############################################################

type fmpz_poly_factor
  c::fmpz
  poly::Ptr{fmpz_poly}
  exp::Ptr{Int} 
  _num::Int
  _alloc::Int
    
  function fmpz_poly_factor()
    z = new()
    ccall((:fmpz_poly_factor_init, :libflint), Void,
            (Ptr{fmpz_poly_factor}, ), &z)
    finalizer(z, _fmpz_poly_factor_clear_fn)
    return z
  end
end

function _fmpz_poly_factor_clear_fn(a::fmpz_poly_factor)
  ccall((:fmpz_poly_factor_clear, :libflint), Void,
          (Ptr{fmpz_poly_factor}, ), &a)
end
 
function factor_to_dict(a::fmpz_poly_factor)
  res = Dict{fmpz_poly,Int}()
  Zx,x = PolynomialRing(FlintZZ, "x")
  for i in 1:fac._num
    f = Zx()
    ccall((:fmpz_poly_set, :libflint), Void, (Ptr{fmpz_poly}, Ptr{fmpz_poly}), &f, &(unsafe_load(a.poly, i)))
    res[f] = unsafe_load(a.exp, i-1)
  end  
  return res
end

function show(io::IO, a::fmpz_poly_factor)
  ccall((:fmpz_poly_factor_print, :libflint), Void, (Ptr{fmpz_poly_factor}, ), &a)
end

type HenselCtx
  f::fmpz_poly
  p::UInt

  LF :: fmpz_poly_factor
  link::Ptr{Int}
  v::Ptr{Ptr{fmpz_poly}}
  w::Ptr{Ptr{fmpz_poly}}
  V::Array{fmpz_poly, 1}
  W::Array{fmpz_poly, 1}
  N::UInt
  lf:: Nemo.nmod_poly_factor

  function HenselCtx(f::fmpz_poly, p::fmpz)
    a = new()
    a.f = f
    a.p = UInt(p)
    Zx,x = PolynomialRing(FlintZZ, "x")
    Rx,x = PolynomialRing(ResidueRing(FlintZZ, p), "x")
    a.lf = Nemo.nmod_poly_factor(UInt(p))
    ccall((:nmod_poly_factor, :libflint), UInt,
          (Ptr{Nemo.nmod_poly_factor}, Ptr{nmod_poly}), &(a.lf), &Rx(f))
    r = a.lf._num
    a.LF = fmpz_poly_factor()
    @assert r > 1  #flint restriction
    a.v = ccall((:flint_malloc, :libflint), Ptr{Ptr{fmpz_poly}}, (Int, ), (2*r-2)*sizeof(Ptr{fmpz_poly}))
    a.w = ccall((:flint_malloc, :libflint), Ptr{Ptr{fmpz_poly}}, (Int, ), (2*r-2)*sizeof(Ptr{fmpz_poly}))
    a.V = Array(fmpz_poly, 2*r-2)
    a.W = Array(fmpz_poly, 2*r-2)
    for i=1:(2*r-2)
      a.V[i] = Zx()
      a.W[i] = Zx()
      unsafe_store!(a.v, pointer_from_objref(a.V[i]), i)
      unsafe_store!(a.w, pointer_from_objref(a.W[i]), i)
    end
    a.link = ccall((:flint_calloc, :libflint), Ptr{Int}, (Int, Int), r, sizeof(Int))
    a.N = 0
    return a
  end
end

function show(io::IO, a::HenselCtx)
  println("factorisation of $(a.f) modulo $(a.p)^$(a.N)")
end

function start_lift(a::HenselCtx, N::Int)
  ccall((:_fmpz_poly_hensel_start_lift, :libflint), Int, 
       (Ptr{fmpz_poly_factor}, Ptr{Int}, Ptr{Ptr{fmpz_poly}}, Ptr{Ptr{fmpz_poly}}, Ptr{fmpz_poly}, Ptr{nmod_poly_factor}, Int),
       &a.LF, a.p, &a.v, &a.w, &a.lf, N)
  return N
end

