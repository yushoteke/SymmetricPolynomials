module SymmetricPolynomials

using TupleTools
using DataStructures

export semi_elementary_monomial, decompose

include("utilities.jl")
include("semi elementary monomial.jl")
include("semi elementary polynomial.jl")
include("evaluate.jl")

end

x = SymmetricPolynomials.semi_elementary_monomial((0, 0, 0, 0, 10), (0, 0, 0, 0, 0))

p1 = SymmetricPolynomials.decompose(x);
p2 = SymmetricPolynomials.decompose2(x);
p3 = SymmetricPolynomials.decompose3(x);
p4 = SymmetricPolynomials.decompose4(x);
p5 = SymmetricPolynomials.decompose5(x);
p6 = SymmetricPolynomials.decompose6(x);
p7 = SymmetricPolynomials.decompose7(x);
p8 = SymmetricPolynomials.decompose8(x);
p9 = SymmetricPolynomials.decompose9(x);

println(p1 == p2 == p3 == p4 == p5 == p6 == p7 == p8 == p9)

for i = 20:10:60
  println("i = $(i)")
  tmp = SymmetricPolynomials.semi_elementary_monomial((0, 0, 0, 0, 0, 0, i), (0, 0, 0, 0, 0, 0, 0))
  @time res1 = SymmetricPolynomials.decompose(tmp)
  @time res2 = SymmetricPolynomials.decompose2(tmp)
  @time res3 = SymmetricPolynomials.decompose3(tmp)
  @time res4 = SymmetricPolynomials.decompose4(tmp)
  @time res5 = SymmetricPolynomials.decompose5(tmp)
  @time res6 = SymmetricPolynomials.decompose6(tmp)
  @time res7 = SymmetricPolynomials.decompose7(tmp)
  @time res8 = SymmetricPolynomials.decompose8(tmp)
  @time res9 = SymmetricPolynomials.decompose9(tmp)
  println("results are same:", res1 == res2 == res3 == res4 == res5 == res6 == res8 == res9)
end

for i = 50:10:80
  println("i = $(i)")
  tmp = SymmetricPolynomials.semi_elementary_monomial((0, 0, 0, 0, 0, 0, i), (0, 0, 0, 0, 0, 0, 0))
  Base.GC.gc()
  @time SymmetricPolynomials.decompose(tmp)
  Base.GC.gc()
  @time SymmetricPolynomials.decompose2(tmp)
  Base.GC.gc()
  @time SymmetricPolynomials.decompose3(tmp)
  Base.GC.gc()
  @time SymmetricPolynomials.decompose4(tmp)
  Base.GC.gc()
  @time SymmetricPolynomials.decompose5(tmp)
  Base.GC.gc()
  @time SymmetricPolynomials.decompose6(tmp)
  Base.GC.gc()
  @time SymmetricPolynomials.decompose7(tmp)
  Base.GC.gc()
  @time SymmetricPolynomials.decompose8(tmp)
  Base.GC.gc()
  @time SymmetricPolynomials.decompose9(tmp)
end

x = SymmetricPolynomials.semi_elementary_monomial((0, 0, 0, 0, 0, 0, 0, 50), (0, 0, 0, 0, 0, 0, 0, 0))
@timev SymmetricPolynomials.decompose3(x);

tmp = SymmetricPolynomials.decompose3(x);

@profview SymmetricPolynomials.decompose3(x)
@profview SymmetricPolynomials.decompose7(x)
@profview SymmetricPolynomials.decompose9(x)

iter = SymmetricPolynomials.ways_place_containers_iterator([2, 1, 3, 4, 5], 11)

cnt = 0
for way in iter
  println(way)
  cnt += 1
end
println(cnt)

x = SymmetricPolynomials.semi_elementary_monomial((0, 0, 0, 0, 10), (0, 0, 0, 0, 0))
p = SymmetricPolynomials.semi_elementary_polynomial(5)
p.terms[x] = 1

using DataStructures

highest_terms = [DataStructures.last(p.terms).first]
cur_semi_token = DataStructures.lastindex(p.terms)
while cur_semi_token != DataStructures.startof(p.terms)
  cur_semi_token = DataStructures.regress((p.terms, cur_semi_token))
  key = DataStructures.deref_key((p.terms, cur_semi_token))
  key.sp_term == highest_terms[1].sp_term ? push!(highest_terms, key) : break
end
SymmetricPolynomials.lower_order_representation3(p, highest_terms)
p.terms

x = SymmetricPolynomials.semi_elementary_monomial((0, 0, 0, 0, 10), (0, 0, 0, 0, 0))
D = SymmetricPolynomials.polynomial_decomposer(5)
D.polynomial.terms[x] = 1

empty!(D.highest_terms)
push!(D.highest_terms, DataStructures.last(D.polynomial.terms).first)
cur_semi_token = DataStructures.lastindex(D.polynomial.terms)
while cur_semi_token != DataStructures.startof(D.polynomial.terms)
  cur_semi_token = DataStructures.regress((D.polynomial.terms, cur_semi_token))
  key = DataStructures.deref_key((D.polynomial.terms, cur_semi_token))
  key.sp_term == D.highest_terms[1].sp_term ? push!(D.highest_terms, key) : break
end

N = 5
sp = D.highest_terms[1].sp_term
empty!(D.original_coefficients)
for x in D.highest_terms
  push!(D.original_coefficients, D.polynomial.terms[x])
end
num_highest = N - findfirst(i -> i == sp[end], sp) + 1
empty!(D.factor)
for i = 1:N
  push!(D.factor, i > N - num_highest ? sp[i] - 1 : sp[i])
end


SymmetricPolynomials.lower_order_representation5(D)
D.polynomial.terms