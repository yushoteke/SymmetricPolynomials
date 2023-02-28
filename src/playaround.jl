using TupleTools
using DataStructures
using TimerOutputs

include("utilities.jl")
include("semi elementary monomial.jl")
include("semi elementary polynomial.jl")
include("evaluate.jl")


x = semi_elementary_monomial((0, 0, 0, 10), (0, 0, 0, 0))

p13_tmp = decompose13(x)
p13 = poly4_to_poly(p13_tmp)
p12_tmp = decompose12(x)
p12 = poly3_to_poly(p12_tmp)
p11_tmp = decompose11(x)
p11 = poly3_to_poly(p11_tmp)
p10_tmp = decompose10(x);
p10 = poly2_to_poly(p10_tmp);
p9 = decompose9(x);
p8 = decompose8(x);
p7 = decompose7(x);
p6 = decompose6(x);
p5 = decompose5(x);
p4 = decompose4(x);
p3 = decompose3(x);
p2 = decompose2(x);
p1 = decompose(x);

println(p1 == p2 == p3 == p4 == p5 == p6 == p7 == p8 == p9 == p10 == p11 == p12 == p13)

for i = 20:10:60
  println("i = $(i)")
  tmp = semi_elementary_monomial((0, 0, 0, 0, 0, 0, i), (0, 0, 0, 0, 0, 0, 0))
  @time res1 = decompose(tmp)
  @time res2 = decompose2(tmp)
  @time res3 = decompose3(tmp)
  @time res4 = decompose4(tmp)
  @time res5 = decompose5(tmp)
  @time res6 = decompose6(tmp)
  @time res7 = decompose7(tmp)
  @time res8 = decompose8(tmp)
  @time res9 = decompose9(tmp)
  @time res10_tmp = decompose10(tmp)
  res10 = poly2_to_poly(res10_tmp)
  @time res11_tmp = decompose11(tmp)
  res11 = poly3_to_poly(res11_tmp)
  @time res12_tmp = decompose12(tmp)
  res12 = poly3_to_poly(res12_tmp)
  @time res13_tmp = decompose13(tmp)
  res13 = poly4_to_poly(res13_tmp)
  println("results are same:", res1 == res2 == res3 == res4 == res5 == res6 ==
                               res8 == res9 == res10 == res11 == res12 == res13)
end

for i = 10:10:80
  println("i=", i)
  tmp = semi_elementary_monomial((0, 0, 0, 0, 0, 0, i), (0, 0, 0, 0, 0, 0, 0))
  decompose10(tmp)
end


for i = 70:10:70
  println("i = $(i)")
  tmp = semi_elementary_monomial((0, 0, 0, 0, 0, 0, i), (0, 0, 0, 0, 0, 0, 0))
  Base.GC.gc()
  @time decompose(tmp)
  Base.GC.gc()
  @time decompose2(tmp)
  Base.GC.gc()
  @time decompose3(tmp)
  Base.GC.gc()
  @time decompose4(tmp)
  Base.GC.gc()
  @time decompose5(tmp)
  Base.GC.gc()
  @time decompose6(tmp)
  Base.GC.gc()
  @time decompose7(tmp)
  Base.GC.gc()
  @time decompose8(tmp)
  Base.GC.gc()
  @time decompose9(tmp)
  Base.GC.gc()
  @time decompose10(tmp)
  Base.GC.gc()
  @time decompose11(tmp)
  Base.GC.gc()
  @time decompose12(tmp)
end


x = semi_elementary_monomial((0, 0, 0, 0, 0, 0, 22), (0, 0, 0, 0, 0, 0, 0))
@timev decompose3(x);

tmp = decompose3(x);

@profview decompose3(x)
@profview decompose7(x)
@profview decompose9(x)
@profview decompose10(x)
@profview decompose12(x)
@profview decompose13(x)
@code_warntype decompose11(x)
@timev decompose12_timeit(x);
@timev decompose13_timeit(x);

iter = ways_place_containers_iterator([2, 1, 3, 4, 5], 11)

cnt = 0
for way in iter
  println(way)
  cnt += 1
end
println(cnt)

x = semi_elementary_monomial((0, 0, 0, 0, 10), (0, 0, 0, 0, 0))
p = semi_elementary_polynomial(5)
p.terms[x] = 1

using DataStructures

highest_terms = [DataStructures.last(p.terms).first]
cur_semi_token = DataStructures.lastindex(p.terms)
while cur_semi_token != DataStructures.startof(p.terms)
  cur_semi_token = DataStructures.regress((p.terms, cur_semi_token))
  key = DataStructures.deref_key((p.terms, cur_semi_token))
  key.sp_term == highest_terms[1].sp_term ? push!(highest_terms, key) : break
end
lower_order_representation3(p, highest_terms)
p.terms

x = semi_elementary_monomial((0, 0, 0, 0, 10), (0, 0, 0, 0, 0))
D = polynomial_decomposer(5)
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


lower_order_representation5(D)
D.polynomial.terms