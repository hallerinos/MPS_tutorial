using TensorOperations, LinearAlgebra, BenchmarkTools, UnicodePlots

## some basic contractions, implemented with the tensor macro
#############################################################
# create some random matrices
A, B, C, D, E, F = [rand(5,5) for i=1:6]

# matrix-matrix multiplication
@tensor C′[i,j] := A[i,k]*B[k,j]
@show C′ ≈ A*B

# trace of 6 matrices
trace = @tensor A[a,b]*B[b,c]*C[c,d]*D[d,e]*E[e,f]*F[f,a]
@show trace ≈ tr(A*B*C*D*E*F)

## many pitfalls -- contraction order matters a lot!
####################################################
m = 64
A, B, C = rand(m,m,m), rand(m,m,m), rand(m,m)

## contraction example of Fig. 1
################################
# the contraction over a given index (here: β₁) reads explicitly
contr_β₁(A,C) = reshape(reshape(A, (m^2, m)) * C, (m,m,m))
# the tensor macro performs reshaping operations automatically
@tensor A′[γ₁,α₁,β₂] := A[γ₁,α₁,β₁]*C[β₁,β₂]
@show A′ ≈ contr_β₁(A,C)

## comparison between contraction sequence 1 and 2
cs_1(A,B,C) = @tensor D1[:] := A[-1,3,2]*B[-2,1,3]*C[2,1]  # this one costs 𝒪(m⁴)
cs_2(A,B,C) = @tensor D2[:] := A[-1,1,2]*B[-2,3,1]*C[2,3]  # this one costs 𝒪(m⁵)
cs_auto(A,B,C) = @tensor D[γ₁,γ₂] := A[γ₁,α₁,β₁]*B[γ₂,β₂,α₁]*C[β₁,β₂]  # this will be evaluated after a contraction sequence optimization
@show cs_auto(A,B,C) ≈ cs_1(A,B,C) ≈ cs_2(A,B,C)  # all contractions are equivalent...

## ... but they differ in runtime due to different overall complexity
@time cs_auto(A,B,C)
@time cs_1(A,B,C)
@time cs_2(A,B,C)

## ... more detailed analysis of the overall runtime (takes some time)
b1 = median(@benchmark cs_1(A,B,C))
b2 = median(@benchmark cs_2(A,B,C))
ba = median(@benchmark cs_auto(A,B,C))
@show judge(b1, b2)
@show judge(ba, b2)
