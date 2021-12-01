using ITensors
using ITensorsVisualization
using Plots

## matrix-matrix multiplication
i = [Index(rand(32:64), "index_$i") for i=1:3]
A = randomITensor(i[1], i[2])
B = randomITensor(i[2], i[3])
C = @visualize A * B
# cast tensors to matrices
Cmat = Matrix(C, i[1], i[3])
Amat = Matrix(A, i[1], i[2])
Bmat = Matrix(B, i[2], i[3])
# perform naive matrix-matrix multiplication
C′ = Amat*Bmat
# result will be the same...
@show Cmat ≈ C′

## trace over six matrices
i = [Index(rand(32:64), "index_$i") for i=1:7]
A, B, C, D, E = [randomITensor(i[j], i[j+1]) for j=1:5]
F = randomITensor(i[6], i[1])
@visualize A*B*C*D*E*F

## contraction sequence matters!
times1 = []
times2 = []
times3 = []
times4 = []
χs = [6*i for i in 1:10]
for χ in χs
    i = [Index(χ, "index_$i") for i=1:5]
    A = randomITensor(i[1], i[2], i[3])
    B = randomITensor(i[4], i[5], i[2])
    C = randomITensor(i[3], i[5])
    append!(times1, @elapsed (A*C)*B)  # first contract A with C, then B => 𝒪(χ⁴)
    append!(times2, @elapsed A*(B*C))  # first contract B with C, then A => 𝒪(χ⁴)
    append!(times3, @elapsed (A*B)*C)  # first contract A with B, then C => 𝒪(χ⁵)
    append!(times4, @elapsed A*B*C)    # is equivalent to (A*B)*C
end
@visualize A*B*C
plot(χs, times1, yaxis=:log, label="(A*C)*B")
plot!(χs, times2, label="A*(B*C)")
plot!(χs, times3, label="(A*B)*C")
plot!(χs, times4, label="A*B*C")
