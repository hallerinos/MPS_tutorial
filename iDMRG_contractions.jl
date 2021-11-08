using TensorOperations, LinearAlgebra

# some random tensors to implement the contractions
MPO = rand(2,3,2,3)
LE = ones(7,size(MPO,4),7); LE /= norm(LE)  # initialize left environment
RE = ones(10,size(MPO,2),10); RE /= norm(RE)  # initialize right environment

AC = rand(size(LE,1), size(MPO,1), size(MPO,1), size(RE,1)); AC /= norm(AC)  # construct random AC
ACmat = reshape(AC, (prod(size(AC)[1:2]), prod(size(AC)[3:end])))  # reshape to matrix
A, S, B = svd(ACmat) # extact the left and right isometries
@show A*Diagonal(S)*B' ≈ ACmat
A = reshape(A, (size(AC,1), size(AC,2), length(S)))
B = reshape(B, (length(S), size(AC,3), size(AC,4)))

# write the following contraction to compute the energy
#######################################################
#
#   ----    ----------    ---- 
#  |    |--|    AC    |  |    |
#  |    |   ----------   |    |
#  |    |    |      |    |    |
#  |    |   ---    ---   |    |
#  | LW |--| W |--| W |--| RW |
#  |    |   ---    ---   |    |
#  |    |    |      |    |    |
#  |    |   ----------   |    |
#  |    |--|    A̅C̅    |  |    |
#   ----    ----------    ---- 
#
function energy(AC,W,LE,RE)
    @tensor LE[a,b,c]*AC[a,d,e,f]*W[g,h,d,b]*W[i,j,e,h]*conj(AC)[c,g,i,k]*RE[f,j,k]
end
energy(AC,MPO,LE,RE);

# write the following contraction to compute the gradient of the energy wrt conj(AC)
####################################################################################
#
#   ----    ----------    ---- 
#  |    |--|    AC    |  |    |
#  |    |   ----------   |    |
#  |    |    |      |    |    |
#  |    |   ---    ---   |    |
#  | LW |--| W |--| W |--| RW |
#  |    |   ---    ---   |    |
#  |    |    |      |    |    |
#  |    |                |    |
#  |    |--            --|    |
#   ----                  ---- 
#
function grad(AC,W,LE,RE)
    @tensor grad[c, g, i, k] := LE[a,b,c]*AC[a,d,e,f]*W[g,h,d,b]*W[i,j,e,h]*RE[f,j,k]
end
grad(AC,MPO,LE,RE);

# write the following contraction to absorb the left isometry A
###############################################################
#
#   ----    ---   
#  |    |--| A |--
#  |    |   ---   
#  |    |    |    
#  |    |   ---   
#  | LW |--| W |--
#  |    |   ---   
#  |    |    |    
#  |    |   ---   
#  |    |--| A̅ |--
#   ----    ---   
#
function absorbA(A,W,LW)
    @tensor newLW[e,g,h] := LW[a,b,c]*A[a,d,e]*W[f,g,d,b]*conj(A)[c,f,h]
end
absorbA(A,MPO,LE);

# write the following contraction to absorb the right isometry B
################################################################
#
#     ---    ---- 
#  --| B |--|    |
#     ---   |    |
#      |    |    |
#     ---   |    |
#  --| W |--| RW |
#     ---   |    |
#      |    |    |
#     ---   |    |
#  --| B̅ |--|    |
#     ---    ---- 
#
function absorbB(B,W,RW)
    @tensor newRW[a,g,h] := B[a,b,c]*RW[c,d,e]*W[f,d,b,g]*conj(B)[h,f,e]
end
absorbB(B,MPO,RE);
