using TensorOperations, LinearAlgebra
include("isingMPO.jl")

# write the following contraction to compute the energy
#######################################################
#
#   ----     -----------      ---- 
#  |    |-1-|    AC     |- 9-|    |
#  |    |    -----------     |    |
#  |    |     |       |      |    |
#  |    |     4       5      |    |
#  |    |     |       |      |    |
#  |    |    ---     ---     |    |
#  | LW |-2-| W |-8-| W |-10-| RW |
#  |    |    ---     ---     |    |
#  |    |     |       |      |    |
#  |    |     6       7      |    |
#  |    |     |       |      |    |
#  |    |    -----------     |    |
#  |    |-3-|    A̅C̅     |-11-|    |
#   ----     -----------      ---- 
#
function energy(AC, W, LW, RW)
    @tensor LW[1,2,3]*AC[1,4,5,9]*W[2,6,8,4]*W[8,7,10,5]*conj(AC)[3,6,7,11]*RW[9,10,11]
end

# write the following contraction to compute the gradient of the energy wrt conj(AC)
####################################################################################
#
#   ----     -----------     ---- 
#  |    |-1-|    AC     |-6-|    |
#  |    |    -----------    |    |
#  |    |     |       |     |    |
#  |    |     3       5     |    |
#  |    |     |       |     |    |
#  |    |    ---     ---    |    |
#  | LW |-2-| W |-4-| W |-7-| RW |
#  |    |    ---     ---    |    |
#  |    |     |       |     |    |
#  |    |    (-2)    (-3)   |    |
#  |    |     |       |     |    |
#  |    |                   |    |
#  |    |-(-1)         (-4)-|    |
#   ----                     ---- 
#
function grad(AC, W, LW, RW)
    @tensor grad[:] := LW[1,2,-1]*AC[1,3,5,6]*W[2,-2,4,3]*W[4,-3,7,5]*RW[6,7,-4]
end

# write the following contraction to absorb the LWft isometry A
###############################################################
#
#   ----     ---       
#  |    |-1-| U |--(-1)
#  |    |    ---       
#  |    |     |        
#  |    |     4        
#  |    |     |        
#  |    |    ---       
#  | LW |-2-| W |--(-2)
#  |    |    ---       
#  |    |     |        
#  |    |     5        
#  |    |     |        
#  |    |    ---       
#  |    |-3-| U̅ |--(-3)
#   ----     ---       
#
function absorbU(U, W, LW)
    @tensor newLW[:] := LW[1,2,3]*U[1,4,-1]*W[2,5,-2,4]*conj(U)[3,5,-3]
end

# write the following contraction to absorb the right isometry B
################################################################
#
#         ----     ---- 
#  (-1)--| V† |-3-|    |
#         ----    |    |
#           |     |    |
#           1     |    |
#           |     |    |
#          ---    |    |
#   (-2)--| W |-4-| RW |
#          ---    |    |
#           |     |    |
#           2     |    |
#           |     |    |
#         ----    |    |
#  (-3)--| V̅†̅ |-5-|    |
#         ----     ---- 
#
function absorbVdag(Vdag, W, RW)
    @tensor newRW[:] := Vdag[-1,1,3]*W[-2,2,4,1]*conj(Vdag)[-3,2,5]*RW[3,4,5]
end

# test if your contractions are working
function check_contractions()
    # the parameters of the Ising model
    J = 1.0
    h = 1.0

    L, W, R = isingMPO(J=J, h=h)

    LW = ones(16,size(W,1),16); LW /= norm(LW)  # initialize Left environment
    RW = ones(16,size(W,3),16); RW /= norm(RW)  # initialize right environment

    AC = rand(size(LW,1), size(W,4), size(W,4), size(RW,1)); AC /= norm(AC)  # construct random AC
    ACmat = reshape(AC, (prod(size(AC)[1:2]), prod(size(AC)[3:end])))  # reshape to matrix
    U, S, V = svd(ACmat) # extact the left and right isometries
    S = Diagonal(S)
    U = reshape(U, (size(AC,1), size(AC,2), size(S,1)))
    Vdag = reshape(V', (size(S,2), size(AC,3), size(AC,4)))
    @tensor ACnew[:] := U[-1,-2,1]*S[1,2]*Vdag[2,-3,-4]
    @show ACnew ≈ AC


    # the following block contains sanity checks
    # if you implemented the contractions correctly, the last line will evaluate true
    ################################################################################
    # this computes energy
    ene1 = energy(AC,W,LW,RW)

    # this is a block which checks absorbU and absorbVdag
    tpL = absorbU(U,W,LW)
    tpR = absorbVdag(Vdag,W,RW)
    ene2 = @tensor tpL[1,2,3]*S[1,4]*conj(S)[3,5]*tpR[4,2,5]

    # this is a block which checks grad
    tpgrad = grad(AC,W,LW,RW)
    ene3 = @tensor tpgrad[1,2,3,4]*conj(AC)[1,2,3,4]
    @show ene1 ≈ ene2 ≈ ene3
    return
end

check_contractions()  # comment after you finished this task