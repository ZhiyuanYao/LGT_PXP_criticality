#==============================================================================#
# Description
# -----------
#   This program performs the ED study of the PXP - mZ model under the open
#   boundary condition.
#
#   Note we use -1 and 1 as spin basis here (instead of -1/2 and 1/2).
#   Therefore, the scaling variables can differ from those in the paper by a
#   constant factor. For example, m_s here is two times that in the paper.
#
#   Version:  1.0
#   Created:  2022-03-20 10:49
#    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#==============================================================================#
using Arpack
using OffsetArrays
using SparseArrays
using LinearAlgebra
using DelimitedFiles

#------------------------------------------------------------------------------
# settings and utility functions
#------------------------------------------------------------------------------
const ⊗ = kron

# matlab style sparse identity matrix
function speye(n::Integer, ::Type{T}=Float64) where T <: Union{Float64, ComplexF64}
    if n < 1
        error("speye(): argument error")
    end
    return sparse(T(1)*I, n, n)
end

#------------------------------------------------------------------------------
# system parameters
#------------------------------------------------------------------------------
const FibonacciTable = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610,
987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393]

#=============================================================================#
# Filter State and prepare state_istate
#=============================================================================#
function getStateMapTable(L::Int)
    stateN = FibonacciTable[L+2]
    state_istate = zeros(Int, stateN)
    # istate_state = OffsetVector(fill(-1, 2^L), 0:2^L-1)
    invalid_stateN = 0
    istate = 0
    for state in 0:2^L-1
        bits = digits(state, base = 2, pad=L) |> reverse
        valid_state = true
        for n in 1:L-1
            if bits[n]*bits[n+1] == 1
                valid_state = false
                break
            end
        end
        if valid_state
            istate += 1
            state_istate[istate] = state
            # istate_state[state]  = istate
        else
            invalid_stateN += 1
        end
    end
    if stateN + invalid_stateN != 2^L
        error("stateN + invalid_stateN != 2^L")
    else
        println("L = $(L)   stateN = $(stateN)")
    end
    return stateN, state_istate #, istate_state
end

#=============================================================================#
# return the value of scaling variable in a given basis
#=============================================================================#
function get_m_s_Table(L::Int, state_istate::Vector{Int})
    stateN = length(state_istate)
    m_s_Table = zeros(stateN)
    for i in 1:stateN
        bits = (digits(state_istate[i], base=2, pad=L) |> reverse)
        m = 0
        for n in 1:L
            m += (-1)^(n-1)*(2*bits[n]-1)
        end
        m_s_Table[i] = m/L
    end
    return m_s_Table
end

# corrlation table (∑ s_i^z*s_{i+1}^z)/(L-1)
function get_corr_Table(L::Int, state_istate::Vector{Int})
    stateN = length(state_istate)
    corr_Table = zeros(stateN)
    for i in 1:stateN
        bits = (digits(state_istate[i], base=2, pad=L) |> reverse)
        corr = 0
        for n in 1:L-1
            corr += (2*bits[n]-1)*(2*bits[n+1]-1)
        end
        corr_Table[i] = corr/(L-1)
    end
    return corr_Table
end

#=============================================================================#
# get PXP -mZ Hamiltonian: for the experiment system with open boundary
# condition, the boundary terms are XP and PX, respectively.
#=============================================================================#
function getPXPHam(L::Int, mass::Real, state_istate::Vector{Int})
    # Note non-standard σ_z has been used here to conform to my Fortran program
    P   = spzeros(2, 2); P[1, 1] = 1                # P = (1 - σ_z)/2
    X   = spzeros(2, 2); X[1, 2] = 1; X[2, 1] = 1
    Z   = spzeros(2, 2); Z[1, 1] =-1; Z[2, 2] = 1
    XP  = kron(X, P); PX  = kron(P, X)
    PXP = kron(kron(P, X), P)
    Ham = spzeros(Float64, 2^L, 2^L)
    for i = 1:L
        Ham .-= speye(2^(i-1)) ⊗ (mass*Z) ⊗ speye(2^(L-i))
        if i == 1 || i == L
            continue
        end
        Ham .+= speye(2^(i-2)) ⊗ PXP ⊗ speye(2^(L-i-1))
    end
    Ham .+= XP ⊗ speye(2^(L-2))
    Ham .+= speye(2^(L-2)) ⊗ PX
    index_istate = broadcast(+, state_istate, 1)
    Ham = Ham[index_istate, :]
    Ham = Ham[:, index_istate]
    return Ham
end


#=============================================================================#
# calculate the expectation values of various scaling variables in the ground
# state (not all of these variables are used in extracting g_c, however)
#=============================================================================#
function getScalingVariables(psi::Vector{T}, corr_Table::Vector{Float64}, m_s_Table::Vector{Float64}) where T <: Union{Float64, ComplexF64}
    corr, m_s, m_s_abs, m_s2, m_s4 = 0.0, 0.0, 0.0, 0.0, 0.0
    for i in 1:length(psi)
        p = abs(psi[i])^2
        corr    += p*corr_Table[i]
        m_s     += p*m_s_Table[i]
        m_s_abs += p*abs(m_s_Table[i])
        m_s2    += p*(m_s_Table[i])^2
        m_s4    += p*(m_s_Table[i])^4
    end
    return corr, m_s, m_s_abs, m_s2, m_s4
end

# calculate the scaling variables for various mass and system size
function getScalingVariableArray(massList::Vector{Float64}, LList::Vector{Int})
    for L in LList
        #----------------------------------------------------------------------
        stateN, state_istate = getStateMapTable(L)
        corr_Table  = get_corr_Table(L, state_istate)
        m_s_Table   = get_m_s_Table(L, state_istate)
        #----------------------------------------------------------------------
        results = zeros(length(massList), 5)
        for (imass, mass) in enumerate(massList)
            # for size L > 20, consider exporting the sparse matrix to Matlab
            # to calculate psi: the eigs() in Arpack is too slow.
            eigvals, eigvects = eigs(Symmetric(Matrix(getPXPHam(L, mass, state_istate))), nev=6, which=:SR)
            psi = eigvects[:, 1]
            results[imass, :] = collect(getScalingVariables(psi, corr_Table, m_s_Table))
        end
        writedlm("./data/FSS_L$(L).dat", [massList results])
    end
end

const massList = collect(0.2:0.02:1.0)
getScalingVariableArray(massList, [5, 7, 9])
