#==============================================================================#
# Description
# -----------
#   This program is to perform the ED study of the time evolution of the PXP - mZ
#   realized under the open boundary condition with fixed system size L=7.
#
#   Note we use -1 and 1 as spin basis here (instead of -1/2 and 1/2).
#   Therefore, the scaling variables can differ from those in the paper by a
#   constant factor. For example, m_s here is two times that in the paper.
#
#   Global variables that need to be set by users are: mass, eigvals, eigvects
#
#   ## Usage:
#      1. run check_relaxation_time(tmax) to extract the relaxation time
#      2. set t0 > relaxation time in getZ2SteadyValueList() and then run
#      getZ2SteadyValueList4mList(mList) to get the desired steady value for a set
#      of mass in mList
#
#   Version:  1.0
#   Created:  2022-03-20 14:49
#    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#==============================================================================#
using Random
using OffsetArrays
using SparseArrays
using LinearAlgebra
using DelimitedFiles

#=============================================================================#
# settings and utility functions
#=============================================================================#
Random.seed!(1234);
const ⊗ = kron
#------------------------------------------------------------------------------
# matlab style sparse identity matrix
#------------------------------------------------------------------------------
function speye(n::Integer, ::Type{T}=Float64) where T <: Union{Float64, ComplexF64}
    if n < 1
        error("speye(): argument error")
    end
    return sparse(T(1)*I, n, n)
end
#------------------------------------------------------------------------------
# name explains
#------------------------------------------------------------------------------
function mean_std(xList::Vector{Float64})
    sum, sum2 = 0.0, 0.0
    N = length(xList)
    for x in xList
        sum  += x
        sum2 += x^2
    end
    return [sum/N sqrt((sum2/N - (sum/N)^2)/(N-1))]
end
#------------------------------------------------------------------------------
# time evolve a state for time dt
#------------------------------------------------------------------------------
function timeEvolve!(psi::Vector{ComplexF64}, dt::Real, eigvals::Vector{Float64}, eigvects::Matrix{T}) where T <: Union{Float64, ComplexF64}
    if (length(psi) != length(eigvals)) || (length(psi) != size(eigvects, 1)) || (size(eigvects, 1) != size(eigvects, 2))
        error("timeEvolve!(): input argument size error")
    end
    basisN = length(psi)
    cList = [eigvects[:, n]'*psi for n in 1:basisN]
    psi .= 0.0
    for n in 1:basisN
        psi[:] += cList[n]*exp(-im*eigvals[n]*dt)*eigvects[:, n]
    end
end
#------------------------------------------------------------------------------
# get the time evolved values of a list of observables
# t0: initial time; dt: time step; Nt: total # of time steps;
# psi: initial state at t=0; funcs... : function list for observable list
#------------------------------------------------------------------------------
function getTimeEvolutionList!(psi::Vector{ComplexF64}, dt::Real, Nt::Int, eigvals::Vector{Float64}, eigvects::Matrix{T}, funs...; t0=0.0::Real) where T <: Union{Float64, ComplexF64}
    obsN  = length(funs)
    tList = zeros(Nt+1)
    obsLists = zeros(Nt+1, obsN)
    for n in 0:Nt
        if n == 0
            timeEvolve!(psi, t0, eigvals, eigvects)
        else
            timeEvolve!(psi, dt, eigvals, eigvects)
        end
        tList[n+1] = t0 + dt*n
        for io in 1:obsN
            obsLists[n+1, io] = funs[io](psi)
        end
    end
    return [tList obsLists]
end
#------------------------------------------------------------------------------
# same as before, except the time list takes a geometric series t = t0*a^it
#------------------------------------------------------------------------------
function getPowerTimeEvolutionList!(psi::Vector{ComplexF64}, tMax::Real, eigvals::Vector{Float64}, eigvects::Matrix{T}, funs...; t0=1.0::Real, a=1.02::Real) where T <: Union{Float64, ComplexF64}
    obsN  = length(funs)
    Nt = ceil(Int, log(tMax/t0)/log(a))
    tList = zeros(Nt+1)
    obsLists = zeros(Nt+1, obsN)
    for n in 0:Nt
        if n == 0
            timeEvolve!(psi, t0, eigvals, eigvects)
        else
            dt = tList[n]*(a - 1)
            timeEvolve!(psi, dt, eigvals, eigvects)
        end
        tList[n+1] = t0*(a^n)
        for io in 1:obsN
            obsLists[n+1, io] = funs[io](psi)
        end
    end
    return [tList obsLists]
end

#=============================================================================#
# system parameters
#=============================================================================#
const L = 7
const FibonacciTable = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377,
610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393]

#=============================================================================#
# Filter State and prepare state_istate
#=============================================================================#
function getStateMapTable()
    stateN = FibonacciTable[L+2]
    state_istate = Vector{Int64}(undef, stateN)
    istate_state = OffsetVector(fill(-1, 2^L), 0:2^L-1)
    invalid_stateN = 0
    istate = 0
    for state in 0:2^L-1
        stateBits = digits(state, base = 2, pad=L) |> reverse
        #---------------------------------------------------------------------
        # site :   1 2 3    ... L-1 L
        #          _________________
        #         |_|_|_|_| ... |_|_|
        #
        # bit :  high weight     low weight
        #---------------------------------------------------------------------
        valid_state = true
        for i in 1:L-1
            if stateBits[i] == 1 && stateBits[i+1] == 1
                valid_state = false
                break
            end
        end
        if valid_state
            istate += 1
            state_istate[istate] = state
            istate_state[state]  = istate
        else
            invalid_stateN += 1
        end
    end
    if stateN + invalid_stateN != 2^L
        error("stateN + invalid_stateN != 2^L")
    end
    return stateN, state_istate, istate_state
end
const stateN, state_istate, istate_state = getStateMapTable()

#=============================================================================#
# get PXP -mZ Hamiltonian: for the experiment system, the boundary term should
# be treated as XP and PX.
#=============================================================================#
function getPXPHam(mass::Float64)
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
# return the value of m, m_s, and Ebar in a given basis conf
#=============================================================================#
function getmTable(state_istate::Vector{Int})
    mTable = zeros(Float64, stateN)
    for i in 1:stateN
        mTable[i] = 2*sum(digits(state_istate[i], base=2, pad=L))/L - 1
    end
    return mTable
end
const mTable = getmTable(state_istate)
#
function get_m_s_Table(state_istate::Vector{Int})
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
const m_s_Table = get_m_s_Table(state_istate)
#------------------------------------------------------------------------------
# not to be confused with electric field (not used, can be neglected)
#------------------------------------------------------------------------------
function getEbarTable(state_istate::Vector{Int})
    EbarTable = zeros(stateN)
    for i in 1:stateN
        EA = 0; EB = 0
        bits = (digits(state_istate[i], base=2, pad=L) |> reverse)
        LB = L ÷ 2; LA = L - LB
        for n in 1:L
            if isodd(n)
                EA += bits[n]
            else
                EB += bits[n]
            end
        end
        EbarTable[i] = EA/LA - EB/LB
    end
    return EbarTable
end
const EbarTable = getEbarTable(state_istate)

#=============================================================================#
# calculate the magnetization, staggered magnetization, and Ebar
#=============================================================================#
function get_m(psi::Vector{T}) where T <: Union{Float64, ComplexF64}
    m= 0.0
    for i in 1:length(psi)
        m+= mTable[i]*abs(psi[i])^2
    end
    return m
end
function get_m_s(psi::Vector{T}) where T <: Union{Float64, ComplexF64}
    m_s = 0.0
    for i in 1:length(psi)
        m_s += m_s_Table[i]*abs(psi[i])^2
    end
    return m_s
end
function get_Ebar(psi::Vector{T}) where T <: Union{Float64, ComplexF64}
    Ebar = 0.0
    for i in 1:length(psi)
        Ebar += EbarTable[i]*abs(psi[i])^2
    end
    return Ebar
end

#=============================================================================#
# time evolution part to benchmark experiment at early stage; This function is
# also need for get the time evolution data in [0, 1] in Fig S8
#=============================================================================#
function doZ2Dynamics(t0::Real, dt::Real, tax::Real)
    psi = zeros(ComplexF64, stateN)
    psi[end] = 1.0   # prepare Z2 state
    Nt = ceil(Int, tax/dt)
    res = getTimeEvolutionList!(psi, dt, Nt, eigvals, eigvects, get_m, get_m_s, get_Ebar; t0=t0)
    run(`mkdir -p ./data`)
    writedlm("./data/Z2dynamics_L$(L)_m=$(mass).dat", res)
end
# global mass = 0.6
# global eigvals, eigvects = eigen(Symmetric(Matrix(getPXPHam(mass))))
# doZ2Dynamics(0, 0.1, 1)

#=============================================================================#
# power time evolve from [0, tmax] to extract the relaxation time
#=============================================================================#
function doZ2PowerDynamics(tax::Real)
    psi = zeros(ComplexF64, stateN)
    psi[end] = 1.0   # prepare Z2 state
    res = getPowerTimeEvolutionList!(psi, tax, eigvals, eigvects, get_m, get_m_s, get_Ebar)
    run(`mkdir -p ./data/Z2Powerdynamics`)
    writedlm("./data/Z2Powerdynamics/Z2Powerdynamics_L$(L)_m=$(mass)_tmax=$(Int(tax)).dat", res)
end
#
function check_relaxation_time(tmax)
    #--------------------------------------------------------------------------
    # this function is used to determine the relaxation time since it is
    # economical: few number of point/computation to reach big time e.g. 1E8
    #--------------------------------------------------------------------------
    for m in 0.0:0.1:1.0
        global mass = m
        global eigvals, eigvects = eigen(Symmetric(Matrix(getPXPHam(mass))))
        doZ2PowerDynamics(tmax)
    end
end

# check_relaxation_time(1E5)

#=============================================================================#
# get long time steady value by sampling values at random time in [t0, 2*t0]
# here t0 is determined from data produced in check_relaxation_time(),
# refer to ./plots/relaxation_time.py for furthur details
#=============================================================================#
function getZ2SteadyValues(t0::Real, dt::Real, Nt::Int)
    psi = zeros(ComplexF64, stateN)
    psi[end] = 1.0   # prepare Z2 state
    funs = [get_m, get_m_s, get_Ebar]; obsN = length(funs)
    obsLists = zeros(Nt+1, obsN)
    for it in 0:Nt
        if it == 0
            timeEvolve!(psi, t0, eigvals, eigvects)
        else
            timeEvolve!(psi, dt*rand(), eigvals, eigvects)
        end
        for io in 1:obsN
            obsLists[it+1, io] = funs[io](psi)
        end
    end
    return [mean_std(obsLists[:, io]) for io in 1:obsN]
end

#------------------------------------------------------------------------------
# start from t0, time evolve for a time of t0; keep increasing t0 until see the
# steady pattern appear
#------------------------------------------------------------------------------
function getZ2SteadyValueList()
    mList, m_sList, EbarList = [], [], []
    std1, std2, std3 = [], [], []
    #--------------------------------------------------------------------------
    t0List = [1000, 2000, 4000, 8000]; dt = 1
    # t0List = [1000, 2000, 4000, 8000, 16000, 32000, 64000]; dt = 1
    #--------------------------------------------------------------------------
    for t0 in t0List
        Nt = ceil(Int, t0/dt)
        res = getZ2SteadyValues(t0, dt, Nt)
        push!(mList, res[1][1]);      push!(std1, res[1][2])
        push!(m_sList, res[2][1]);    push!(std2, res[2][2])
        push!(EbarList, res[3][1]);   push!(std3, res[3][2])
    end
    run(`mkdir -p ./data/steady`)
    open("./data/steady/Z2_steady_L$(L)_m=$(mass).dat", "w") do io
        writedlm(io, [t0List mList std1 m_sList std2 EbarList std3])
    end
end
function getZ2SteadyValueList4mList(mList::Vector{Float64})
    for m in mList
        global mass = m
        global eigvals, eigvects = eigen(Symmetric(Matrix(getPXPHam(mass))))
        getZ2SteadyValueList()
    end
end

# getZ2SteadyValueList4mList(collect(0:0.1:1.0))
# getZ2SteadyValueList4mList([0.26, 0.39, 0.52, 0.65, 0.78, 0.91])
