#==============================================================================#
# Description
# -----------
#   This program use ED to calculate the thermal values of various quantities of
#   the PXP - mZ under the open boundary condition with fixed system size L=7.
#
#   Note we use -1 and 1 as spin basis here (instead of -1/2 and 1/2).
#   Therefore, the scaling variables can differ from those in the paper by a
#   constant factor. For example, m_s here is two times that in the paper.
#
#   Version:  1.0
#   Created:  2022-03-20 14:51
#    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#==============================================================================#
using OffsetArrays
using SparseArrays
using LinearAlgebra
using DelimitedFiles

#=============================================================================#
# settings and utility functions
#=============================================================================#
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
# personal implementation of Newton bisection method to get f(x) = y* to an
# accuracy or std_x for xm and std_y for y; x0: initial guess; step: initial
# step to locate a region; Nmax: maximum number of iterations.
# The function f(x) has to be monotonic for it to work properly
#------------------------------------------------------------------------------
function doNewtonBisect(x0::T, ystar::T, f::Function; step=0.1, std_x=1E-6, std_y=Inf, Nmax=1000) where T <: Real
    if step < 0
        error("doNewtonBisect(): step must be bigger than zero!")
    end
    y0 = f(x0)
    x1 = x0 + step
    y1 = f(x1)
    slope = (y1 > y0 ? 1 : -1)
    if slope == 1
        if ystar > y0 && ystar > y1
            while ystar > y1
                step *= 2
                x0, y0 = x1, y1
                x1 = x1 + step
                y1 = f(x1)
            end
        elseif ystar < y0 && ystar < y1
            while ystar < y0
                step *= 2
                x1, y1 = x0, y0
                x0 = x0 - step
                y0 = f(x0)
            end
        end
        xa, ya, xb, yb = x0, y0, x1, y1
    else
        if ystar > y0 && ystar > y1
            while ystar > y0
                step *= 2
                x1, y1 = x0, y0
                x0 = x0 - step
                y0 = f(x0)
            end
        elseif ystar < y0 && ystar < y1
            while ystar < y1
                step *= 2
                x0, y0 = x1, y1
                x1 = x1 + step
                y1 = f(x1)
            end
        end
        xa, ya, xb, yb = x1, y1, x0, y0
    end
    itime = 0
    #------------------------------------------------------------------------------
    # make sure ya < yb to facillate robust programming; although this should
    # be guranteed from the above code
    if ya > yb
        xa, xb = xb, xa
        ya, yb = yb, ya
    end
    #------------------------------------------------------------------------------
    while true
        if itime > Nmax
            error("doNewtonBisect(): bisection performed $(Nmax) times without convergence!")
        end
        x  = (xa + xb)/2
        dx = abs(xb - xa)/2
        y  = f(x)
        dy = abs(y - ystar)
        if dx < std_x && dy < std_y
            return x, dx, dy
        elseif y > ystar
            xb = x
        else
            xa = x
        end
        itime += 1
    end
end
#-------------------------------------------------------------------------#
# calcualte the temperature corresponding to the energy
#-------------------------------------------------------------------------#
function getThermalBeta(E::Real, eigvals::Vector{Float64}, beta0=1.0::Real; std=1E-6)
    function func(beta::Real)
        Elist = (issorted(eigvals) ? eigvals : sort(eigvals))
        if beta < 0
            # to make sure in the for sum below p is decreasing from 1
            Elist = reverse(Elist)
        end
        E_0 = Elist[1]
        Z, Esum = 0.0, 0.0
        for E in Elist
            p = exp(-beta*(E-E_0))
            Z += p
            Esum += E*p
        end
        E = Esum/Z
        return E
    end
    beta, std_x, std_y = doNewtonBisect(beta0, E, func; step=0.1, std_x=std, std_y=std)
    if beta < 0
        println("!!!==== getThermalBeta(): beta < 0 ====!!!")
    end
    return beta
end
#-------------------------------------------------------------------------#
# calcualte the thermal average: OList is the list of ⟨ψ_n|Ô|ψ_n⟩
#-------------------------------------------------------------------------#
function getThermalAverage(beta::Real, eigvals::Vector{Float64}, OList::Vector{Float64})
    if beta < 0
        println("!!!==== getThermalAverage(): input beta < 0 ====!!!")
    end
    if length(eigvals) != length(OList)
        error("getThermalAverage(): input argument size error!")
    end
    Elist = (issorted(eigvals) ? eigvals : sort(eigvals))
    if beta < 0
        # to make sure in the for sum below p is decreasing from 1
        # do NOT use reverse!() to modify data
        Elist = reverse(Elist)
        OList = reverse(OList)
    end
    E_0 = Elist[1]
    Z, sum = 0.0, 0.0
    for (E, O) in zip(Elist, OList)
        p = exp(-beta*(E-E_0))
        Z += p
        sum += p*O
    end
    return sum/Z
end

#=============================================================================#
# system parameters
#=============================================================================#
const L = 7
const FibonacciTable = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610,
987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393]

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
function getPXPHam(mass::Real)
    # Note non-standard σ_z has been used here to conform to my Fortran program
    P   = spzeros(2, 2); P[1, 1] = 1                # P = (1 - σ_z)/2
    X   = spzeros(2, 2); X[1, 2] = 1; X[2, 1] = 1
    Z   = spzeros(2, 2); Z[1, 1] =-1; Z[2, 2] = 1
    XP  = kron(X, P); PX  = kron(P, X)
    PXP = kron(kron(P, X), P)
    Ham = spzeros(Float64, 2^L, 2^L)
    for i = 1:L
        Ham .+= speye(2^(i-1)) ⊗ (-mass*Z) ⊗ speye(2^(L-i))
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
# const Ham = getPXPHam(mass)
# const eigvals, eigvects = eigen(Symmetric(Matrix(Ham)))

#=============================================================================#
# return the number of N↑ spins in a given basis conf.
#=============================================================================#
function getmTable(state_istate::Vector{Int})
    mTable = zeros(Float64, stateN)
    for i in 1:stateN
        mTable[i] = 2*sum(digits(state_istate[i], base=2, pad=L))/L - 1
    end
    return mTable
end
const mTable = getmTable(state_istate)
#------------------------------------------------------------------------------
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
        bits = (digits(state_istate[i], base=2, pad=L) |> reverse)
        LB = L ÷ 2; LA = L - LB
        EA = 0; EB = 0
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
# calculate expectation of various quantities
#=============================================================================#
function get_m(psi::Vector{T}) where T <: Union{Float64, ComplexF64}
    if abs(norm(psi) - 1.0) > 1E-4
        error("get_Ebar(): psi not normalized")
    end
    mag = 0.0
    for i in 1:length(psi)
        mag += mTable[i]*abs(psi[i])^2
    end
    return mag
end
#------------------------------------------------------------------------------
function get_m_s(psi::Vector{T}) where T <: Union{Float64, ComplexF64}
    if abs(norm(psi) - 1.0) > 1E-4
        error("get_Ebar(): psi not normalized")
    end
    m_s = 0.0
    for i in 1:length(psi)
        m_s += abs(psi[i])^2*m_s_Table[i]
    end
    return m_s
end
#------------------------------------------------------------------------------
function get_Ebar(psi::Vector{T}) where T <: Union{Float64, ComplexF64}
    if abs(norm(psi) - 1.0) > 1E-4
        error("get_Ebar(): psi not normalized")
    end
    Ebar = 0.0
    for i in 1:length(psi)
        Ebar += EbarTable[i]*abs(psi[i])^2
    end
    return Ebar
end

#=============================================================================#
# calcualte thermal values of various quantities
#=============================================================================#
function getZ2_thermal_obs(mass::Real)
    Ham = getPXPHam(mass)
    eigvals, eigvects = eigen(Symmetric(Matrix(Ham)))
    psi = zeros(stateN); psi[end] = 1.0
    beta = getThermalBeta(psi'*Ham*psi, eigvals, 0.0001)
    println("beta = $(beta)")
    OList  = [get_Ebar(eigvects[:, n]) for n in 1:stateN]
    Ebar   = getThermalAverage(beta, eigvals, OList)
    OList  = [get_m(eigvects[:, n]) for n in 1:stateN]
    m      = getThermalAverage(beta, eigvals, OList)
    OList  = [get_m_s(eigvects[:, n]) for n in 1:stateN]
    m_s    = getThermalAverage(beta, eigvals, OList)
    return round(Ebar, digits=6), round(m, digits=6), round(m_s, digits=6)
end
#------------------------------------------------------------------------------
function getZ2_thermal_obs_list()
    massList = collect(0.0:0.1:1.0)
    mList    = zeros(length(massList))
    m_sList  = zeros(length(massList))
    EbarList = zeros(length(massList))
    for (imass, mass) in enumerate(massList)
        EbarList[imass], mList[imass], m_sList[imass] = getZ2_thermal_obs(mass)
    end
    run(`mkdir -p ./data/Z2Thermal`)
    writedlm("./data/Z2Thermal/Z2_thermal_obs_L$(L).dat", [massList EbarList mList m_sList])
end
#------------------------------------------------------------------------------
getZ2_thermal_obs_list()
