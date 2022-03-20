#==============================================================================#
# Description
# -----------
#   This program extrapolate to infinite system size to get an estimate g_c for
#   the infinite system size. To view all the plots, run in interactive mode.
#
#   Version:  1.0
#   Created:  2022-03-20 11:56
#    Author:  Zhiyuan Yao, zhiyuan.yao@icloud.com
# Institute:  Insitute for Advanced Study, Tsinghua University
#==============================================================================#
using Plots
using EasyFit
using DelimitedFiles

#------------------------------------------------------------------------------
# use system size L in LList and values of mass (setting t̃=1) in [m_min, m_max]
# to perform quadratic fit to determine the crossing point g_c(L)
#------------------------------------------------------------------------------
const LList = [11, 13, 15, 17, 19, 21, 23, 25]
const m_min = 0.56
const m_max = 0.7

const scaling_choice = 2  # 1: m_s_abs*L^{1/8}  2: m_s2*L^{1/4}  3: Binder ratio R

#------------------------------------------------------------------------------
# plot scaled variable against mass to make sure everything goes smoothly
#------------------------------------------------------------------------------
p = plot()
fitList = []
for L in LList
    filename = "./data/FSS_L$(L).dat"
    datas = readdlm(filename)
    xdat, ydat = Float64[], Float64[]
    for (i, x) in enumerate(datas[:, 1])
        if x >= m_min && x <= m_max
            push!(xdat, x)
            if scaling_choice == 1
                # 4 here corresponds to m_s_abs
                push!(ydat, datas[i, 4]*(L^0.125))
            elseif scaling_choice == 2
                # 5 corresponds to m_s2
                push!(ydat, datas[i, 5]*(L^0.25))
            elseif scaling_choice == 3
                # 6 corresponds to m_s4, and Binder ratio R = 1 - m_s4/(3*m_s2)
                push!(ydat, 1-datas[i, 6]/(3*datas[i, 5]^2))
            end
        end
    end
    plot!(p, xdat, ydat, markershape=:circle, label="L=$L")
    fit = fitquad(xdat, ydat)
    println(" goodness of fit = $(fit.R)")
    a, b, c, xfit, yfit = fit.a, fit.b, fit.c, fit.x, fit.y
    push!(fitList, (a, b, c, xfit, yfit))
end
p

#------------------------------------------------------------------------------
# find the intersection of different lines: ϵ stands for 1/L
#------------------------------------------------------------------------------
ϵList = Float64[]
gcList = Float64[]
for (i, L) in enumerate(LList[1:end-1])
    a1, b1, c1, xfit1, yfit1 = fitList[i]
    a2, b2, c2, xfit2, yfit2 = fitList[i+1]
    a = a1 -a2; b = b1 - b2; c = c1 - c2;
    g1 = (-b + sqrt(b^2 - 4*a*c))/(2*a)
    g2 = (-b - sqrt(b^2 - 4*a*c))/(2*a)
    gc = (abs(g1 - 0.6) < abs(g2-0.6) ? g1 : g2)
    push!(ϵList, 1/L)
    push!(gcList, gc)
end

plot(ϵList, gcList, markershape=:circle, xlim=(0, maximum(ϵList)))

#------------------------------------------------------------------------------
# extrapolate for infinite system size: note ϵList[1] is the biggest value
# the fitting result for various number of data give consistent value of g_c
#------------------------------------------------------------------------------
for n in 1:min(5, length(ϵList) - 2)   # so that we have at least 3 data to fit
    fit = fitlinear(ϵList[n:end], gcList[n:end]);
    println(round(fit.b, digits=4), "   goodness of fit = ", fit.R)
end
