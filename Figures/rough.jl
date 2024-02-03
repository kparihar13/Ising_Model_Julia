using PyPlot
using StatsBase
using JLD2
using LaTeXStrings

#=taum = []
m = []
cv = []
x = []

for L in 10:10:20
    data = jldopen("data_3D_updated_$L.jld2")

    corrm = []
    m_L = []
    cv_L = []
    x_L = []

    for i in 1:size(data["M"])[1]
        if data["M"][i][1] <= 4.2 || data["M"][i][1] >= 4.8
            tau = minimum(findall(<=(exp(-1)),autocor(abs.(data["M"][i][2][1000:9999]/L^3),1:100)))        
            samplet = maximum([100,3*tau])
            Mavg = mean(abs.(data["M"][i][2][1000:samplet:9999]))
            M2avg = mean(data["M2"][i][2][1000:samplet:9999])
            Eavg = mean(data["E"][i][2][1000:samplet:9999])
            E2avg = mean(data["E2"][i][2][1000:samplet:9999])

        elseif data["M"][i][1]>4.2 && data["M"][i][1]<4.8
            tau = minimum(findall(<=(exp(-1)),autocor(abs.(data["M"][i][2][5000:49999]/L^2),1:1000)))
            samplet = maximum([100,3*tau])
            Mavg = mean(abs.(data["M"][i][2][5000:samplet:49999]))
            M2avg = mean(data["M2"][i][2][5000:samplet:49999])
            Eavg = mean(data["E"][i][2][5000:samplet:49999])
            E2avg = mean(data["E2"][i][2][5000:samplet:49999])
        end

        push!(corrm, [data["M"][i][1],tau])
        push!(m_L,[data["M"][i][1],Mavg/L^3])
        push!(x_L,[data["M"][i][1],(1/data["M"][i][1])*(M2avg - Mavg^2)/L^3])
        push!(cv_L,[data["M"][i][1],(1/data["M"][i][1]^2)*(E2avg - Eavg^2)/L^3])
    end

    push!(taum,[L,corrm])
    push!(m,[L,m_L])
    push!(x,[L,x_L])
    push!(cv,[L,cv_L])
end


fig, ax = plt.subplots()
for i in 1:2
    xdata = []
    ydata = []
    for j in 1:size(taum[i][2])[1]
        push!(xdata,taum[i][2][j][1])
        push!(ydata,taum[i][2][j][2])
    end
    ax.scatter(xdata,ydata,label="L = $(taum[i][1])")
end
ax.set_ylabel(L"Autocorrelation time, $\tau$ ")
ax.set_xlabel(L"Temperature, T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_autocorr_3D.png")


fig, ax = plt.subplots()
for i in 1:2
    xdata = []
    ydata = []
    for j in 1:size(m[i][2])[1]
        push!(xdata,m[i][2][j][1])
        push!(ydata,m[i][2][j][2])
    end
    ax.scatter(xdata,ydata,label="L = $(m[i][1])")
end
ax.set_ylabel(L"Magnetization per spin, $\langle |m| \rangle$ ")
ax.set_xlabel(L"Temperature, T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_m_3D.png")


fig, ax = plt.subplots()
for i in 1:2
    xdata = []
    ydata = []
    for j in 1:size(cv[i][2])[1]
        push!(xdata,cv[i][2][j][1])
        push!(ydata,cv[i][2][j][2])
    end
    ax.scatter(xdata,ydata,label="L = $(cv[i][1])")
end
ax.set_ylabel(L"Specific Heat, $C_v$ ")
ax.set_xlabel(L"Temperature, T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_cv_3D.png")


fig, ax = plt.subplots()
for i in 1:2
    xdata = []
    ydata = []
    for j in 1:size(x[i][2])[1]
        push!(xdata,x[i][2][j][1])
        push!(ydata,x[i][2][j][2])
    end
    ax.scatter(xdata,ydata,label="L = $(x[i][1])")
end
ax.set_ylabel(L"Susceptibility $\mathcal{X}$ ")
ax.set_xlabel(L"Temperature, T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_x_3D.png")=#



#=taum = []
m = []
cv = []
x = []

for L in 25:25:100
    data = jldopen("data_updated_$L.jld2")

    corrm = []
    m_L = []
    cv_L = []
    x_L = []

    for i in 1:size(data["M"])[1]
        if data["M"][i][1] <= 2.0 || data["M"][i][1] >= 2.6
            tau = minimum(findall(<=(exp(-1)),autocor(abs.(data["M"][i][2][500:5000]/L^2),1:100)))        
            samplet = maximum([100,3*tau])
            Mavg = mean(abs.(data["M"][i][2][500:samplet:5000]))
            M2avg = mean(data["M2"][i][2][500:samplet:5000])
            Eavg = mean(data["E"][i][2][500:samplet:5000])
            E2avg = mean(data["E2"][i][2][500:samplet:5000])

        elseif (data["M"][i][1]>2.0 && data["M"][i][1]<2.25) || (data["M"][i][1]>2.34 && data["M"][i][1]<=2.55)
            tau = minimum(findall(<=(exp(-1)),autocor(abs.(data["M"][i][2][5000:25000]/L^2),1:1000)))
            samplet = maximum([100,3*tau])
            Mavg = mean(abs.(data["M"][i][2][5000:samplet:25000]))
            M2avg = mean(data["M2"][i][2][5000:samplet:25000])
            Eavg = mean(data["E"][i][2][5000:samplet:25000])
            E2avg = mean(data["E2"][i][2][5000:samplet:25000])

        elseif data["M"][i][1]>=2.25 && data["M"][i][1]<=2.34
            tau = minimum(findall(<=(exp(-1)),autocor(abs.(data["M"][i][2][10000:150000]/L^2),1:4000)))
            samplet = maximum([100,3*tau])
            Mavg = mean(abs.(data["M"][i][2][10000:samplet:150000]))
            M2avg = mean(data["M2"][i][2][10000:samplet:150000])
            Eavg = mean(data["E"][i][2][10000:samplet:150000])
            E2avg = mean(data["E2"][i][2][10000:samplet:150000])
        end

        push!(corrm, [data["M"][i][1],tau])
        push!(m_L,[data["M"][i][1],Mavg/L^2])
        push!(x_L,[data["M"][i][1],(1/data["M"][i][1])*(M2avg - Mavg^2)/L^2])
        push!(cv_L,[data["M"][i][1],(1/data["M"][i][1]^2)*(E2avg - Eavg^2)/L^2])
    end

    push!(taum,[L,corrm])
    push!(m,[L,m_L])
    push!(x,[L,x_L])
    push!(cv,[L,cv_L])
end

for i in 1:4
    xdata = []
    ydata = []
    for j in 1:size(taum[i][2])[1]
        if taum[i][2][j][1] == 2.27
            println(taum[i][2][j][2])
        end
    end
end

fig, ax = plt.subplots()
for i in 1:4
    xdata = []
    ydata = []
    for j in 1:size(taum[i][2])[1]
        push!(xdata,taum[i][2][j][1])
        push!(ydata,taum[i][2][j][2])
    end
    ax.scatter(xdata,ydata,label="L = $(taum[i][1])")
end
ax.set_ylabel(L"Autocorrelation time, $\tau$ ")
ax.set_xlabel(L"Temperature, T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_autocorr_2D.png")


fig, ax = plt.subplots()
for i in 1:4
    xdata = []
    ydata = []
    for j in 1:size(m[i][2])[1]
        push!(xdata,m[i][2][j][1])
        push!(ydata,m[i][2][j][2])
    end
    ax.scatter(xdata,ydata,label="L = $(m[i][1])")
end
ax.set_ylabel(L"Magnetization per spin, $\langle |m| \rangle$ ")
ax.set_xlabel(L"Temperature, T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_m_2D.png")


fig, ax = plt.subplots()
for i in 1:4
    xdata = []
    ydata = []
    for j in 1:size(cv[i][2])[1]
        push!(xdata,cv[i][2][j][1])
        push!(ydata,cv[i][2][j][2])
    end
    ax.scatter(xdata,ydata,label="L = $(cv[i][1])")
end
ax.set_ylabel(L"Specific Heat, $\langle C_v \rangle$ ")
ax.set_xlabel(L"Temperature, T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_cv_2D.png")


fig, ax = plt.subplots()
for i in 1:4
    xdata = []
    ydata = []
    for j in 1:size(x[i][2])[1]
        push!(xdata,x[i][2][j][1])
        push!(ydata,x[i][2][j][2])
    end
    ax.scatter(xdata,ydata,label="L = $(x[i][1])")
end
ax.set_ylabel(L"Susceptibility $\langle \mathcal{X} \rangle$ ")
ax.set_xlabel(L"Temperature, T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_x_2D.png")=#



#=data = jldopen("data_corr_func_100.jld2")
fig, ax = plt.subplots()
for i in 1:size(data["M"])[1]
    if data["M"][i][1] != 2.4
        ax.plot(0:50,data["M"][i][2],label="T = $(data["M"][i][1])")
    end
end

ax.set_ylabel(L"$g(r) = \langle s(k)s(k+r) \rangle$")
ax.set_xlabel(L"Distance, r")
ax.set_xlim(0,50)
ax.set_ylim(-0.03,1.03)
ax.grid()
ax.legend(bbox_to_anchor=(0.17,1.02,1,0.2), loc="lower left",frameon=false, borderaxespad=0, ncol=3)
fig.tight_layout()
plt.savefig("Figure_corr_func.png")



temperature = [2.27,2.3,2.35]
data = jldopen("data_wolff_autocorr.jld2")
fig, ax = plt.subplots()
s = 0:200
for T in temperature
    for i in 1:size(data["M"])[1]
        if data["M"][i][1] == T 
            ax.plot(s,autocor(abs.(data["M"][i][2][100:15000]/L^2),0:200),label="T = $T")
        end
    end
end

ax.plot([0,200],[exp(-1),exp(-1)],"k--",label=L"e^{-1}")
ax.set_xlim([0,201])
ax.set_ylim([-0.1,1.0])
ax.set_ylabel("Autocorrelation for |m|")
ax.set_xlabel("Time (# of sweeps)")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_autocorr_wolff.png")=#

using GLM
using DataFrames

L = 20:256
lattice = jldopen("data_finite_size_scaling.jld2")["L"]
data = jldopen("data_finite_size_scaling.jld2")["M"]
m = []; x = []; cv = [];

for i in 1:3:size(data)[1]
    push!(m,data[i])
    push!(x,data[i+1])
    push!(cv,data[i+2])
end

dataset = DataFrame(X = log10.(lattice), Y = log10.(m))
ols_m = lm(@formula(Y ~ 1 + X), dataset)
a, b = coef(ols_m)
mhat_log = a .+ b.*log10.(L)
mhat = []
for i in 1:size(mhat_log)[1]
    push!(mhat,10^mhat_log[i])
end
fig, ax = plt.subplots()
ax.plot(L,mhat,"k--",label="Fit: slope = $(round(b,digits=3))")
ax.scatter(lattice,m,label="Wolff algorithm based Monte Carlo")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylabel(L"Magnetization per spin, $\langle |m| \rangle$ ")
ax.set_xlabel("Lattice size, L")
ax.legend()
fig.tight_layout()
plt.savefig("Figure_fss_m.png")


dataset = DataFrame(X = log10.(lattice), Y = log10.(x))
ols_x = lm(@formula(Y ~ 1 + X), dataset)
a, b = coef(ols_x)
xhat_log = a .+ b.*log10.(L)
xhat = []
for i in 1:size(xhat_log)[1]
    push!(xhat,10^xhat_log[i])
end
fig, ax = plt.subplots()
ax.plot(L,xhat,"k--",label="Fit: slope = $(round(b,digits=3))")
ax.scatter(lattice,x,label="Wolff algorithm based Monte Carlo")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylabel(L"Susceptibility $\mathcal{X}$ ")
ax.set_xlabel("Lattice size, L")
ax.legend()
fig.tight_layout()
plt.savefig("Figure_fss_x.png")