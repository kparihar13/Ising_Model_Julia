using PyPlot
using StatsBase
using JLD2
using LaTeXStrings

##########################Thermalization#######################
## Loading data
#L = 100
#data = jldopen("data_$L.jld2")
#=data_random = jldopen("data_random_$L.jld2")
data_m = []
data_e = []

for i in 1:3
    for j in 1:size(data["M"])[1]
        if data["M"][j][1] == data_random["M"][i][1]
            push!(data_m,[data["M"][j][1],data["M"][j][2],data_random["M"][i][2]])
            push!(data_e,[data["E"][j][1],data["E"][j][2],data_random["E"][i][2]])
        end
    end
end

## Ploting
for i in 1:3
    fig, ax = plt.subplots(2)
    s = 1:size(data_m[i][2])[1]
    ax[1].plot(s, data_m[i][2]/L^2,label="Uniform initialization")
    ax[1].set_ylim([-1.05,1.05])
    ax[1].legend(frameon=false)
    ax[1].grid()
    ax[2].plot(s, data_m[i][3]/L^2,label="Random initialization")
    ax[2].set_ylim([-1.05,1.05])
    ax[2].legend(frameon=false)
    ax[2].grid()
    fig.text(0.03, 0.5, "Magnetization per spin, m", va="center", rotation="vertical")
    fig.text(0.5, 0.03, "Number of sweeps", ha="center")
    plt.savefig("Figure_m$i.png")
end

for i in 1:3
    fig, ax = plt.subplots(2)
    s = 1:size(data_e[i][2])[1]
    ax[1].plot(s, data_e[i][2]/L^2,label="Uniform initialization")
    ax[1].set_ylim([-2.0,0.0])
    ax[1].legend(frameon=false)
    ax[1].grid()
    ax[2].plot(s, data_e[i][3]/L^2,label="Random initialization")
    ax[2].set_ylim([-2.0,0.0])
    ax[2].legend(frameon=false)
    ax[2].grid()
    fig.text(0.03, 0.5, "Energy per spin", va="center", rotation="vertical")
    fig.text(0.5, 0.03, "Number of sweeps", ha="center")
    plt.savefig("Figure_e$i.png")
end=#



##########################Autocorrelation#####################
Temp = [2.0,2.1,2.2,2.3,2.4,2.6,3.0]

fig, ax = plt.subplots()
s = 0:700
for T in Temp
    for i in 1:size(data["M"])[1]
        if data["M"][i][1] == T 
            if (T<=2.0 || T>2.5)
                ax.plot(s,autocor(abs.(data["M"][i][2][500:5000]/L^2),0:700),label="T = $T")
            elseif T>2 && T<2.5
                ax.plot(s,autocor(abs.(data["M"][i][2][5000:25000]/L^2),0:700),label="T = $T")
            end
        end
    end
end

ax.plot([0,700],[exp(-1),exp(-1)],"k--",label=L"e^{-1}")
ax.set_xlim([0,701])
ax.set_ylim([-0.2,1.0])
ax.set_ylabel("Autocorrelation for |m|")
ax.set_xlabel("Time (# of sweeps)")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_autocorr_100.png")

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
ax.set_xlabel(L"T")
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
ax.set_xlabel(L"T")
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
ax.set_xlabel(L"T")
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
ax.set_xlabel(L"T")
ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig("Figure_x_2D.png")=#




