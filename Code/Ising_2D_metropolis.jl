module Ising_2D_metropolis

    using Random
    using StatsBase

    struct run_parameters
        ## Paramters for running Monte Carlo
        N_sweeps::Int64     #Total MCS
        N_skip::Int64       #Sampling time
        N_thermalize::Int64 #Steps to discard initially
    end

    mutable struct lattice
        ## Lattice parameters
        L::Int64            #Lattice size
        beta::Float64       #1/(kB*T)
        config              #storing -1/+1 for lattice sites
    end


    function initialize(init::String,L::Int)
        #=
            Initial spin configuration

            Parameter:
            lattice - struct (ising_pars)

            Return: 
            config - initial lattice configuration
        =#

        if init == "uniform"
            config = ones(L,L)
        elseif init == "random"
            config = 2*bitrand((L,L)) - ones(L,L)
        end

        return config
    end


    function periodic_boundary(index::Int64,L::Int64)
        #=
            Check for boundary index
            
            Parameter:
            index - Int (index to be checked for boundary condition)
            
            Return:
            integer
        =#
        
        if index == 0
            return L 
        elseif index == L + 1
            return 1
        else
            return index
        end
    end


    function del_energy(lat::lattice,x::Int64,y::Int64)
        #=
            Calculate the change in energy after spin flip
            
            Parameters:
            lat - struct lattice 
            x, y - Int (indices of the spin flip)
            
            Return:
            e - Float  (change in energy)
        =#
            
        e = begin 2*lat.config[x,y]*(lat.config[periodic_boundary(x-1,lat.L),y] + 
                                lat.config[periodic_boundary(x+1,lat.L),y] +
                                lat.config[x,periodic_boundary(y-1,lat.L)] +
                                lat.config[x,periodic_boundary(y+1,lat.L)]) end
        
        return e
    end


    function magnetization(lat::lattice)
        #= 
            To calculate magnetizaton (order parameter)

            Parameters:
            lat - struct lattice
        =#
        M = sum(sum(lat.config))
        return M, M^2 
    end


    function energy(lat::lattice)
        #=
            To calculate energy  

            Parameters:
            lat - struct lattice
        =#
        E = 0
        for i in 1:lat.L
            for j in 1:lat.L
                temp = -0.5*del_energy(lat,i,j)
                E += temp
            end
        end
        E /=2
        return E, E^2
    end
        

    function run_MC(lat::lattice,p::run_parameters,corr_check::Bool)
        #=
            To run Monte Carlo

            Parameters:
            lat - struct lattice    
            p - parameters for running MC
            corr_check - whether to calculate correlation function or not (true or false)
        =#

        m = []      #to store magnetization
        m2 = []     #to store (magnetization)^2
        e = []      #to store energy
        e2 = []     #to store (energy)^2

        if corr_check
            corr = zeros(Int64(floor(lat.L/2))+1)
        end

        ## Start Monte Carlo
        for i in 1:p.N_sweeps
            for _ in 1:lat.L^2
                ## Choose some random spin to flip
                x, y = rand(1:lat.L,2)

                ## Change in energy due to flip
                delE = del_energy(lat,x,y)

                ## Perform Metropolis check and flip
                if delE <= 0
                    lat.config[x,y] = - lat.config[x,y]
                else
                    check = rand(Float64)
                    if check < exp(-lat.beta*delE)
                        lat.config[x,y] = - lat.config[x,y]
                    end
                end
            end

            ## Check if thermalization time is over
            ## if yes, check for sampling time 
            if i >= p.N_thermalize && i%p.N_skip == 0
                temp, temp2 = magnetization(lat)
                push!(m,temp)
                push!(m2,temp2)
                temp, temp2 = energy(lat)
                push!(e,temp)
                push!(e2,temp2)
                if corr_check
                    corr += correlation(lat)
                end
            end

            ## Print to keep track of Monte Carlo steps
            if i%1000 == 0
                println("MC step: $i")
            end
        end

        if corr_check
            ## Returning magnetization per spin, susceptibility, specific heat and correlation function
            return  mean(abs.(m))/lat.L^2, lat.beta*(mean(m2) - mean(abs.(m))^2)/L^2, lat.beta^2*(mean(e2) - mean(e)^2)/L^2, corr/size(m)[1]
        else
            ## Returning magnetization per spin, susceptibility and specific heat
            return  mean(abs.(m))/lat.L^2, lat.beta*(mean(m2) - mean(abs.(m))^2)/L^2, lat.beta^2*(mean(e2) - mean(e)^2)/L^2
        end
    end


    function correlation(lat::lattice)
        #=
            to calculate correlation function for a Ising 2D square lattice
            
                Parameters:
                lat: Struct lattice

                Returns an array with ith element corresponding to correlation function value at g(i-1)
        =#

        corr = zeros(Int64(floor(lat.L/2))+1)
        sample = zeros(Int64(floor(lat.L/2))+1)

        for i in 1:lat.L
            for j in 1:lat.L
                l = min(Int64(floor(lat.L/2)),lat.L-j)
                for k in j:j+l
                    corr[k-j+1] += lat.config[i,j]*lat.config[i,k]
                    sample[k-j+1] += 1
                end
            end
        end

        return corr./sample

    end


    function main(N_sweeps::Int64,N_skip::Int64,N_thermalize::Int64,L::Int64,init::String,T::Float64,corr_check::Bool)
        #=
            main function to perform Metropolis Monte Carlo
            
            Parameters:
            N_sweeps: Total number of MCS
            N_skip: Sampling time
            N_thermalize: Number of MCS to skip initially
            L: Lattice size
            init: Initial configuration to use ("uniform" or "random")
            T: Temperature
            corr_check: whether to perform correlation function or not (true or false)
        =#

        p = run_parameters(N_sweeps,N_skip,N_thermalize)
        lat =  lattice(L,1/T,initialize(init,L))
        temp = run_MC(lat,p,corr_check)

        ## Return magnetization per spin, susceptibility, specific heat (and correlation function)
        return temp
    end

end
