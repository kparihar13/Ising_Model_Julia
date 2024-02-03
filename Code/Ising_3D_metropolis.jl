module Ising_3D_metropolis

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
        L::Int64            #cibe size
        beta::Float64       #1/(kB*T)
        config              #storing -1/+1 for cube sites
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
            config = ones(L,L,L)
        elseif init == "random"
            config = 2*bitrand((L,L,L)) - ones(L,L,L)
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


    function del_energy(lat::lattice,x::Int64,y::Int64,z::Int64)
        #=
            Calculate the change in energy after spin flip
            
            Parameters:
            lat - struct lattice 
            x, y - Int (indices of the spin flip)
            
            Return:
            e - Float  (change in energy)
        =#
            
        e = begin 2*lat.config[x,y,z]*
            (lat.config[periodic_boundary(x-1,lat.L),y,z] + 
             lat.config[periodic_boundary(x+1,lat.L),y,z] +
             lat.config[x,periodic_boundary(y-1,lat.L),z] +
             lat.config[x,periodic_boundary(y+1,lat.L),z] +
             lat.config[x,y,periodic_boundary(z-1,lat.L)] +
             lat.config[x,y,periodic_boundary(z+1,lat.L)]) end
        
        return e
    end


    function magnetization(lat::lattice)
        #= 
            To calculate magnetizaton (order parameter)
        =#
        M = sum(sum(sum(lat.config)))
        return M, M^2 
    end


    function energy(lat::lattice)
        #=
            To calculate energy  
        =#
        E = 0
        for i in 1:lat.L
            for j in 1:lat.L
                for k in 1:lat.L
                    temp = -0.5*del_energy(lat,i,j,k)
                    E += temp
                end
            end
        end
        E /=2
        return E, E^2
    end
        

    function run_MC(lat::lattice,p::run_parameters)
        #=
            To run Monte Carlo

            Parameters:
            lat - struct lattice    
            p - parameters for running MC
        =#

        m = []      #to store magnetization
        m2 = []     #to store (magnetization)^2
        e = []      #to store energy
        e2 = []     #to store (energy)^2
        
        ## Start Monte Carlo
        for i in 1:p.N_sweeps
            for _ in 1:lat.L^3
                ## Choose some random spin to flip
                x, y, z = rand(1:lat.L,3)

                ## Change in energy due to flip
                delE = del_energy(lat,x,y,z)

                ## Perform Metropolis check and flip
                if delE <= 0
                    lat.config[x,y,z] = - lat.config[x,y,z]
                else
                    check = rand(Float64)
                    if check < exp(-lat.beta*delE)
                        lat.config[x,y,z] = - lat.config[x,y,z]
                    end
                end
            end

            ## Check if thermalization time is over
            ## if yes, check for sampling time
            if i > p.N_thermalize && i%p.N_skip == 0
                temp, temp2 = magnetization(lat)
                push!(m,temp)
                push!(m2,temp2)
                temp, temp2 = energy(lat)
                push!(e,temp)
                push!(e2,temp2)
            end

            ## Print to keep track of Monte Carlo steps
            if i%1000 == 0
                println("MC step: $i")
            end
        end

        ## Returning magnetization per spin, susceptibility 
        ## and specific heat
        return  mean(abs.(m))/lat.L^3, lat.beta*(mean(m2) - mean(abs.(m))^2)/L^3, lat.beta^2*(mean(e2) - mean(e)^2)/L^3
    end


    function main(N_sweeps::Int64,N_skip::Int64,N_thermalize::Int64,L::Int64,init::String,T::Float64)
        #=
            main function to perform Metropolis Monte Carlo
            
            Parameters:
            N_sweeps: Total number of MCS
            N_skip: Sampling time
            N_thermalize: Number of MCS to skip initially
            L: Lattice size
            init: Initial configuration to use ("uniform" or "random")
            T: Temperature
        =#

        p = run_parameters(N_sweeps,N_skip,N_thermalize)
        lat =  lattice(L,1/T,initialize(init,L))
        temp = run_MC(lat,p)

        ## Return magnetization per spin, susceptibility, specific heat (and correlation function)
        return temp
    end
end
