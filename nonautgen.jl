# PART 0: Q Assembly


# PART 1: Right eigenvectors of Q and PCCA+
#     ?    The second largest eigenvalue/vector is used to compute
#         the membership function of the "right" basin of the potential
#         This function can be time-dependent (if the potential is so)
#         The membership function is plotted (per time-interval)

# PART 2: Left largest eigenvector of Q and "invariant density"
#         The "invariant density" is decomposed into parts
#         which belong to the different time-steps
#         If sc2==0 then the invariant density is only non-zero for the
#         last time-step
#         The density is plotted (per time-interval, blue) and summed up (red)

# PART 3: Computation of the exit rate of the "membership fucntion"
#         according to new ideas in archive paper
#         Computation of mean first exit times according to "old school"
#         with problems of interpretation (moving time)
#         Plotted: The comparison of the two concepts
#         only makes sense for sc2==0 or sc==2


## SETTINGS
include("movingpotential.jl")
include("bickleyjet.jl")
using LinearAlgebra, SparseArrays, Arpack, Plots, Random

function nonautgen(gen)
    #scenario of Q-assembly
    sc2 = 0 # 0=directed time (1=backward, 2=rows&columns, 3=doubly stochastic)


    if (sc2 == 3)
        # for the doubly stochastic approach
        nr = 2 # 0=sytematically >0=number of randperm-matrices to be added for doubly stochastic
    end

    #discretization depth which should not influence eigenvalues of Q
    timesteps=10  #number of intervals in [0,1]-time (e.g. 10)
    t=1/timesteps:1/timesteps:1

    # threshold for turning fuzzy sets into hard sets
    thresh=0.3

    ##

    eye(n) = Diagonal(ones(n))

    #sc2 ={0,1,2} assemble big Q-matrix starting from last time step
    if (sc2<3)
        Q=[]
        for s=timesteps:-1:1
            Qs=generator(t[s], gen)
            grids=size(Qs,1)
            ds=sum(Qs, dims = 2)       #and entropy part of Q

            if (s==1)
                delta_t=t[s]
            else
                delta_t=t[s]-t[s-1]
            end
            H=Diagonal(exp.(-delta_t*ds)|>vec)  #holding probabilities

            #new row of Q
            s1=size(Q,1)
            if (s1==0)
                Q=Qs #first time step = corresponding rate matrix
            else
                #EITHER: with directed time
                if (sc2==0)
                    Q=[(I-H)*Qs  H*Q[1:grids,:]
                    zeros(s1,grids)    Q]
                end

                #OR: with backward steps based columnwise on same rate matrix
                if (sc2==1)
                    block=[]
                    fact=1/2
                    for j=s:timesteps-1
                        block=[block; fact*Qs]
                        fact=fact*1/2
                    end
                    Q=[(eye(grids)-H)*Qs       H*Q[1:grids,:]
                    block  1/2*Q[:,1:grids]  Q[:,grids+1:end]]
                end

                #OR: with last row making columns sums equal to 1
                if (sc2==2)
                    Q=[(eye(grids)-H)*Qs  H*Q[1:grids,:]
                    zeros(s1,grids)    Q]
                    Q[end-grids+1:end,1:grids]=Qs
                    Q[end-grids+1:end,:] = Q[end-grids+1:end,:]-Q[1:grids,:]
                end
            end
        end
    end


    # sc2==3: Assemble big Q-matrix according to doubly stochastic matrix
    if (sc2==3)
        #generate doubly stochastic matrix R
        R=zeros(timesteps, timesteps)
        weights=rand(1,nr)
        weights=weights/sum(weights)
        for i=1:nr
            E=eye(timesteps)
            R=R+weights[i]*E[1:timesteps,randperm(timesteps)]
        end
        if (nr==0)
            R=1/timesteps*ones(timesteps,timesteps)
        end

        #Fill columns of Q with weighted Qs
        Qs=infinit_generator(t[1])
        grids=size(Qs,1)
        Q=zeros(timesteps*grids, timesteps*grids)
        for i=1:timesteps
            Qs=infinit_generator(t[i])
            grids=size(Qs,1)
            ds=sum(Qs,dims=2)       #and entropy part of Q

            for j=1:timesteps
                Q[(j-1)*grids+1:j*grids,(i-1)*grids+1:i*grids]=R[j,i]*Qs
            end
        end
    end

    Q=Q-Diagonal(sum(Q,dims=2)|>vec) #negative elements are only on the diagonal
    Q=sparse(Q)

    spectrum, _ = eigs(Q-0.1*I,nev=6,which=:SM)
    spectrum = spectrum .+ 0.1

    #PCCA+ based on second largest eigenvalue/vector
    v, u = eigs(Q-0.1*I, nev=2, which=:SM)
    u = u*sign(u[grids,2]) # such that chi is the "right part"
    u2 = real(u[:,2]) # Alex: the real is adhoc, correct?)
    chi=1/(maximum(u2)-minimum(u2))*(u2.-minimum(u2))

    #plot membership furnction of the right cluster per time step
    plot(title="membership of right cluster per time")
    for i=1:timesteps
        plot!(chi[(i-1)*grids+1:i*grids])
    end
    plot!() |> display

    # "invariant density"
    v, u = eigs(Q'-0.1*I, nev=2,which=:SM)
    density = u[:,1] / sum(u[:,1])

    #plot "invariant density" density per time step
    plot(title="invariant density per time")
    sumdens=zeros(grids,1)
    for i=1:timesteps
        dd = density[(i-1)*grids+1:i*grids] |> real
        plot!(dd)
        sumdens=sumdens+dd
    end
    plot!(sumdens,color="red") |> display

    # Exit rate estimation according to Natalia Ernst

    u2 = real(u[:,2])
    rate = -spectrum[2]*maximum(u2)/(maximum(u2)-minimum(u2)) |> real
    epsilon2 = spectrum[2]*minimum(u2)/(maximum(u2)-minimum(u2))

    # Mean first passage time computation
    indicator = findall(chi.>thresh)
    Qred = Q[indicator, indicator]
    mfpts = Qred\(-1*ones(length(indicator),1))
    mfpt = zeros(size(Q,1),1)
    mfpt[indicator,1]=mfpts

    #plot mean first exit times from both concepts (comparison)
    plot(title="mean first exit time comparison")
    plot!(chi[grids+1:end]/rate,mfpt[grids+1:end]) #,"y.")
    plot!(chi[1:grids]/rate,mfpt[1:grids]) |> display #,"ko")

end
