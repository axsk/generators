using LinearAlgebra, SparseArrays, Arpack, Plots, Random, Parameters

@with_kw struct AugmentedGenerator
    gen = MovingPotential() # the underlying generator
    ntime = 10     # discretization depth which should not influence eigenvalues of Q
    sc2 = 0            # 0=directed time (1=backward, 2=rows&columns, 3=doubly stochastic)
end

length(g::AugmentedGenerator) = length(g.gen)


""" assemble the Q-matrix for a nonautonomous system from their fixed-time generators """
function generator(c::AugmentedGenerator)
    @unpack sc2, gen, ntime = c

    t=1/ntime:1/ntime:1
    n = length(gen)

    #assemble big Q-matrix starting from last time step
    Q=[]
    for s=ntime:-1:1
        Qs=generator(t[s], gen)
        # remove the diagnal entries since we dont want negatives (?)
        Qs = Qs - Diagonal(Qs)

        # entropy part of Q (?)
        ds=sum(Qs, dims = 2)

        if (s==1)
            delta_t=t[s]
        else
            delta_t=t[s]-t[s-1]
        end

        # holding probabilities
        H = Diagonal(exp.(-delta_t*ds)|>vec)

        s1 = size(Q,1)

        #new row of Q
        if (s1==0)
            #first time step = corresponding rate matrix
            Q=Qs
        else
            # with directed time
            if (sc2==0)
                Q = [(I-H)*Qs       H*Q[1:n,:]
                     zeros(s1,n)    Q]
            end

            # with last row making columns sums equal to 1
            if (sc2==2)
                Q = [(I-H)*Qs       H*Q[1:n,:]
                     zeros(s1,n)    Q]
                Q[end-n+1:end,1:n]=Qs
                Q[end-n+1:end,:] = Q[end-n+1:end,:]-Q[1:n,:]
            end
        end
    end

    #negative elements are only on the diagonal
    Q=Q-Diagonal(sum(Q,dims=2)|>vec)
    Q=sparse(Q)
end

function analyse(ag::AugmentedGenerator = AugmentedGenerator(), Q=generator(ag))
    # threshold for turning fuzzy sets into hard sets
    thresh=0.3
    ntime = ag.ntime
    n = length(ag)

    spectrum, _ = eigs(Q-0.1*I,nev=6,which=:SM)
    spectrum = spectrum .+ 0.1

    # PART 1: PCCA+ based on second largest eigenvalue/vector
    #         The second largest eigenvalue/vector is used to compute
    #         the membership function of the "right" basin of the potential
    #         This function can be time-dependent (if the potential is so)
    #         The membership function is plotted (per time-interval)

    v, u = eigs(Q-0.1*I, nev=2, which=:SM)
    u = u*sign(u[n,2]) # such that chi is the "right part"
    u2 = real(u[:,2]) # Alex: the real is adhoc, correct?)
    chi=1/(maximum(u2)-minimum(u2))*(u2.-minimum(u2))

    #plot membership furnction of the right cluster per time step
    plot(title="membership of right cluster per time")
    for i=1:ntime
        plot!(chi[(i-1)*n+1:i*n])
    end
    plot!() |> display

    # PART 2: Left largest eigenvector of Q and "invariant density"
    #         The "invariant density" is decomposed into parts
    #         which belong to the different time-steps
    #         If sc2==0 then the invariant density is only non-zero for the
    #         last time-step
    #         The density is plotted (per time-interval, blue) and summed up (red)

    v, u = eigs(Q'-0.1*I, nev=2,which=:SM)
    density = u[:,1] / sum(u[:,1])

    #plot "invariant density" density per time step
    plot(title="invariant density per time")
    sumdens=zeros(n,1)
    for i=1:ntime
        dd = density[(i-1)*n+1:i*n] |> real
        plot!(dd)
        sumdens=sumdens+dd
    end
    plot!(sumdens,color="red") |> display

    # PART 3: Computation of the exit rate of the "membership fucntion"
    #         according to new ideas in archive paper
    #         Computation of mean first exit times according to "old school"
    #         with problems of interpretation (moving time)
    #         Plotted: The comparison of the two concepts
    #         only makes sense for sc2==0 or sc==2

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
    plot!(chi[n+1:end]/rate,mfpt[n+1:end]) #,"y.")
    plot!(chi[1:n]/rate,mfpt[1:n]) |> display #,"ko")
end
