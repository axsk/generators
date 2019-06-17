
using Parameters
@with_kw struct MovingPotential
    sc1 = 2
    grids = 30
    beta = 10
end

length(m::MovingPotential) = m.grids

""" generate the infinitesimal generator for time t

PART 1: Defining the scenario
    sc1 defines the potential that is used for the example {0,1,2}

PART 2: Defining the physical constants
    The inverse temperature beta defines the Boltzmann density
    The computations can be based on finer or coarser
    discretizations which should lead to an adjusted flux (LUCA!)
    The potential energy function is plotted

PART 3: The Q-matrix assembly
    If the potential is chosen to be time-independent
    then the results should be the same as for
    "ordinary" SQRA => the Q-assembly is done like here

Arguments:
sc1 = 2: scenario of potential
    (0=time-indep. potential, 1= vanishing potential, 2 = moving potential)

grids=30: number of intervals in [0,1]-space (e.g. 30)

beta=10: inverse temperature for Boltzmann (e.g. 10)
    (physical constant which has an effect on eigenvalues of Q)
"""
function generator(t, c::MovingPotential)
    @unpack sc1, grids, beta = c

    # grid & flux computation according to Luca Donati
    grid=1/grids:1/grids:1;
    flux=grids^2/beta*1/9;

    #evaluate potential for discretization points. The potential is...

    #EITHER: a Gaussian that vanishes in time
    if (sc1==1)
        potential=t*exp.(-50*(grid.-0.6).^2);
    end

    #OR: a Gaussian that moves in time
    if (sc1==2)
        potential=exp.(-50*(grid .- t*0.4 .- (1-t)*0.6).^2);
    end

    #OR: a time-independent potential
    if (sc1==0)
        potential=exp.(-50*(grid.-0.6).^2);
    end

    # adjacency matrix of the intervals for Q-assembly
    # A=[[zeros(grids-1,1),eye(grids-1)];zeros(1,grids)];
    A = diagm(1=>ones(grids-1))
    A = A+A'

    # Q assembly for time-step
    sr=sqrt.(exp.(-beta*potential)); #Boltzmann distribution
    sr=sr/sum(sr);

    D  = diagm(0=>sr);
    D1 = diagm(0=>1 ./sr);         #SQRA: potential dependent part
    Q=flux*D1*A*D;
end
