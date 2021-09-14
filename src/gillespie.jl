
using Distributions

# simulate a single trajectory given rates and and their switching times
function gillespie(qs, qts, x0, tmax)
    # begin if first frame
    i = 1
    t = qts[i]
    q = qs[i]
    x = x0
    xs = [x]
    ts = [t]
    while t < tmax
        qout = -q[x,x]
        dt = rand(Distributions.Exponential(1/qout))
        if dt + t > qts[i+1]
            i += 1
            q = qs[i]
            t = qts[i]
        else
            x = drawnext(q, x)
            t = t + dt
            push!(xs, x)
            push!(ts, t)
        end
    end
    return ts, xs
end

function drawnext(q, x)
    qq = q[x, :]
    qq[x] = 0
    qq ./= sum(qq)
    rand(Categorical(qq))
end