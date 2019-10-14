#==
Random tests to compare code speeds.
==#

#==
Profiling of `solve()` shows that most time is consumed by occurrences of
`norm(a-b) < thresh`, e.g. in `vote()`. However, there is very little gain by
writing this in a different way.
==#

# Comparing these three, they get faster at each step, but only by a fraction of a millisecond
# when performed on 100_000_000 `a`s and `b`s
jnk = [@SVector rand(Float64, 2) for i = 1:100_000_000]
jnk2 = [@SVector rand(Float64, 2) for i = 1:100_000_000]

const thresh = 4
function normdiff(a, b)
    return norm(a-b) < thresh
end

function longdiff(a, b)
    return abs(a[1]-b[1]) < thresh && abs(a[2]-b[2]) < thresh
end

function longerdiff(a,b)
    return ((a[1]-thresh) < b[1] < (a[1]+thresh)) && ((a[2]-thresh) < b[2] < (a[2]+thresh))
end

@btime begin
    for i in 1:100_000_000
        normdiff($jnk[i],$jnk2[i])
    end
end
