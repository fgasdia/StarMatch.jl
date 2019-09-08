# Test hashing algorithm

using BenchmarkTools

# Dictionary for approach 1
dict1 = Dict((1:10^6) .=> zip(rand(1:10, 1_000_000),
                           rand(1:1_000_000, 1_000_000),
                           rand(1:1_000_000, 1_000_000),
                           rand(1:1_000_000, 1_000_000),
                           rand(1:1_000_000, 1_000_000)))

@benchmark $dict1[1000]
@benchmark $dict1[100]
@benchmark $dict1[14203]
# This is consistently about 25 ns per lookup
25e-9*1_000_000 # -> about 25 ms to lookup entire dictionary


# Customish solution
# https://discourse.julialang.org/t/compiling-to-branch-table/16599
const CODES = 1:10^6

function findcode(x)
    index = findfirst(isequal(x), CODES)
    index = index == nothing ? 0 : index
end
# This does _not_ scale well.

# Here's a lookup table version (from same url)
lookuptable(v) = (a = zeros(Int, maximum(v)); for (i,x) in enumerate(v); a[x]=i; end; a)
const LUT = lookuptable(CODES)

function findcode_lookuptable(x)
    x < 1 || x > length(LUT) ? 0 : LUT[x]
end

@benchmark findcode_lookuptable(1323)

singledict = Dict(1:10^6 .=> 1:10^6)
@benchmark $singledict[1323]
# The lookuptable approache averages about 10 ns and is faster than a regular dictionary

# With 5 indices...
function lookuptable1(v)
    a = zeros(NTuple{5,Int}, maximum(v))
    for x in eachindex(v)
        a[x] = (rand(1:10), rand(1:1_000_000), rand(1:1_000_000), rand(1:1_000_000), rand(1:1_000_000))
    end
    return a
end
const LUT1 = lookuptable(CODES)

function findcode_lookuptable1(x)
    x < 1 || x > length(LUT1) ? 0 : LUT1[x]
end

@benchmark findcode_lookuptable1(10232)
# This is still faster than a dictionary

# What if the dict is constant?
const dict1c = Dict((1:10^6) .=> zip(rand(1:10, 1_000_000),
                                    rand(1:1_000_000, 1_000_000),
                                    rand(1:1_000_000, 1_000_000),
                                    rand(1:1_000_000, 1_000_000),
                                    rand(1:1_000_000, 1_000_000)))

@benchmark $dict1c[1000]  # why is this faster without $?
# without, this is about 15ns



# TODO: How does size of tuple scale wrt time?
