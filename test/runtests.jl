using Test

println("=== Starting tests ===")
println("NOTE: SPD generation failures may cascade into image solve.")

@time include("spr_riav_comparison.jl")
# @time include("imagesolve.jl")
