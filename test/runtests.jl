using Test

println("\n=== Starting tests ===")
println("NOTE: SPD generation failures may cascade into image solve.\n")

@time include("spr_riav_comparison.jl")
@time include("imagesolve.jl")
