using Test
using CSV
using JLD2
using StaticArrays

using StarMatch

img_width = 1024
img_height = 1024
pixelsize = 13.3e-6
fov = 15.0
fl = img_height*pixelsize/2/tand(fov/2)
camera = StarMatch.Camera(img_width, img_height, pixelsize, fl, fov)

f = CSV.file("SKY2000_Magnitude6_doublestars_0.12.txt")

catalog = Array{StarMatch.CatalogStar}(undef, f.rows)
for (i, s) in enumerate(f)
    catalog[i] = StarMatch.CatalogStar(s.ID, s.RAJ2000, s.DEJ2000, s.mag)
end

function readspdfile(f::IOStream)
    vectorpattern = []
    for line in eachline(f)
        numstrs = split(line, "\t")

        nums = Int[]
        for str in numstrs
            if str != ""
                push!(nums, parse(Int, str))
            end
        end

        push!(vectorpattern, SVector.(nums[1:2:end], nums[2:2:end]))
    end

    return vectorpattern
end

# spd = StarMatch.generatespd(camera, catalog)

# @save "SKY2000_Magnitude6_doublestars_0.12_spd.jd2" spd
# @load "SKY2000_Magnitude6_doublestars_0.12_spd.jd2" spd

@testset "SPD creation" begin
    spd = StarMatch.generatespd(camera, catalog)

    d1 = open(readspdfile, "SPD_vect_patt_Mv_6_dist_1.txt")
    d2 = open(readspdfile, "SPD_vect_patt_Mv_6_dist_2.txt")
    d3 = open(readspdfile, "SPD_vect_patt_Mv_6_dist_3.txt")
    d4 = open(readspdfile, "SPD_vect_patt_Mv_6_dist_4.txt")

    # Choose a couple random stars as center star (`starO`) and then compare
    # paths for each `SP` star

    # starid = 14
    spd_14 = findall(x->x.Oid == 14, spd)

    d1_diff = d1[14] - spd[spd_14[1]].path
    d1_diffr = [d[1] for d in d1_diff]
    d1_diffc = [d[2] for d in d1_diff]
    @test all(abs.(d1_diffc) .< 1)
    @test all(abs.(d1_diffr) .< 1)

    d2_diff = d2[14] - spd[spd_14[2]].path
    d2_diffr = [d[1] for d in d2_diff]
    d2_diffc = [d[2] for d in d2_diff]
    @test all(abs.(d2_diffc) .< 1)
    @test all(abs.(d2_diffr) .< 1)

    d3_diff = d3[14] - spd[spd_14[3]].path
    d3_diffr = [d[1] for d in d3_diff]
    d3_diffc = [d[2] for d in d3_diff]
    @test all(abs.(d3_diffc) .< 1)
    @test all(abs.(d3_diffr) .< 1)

    d4_diff = d4[14] - spd[spd_14[4]].path
    d4_diffr = [d[1] for d in d4_diff]
    d4_diffc = [d[2] for d in d4_diff]
    @test all(abs.(d4_diffc) .< 1)
    @test all(abs.(d4_diffr) .< 1)

    # starid = 3864
    spd_3864 = findall(x->x.Oid == 3864, spd)

    d1_diff = d1[3864] - spd[spd_3864[1]].path
    d1_diffr = [d[1] for d in d1_diff]
    d1_diffc = [d[2] for d in d1_diff]
    @test all(abs.(d1_diffc) .< 1)
    @test all(abs.(d1_diffr) .< 1)

    d2_diff = d2[3864] - spd[spd_3864[2]].path
    d2_diffr = [d[1] for d in d2_diff]
    d2_diffc = [d[2] for d in d2_diff]
    @test all(abs.(d2_diffc) .< 1)
    @test all(abs.(d2_diffr) .< 1)

    d3_diff = d3[3864] - spd[spd_3864[3]].path
    d3_diffr = [d[1] for d in d3_diff]
    d3_diffc = [d[2] for d in d3_diff]
    @test all(abs.(d3_diffc) .< 1)
    @test all(abs.(d3_diffr) .< 1)

    d4_diff = d4[3864] - spd[spd_3864[4]].path
    d4_diffr = [d[1] for d in d4_diff]
    d4_diffc = [d[2] for d in d4_diff]
    @test all(abs.(d4_diffc) .< 1)
    @test all(abs.(d4_diffr) .< 1)
end
