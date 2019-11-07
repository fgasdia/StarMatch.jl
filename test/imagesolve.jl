using Test
using CSV
using FITSIO
using FileIO
using JLD2
using StaticArrays

using StarMatch

const gaiaspd = load("gaia_dr2_+11.jld2", "spd")

function hms2deg(h,m,s)
    dechrs = h + m/60 + s/3600
    return dechrs/24*360
end

function dms2deg(d,m,s)
    return sign(d)*(abs(d) + m/60 + 10/3600)
end

"""
Calculate pixel scale (i.e. resolution) in degrees
"""
function resolution(camera::StarMatch.Camera)
    return 2*atand(camera.pixelsize/2/camera.focallength)
end

function buildgaiaspd()
    camera = StarMatch.Camera(1936, 1216, 5.86e-6, 620e-3)

    f = CSV.file("gaia_dr2_+11.csv")

    # TODO: Remove double stars
    catalog = [StarMatch.CatalogStar(s.source_id, s.ra, s.dec, s.phot_g_mean_mag) for s in f]

    spd = StarMatch.generatespd(camera, catalog)

    save("gaia_dr2_+11.jld2", Dict("spd"=>spd))
end

function solvenarrowimage()
    camera = StarMatch.Camera(1936, 1216, 5.86e-6, 620e-3)

    # TODO: Filtering for latitude; if s.dec > (29-90)  # Daytona, FL is 29° latitude
    # TODO: and what stars are in the sky at observation time
    # XXX: We really need to do this for a significant speedup

    # Star positions in simulated image
    imagedata = CSV.File("synthetic_25deg_+10.5.txt"; header=3)

    imagestars = StarMatch.CoordinateVector{Float64}([SVector(x, y) for (x, y) in zip(imagedata.PixelX,
        imagedata.PixelY)])

    @time matches = StarMatch.solve(camera, imagestars, gaiaspd, 2, 2)

    return imagedata, matches
end

@testset "Image solve" begin
    #==
    Simulated wide, shallow field (mag +6)

    Use Tycho-2 catalog
    ==#
    camera = StarMatch.Camera(2048, 2048, 5e-6, 100e-3)

    f = CSV.file("SKY2000_Magnitude6_doublestars_0.12.txt")
    catalog = [StarMatch.CatalogStar(s.ID, s.RAJ2000, s.DEJ2000, s.mag) for s in f]

    spd = StarMatch.generatespd(camera, catalog)  # only takes a second

    # Star positions in simulated image
    f = CSV.File("synthetic_254deg_+6.txt"; header=3)

    imagestars = StarMatch.CoordinateVector{Float64}([SVector(x, y) for (x, y) in zip(f.PixelX, f.PixelY)])

    matches = StarMatch.solve(camera, imagestars, spd, 3, 4)

    trueRAs = getproperty(f, Symbol("RA(deg)"))
    trueDECs = getproperty(f, Symbol("Dec(deg)"))
    for m in matches
        truestaridx = findfirst((f.PixelX .== m.xy[1]) .& (f.PixelY .== m.xy[2]))
        # println(truestaridx)
        @test isapprox(catalog[m.catalogidx].ra, trueRAs[truestaridx], atol=10*resolution(camera))
        @test isapprox(catalog[m.catalogidx].dec, trueDECs[truestaridx], atol=10*resolution(camera))
    end

    #==
    Simulated narrow, deep field (mag +10)

    Use Gaia DR2
    ==#
    camera = StarMatch.Camera(1936, 1216, 5.86e-6, 620e-3)

    f = CSV.file("gaia_dr2_+11.csv")
    # TODO: Remove double stars
    catalog = [StarMatch.CatalogStar(s.source_id, s.ra, s.dec, s.phot_g_mean_mag) for s in f]

    imagedata, matches = solvenarrowimage()

    trueRAs = getproperty(imagedata, Symbol("RA(deg)"))
    trueDECs = getproperty(imagedata, Symbol("Dec(deg)"))

    for m in matches
        truestaridx = findfirst((imagedata.PixelX .== m.xy[1]) .& (imagedata.PixelY .== m.xy[2]))
        @test isapprox(catalog[m.catalogidx].ra, trueRAs[truestaridx], atol=10*resolution(camera))
        @test isapprox(catalog[m.catalogidx].dec, trueDECs[truestaridx], atol=10*resolution(camera))
    end

    using Plots
    scatter(trueRAs, trueDECs, xlims=(71.2,72.4), ylims=(26.9,27.6))
    scatter!(getfield.(catalog[getfield.(matches, :catalogidx)], :ra), getfield.(catalog[getfield.(matches, :catalogidx)], :dec), markeralpha=0.4)
    scatter!(getfield.(catalog, :ra), getfield.(catalog,:dec), markersize=2)
end
