using Test
using CSV
using FITSIO
using FileIO
using JLD2
using StaticArrays

using StarMatch

println("\nNOTE: Loading `gaia_dr2_+11.jld2` may take a minute.\n")
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

"""
    buildgaiaspd()

Build and save the star pattern database `gaia_dr2_+11.jld2`.

!!! note

    This function can take a long time (~ 1 hour) to run.
    A progress meter will display progress.
"""
function buildgaiaspd()
    camera = StarMatch.Camera(1936, 1216, 5.86e-6, 620e-3)

    f = CSV.file("gaia_dr2_+11.csv")

    # TODO: Remove double stars, see Distances.jl
    catalog = [StarMatch.CatalogStar(s.source_id, s.ra, s.dec, s.phot_g_mean_mag) for s in f]

    spd = StarMatch.generatespd(camera, catalog)

    save("gaia_dr2_+11.jld2", Dict("spd"=>spd))
end

function solvewideimage()
    println("\nVisual magnitude, wide field image")
    camera = StarMatch.Camera(2048, 2048, 5e-6, 100e-3)

    f = CSV.file("SKY2000_Magnitude6_doublestars_0.12.txt")
    catalog = [StarMatch.CatalogStar(s.ID, s.RAJ2000, s.DEJ2000, s.mag) for s in f]

    # Generate the spd (usually this is a prebuilt file)
    spd = StarMatch.generatespd(camera, catalog)  # only takes a second

    # Star positions in simulated image
    imagedata = CSV.File("synthetic_254deg_+6.txt"; header=3)

    imagestars = StarMatch.CoordinateVector{Float64}([SVector(x, y) for (x, y) in
        zip(imagedata.PixelX, imagedata.PixelY)])

    @time matches = StarMatch.solve(camera, imagestars, spd, 2, 3)

    return camera, catalog, imagedata, matches
end

function solvenarrowimage()
    println("\nDeep, narrow field image")
    camera = StarMatch.Camera(1936, 1216, 5.86e-6, 620e-3)

    # TODO: Filtering for latitude; if s.dec > (29-90)  # Daytona, FL is 29° latitude
    # TODO: and what stars are in the sky at observation time
    # TODO: Remove double stars

    f = CSV.file("gaia_dr2_+11.csv")
    catalog = [StarMatch.CatalogStar(s.source_id, s.ra, s.dec, s.phot_g_mean_mag) for s in f]

    # Star positions in simulated image
    imagedata = CSV.File("synthetic_25deg_+10.5.txt"; header=3)

    imagestars = StarMatch.CoordinateVector{Float64}([SVector(x, y) for (x, y) in
        zip(imagedata.PixelX, imagedata.PixelY)])

    @time matches = StarMatch.solve(camera, imagestars, gaiaspd, 2, 3)

    return camera, catalog, imagedata, matches
end

@testset "Image solve" begin
    #==
    Simulated wide, shallow field (mag limit +6)

    Use Tycho-2 catalog
    ==#
    camera, catalog, imagedata, matches = solvewideimage()

    trueRAs = getproperty(imagedata, Symbol("RA(deg)"))
    trueDECs = getproperty(imagedata, Symbol("Dec(deg)"))

    for (i, m) in enumerate(matches)
        truestaridx = findfirst((imagedata.PixelX .== m.xy[1]) .& (imagedata.PixelY .== m.xy[2]))
        @test isapprox(catalog[m.catalogidx].ra, trueRAs[truestaridx], atol=3*resolution(camera))
        @test isapprox(catalog[m.catalogidx].dec, trueDECs[truestaridx], atol=3*resolution(camera))
    end

    #==
    Simulated narrow, deep field (mag limit +10.5)

    Use Gaia DR2
    ==#
    camera, catalog, imagedata, matches = solvenarrowimage()

    trueRAs = getproperty(imagedata, Symbol("RA(deg)"))
    trueDECs = getproperty(imagedata, Symbol("Dec(deg)"))

    for m in matches
        truestaridx = findfirst((imagedata.PixelX .== m.xy[1]) .& (imagedata.PixelY .== m.xy[2]))
        @test isapprox(catalog[m.catalogidx].ra, trueRAs[truestaridx], atol=10*resolution(camera))
        @test isapprox(catalog[m.catalogidx].dec, trueDECs[truestaridx], atol=10*resolution(camera))
    end
end
