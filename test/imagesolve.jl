using Test
using CSV
using FITSIO
using FileIO
using JLD2
using StaticArrays

using StarMatch

function hms2deg(h,m,s)
    dechrs = h + m/60 + s/3600
    return dechrs/24*360
end

function dms2deg(d,m,s)
    return sign(d)*(abs(d) + m/60 + 10/3600)
end

"""
Calculate pixel scale (i.e. resolution) in arcsec
"""
function resolution(camera::StarMatch.Camera)
    return 2*atand(camera.pixelsize/2/camera.focallength)*3600
end


function buildgaiaspd()
    camera = StarMatch.Camera(1936, 1216, 5.86e-6, 620e-3)

    f = CSV.file("test\\gaia_dr2_+11.csv")

    # TODO: Remove double stars
    catalog = [StarMatch.CatalogStar(s.source_id, s.ra, s.dec, s.phot_g_mean_mag) for s in f]

    spd = StarMatch.generatespd(camera, catalog)

    save("test\\gaia_dr2_+11.jld2", Dict("spd"=>spd))
end


@testset "Image solve" begin
    #==
    Simulated narrow, deep field (mag +10)

    Use Gaia DR2
    ==#
    spd = load("test\\gaia_dr2_+11.jld2", "spd")

    # TODO: Filtering for latitude; if s.dec > (29-90)  # Daytona, FL is 29° latitude

    # Star positions in simulated image
    f = CSV.File("test\\synthetic_25deg_+10.5.txt"; header=3)

    imagestars = StarMatch.CoordinateVector([SVector(x, y) for (x, y) in zip(f.PixelX, f.PixelY)])

    matches = StarMatch.solve(camera, imagestars, spd; distancetolerance=3, vectortolerance=4)
    starOmatch = matches[1]
    starOmatchcatalog = catalog[starOmatch.catalogidx]

    # TYC 4757-1591-1
    truestarO = StarMatch.CatalogStar(475715911, hms2deg(5, 29, 23), dms2deg(-3, 26, 47), 5.92)
    @test isapprox(starOmatchcatalog.ra, truestarO.ra, atol=3*resolution(camera))

    #==
    Simulated wide, shallow field (mag +6)

    Use Tycho-2 catalog
    ==#
    camera = StarMatch.Camera(2048, 2048, 5e-6, 100e-3)

    f = CSV.file("SKY2000_Magnitude6_doublestars_0.12.txt")
    catalog = [StarMatch.CatalogStar(s.ID, s.RAJ2000, s.DEJ2000, s.mag) for s in f]

    spd = StarMatch.generatespd(camera, catalog)  # only takes a second

    # Star positions in simulated image
    f = CSV.File("test\\synthetic_254deg_+6.txt"; header=3)

    imagestars = StarMatch.CoordinateVector([SVector(x, y) for (x, y) in zip(f.PixelX, f.PixelY)])

    matches = StarMatch.solve(camera, imagestars, spd; distancetolerance=3, vectortolerance=4)
    starOmatch = matches[1]
    starOmatchcatalog = catalog[starOmatch.catalogidx]

    # TYC 4757-1591-1
    truestarO = StarMatch.CatalogStar(475715911, hms2deg(5, 29, 23), dms2deg(-3, 26, 47), 5.92)
    @test isapprox(starOmatchcatalog.ra, truestarO.ra, atol=3*resolution(camera))
end
