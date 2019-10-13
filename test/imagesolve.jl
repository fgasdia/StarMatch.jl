using Test
using CSV
using FITSIO
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

@testset "Image solve" begin
    # Simulated image
    camera = StarMatch.Camera(2048, 2048, 5e-6, 100e-3)

    f = CSV.file("SKY2000_Magnitude6_doublestars_0.12.txt")
    catalog = [StarMatch.CatalogStar(s.ID, s.RAJ2000, s.DEJ2000, s.mag) for s in f]
    spd = StarMatch.generatespd(camera, catalog)

    # Star positions in simulated image without rotation
    using FITSIO
    f = FITS("axy_syntheticim_0deg_+6.fits")
    X = read(f[2], "X")
    Y = read(f[2], "Y")

    imagestars = StarMatch.CoordinateVector([SVector(x, y) for (x, y) in zip(X, Y)])

    matches = StarMatch.solve(camera, imagestars, spd; distancetolerance=3, vectortolerance=4)
    starOmatch = matches[1]
    starOmatchcatalog = catalog[starOmatch.catalogidx]

    # TYC 4757-1591-1
    truestarO = StarMatch.CatalogStar(475715911, hms2deg(5, 29, 23), dms2deg(-3, 26, 47), 5.92)
    @test isapprox(starOmatchcatalog.ra, truestarO.ra, atol=3*resolution(camera))

    # Star positions in simulated image with 118° rotation
    f = FITS("axy_syntheticim_118deg_+6.fits")
    X = read(f[2], "X")
    Y = read(f[2], "Y")

    imagestars = StarMatch.CoordinateVector([SVector(x, y) for (x, y) in zip(X, Y)])

    matches = StarMatch.solve(camera, imagestars, spd; distancetolerance=3, vectortolerance=4)
    starOmatch = matches[1]
    starOmatchcatalog = catalog[starOmatch.catalogidx]

    # TYC 4757-1591-1
    truestarO = StarMatch.CatalogStar(475715911, hms2deg(5, 29, 23), dms2deg(-3, 26, 47), 5.92)

    @test isapprox(starOmatchcatalog.ra, truestarO.ra, atol=3*resolution(camera))
end
