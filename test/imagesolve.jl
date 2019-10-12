using Test
using CSV
using FITSIO
using JLD2
using StaticArrays

using StarMatch


f = FITS("axy.fits")

X = read(f[2], "X")
Y = read(f[2], "Y")

# RASA/Manta
camera = Camera(1936, 1216, 5.86e-6, 620e-3)


f = FITS("test\\axy.fits")

X = read(f[2], "X")
Y = read(f[2], "Y")

imagestars = [SVector(x, y) for (x, y) in zip(X, Y)]

f = CSV.file("URAT1-10mag.csv")

catalog = Array{CatalogStar}(undef, f.rows)
for (i, s) in enumerate(f)
  a, b = split(s.URAT1, '-')
  starid = parse(Int, a*b)
  catalog[i] = CatalogStar(starid, s.RAJ2000, s.DEJ2000, s.f_mag)
end

spd = generatespd(camera, catalog)


@testset "Image solve" begin
  # Simulated image
  camera = StarMatch.Camera(2048, 2048, 5e-6, 100e-3)

  f = CSV.file("test\\SKY2000_Magnitude6_doublestars_0.12.txt")
  catalog = [StarMatch.CatalogStar(s.ID, s.RAJ2000, s.DEJ2000, s.mag) for s in f]
  spd = StarMatch.generatespd(camera, catalog)

  # Star positions in simulated image
  using FITSIO
  f = FITS("test\\axy_syntheticim_0deg.fits")
  X = read(f[2], "X")
  Y = read(f[2], "Y")

  imagestars = [SVector(x, y) for (x, y) in zip(X, Y)]

  winner = StarMatch.solve(camera, imagestars, spd; distancetolerance=4, vectortolerance=6)
  # dec: 25.897091
  # -5 h, 23.38176 m
  # -5 h, 23 m, 22.9056 s
  # ra: 84.934162
  # 5 h, 39.7366480 m
  # 5 h, 39 m, 44.19888 s


  @test winner.starO.ra
end
