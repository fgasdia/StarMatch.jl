using StarMatch

using FITSIO

f = FITS("axy.fits")

X = read(f[2], "X")
Y = read(f[2], "Y")

# RASA/Manta
camera = Camera(1936, 1216, 5.86e-6, 620e-3)

f = CSV.file("URAT1-10mag.csv")

catalog = Array{CatalogStar}(undef, f.rows)
for (i, s) in enumerate(f)
  a, b = split(s.URAT1, '-')
  starid = parse(Int, a*b)
  catalog[i] = CatalogStar(starid, s.RAJ2000, s.DEJ2000, s.f_mag)
end

imagestars = [SVector(x, y) for (x, y) in zip(X, Y)]



# Simulated image
camera = Camera(2048, 2048, 5e-6, 100e-3)
@load "SKY2000_Magnitude6_doublestars_0.12_spd.jd2" spd

f = FITS("test\\axy_m42.fits")
X = read(f[2], "X")
Y = read(f[2], "Y")

imagestars = [SVector(x, y) for (x, y) in zip(X, Y)]
