using CSV
using JLD2


f = CSV.file("URAT1-10mag.csv")

catalog = Array{CatalogStar}(undef, f.rows)
for (i, s) in enumerate(f)
  a, b = split(s.URAT1, '-')
  starid = parse(Int, a*b)
  catalog[i] = CatalogStar(starid, s.RAJ2000, s.DEJ2000, s.f_mag)
end

@save "URAT1-10mag_spd.jld2" catalog
