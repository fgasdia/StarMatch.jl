using StarMatch

camera = Camera(1024, 1024, 13.3e-6, 0.0, 15.)

f = CSV.file("test\\SKY2000_Magnitude6_doublestars_0.12.txt")

catalog = Array{CatalogStar}(undef, f.rows)
for (i, s) in enumerate(f)
  catalog[i] = CatalogStar(s.ID, s.RAJ2000, s.DEJ2000, s.mag)
end

spd = generatespd(camera, catalog)
