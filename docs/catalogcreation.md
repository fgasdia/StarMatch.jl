# URAT-1 Catalog

URAT-1 is a USNO astrometric catalog from magnitude ~3 to 18.5 for all stars in
the northern hemisphere down to declination -15°. Completed 2015.

Download data from the TAPVizieR service using ADQL (a SQL extension for astronomy).

Website: http://tapvizier.u-strasbg.fr/adql/?%20I/329/urat1

- Select the `I/329/urat1` checkbox for URAT1 catalog
- Choose Max records `all`
- Output format: `csv`
- Under `Favorite tables available to construct queries`, click `Construct your query`
- Use the query
```
-- output format : csv
SELECT "I/329/urat1".URAT1,  "I/329/urat1".RAJ2000,  "I/329/urat1".DEJ2000,  "I/329/urat1".Epoch,  "I/329/urat1"."f.mag",  "I/329/urat1"."e_f.mag",  "I/329/urat1".pmRA,
"I/329/urat1".pmDE
FROM "I/329/urat1"
WHERE "I/329/urat1"."f.mag"<10
```
- Click `Run`
- When phase is `COMPLETED`, download csv

# Gaia Data Release 2 (Gaia DR2)

Gaia is an ESA spacecraft to map positions of over 1 billion stars from magnitude
~3 to 21. This is more complete and should generally be preferred to URAT.
Completed 2018.

Download data from the ESA Gaia Archive. Reference epoch is Julian Year in TCB.
For DR2 the epoch is always J2015.5. RA and DEC are in decimal degrees ICRS
at the reference epoch. Proper motions are in mas/year.

Website: https://gea.esac.esa.int/archive/

- Click Search
- Select Advanced (ADQL)
- Choose Download format: CSV
- Use the query
```
SELECT source_id, ref_epoch, ra, dec, pmra, pmdec, phot_g_mean_mag
FROM gaiadr2.gaia_source
WHERE phot_g_mean_mag<11 AND dec>-70
```
where we restrict declinations to > -70°. The lowest declination that can be
observed from a site of latitude ϕ is ϕ - 90°.
- Click `Submit Query`
- When Status is complete, click `Download results`
