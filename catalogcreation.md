To download catalog

Website: http://tapvizier.u-strasbg.fr/adql/?%20I/329/urat1

Which is the TAPVizieR service, providing VizieR tables using ADQL (a SQL extension)
in Astronomy.

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
