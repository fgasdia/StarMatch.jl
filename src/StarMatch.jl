"""
Solve the ``lost in space'' problem.

An implementation of Approach 3 (_Rotation Invariant Vector Frame_) from:
``A reliable and fast lost-in-space mode star tracker'' by Metha Deval Samirbhai.
"""
module StarMatch

using LinearAlgebra
using StaticArrays
using CSV
using ProgressMeter

const N_NEAREST = 4
const VOTE_THRESHOLD = 3

struct CatalogStar
  id::Int
  ra::Float64
  dec::Float64
  mag::Float64
end

struct UnknownStar
    id::Int
    xy::SVector{2, Float64}
    vec::SVector{2, Float64}
    uvec::SVector{2, Float64}
    d::Float64
end

struct KnownStar
    id::Int
    ra::Float64
    dec::Float64
    xy::SVector{2, Float64}
    vec::SVector{2, Float64}
    uvec::SVector{2, Float64}
    d::Float64
end

struct SPDEntry
    starO::CatalogStar
    SPid::Int
    d::Float64
    path::Array{SVector{2, Float64}}
end

struct Camera
    img_width::Int
    img_height::Int
    pixelsize::Float64
    focallength::Float64
    fov::Float64
end
function Camera(w, h, pixelsize, focallength)
    imdim = max(h, w)
    fov = atand(imdim*pixelsize/2/focallength)*2
    Camera(w, h, pixelsize, focallength, fov)
end

"""
Transform RA/DEC to ECI unit vector.

RA and DEC in degrees.

See Vallado 4-1.
"""
function radec2eci(ra, dec)
    Cδ = cosd(dec)
    return SVector(cosd(ra)*Cδ, sind(ra)*Cδ, sind(dec))
end

"""
Transform from RA/DEC to image plane XY.

See Vallado chapter 4 and the SPD code from Samirbhai.
"""
function cameraattitude(Veci)
    tmp = sqrt(Veci[1]^2 + Veci[2]^2)
    c1x = Veci[2]/tmp
    c1y = -Veci[1]/tmp
    c1z = 0
    c1 = SVector(c1x, c1y, c1z)

    c2 = cross(Veci, c1)

    C = @SMatrix [c1[1] c1[2] c1[3];
                  c2[1] c2[2] c2[3];
                  Veci[1] Veci[2] Veci[3]]

    return C
end

function camera2image(Vcam, camera)
    f = camera.focallength
    ρ = camera.pixelsize
    w = camera.img_width
    h = camera.img_height

    x = f/ρ*Vcam[1]/Vcam[3] + w/2
    y = f/ρ*Vcam[2]/Vcam[3] + h/2

    return SVector(x, y)
end

function buildpath(neighbors::AbstractArray, V1, V1u)
    path = Array{SVector{2, Float64}}(undef, length(neighbors)-1)
    Vbsum = zero(V1)
    @inbounds for j in 1:length(neighbors)-1
        k = j + 1

        Vj = neighbors[j].vec
        Vk = neighbors[k].vec
        Vku = neighbors[k].uvec

        Vbsum += (Vk - Vj)

        path[j] = SVector(dot(Vbsum, V1u), dot(Vbsum, Vku))
    end
    return path
end

function vote(
    candidateindices,
    neighborpaths,
    spd::AbstractArray{SPDEntry};
    vectortolerance=5)

    votes = zeros(Int, length(candidateindices))
    for (c, (i, j)) in enumerate(candidateindices)
        for sp in spd[j].path
            for np in neighborpaths[i]
                if norm(sp - np) < vectortolerance
                    votes[c] += 1
                end
            end
        end
    end

    return votes
end

"""
`f` needs to use Julia's Tables.jl interface with the following fields:
    - `id`
    - `RA`
    - `DEC`
"""
function generatespd(camera::Camera, catalog::Array{CatalogStar})
    # Common values
    halffov = camera.fov/2
    coshalffov = cosd(halffov)

    spd = SPDEntry[]
    @showprogress 5 "Building SPD..." for starO in catalog
        # Transform from RA/DEC to ECI XYZ
        VeciO = radec2eci(starO.ra, starO.dec)
        CO = cameraattitude(VeciO)
        xyO = camera2image(CO*VeciO, camera)  # should be the center pixel

        # Identify neighboring stars in FOV
        neighbors = KnownStar[]
        @inbounds for s in catalog
            Veci = radec2eci(s.ra, s.dec)

            # Is star `s` within fov/2 of `starO`?
            # based on dot product / cosine similarity
            if dot(Veci, VeciO) > coshalffov && Veci != VeciO
                # Star coordinates in camera frame
                Vcam = CO*Veci

                # Star coordinates on image plane
                xy = camera2image(Vcam, camera)

                # Vector from `starO` (center of image) to star `s`
                svec = xy - xyO

                # Distance from `starO`
                d = norm(svec)

                # Unit vector from `starO` to star `s`
                usvec = svec/d

                push!(neighbors, KnownStar(s.id, s.ra, s.dec, xy, svec, usvec, d))
            end
        end

        # Order by distance
        sort!(neighbors, by = x -> x.d)

        #==
        NOTE See line 117 in spr_riav-master SPD_generate.m

        He removes the first 4-1 entries of his bound_vector. Is this necessary,  or can I
        include them as I make the path? (obviously not including the stars that are closer
        than the current SP star)

        We want to include them if possible because we have relatively few stars
        ==#
        length(neighbors) < N_NEAREST && continue

        for i in 1:N_NEAREST
            # Save these outside the inner loop because it's used repeatedly
            V1 = neighbors[i].vec
            V1u = neighbors[i].uvec

            # Generate path for current SP
            sn = @view neighbors[i:end]
            path = buildpath(sn, V1, V1u)

            push!(spd, SPDEntry(starO, neighbors[i].id, neighbors[i].d, path))
        end
    end
    return spd
end

function solve(
    camera::Camera,
    imagestars::AbstractArray{SVector{2, T}},
    spd::AbstractArray{SPDEntry};
    distancetolerance=1) where T <: Real

    @assert length(imagestars) > N_NEAREST "A minimum of 5 stars are required to solve image."

    imagecenter = SVector(camera.img_height/2, camera.img_width/2)

    # Determine star closest to center. This will be `starO`
    mindist = typemax(T)
    Oidx = 0
    for (i, star) in enumerate(imagestars)
        # Vector from star to image center
        svec = star - imagecenter

        # Distance from star to image center
        d = norm(svec)

        # Check if this star is closer to center than other stars
        if d < mindist
            mindist = d
            Oidx = i
        end
    end
    starO = imagestars[Oidx]

    # Calculate vectors from each star to `starO`
    neighbors = Array{UnknownStar}(undef, length(imagestars)-1)
    ni = 1
    for (i, star) in enumerate(imagestars)
        if star != starO
            svec = star - starO
            d = norm(svec)
            usvec = svec/d
            neighbors[ni] = UnknownStar(i, SVector(star), svec, usvec, d)
            ni += 1
        # else starOid = i  # just fyi
        end
    end
    sort!(neighbors, by = x -> x.d)

    # Check for distance matches with SPD
    candidateindices = Tuple{Int,Int}[]
    for (j, entry) in enumerate(spd)
        for i in 1:N_NEAREST
            if norm(entry.d - neighbors[i].d) < distancetolerance
                push!(candidateindices, (i, j))
            end
        end
    end

    neighborpaths = Array{Array{SVector{2, Float64}}}(undef, 4)
    for i in 1:N_NEAREST
        V1 = neighbors[i].vec
        V1u = neighbors[i].uvec

        sn = @view neighbors[i:end]
        neighborpaths[i] = buildpath(sn, V1, V1u)
    end

    votecounts = vote(candidateindices, neighborpaths, spd)

    winnerid = argmax(votecounts)
    winnervotecount = votecounts[winnerid]
    winnervotecount > VOTE_THRESHOLD && println("Not enough votes!")

    winner = spd[candidateindices[winnerid][2]]
end





# f = FITS("test\\axy.fits")
#
# X = read(f[2], "X")
# Y = read(f[2], "Y")
#
# imagestars = [SVector(x, y) for (x, y) in zip(X, Y)]
#
#
# # RASA/Manta
# camera = Camera(1936, 1216, 5.86e-6, 620e-3)
#
# f = CSV.file("URAT1-10mag.csv")
#
# catalog = Array{CatalogStar}(undef, f.rows)
# for (i, s) in enumerate(f)
#   a, b = split(s.URAT1, '-')
#   starid = parse(Int, a*b)
#   catalog[i] = CatalogStar(starid, s.RAJ2000, s.DEJ2000, s.f_mag)
# end
#
# spd = generatespd(camera, catalog)


using JLD2
using FITSIO

# Simulated image
camera = Camera(2048, 2048, 5e-6, 100e-3)

# @load "SKY2000_Magnitude6_doublestars_0.12_spd.jd2" spd

f = CSV.file("test\\SKY2000_Magnitude6_doublestars_0.12.txt")

catalog = Array{CatalogStar}(undef, f.rows)
for (i, s) in enumerate(f)
    catalog[i] = CatalogStar(s.ID, s.RAJ2000, s.DEJ2000, s.mag)
end
spd = generatespd(camera, catalog)


f = FITS("test\\axy_m42.fits")
X = read(f[2], "X")
Y = read(f[2], "Y")

imagestars = [SVector(x, y) for (x, y) in zip(X, Y)]



imagecenter = SVector(camera.img_height/2, camera.img_width/2)

# Determine star closest to center. This will be `starO`
mindist = typemax(T)
Oidx = 0
for (i, star) in enumerate(imagestars)
    # Vector from star to image center
    svec = star - imagecenter

    global mindist
    global Oidx

    # Distance from star to image center
    d = norm(svec)

    # Check if this star is closer to center than other stars
    if d < mindist
        mindist = d
        Oidx = i
    end
end
starO = imagestars[Oidx]

# Calculate vectors from each star to `starO`
neighbors = Array{UnknownStar}(undef, length(imagestars)-1)
ni = 1
for (i, star) in enumerate(imagestars)
    global ni
    if star != starO
        svec = star - starO
        d = norm(svec)
        usvec = svec/d
        neighbors[ni] = UnknownStar(i, SVector(star), svec, usvec, d)
        ni += 1
    # else starOid = i  # just fyi
    end
end
sort!(neighbors, by = x -> x.d)

# Check for distance matches with SPD
candidateindices = Tuple{Int,Int}[]
for (j, entry) in enumerate(spd)
    for i in 1:N_NEAREST
        if norm(entry.d - neighbors[i].d) < distancetolerance
            push!(candidateindices, (i, j))
        end
    end
end

neighborpaths = Array{Array{SVector{2, Float64}}}(undef, 4)
for i in 1:N_NEAREST
    V1 = neighbors[i].vec
    V1u = neighbors[i].uvec

    sn = @view neighbors[i:end]
    neighborpaths[i] = buildpath(sn, V1, V1u)
end

votecounts = vote(candidateindices, neighborpaths, spd)

winnerid = argmax(votecounts)
winnervotecount = votecounts[winnerid]
winnervotecount > VOTE_THRESHOLD && println("Not enough votes!")

winner = spd[candidateindices[winnerid][2]]





end # module
