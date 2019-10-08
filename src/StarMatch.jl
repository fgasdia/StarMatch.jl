"""
Solve the ``lost in space'' problem.

An implementation of Approach 3 (_Rotation Invariant Vector Frame_) from:
``A reliable and fast lost-in-space mode star tracker'' by Metha Deval Samirbhai.
"""
module StarMatch

using LinearAlgebra
using StaticArrays
using CSV

struct CatalogStar
  id::Int
  ra::Float64
  dec::Float64
  mag::Float64
end

struct Star
    id::Int
    ra::Float64
    dec::Float64
    xy::SVector{2, Float64}
    vec::SVector{2, Float64}
    uvec::SVector{2, Float64}
    d::Float64
end

struct SPDEntry{N}
    Oid::Int
    SPid::Int
    d::Float64
    path::Array{SVector{2, Float64}, N}
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
    for starO in catalog
        # Transform from RA/DEC to ECI XYZ
        VeciO = radec2eci(starO.ra, starO.dec)
        CO = cameraattitude(VeciO)
        xyO = camera2image(CO*VeciO, camera)  # should be the center pixel

        # Identify neighboring stars in FOV
        neighbors = Star[]
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

                push!(neighbors, Star(s.id, s.ra, s.dec, xy, svec, usvec, d))
            end
        end

        # Order by distance
        sort!(neighbors, by = x -> x.d)

        #==
        NOTE See line 117 in spr_riav-master SPD_generate.m

        He removes the first 4-1 entries of his bound_vector. Is this necessary,
        or can I include them as I make the path? (obviously not including the
        stars that are closer than the current SP star)

        We want to include them if possible because we have relatively few stars
        ==#

        # We try 4 nearest stars to star O as starting points
        # We only loop to 1 less than length(neighbors)-1 because at the last
        # neighbor, there are not enough stars for a path
        # NOTE Samirbhai requires 5, we may need more than 3 for a robust solution
        length(neighbors) < 3 && continue

        for i in 1:min(4, length(neighbors)-1)
            # Save these outside the inner loop because it's used repeatedly
            V1 = neighbors[i].xy - xyO
            V1u = V1/norm(V1)

            # Generate path for current SP
            sn = @view neighbors[i:end]
            path = Array{SVector{2, Float64}}(undef, length(sn)-1)
            Vbsum = zero(V1)
            for j in 1:length(sn)-1
                k = j + 1

                Vj = sn[j].vec
                Vk = sn[k].vec
                Vku = sn[k].uvec

                Vbsum += (Vk - Vj)

                path[j] = SVector(dot(Vbsum, V1u), dot(Vbsum, Vku))
            end
            push!(spd, SPDEntry(starO.id, neighbors[i].id, neighbors[i].d, path))
        end
    end
    return spd
end

end # module
