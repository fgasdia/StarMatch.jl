"""
Solve the ``lost in space'' problem.

An implementation of Approach 3 (_Rotation Invariant Vector Frame_) from:
``A reliable and fast lost-in-space mode star tracker'' by Metha Deval Samirbhai.
"""
module StarMatch

using LinearAlgebra
using StaticArrays
using CSV


struct Camera
    img_width::Int
    img_height::Int
    pixelsize::Float64
    focallength::Float64
    fov::Float64
    function Camera(w, h, pixelsize, focallength)
        imdim = max(h, w)
        fov = atand(imdim*pixelsize/2/focallength)*2
        new(w, h, pixelsize, focallength, fov)
    end
end

struct Star
    ida::Int
    idb::Int
    mag::Float64
    ra::Float64
    dec::Float64
    xy::SVector
    d::Float64
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

    return (x, y)
end


function read_database(file="URAT1.csv")
    return CSV.file(file)
end

function distance(xy, xyo)
    x, y = xy
    xo, yo = xyo
    return sqrt((x-xo)^2 + (y-yo)^2)
end

function generatespd(camera)
    # Common values
    halffov = camera.fov/2
    coshalffov = cosd(halffov)

    f = read_database()

    spd = []

    # for starO in f
        starO = iterate(f, 1)[1]

        # Transform from RA/DEC to ECI XYZ
        VeciO = radec2eci(starO.RAJ2000, starO.DEJ2000)
        CO = cameraattitude(VeciO)
        xyO = camera2image(CO*VeciO, camera)  # should be the center pixel

        # Identify neighboring stars in FOV
        neighbors = Array{Star}(undef,0)
        for s in f
            Veci = radec2eci(s.RAJ2000, s.DEJ2000)

            # Is star `s` within fov/2 of `starO`?
            # based on dot product / cosine similarity
            if dot(Veci, VeciO) > coshalffov && Veci != VeciO
                # Star coordinates in camera frame
                Vcam = CO*Veci

                # Star coordinates on image plane
                xy = camera2image(Vcam, camera)

                # Distance from `starO`
                d = distance(xy, xyO)

                # TEMP
                ida, idb = parse.(Int, split(s.URAT1, '-'))
                push!(neighbors, Star(ida, idb, s.f_mag, s.RAJ2000, s.DEJ2000, xy, d))
            end
        end

        # Order by distance
        sort!(neighbors, by = x -> x.d)

        # We use 4 nearest stars
        @assert length(neighbors) >= 3 "Fewer than 3 neighboring stars!"
        for i in 1:min(4, length(neighbors))
            # Save these outside the inner loop because it's used repeatedly
            V₁ = neighbors[i].xy
            V̂₁ = V₁/norm(V₁)

            sn = @view neighbors[i:end]
            path = MVector{length(sn)-1, Tuple{Float64, Float64}}(undef)
            Vbsum = zero(V₁)
            for j in 1:length(sn)-1
                k = j + 1

                Vⱼ = sn[j].xy
                Vₖ = sn[k].xy
                V̂ₖ = Vₖ/norm(Vₖ)

                Vbsum += (Vₖ - Vⱼ)

                path[j] = (dot(Vbsum, V̂₁), dot(Vbsum, V̂ₖ))
            end
            push!(spd, (neighbors[i].ida, neighbors[i].idb, neighbors[i].d, path))
        end

        return neighbors
    # end
end

# RASA 11 / Manta G235
camera = Camera(1936, 1216, 5.86e-6, 620e-3)


end # module
