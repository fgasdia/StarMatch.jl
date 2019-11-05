__precompile__()

module StarMatch

using Logging
using LinearAlgebra
using StaticArrays
using ProgressMeter
using PushVectors

export generatespd, solve
export CoordinateVector, Camera, CatalogStar

const N_NEAREST = 5
const VOTE_THRESHOLD = 3
const VECTORTOLERANCE = 4

struct CoordinateVector{T<:Real} <: AbstractVector{T}
    p::Vector{SVector{2,T}}
end
CoordinateVector(::Type{T}, dims::Int) where {T} = CoordinateVector{T}(Vector{SVector{2,T}}(undef, dims));
Base.size(A::CoordinateVector) = size(A.p)
Base.IndexStyle(::Type{<:CoordinateVector}) = IndexLinear()
Base.getindex(A::CoordinateVector, i::Int) = A.p[i]
Base.setindex!(A::CoordinateVector{T}, v, i::Int) where {T} = (A.p[i] = v)
Base.similar(A::CoordinateVector{T}) where {T} = CoordinateVector(T, size(A))
Base.similar(A::CoordinateVector, ::Type{T}, dims::Dims) where {T} = CoordinateVector(T, dims)
Base.eltype(A::CoordinateVector{T}) where {T} = T

struct CatalogStar
  id::Int
  ra::Float64
  dec::Float64
  mag::Float64
end

abstract type ImageStar end

struct UnknownStar <: ImageStar
    xy::SVector{2, Float64}
    vec::SVector{2, Float64}
    uvec::SVector{2, Float64}
    d::Float64
end

struct KnownStar <: ImageStar
    catalogidx::Int
    xy::SVector{2, Float64}
    vec::SVector{2, Float64}
    uvec::SVector{2, Float64}
    d::Float64
end
Base.zero(::Type{KnownStar}) = KnownStar(0, SVector(0.0,0.0), SVector(0.0,0.0), SVector(0.0,0.0), 0.0)

struct SPDEntry
    starOidx::Int  # catalog idx
    SPidx::Int  # catalog idx
    d::Float64
    path::CoordinateVector{Float64}
    pathidxs::Vector{Int}  # catalog idx
end

struct Camera
    img_width::Int
    img_height::Int
    img_center::SVector{2,Float64}
    pixelsize::Float64
    focallength::Float64
    fov::Float64
end
function Camera(w, h, pixelsize, focallength)
    imdim = max(h, w)
    fov = atand(imdim*pixelsize/2/focallength)*2
    center = SVector(w/2, h/2)
    Camera(w, h, center, pixelsize, focallength, fov)
end

"""
Transform RA/DEC to ECI unit vector.

RA and DEC in degrees.

See Vallado 4-1.
"""
function radec2eci(ra, dec)
    sinδ, cosδ = sincos(deg2rad(dec))
    sinα, cosα = sincos(deg2rad(ra))
    return SVector(cosα*cosδ, sinα*cosδ, sinδ)
end

"""
Transform from RA/DEC to image plane XY.

See Vallado chapter 4 and the SPD code from Samirbhai.
"""
function cameraattitude(Veci)
    tmp = 1/sqrt(Veci[1]^2 + Veci[2]^2)
    c1x = Veci[2]*tmp
    c1y = -Veci[1]*tmp
    c1z = 0
    c1 = SVector(c1x, c1y, c1z)

    c2 = cross(Veci, c1)

    C = @SMatrix [c1[1] c1[2] c1[3];
                  c2[1] c2[2] c2[3];
                  Veci[1] Veci[2] Veci[3]]

    return C
end

function camera2image(Vcam, camera::Camera)
    f = camera.focallength
    ρ = camera.pixelsize
    w = camera.img_width
    h = camera.img_height

    x = f/ρ*Vcam[1]/Vcam[3] + w/2
    y = f/ρ*Vcam[2]/Vcam[3] + h/2

    return SVector(x, y)
end

function buildpath(neighbors::AbstractVector{<:ImageStar}, V1, V1u)
    path = CoordinateVector(Float64, length(neighbors)-1)
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



# TODO
#XXX: This function completely dominates the runtime of `solve`
# There is a huge memory allocation too, although I do not know
# why. The norm(sp-np) is allocating 112 bytes and if done hundreds of thousands of times
# this might add up?
function vote(candidateindices::AbstractVector{Tuple{Int,Int}}, neighborpaths::AbstractVector{CoordinateVector},
    spd::AbstractVector{SPDEntry}, vectortolerance=VECTORTOLERANCE)

    # TODO: Reuse `votes` when there are multiple images
    votes = zeros(UInt32, length(candidateindices))
    @inbounds for (c, (i, j)) in enumerate(candidateindices)
        votecount = 0
        @inbounds for sp in spd[j].path
            @inbounds for np in neighborpaths[i]
                if norm(sp - np) < vectortolerance
                    votecount += 1
                end
            end
        end
        votes[c] += votecount
    end
    return votes
end

function closest(X::CoordinateVector{T}, a::SVector{2,T}) where T <: Real
    mindist = typemax(T)
    mindistidx = 0
    for (i, x) in enumerate(X)
        dist = norm(a - x)
        if dist < mindist
            mindist = dist
            mindistidx = i
        end
    end
    return mindistidx
end

function match(starO::SVector{2,T}, winner::SPDEntry, neighbors::AbstractVector{UnknownStar},
    neighborpath::CoordinateVector, vectortolerance=VECTORTOLERANCE) where T <: Real

    matches = PushVector{KnownStar}(length(neighbors))

    # conversion needed of `starO` in case its type is not Float64
    push!(matches, KnownStar(winner.starOidx, SVector{2,Float64}(starO), SVector(0.0, 0.0), SVector(0.0, 0.0), 0.0))

    for i in eachindex(winner.path)
        p = winner.path[i]

        idx = closest(neighborpath, p)
        # idx = findfirst(x -> norm(x - p) < vectortolerance, neighborpath)

        if norm(neighborpath[idx] - p) < vectortolerance
            push!(matches, KnownStar(winner.pathidxs[i],
                neighbors[idx].xy, neighbors[idx].vec, neighbors[idx].uvec, neighbors[idx].d))
        end
    end
    # TODO: handle the last neighbor (neighborpath is 1 shorter than neighbors)

    return matches
end

"""
    generatespd(camera, catalog)

Generate a Star Pattern Database for `camera` and a `catalog` of `CatalogStar`s.
"""
function generatespd(camera::Camera, catalog::Vector{CatalogStar})
    # Common values
    halffov = camera.fov/2
    coshalffov = cosd(halffov)
    imagecenter = camera.img_center

    # Calculating star position in ECI is surprisingly expensive, so we do it once
    catalogeci = [radec2eci(s.ra, s.dec) for s in catalog]

    spd = PushVector{SPDEntry}(length(catalog)*N_NEAREST)
    neighbors = PushVector{KnownStar}()

    @showprogress 5 "Building SPD..." for (oi, starO) in enumerate(catalog)
        VeciO = catalogeci[oi]
        CO = cameraattitude(VeciO)
        xyO = imagecenter  # = camera2image(CO*VeciO, camera)

        # Identify neighboring stars in FOV
        @inbounds for (si, s) in enumerate(catalog)
            Veci = catalogeci[si]

            # Is star `s` within fov/2 of `starO`?
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

                push!(neighbors, KnownStar(si, xy, svec, usvec, d))
            end
        end

        sort!(neighbors, by = x -> x.d)

        #==
        NOTE Samirbhai removes the first N_NEAREST-1 entries of his boundary vector (path).

        They are included here because we tend to have few stars in our FOV and want to use
        as many of them as possible.

        By enforcing there be N_NEAREST+2 neighbors, we ensure the path is at least 2 long
        ==#
        length(neighbors) < (N_NEAREST+2) && continue

        @inbounds for i in 1:N_NEAREST
            V1 = neighbors[i].vec
            V1u = neighbors[i].uvec

            sn = @view neighbors[i:end]
            path = buildpath(sn, V1, V1u)
            pathidxs = [n.catalogidx for n in sn]

            push!(spd, SPDEntry(oi, neighbors[i].catalogidx, neighbors[i].d, path, pathidxs))
        end

        empty!(neighbors)
    end
    PushVectors.finish!(spd)
    return spd
end

"""
    solve(camera, imagestars, spd; distancetolerance=3, vectortolerance=VECTORTOLERANCE)

Return catalog matches from `spd` for `imagestars` extracted from an image taken with
`camera`. `distancetolerance` is the number of pixels within which the distance from the
nearest neighbor stars must be to the center star compared to the distances stored in the
`spd`. `vectortolerance` is the tolerance used in `vote` and `match` to check the `norm` of
the difference of path vectors.
"""
function solve(camera::Camera, imagestars::CoordinateVector{T},
    spd::AbstractVector{SPDEntry}, distancetolerance=3, vectortolerance=VECTORTOLERANCE) where T <: Real

    @assert length(imagestars) > N_NEAREST "A minimum of $(N_NEAREST+1) stars required to solve image."

    imagecenter = camera.img_center

    # Determine star closest to center. This will be `starO`
    starOidx = closest(imagestars, imagecenter)
    starO = imagestars[starOidx]

    # Calculate vectors from each star to `starO`
    neighbors = Vector{UnknownStar}(undef, length(imagestars)-1)
    ni = 1
    for star in imagestars
        if star != starO
            svec = star - starO
            d = norm(svec)
            usvec = svec/d

            neighbors[ni] = UnknownStar(star, svec, usvec, d)
            ni += 1
        # else starOid = i  # just fyi
        end
    end
    sort!(neighbors, by = x -> x.d)

    # TODO: Reuse `candidateindices` when we handle multiple images
    # Check for distance matches with SPD
    candidateindices = PushVector{Tuple{Int,Int}}()
    for j in eachindex(spd)
        spddist = spd[j].d
        @inbounds for i in 1:N_NEAREST
            if norm(spddist - neighbors[i].d) < distancetolerance
                push!(candidateindices, (i, j))
            end
        end
    end
    PushVectors.finish!(candidateindices)
    # @info "$(length(candidateindices)) candidate indices"

    neighborpaths = Vector{CoordinateVector}(undef, N_NEAREST)
    @inbounds for i in 1:N_NEAREST
        V1 = neighbors[i].vec
        V1u = neighbors[i].uvec

        sn = @view neighbors[i:end]
        neighborpaths[i] = buildpath(sn, V1, V1u)
    end

    votecounts = vote(candidateindices, neighborpaths, spd, vectortolerance)

    winnerid = argmax(votecounts)
    winnervotecount = votecounts[winnerid]
    # TODO: Warn if there is a tie

    if winnervotecount < VOTE_THRESHOLD
        @warn "Winner vote count is below threshold! ($winnervotecount/$VOTE_THRESHOLD)"
    end

    winningindices = candidateindices[winnerid]
    winner = spd[winningindices[2]]

    sn = @view neighbors[winningindices[1]:end]
    matches = match(starO, winner, sn, neighborpaths[winningindices[1]], vectortolerance)

    return matches
end


#
# #==
# Utility functions
# ==#
#
end # module
