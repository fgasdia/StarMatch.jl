"""
Solve the ``lost in space'' problem.

An implementation of Approach 3 (_Rotation Invariant Vector Frame_) from:
``A reliable and fast lost-in-space mode star tracker'' by Metha Deval Samirbhai.
"""
module StarMatch

using LinearAlgebra
using StaticArrays
using CSV

# Construct a rotation-invariant vector frame for 2D vectors connecting the stars
# on a normal X-Y plane.

const R = MMatrix{2,2}(rand(2,2))

"""
θ is a rotation angle in radians measured in an anti-clockwise direction.
"""
function rotate!(R, θ)
    S, C = sin(θ), cos(θ)
    R[1] = C
    R[2] = S
    R[3] = -S
    R[4] = C
    return nothing
end

function read_database(file="URAT1.csv")
    return CSV.file(file)
end

f = read_database()

end # module
