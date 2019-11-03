# StarMatch

Solve the ``lost in space'' problem.

An implementation of Approach 3 (_Rotation Invariant Vector Frame_) from:
``A reliable and fast lost-in-space mode star tracker'' by Metha Deval Samirbhai.

# TODO

- Proper motion
- Usage interface
- Documentation
- Implement this so that I can keep track of whether or not SPD entry is sufficiently close or whatever else
  https://discourse.julialang.org/t/how-best-to-maintain-collection-of-indices-into-2-d-array/29921/2

- Most importantly, really need to address the voting function which appears to account for the vast majority of the solve time and allocates a huge amount of memory although I don't know why.
