nearest(A::AbstractArray,t) = A[findmin(abs.(A.-t))[2]]
findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]
