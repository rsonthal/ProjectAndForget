# """
# Functions used by metric optimization code for the paper:
#
# A Projection Method for Metric Constrained Optimization
# Nate Veldt, David Gleich, Anthony Wirth, and James Saunderson
# https://arxiv.org/abs/1806.01678
#
# """
# This version of the code is specifically designed for Julia 1.0

using SparseArrays

# """
# LamCC_DandW:
# Builds the weights matrix W and dissimilarity matrix D
# for a LambdaCC LP relaxation problem
# """
function LamCC_DandW(A,lam::Float64)
  n = size(A,1)
  D = zeros(Float64,n,n)
  W = (1-lam)*ones(n,n)
  for i = 1:n-1
    for j = i+1:n
      if A[i,j] < .1
        D[j,i] = 1
        W[j,i] = lam
      end
    end
  end
  D, W
end

# """
# LeightonRaoQP_Initialize:
# Initialize solution vector and weights matrix for the Leighton Rao quadratic
# programming relaxation.
# """
function LeightonRaoQP_Initialize(A,lam::Float64,gam::Float64)

  n = size(A,1)
  # X is a dense matrix
  X = zeros(Float64,n,n)
  W = lam*ones(n,n)
  @inbounds for i = 1:n-1
    for j = i+1:n
      if A[i,j] == 1
        X[j,i] = -gam
        W[j,i] = 1
      end
    end
  end
  for i = 1:n
    W[i,i] = 0.0
  end
  X, tril(W)
end


# """
# TriangleCheck:
# Returns whether or not all triangle inequality constraints are satisfied
# to within the desired tolerance. Returns as soon as it finds a constraint
# violation.
# """
function TriangleCheck(D::Matrix{Float64},tol::Float64)
  n = size(D,1)
  @inbounds for i = 1:n-2
    for j = i+1:n-1
      for k = j+1:n
        a = D[j,i]
        b = D[k,i]
        c = D[k,j]
        if a - b - c > tol || b - c - a > tol || c - a - b > tol
          return false
        end
      end
    end
  end
  return true
end

# FullTriangleCheck
# Returns the worst triangle violation in the whole matrix
function FullTriangleCheck(D::Matrix{Float64})
  n = size(D,1)
  maxi = 0.0
  @inbounds for i = 1:n-2
    for j = i+1:n-1
      a = D[j,i]
      for k = j+1:n
        b = D[k,i]
        c = D[k,j]

        #vio = maximum([a-b-c,b-a-c,c-a-b])
        if a-b > maxi && a-c > maxi && a-b-c > maxi
          maxi = a-b-c
        end
        if b - a > maxi && b-c > maxi && b-c-a > maxi
          maxi = b-c-a
        end
        if c-a > maxi && c-b > maxi && c-a-b > maxi
          maxi = c-a-b
        end

      end
    end
  end
  maxi
end

# """
# DoubleCheck:
# The L1 metric nearness problem (which the correlation clustering LP relaxation
# is a special case of) contains constraints eij - fij <= 0, -eij - fij <= 0.
# This function checks if they are satisfied to within a certain tolerance.
# """
function DoubleCheck(E::Matrix{Float64},F::Matrix{Float64},tol::Float64)
n = size(E,1)
for i = 1:n-1
  for j = i+1:n
    eij = E[j,i]
    fij = F[j,i]
    if eij - fij > tol || -eij - fij > tol
      return false
    end
  end
end
return true
end


# FullTriangleCheck
# Returns the worst triangle violation in the whole matrix
function FullTriangleCheck(D::Matrix{Float64})
  n = size(D,1)
  maxi = 0.0
  @inbounds for i = 1:n-2
    for j = i+1:n-1
      a = D[j,i]
      for k = j+1:n
        b = D[k,i]
        c = D[k,j]

        #vio = maximum([a-b-c,b-a-c,c-a-b])
        if a-b > maxi && a-c > maxi && a-b-c > maxi
          maxi = a-b-c
        end
        if b - a > maxi && b-c > maxi && b-c-a > maxi
          maxi = b-c-a
        end
        if c-a > maxi && c-c > maxi && c-a-b > maxi
          maxi = c-a-b
        end

      end
    end
  end
  maxi
end

# Checking all constraints for the correlation clustering LP problem
function FullConstraintCheck(X::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})

  tri = FullTriangleCheck(X)
  n = size(X,1)
  maxi = 0.0

  for i = 1:n-1
    for j = i+1:n
      eij = E[j,i]
      fij = F[j,i]
      if eij - fij > maxi
        maxi = eij - fij
      end
      if -eij - fij > maxi
        maxi = -eij - fij
      end
    end
  end
  doublecheck = maxi

  return round(max(tri,doublecheck),digits = 4), round(tri,digits = 4)
end


# Evaluate the LambdaCC LP relaxed objective, given a distance matrix D
function LPcc_obj(A,D::Matrix{Float64},lam::Float64)
  n = size(A,1)
  # assert(issymmetric(D))
  # assert(issymmetric(A))
  numedges = countnz(A)/2
  lccBound = sum((A[j,i]-lam)*D[j,i] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
  return lccBound
end

# This is effectively a vector dot product, but for matrices.
# Specifically, this corresponds to the linear program objective score for
# variables F.
function LPobj(F::Matrix{Float64},W::Matrix{Float64})
  n = size(F,1)
  obj = 0.0
  for i = 1:n-1
    for j = i+1:n
      obj += F[j,i]*W[j,i]
    end
  end
  return obj
end

# Computes the norm for the vector of variables for the correlation clustering problem
function Wnorm(W::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})

n = size(W,1)
out = 0.0
for i = 1:n-1
  for j = i+1:n
    out += (E[j,i]^2 + F[j,i]^2)*W[j,i]
  end
end
return out

end

# Computes the norm for the vector of variables for the Leighton-Rao relaxaton LP
function xWnorm(W::Matrix{Float64},X::Matrix{Float64})
  n = size(W,1)
  out = 0.0
  for i = 1:n-1
    for j = i+1:n
      out += (X[j,i]^2)*W[j,i]
    end
  end
  return out
end

# Computes the objective for the Leighton Rao sparsest cut LP relaxation
function LR_obj(A,X::Matrix{Float64})
  n = size(A,1)
  lr = sum((A[j,i]*X[j,i] for i=1:n-1 for j = i+1:n))
  return lr
end
