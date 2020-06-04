include("MetricOpt_helper.jl")

using Base.Threads
using SparseArrays
using LinearAlgebra

function parallel_cyclic_triple_loopE!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},Bookshelf::Vector{Vector{Tuple{Int64, Float64}}}, newBookshelf::Vector{Vector{Tuple{Int64, Float64}}})
   #cyclic_triple_loop(D, E, W, now_corrections, next_corrections)
    n::Int = size(D,1)
    o::Int = 10
    nthds::Int = Threads.nthreads()
    epsi = 0
    # Construct the initial dissimilarity matrix, set X_0 = D essentially
    # We only edit the lower triangular part of the matrix. Inside of loops,
    # whenever possible we iterate through columns first, then rows.

    # This is an array of arrays of size nThreads, each thread gets it's own update array.

    #These arrays are used to store information about block size, e.g. index ranges in blocks of size 10
    A::Vector{Tuple{Int, Int}} = Vector{Tuple{Int, Int}}()
    B::Vector{Tuple{Int, Int}} = Vector{Tuple{Int, Int}}()
    for k in n:-o:1
        a = k-o+1
        b = k
        if a < 1
            a = 1
        end
        pushfirst!(B, (a,b))
    end

        A=B

    # Now perform Dykstra's method and plot intermediate steps
    # Each global iteration visits every constraint once

    indexing = Vector{Int}()
    for i in 1:nthds
        push!(indexing, 1)
    end

    Duals = Vector{Float64}()
    for i in 1:nthds
        push!(Duals, 0.0)
    end



    begin
        r_endpoint::Int = round(Int64, ceil(n/o))
        for m::Int in r_endpoint:-1:1
            #@threads
            Threads.@threads for threadholder in 1:nthds
                c::Int64 = Threads.threadid()

                nowInd = indexing[c]
                BtDual = Duals[c]
                old_triplet_corrections::Vector{Tuple{Int64, Float64}} = Bookshelf[c]
                new_triplet_corrections::Vector{Tuple{Int64, Float64}} = newBookshelf[c]

                numthreads::Int64 = Threads.nthreads()
                threadnum::Int64 = c
                start::Int64 = 1
                endpoint::Int64 = round(Int64, ceil(m/2))
                interval::Int = endpoint-start+1

                threadrange::Int = ceil(Int64, interval/numthreads)
                s::Int = start + threadrange*(threadnum-1)
                e::Int = min(start + threadrange*threadnum -1, endpoint)

                for iIter in s:e    #iIter::Int in 1:round(Int64, ceil(m/2))
                    kIter::Int = m-iIter+1
                    nowInd, BtDual = iterate_E!(D, E, W, BtDual, iIter, kIter, o, n, A, B, old_triplet_corrections, new_triplet_corrections, nowInd)
                end
                indexing[c] = nowInd
                Duals[c] = BtDual
            end
        end
    end

    begin
        r_endpoint = round(Int64, ceil(n/o))
        for m::Int in 2:r_endpoint
            #Threads.@threads
            Threads.@threads for threadholder in 1:nthds
                c::Int64 = Threads.threadid()
                nowInd::Int64 = indexing[c]
                BtDual::Float64 = Duals[c]
                #new_triplet_corrections = nothing#::Array{Tuple, 1} = Tuple[]#Bookshelf[c]
                old_triplet_corrections::Vector{Tuple{Int64, Float64}} = Bookshelf[c]
                new_triplet_corrections::Vector{Tuple{Int64, Float64}} = newBookshelf[c]

                numthreads::Int = Threads.nthreads()
                threadnum::Int = c
                start::Int = m
                endpoint::Int = round(Int64,floor((m+r_endpoint)/2))
                interval::Int = endpoint-start+1

                threadrange::Int = ceil(Int64, interval/numthreads)
                s::Int = start + threadrange*(threadnum-1)
                e::Int = min(start + threadrange*threadnum -1, endpoint)

                for iIter::Int in s:e #m:round(Int64,floor((m+r_endpoint)/2)) #thread shplit
                    kIter::Int = r_endpoint-iIter+m
                    nowInd, BtDual = iterate_E!(D, E, W, BtDual, iIter, kIter, o, n, A, B, old_triplet_corrections, new_triplet_corrections, nowInd)
                end
                indexing[c] = nowInd
                Duals[c] = BtDual
            end
        end
    end


    # We iterate through a triple for loop to sequentially visit constraints of
    # of the form X[j,i] - X[k,i] - X[k,j] <= 0
    DualVar = 0.0
    for dual in Duals
        DualVar += dual
    end


    DualVar

end # function end


function iterate_E!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64}, BtDual::Float64, iIter::Int, kIter::Int, o::Int, n::Int, A::Vector{Tuple{Int, Int}}, B::Vector{Tuple{Int, Int}}, now_corrections::Vector{Tuple{Int64, Float64}}, next_corrections::Vector{Tuple{Int64, Float64}}, nowInd::Int)
    epsi = 0
    iRange::Tuple{Int, Int} = A[iIter]
    kRange::Tuple{Int, Int} = B[kIter]

    i::Int = minimum(iRange)
    k::Int = maximum(kRange)

    C = Vector{Tuple{Int, Int}}()
    for j in k-1:-o:i+1
        a::Int = j-o+1
        b::Int = j
        if a < i+1
            a = i+1
        end
        push!(C, (a,b))
    end
    nextKey = now_corrections[nowInd]
    correctionsLength = length(now_corrections)

    for jIter::Int in 1:length(C)
        jRange::Tuple{Int, Int} = C[jIter]

        for i in iRange[1]:iRange[2]
            for j::Int in max(i+1, jRange[1]):jRange[2]


                 Dij = D[j,i]
                 Wij = W[j,i]

                for k in max(j+1,kRange[1]):kRange[2]
                    # i,j,k here satisfies i < j < k
                    # For each ordered triplet, there are 3 constriants:

                    #  X[j,i] - X[k,i] - X[k,j] <= 0
                    #  -X[j,i] + X[k,i] - X[k,j] <= 0
                    #  -X[j,i] - X[k,i] + X[k,j] <= 0
                        Dik = D[k,i]
                        Djk = D[k,j]
                        Eij = E[j,i]
                        Eik = E[k,i]
                        Ejk = E[k,j]
                        Wik = W[k,i]
                        Wjk = W[k,j]

                        ### Check triangle i,j,k
                        ijkKey = (i-1)*n^2+(j-1)*n+k
                        # First see if this is the next triangle with a nonzero correction variable
                        if ijkKey == nextKey[1]

                                cor = nextKey[2]
                          # We need to scale the correction since we're projecting into the minimum
                          # weighted 2-norm direction

                          denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                          E[j,i] = Eij + cor*(Wik*Wjk/denom)
                          E[k,i] = Eik - cor*(Wij*Wjk/denom)
                          E[k,j] = Ejk - cor*(Wik*Wij/denom)
                          Eij = E[j,i]
                          Eik = E[k,i]
                          Ejk = E[k,j]
                          # Move along in the list of triplets with corrections
                          if nowInd < correctionsLength
                                  nowInd +=1
                                  nextKey = now_corrections[nowInd]
                          end
                        end

                        b = Dik + Djk - Dij
                        mu = (Eij - Ejk - Eik - b)

                        if mu > epsi
                          denom = Wij*Wjk + Wik*Wij + Wjk*Wik
                          E[j,i] = Eij - mu*(Wik*Wjk/denom)
                          E[k,i] = Eik + mu*(Wij*Wjk/denom)
                          E[k,j] = Ejk + mu*(Wik*Wij/denom)
                          # Next time we see this triple we have to correct
                          push!(next_corrections,(ijkKey,mu))
                          BtDual += mu*b*(Wij*Wik*Wjk/denom)
                        end

                        ### Done checking triangle i,j,k


                        ### Check triangle i,k,j
                        ijkKey = (i-1)*n^2+(k-1)*n+j
                        # First see if this is the next triangle with a nonzero correction variable
                        if ijkKey == nextKey[1]

                                cor = nextKey[2]

                          denom = Wij*Wjk + Wik*Wij + Wjk*Wik
                          E[j,i] = E[j,i] - cor*(Wik*Wjk/denom)
                          E[k,i] = E[k,i] + cor*(Wij*Wjk/denom)
                          E[k,j] = E[k,j] - cor*(Wik*Wij/denom)
                          Eij = E[j,i]
                          Eik = E[k,i]
                          Ejk = E[k,j]
                          # Move along in the list of triplets with corrections
                          if nowInd < correctionsLength
                                  nowInd +=1
                                  nextKey = now_corrections[nowInd]
                          end
                        else
                          Eij = E[j,i]
                          Eik = E[k,i]
                          Ejk = E[k,j]
                        end

                        b = -Dik + Djk + Dij
                        mu = (-Eij - Ejk + Eik - b)
                        if mu > epsi

                          denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                          E[j,i] = Eij + mu*(Wik*Wjk/denom)
                          E[k,i] = Eik - mu*(Wij*Wjk/denom)
                          E[k,j] = Ejk + mu*(Wik*Wij/denom)

                          # Next time we see this triple we have to correct
                          push!(next_corrections,(ijkKey,mu))
                          BtDual += mu*b*(Wij*Wik*Wjk/denom)
                        end
                        ### Done checking triangle i,k,j

                        ### Triangle j,k,i
                        ijkKey = (j-1)*n^2+(k-1)*n+i
                        # First see if this is the next triangle with a nonzero correction variable
                        if ijkKey == nextKey[1]

                                cor = nextKey[2]

                          denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                          E[j,i] = E[j,i] - cor*(Wik*Wjk/denom)
                          E[k,i] = E[k,i] - cor*(Wij*Wjk/denom)
                          E[k,j] = E[k,j] + cor*(Wik*Wij/denom)
                          Eij = E[j,i]
                          Eik = E[k,i]
                          Ejk = E[k,j]
                          # Move along in the list of triplets with corrections
                          if nowInd < correctionsLength
                                  nowInd +=1
                                  nextKey = now_corrections[nowInd]
                          end
                        else
                          Eij = E[j,i]
                          Eik = E[k,i]
                          Ejk = E[k,j]
                        end

                        b = Dik - Djk + Dij
                        mu = (-Eij + Ejk - Eik - b)

                        if mu > epsi
                          denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                          E[j,i] = Eij + mu*(Wik*Wjk/denom)
                          E[k,i] = Eik + mu*(Wij*Wjk/denom)
                          E[k,j] = Ejk - mu*(Wik*Wij/denom)

                          # Next time we see this triple we have to correct
                          push!(next_corrections,(ijkKey,mu))
                          BtDual += mu*b*(Wij*Wik*Wjk/denom)
                        end
                        ### Done checking triangle j,k,i


                end
            end
        end
    end
    nowInd, BtDual
end


function Parallel_Dykstra_lamCC_TFA(A,GapTol::Float64=1e-3,ConTol::Float64=1e-3,
                lam::Float64=0.5,filename::String="DykstraLamCCoutput",gam::Float64=10.0,maxits::Int64=1000,statusFrequency::Int64=20)

        D,W = LamCC_DandW(A,lam)
        Parallel_Dykstra_CC_TFA(A,W,D,GapTol,ConTol,filename,gam,maxits,statusFrequency)

end


function Parallel_Dykstra_CC_TFA(A,W::Matrix{Float64},D::Matrix{Float64},GapTol::Float64=1e-5,ConTol::Float64=1e-5,
                        filename::String="DykstraCCoutput",gam::Float64=10.0,maxits::Int64=1000,statusFrequency::Int64=20)
        n = size(A,1)

        nthds::Int = Threads.nthreads()
        open(filename, "w") do f
                write(f, "Output from DykstraCC \n")
                write(f, "Gamma = $gam, Primal/Dual Gap tolerance = $GapTol, Constraint Tolerance = $ConTol \n")
        end

        # E = "error", eij = (xij - dij)
        # F = |E| i.e. fij = |eij| at optimum
        # D is the "desired distance" or just 1-Aij (anti-adjacency)
        # W is the weights matrix for LambdaCC, 1-lam or lam for each entry
        E = zeros(n,n)
        F = -gam*ones(n,n)

        # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)
        Q = zeros(n,n)

        # Allocate space for the output
        LPobjs = Vector{Float64}()
        duals = Vector{Float64}()
        primals = Vector{Float64}()
        gaps = Vector{Float64}()
        ConViolation = Vector{Float64}()

        # Correction term vector for triangle constraints

        current_corrections = Vector{Vector{Tuple{Int64, Float64}}}()
        next_corrections = Vector{Vector{Tuple{Int64, Float64}}}()

        current_corrections = Vector{Vector{Tuple{Int64, Float64}}}()
        for i in 1:nthds
            Ac = Vector{Tuple{Int64, Float64}}()
            push!(Ac, (0, 0.0))
            push!(current_corrections, Ac)
        end


        next_corrections = Vector{Vector{Tuple{Int64, Float64}}}()
        for i in 1:nthds
            Ac = Vector{Tuple{Int64, Float64}}()
            push!(next_corrections, Ac)
        end

        # First time through constraints
        parallel_cyclic_triple_loopE!(D,E,W,current_corrections,next_corrections)
        cyclic_double_loop!(E,F,P,Q,W)

        iter = 0

        # Keep track of worst constraint violations
        lastConCheck = 1.0
        lastTriCheck = 1.0

        # This keeps track of the dual variables times the right hand side.
        # It is used for computing the dual objective function.
        BtDual = 0.0

        # Initialize outside while loop
        X = zeros(n,n)
        FinalCon = 0.0
        Finalits = 0.0
        Finalobj = 0.0
        FinalGap = 0.0
        R = 0.0         # Ratio between the LP part and 2-norm penalty in objective
        Bty = 0.0

        while true
                iter += 1

                # Update which corretions to make in this iteration
                current_corrections = next_corrections
                next_corrections = Vector{Vector{Tuple{Int64, Float64}}}()

                for i in 1:nthds
                    Ac = Vector{Tuple{Int64, Float64}}()
                    push!(next_corrections, Ac)
                end

                for i in 1:nthds
                    Ac = current_corrections[i]
                    push!(Ac, (0, 0.0))
                end

                start = time()
                BtDual = parallel_cyclic_triple_loopE!(D,E,W,current_corrections,next_corrections)
                TriTime = time()-start

                cyclic_double_loop!(E,F,P,Q,W)

                num_corrections = 0
                for i in 1:length(current_corrections)
                    num_corrections += length(current_corrections[i])
                end
                # Check convergence and give a progress update
                ConCheck, lastConCheck, lastTriCheck, objective, LPobjs, primals, duals, gaps, ConViolation,
                gap, roundconverge, roundR, roundedGap, R, Bty =
                report_progress_CCDykstra(A,D,E,F,W,num_corrections,filename,iter,ConTol,GapTol,
                        lastConCheck,lastTriCheck,BtDual,TriTime,gam,LPobjs,
                        primals,duals,gaps,ConViolation,statusFrequency)

                # Stop if a desired constraint tolerance AND relative
                # primal/dual gap is small enough
                if abs(gap) < GapTol && ConCheck
                  println("Converged to within tolerance after $iter iterations")
                  X = E+D
                  Finalobj = LPobj(abs.(X-D),W)
                  Finalits = iter
                  FinalGap = gap
                  FinalCon = FullTriangleCheck(X)
                  break
                end

                # If rounding the current iterate produces a solution, stop and
                # return that rounded solution
                if roundconverge
                  println("\t Converged to within desired tolerance by rounding the solution")
                  r = roundR
                  X = round.(E+D,digits = r)
                  PrimalRound = Wnorm(W,X-D,abs.(X-D))/(2*gam) + LPobj(abs.(X-D),W)
                  FinalCon = FullTriangleCheck(X)
                  FinalGap = roundedGap
                  Finalobj = LPobj(abs.(X-D),W)
                  break
                end

                if iter >= maxits
                  println("Reached maximum number of iterations")
                  X = E+D
                  Finalobj = LPobj(abs.(X-D),W)
                  FinalGap = gap
                  break
                end

        end

        Finalits = iter
        return X, LPobjs, primals, duals, gaps, ConViolation, FinalCon, Finalits, Finalobj, FinalGap, R, Bty, next_corrections
end

# Collect variables and report progress each iteration
function report_progress_CCDykstra(A,D::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},current_corrections::Int,
    filename::String,iter::Int64,ConTol::Float64,GapTol::Float64,lastConCheck::Float64,lastTriCheck::Float64,BtDual::Float64,triTime::Float64,
    gam::Float64,LPobjs::Vector{Float64},primals::Vector{Float64},duals::Vector{Float64},gaps::Vector{Float64},ConViolation::Vector{Float64},statusFrequency::Int64)

        n = size(A,1)
        # Check whether constraint tolerance for metric constraints is satisfied
        tricheck = TriangleCheck(D+E,ConTol)

        # Constraint check for non-metric constraints
        doublecheck = DoubleCheck(E,F,ConTol)

        # True or false--are constraints satisfied to within tolerance?
        ConCheck = tricheck*doublecheck

        # Keep an updated record of the LP objective corresponding to the X variables
        #objective = LPcc_obj(A,D+E,lam)
        objective = LPobj(abs.(E),W)

        push!(LPobjs,objective)

        # Compute the objective with the F vector, which in theory
        # at the optimum will be the same, but isn't exaclty during the process
        fobj = LPobj(F,W)

        # x'*W*x where x = [e, f]
        xWx = Wnorm(W,E,F)

        # Primal objective = 1/(2*gamma)*x'*W*x + c'*x; the c'*x is the fobj function
        PrimalQP = xWx/(2*gam) + fobj
        R = xWx/(2*gam*fobj)
        # Dual objective, -b'*y - 1/(2*gamma)x'*W*x
        DualQP = -BtDual/gam - xWx/(2*gam)
        Bty = -BtDual/gam
        gap = (PrimalQP-DualQP)/DualQP

        # Round solutions for readability in output
        gapround = round(gap,digits = 5)
        PriRound = round(PrimalQP,digits = 5)
        DuRound = round(DualQP,digits = 5)
        ob = round(objective,digits = 3)
        time = round(triTime,digits = 3)

        # Check how dense the dual variable vector is--how many of (n choose 3)
        # triangle constraints have a nonzero correction term?
        nnz = current_corrections
        ysparsity = round(nnz*6/(n*(n-1)*(n-2)),digits = 5)

        # Do a full constraint check every so often, which takes longer
        # than just answering true or false whether constraints are satisfied
        roundconverge = false
        specialPrint = false
        PrimalRound = 1.0
        roundTri = 1.0
        roundedGap = 1.0
        GapRoundTol = 5e-1
        ConRoundTol = 5e-1
        roundR = 0
        roundedGap = 0
        if iter%statusFrequency == statusFrequency-1
          lastConCheck, lastTriCheck = FullConstraintCheck(D+E,E,F)

          # If we are close to optimal, try rounding the solution vector
          # and checking feasibility and the new gap. Stop early if this
          # rounded solution is good enough.

          if abs(gap) < GapRoundTol && lastConCheck < ConRoundTol
                # Round to a few decimals. If any one is good enough, stop
                for r = 5:-1:1
                  Xr = round.(E+D,digits = r)
                  xWxr = Wnorm(W,Xr-D,abs.(Xr-D))
                  objr = LPobj(abs.(Xr-D),W)
                  PrimalRound = xWxr/(2*gam) + objr
                  R = xWxr/(2*gam*objr)
                  roundTri = FullTriangleCheck(Xr)
                  #@show roundTri
                  roundedGap = (PrimalRound-DualQP)/DualQP
                  if roundTri < ConTol && abs(roundedGap) < GapTol
                          E = Xr - D
                          roundconverge = true
                          roundR = r  # save the
                          break
                  end
                end
                specialPrint = true
          end
        end

        # Save current progress
        push!(primals,PrimalQP)
        push!(duals,DualQP)
        push!(gaps,gap)
        push!(ConViolation,lastConCheck)


        println("$iter:\t Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio: $lastConCheck, TriVio = $lastTriCheck, 3Loop: $time, LPobj: $ob")
        if specialPrint
                PR = round(PrimalRound,digits = 5); rgap = round(roundedGap,digits = 5); rTri = round(roundTri,digits = 5)
                println("\t Rounded: \tNewPrimal = $PR, gap = $rgap, ConVio = $rTri")
        end

        open(filename, "a") do f
          write(f,"$iter:\t Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio: $lastConCheck, TriVio = $lastTriCheck, 3Loop: $time, LPobj: $ob, nnz(y) = $nnz, sparsity = $ysparsity\n")
        end

        return ConCheck, lastConCheck, lastTriCheck, objective, LPobjs, primals, duals, gaps, ConViolation, gap, roundconverge, roundR, roundedGap, R, Bty
end


# Visits double loop constraints for the correlation clustering LP
function cyclic_double_loop!(E::Matrix{Float64},F::Matrix{Float64},
    P::Matrix{Float64},Q::Matrix{Float64},W::Matrix{Float64})

 n = size(P,1)
 @inbounds for i = 1:n-1
  for j = i+1:n

    #  -E - F <= 0
    # corrections
    cor = P[j,i]
    Eij = E[j,i] - cor
    Fij = F[j,i] - cor
    delta = -Eij - Fij

    if delta > 0.0
      E[j,i] = Eij + delta/2
      F[j,i] = Fij + delta/2
      P[j,i] = delta/2
    else
      P[j,i] = 0.0
      E[j,i] = Eij
      F[j,i] = Fij
    end
   end
 end

 for i = 1:n-1
        for j = i+1:n
    #  E - F <= 0
    # corrections
    cor = Q[j,i]
    Eij = E[j,i] + cor
    Fij = F[j,i] - cor
    delta = Eij - Fij
    if delta > 0.0
      E[j,i] = Eij - delta/2
      F[j,i] = Fij + delta/2
      Q[j,i] = delta/2
    else
      Q[j,i] = 0.0
      E[j,i] = Eij
      F[j,i] = Fij
    end

  end
 end

end


## DykstraSC code

#
# Code for quadratic programming regulaziation of Leighton-Rao LP relaxation
# for sparsest cut.
#

# DYKSTRA_SC
#
# Dykstra-based projection method for solving a quadratic program which is the
# regularized version of the Leighton-Rao Linear Programming relaxation for
# sparsest cut.
#
# Paramters:
#
# A = adjacency matrix of an undirected, unweighted graph
# GapTol = the desired relative duality gap tolerance we wish to achieve for convergence
# ConTol = tolerance for constraint violations
# lam, gam = paramters controlling the relationship between the original LP
#               and the quadratic program which is solved in practice here.
# statusFreqcuency = controls how often we perform a full convergence check,
#                       which involves a full check for the maximum constraint violations
#                       and includes the "entrywise rounding step" (see paper)
# maxits = maximum number of iterations to run before terminating
# stagnationTol = if the QP objective score doesn't change by at least this much in one pass
#               through the constraints, terminate early. This isn't very useful and can be ignored.
#               It was only useful for catching bugs in early stages of the development of this code, and
#               shouldn't come into play unless perhaps you want to stop the code
#               early if Dykstra's method isn't making progress fast enough for it to
#               be worthwhile.
function ParallelDykstraSC(A::SparseMatrixCSC{Float64,Int64},GapTol::Float64=1e-3,ConTol::Float64=1e-5,
                lam::Float64=0.1,filename::String="DykstraLeightonRaoOutput",gam::Float64=10.0,
                maxits::Int64=1000,statusFrequency::Int64=5)

        n = size(A,1)
        nthds::Int = Threads.nthreads()
        open(filename, "w") do f
                write(f, "Parallel DykstraSC\n")
                write(f, "Lambda = $lam, gamma = $gam, tol = $GapTol, ConTol = $ConTol \n")
        end

        # Initialize X = -gam*A,
        # Wij = 1 if A_ij = 1
        # Wij = lam is A_ij = 0
        X,W = LeightonRaoQP_Initialize(A,lam,gam)

        # Matrices for double-loop correction terms (purposely dense)
        P = zeros(n,n)

        # Correction variable for constraint sum_ij xij = n
        SumCorrection = 0.0

        # Constant used in performing projections
        InvSumWeights = 0.0
        for i = 1:n-1
            for j = i+1:n
                InvSumWeights += 1/W[j,i]
            end
        end

        # Allocate space for the output
        LPobjs = Vector{Float64}()
        duals = Vector{Float64}()
        primals = Vector{Float64}()
        gaps = Vector{Float64}()
        ConViolation = Vector{Float64}()

        # Correction term vector (i.e. dual varialbes) for triangle constraints
        current_corrections = Vector{Vector{Tuple{Int64, Float64}}}()
        next_corrections = Vector{Vector{Tuple{Int64, Float64}}}()

        current_corrections = Vector{Vector{Tuple{Int64, Float64}}}()
        for i in 1:nthds
            Ac = Vector{Tuple{Int64, Float64}}()
            push!(Ac, (0, 0.0))
            push!(current_corrections, Ac)
        end


        next_corrections = Vector{Vector{Tuple{Int64, Float64}}}()
        for i in 1:nthds
            Ac = Vector{Tuple{Int64, Float64}}()
            push!(next_corrections, Ac)
        end

        # First time through constraints
        parallel_cyclic_triple_loop_X!(X,W,current_corrections,next_corrections)
        SumCorrection = the_sum_constraint!(X,W,SumCorrection,InvSumWeights)
        box_constraints!(X,W,P)

        iter = 0
        lastConCheck = 1.0

        # Make these variables global so they exist outside the while loop
        objective = 0.0
        R = 0.0
        FinalGap = 0.0
        FinalCon = 0.0
        Finalobj = 0.0
        Finalits = 0.0
        Bty = 0.0

        while true
                iter += 1

                # An empty vector to fill future correction terms in
                current_corrections = next_corrections
                next_corrections = Vector{Vector{Tuple{Int64, Float64}}}()

                for i in 1:nthds
                    Ac = Vector{Tuple{Int64, Float64}}()
                    push!(next_corrections, Ac)
                end

                for i in 1:nthds
                    Ac = current_corrections[i]
                    push!(Ac, (0, 0.0))
                end

                start = time()
                parallel_cyclic_triple_loop_X!(X,W,current_corrections,next_corrections)
                SumCorrection = the_sum_constraint!(X,W,SumCorrection,InvSumWeights)
                box_constraints!(X,W,P)
                TriTime = time() - start

                num_corrections = 0
                for i in 1:length(current_corrections)
                    num_corrections += length(current_corrections[i])
                end

                tricheck, lastConCheck, objective, LPobjs, primals, duals, gaps,
                ConViolation, gap, roundconverge, roundR, R, roundedGap, roundTri, Bty = report_progress_DykstraLR(A,X,W,P,num_corrections,filename,iter,
                ConTol,GapTol,lam,lastConCheck,TriTime,SumCorrection,gam,LPobjs,primals,duals,gaps,ConViolation,statusFrequency)

                # Return if the gap is less than tolerance and constraints
                # are satisfied to within a given tolerance
                if abs(gap) < GapTol && tricheck
                  open(filename,"a") do f
                     write(f, "Converged without rounding procedure.\n")
                  end
                  FinalCon = FullTriangleCheck(X)
                  Finalobj = LR_obj(A,X)
                  Finalits = iter
                  FinalGap = gap
                  break
                end

                # If rounding the current iterate produces a solution, stop and
                # return that rounded solution
                if roundconverge
                  println("\t Converged to within desired tolerance by rounding the solution to $roundR decimals")
                  open(filename,"a") do f
                      write(f, "\t Converged to within desired tolerance by rounding the solution to $roundR decimals\n")
                  end
                  Xr = round.(X,digits = roundR)
                  adjust = n/sum(tril(Xr))
                  Xr = adjust*Xr
                  X = Xr
                  FinalCon = FullTriangleCheck(X)
                  Finalobj = LR_obj(A,X)
                  Finalits = iter
                  FinalGap = roundedGap
                  break
                end

                if iter >= maxits
                  println("Maximum number of iterations reached")
                  open(filename,"a") do f
                      write(f, "Maximum number of iterations reached\n")
                  end
                  FinalCon = FullTriangleCheck(X)
                  Finalobj = LR_obj(A,X)
                  Finalits = iter
                  FinalGap = gap
                  break
                end
        end

        # Output final statistics and iterates
        return X, FinalCon, FinalGap, Finalobj, Finalits, R, LPobjs, duals, primals, gaps, ConViolation, Bty
end


function report_progress_DykstraLR(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},current_corrections::Int,
    filename::String,iter::Int64,ConTol::Float64,GapTol::Float64,lam::Float64,lastConCheck::Float64,triTime::Float64,SumCorrection::Float64,gam::Float64,
    LPobjs::Vector{Float64},primals::Vector{Float64},duals::Vector{Float64},gaps::Vector{Float64},ConViolation::Vector{Float64},statusFrequency::Int64)

        n = size(X,1)

        # True or false convergence check
        tricheck = TriangleCheck(X,ConTol)

        sumcheck = abs(n-sum(tril(X)))

        tricheck = tricheck*(sumcheck<ConTol)

        # nonzeros in dual vector
        nnzDelta = current_corrections

        # Compute primal and dual objectives
        objective = LR_obj(A,X)
        xWx = xWnorm(W,X)
        BtDual = (n*SumCorrection)
        PrimalQP = xWx/(2*gam) + objective
        R = xWx/(2*gam*objective)
        DualQP = -BtDual/(gam) - xWx/(2*gam)
        gap = (PrimalQP-DualQP)/DualQP

        # Every "statusFrequency" number of iterations, fully check constraints
        # to find the magnitude of the worst violation.
        # Also, try rounding the current solution if we are close to convergence
        roundconverge = false
        specialPrint = false
        PrimalRound = 1.0
        roundTri = 1.0
        roundedGap = 1.0
        GapRoundTol = 1e-1
        ConRoundTol = 1e-1
        roundR = 0
        #R = 0
        if iter%statusFrequency == statusFrequency-1
          lastConCheck = max(sumcheck,FullTriangleCheck(X))

          #tic()
          if abs(gap) < GapRoundTol && lastConCheck < ConRoundTol
                # Round to a few decimals. If any one is good enough, stop
                for r = 6:-1:2
                        Xr = round.(X,digits = r)
                        adjust = n/sum(tril(Xr))
                        Xr = adjust*Xr  # make sure sum of entries equals n
                        objr = LR_obj(A,Xr)
                        xWxr = xWnorm(W,Xr)
                        PrimalRound = xWxr/(2*gam) + objr
                        roundTri = FullTriangleCheck(Xr)
                        roundedGap = (PrimalRound-DualQP)/DualQP
                        R = xWxr/(2*gam*objr)
                        if roundTri < ConTol && abs(roundedGap) < GapTol
                                roundR = r
                                roundconverge = true
                                break
                        end
                end
                specialPrint = true
          end
          #timeRound = toc()
          #@show timeRound
        end

        # Save current progress
        push!(primals,PrimalQP)
        push!(duals,DualQP)
        push!(gaps,gap)
        push!(ConViolation,lastConCheck)
        push!(LPobjs,objective)

        # Round things to print out results nicely
        gapround = round(gap,digits = 5)
        PriRound = round(PrimalQP,digits = 5)
        DuRound = round(DualQP,digits = 5)
        ob = round(objective,digits = 3)
        time = round(triTime,digits = 3)
        tr = round(lastConCheck,digits = 5)
        Bty =  round(-BtDual/(gam),digits = 4)
        println("Iter $iter: Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio = $tr, 3Loop = $time, Obj: $ob")
        open(filename,"a") do f
            write(f, "Iter $iter: Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio = $tr, 3Loop = $time, Obj: $ob \n")
        end

        # Print something extra if you perform the entrywise rounding step
        if specialPrint
                PR = round(PrimalRound,digits = 5); rgap = round(roundedGap,digits = 5); rTri = round(roundTri,digits = 5)
                println("\t Rounded: \tNewPrimal = $PR, gap = $rgap, ConVio = $rTri")
                open(filename,"a") do f
                    write(f, "\t Rounded: \tNewPrimal = $PR, gap = $rgap, ConVio = $rTri \n")
                end
        end
        return tricheck, lastConCheck, objective, LPobjs, primals, duals, gaps, ConViolation, gap, roundconverge, roundR,R, roundedGap, roundTri,Bty
end

# Enforce the constraint \sum_ij x_ij = n
function the_sum_constraint!(X::Matrix{Float64},W::Matrix{Float64},SumCorrection::Float64,InvSumW::Float64)

    n = size(X,1)
    # Correction step
    sumX = 0
    for i = 1:n-1
        for j = i+1:n
            X[j,i] = X[j,i] + SumCorrection/W[j,i]
            sumX += X[j,i]
        end
    end

    constant = (sumX - n)/InvSumW
    for i = 1:n-1
        for j = i+1:n
            X[j,i] = X[j,i] - constant/W[j,i]
        end
    end
    return constant

end


# Enforce constraint X_ij >= 0
function box_constraints!(X::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64})

    n = size(X,1)
    for i = 1:n-1
        for j = i+1:n
            X[j,i] -= P[j,i]/W[j,i]

            thetaIplus = -X[j,i]*W[j,i]
            if thetaIplus > 0
                X[j,i] = 0.0
                P[j,i] = thetaIplus
            else
                P[j,i] = 0.0
            end

        end
    end
end

# Enforce triangle inequality constraints


function parallel_cyclic_triple_loop_X!(X::Matrix{Float64},W::Matrix{Float64},
      Bookshelf::Vector{Vector{Tuple{Int64, Float64}}}, newBookshelf::Vector{Vector{Tuple{Int64, Float64}}})
   #cyclic_triple_loop(D, E, W, now_corrections, next_corrections)
    n::Int = size(X,1)
    o::Int = 10
    nthds::Int = Threads.nthreads()
    epsi = 0
    # Construct the initial dissimilarity matrix, set X_0 = D essentially
    # We only edit the lower triangular part of the matrix. Inside of loops,
    # whenever possible we iterate through columns first, then rows.

    # This is an array of arrays of size nThreads, each thread gets it's own update array.


    #These arrays are used to store information about block size, e.g. index ranges in blocks of size 10
    A::Vector{Tuple{Int, Int}} = Vector{Tuple{Int, Int}}()
    B::Vector{Tuple{Int, Int}} = Vector{Tuple{Int, Int}}()
    for k in n:-o:1
        a = k-o+1
        b = k
        if a < 1
            a = 1
        end
        pushfirst!(B, (a,b))
    end

        A=B

    # Now perform Dykstra's method and plot intermediate steps
    # Each global iteration visits every constraint once



    indexing = Vector{Int}()
    for i in 1:nthds
        push!(indexing, 1)
    end




    begin
        r_endpoint::Int = round(Int64, ceil(n/o))
        for m::Int in r_endpoint:-1:1
            #@threads
            Threads.@threads for threadholder in 1:nthds
                c::Int64 = Threads.threadid()

                nowInd = indexing[c]
                old_triplet_corrections::Vector{Tuple{Int64, Float64}} = Bookshelf[c]
                new_triplet_corrections::Vector{Tuple{Int64, Float64}} = newBookshelf[c]

                numthreads::Int64 = Threads.nthreads()
                threadnum::Int64 = c
                start::Int64 = 1
                endpoint::Int64 = round(Int64, ceil(m/2))
                interval::Int = endpoint-start+1

                threadrange::Int = ceil(Int64, interval/numthreads)
                s::Int = start + threadrange*(threadnum-1)
                e::Int = min(start + threadrange*threadnum -1, endpoint)

                for iIter in s:e    #iIter::Int in 1:round(Int64, ceil(m/2))
                    kIter::Int = m-iIter+1
                    nowInd = iterate_X!(X, W, iIter, kIter, o, n, A, B, old_triplet_corrections, new_triplet_corrections, nowInd)
                end
                indexing[c] = nowInd
            end
        end
    end

    begin
        r_endpoint = round(Int64, ceil(n/o))
        for m::Int in 2:r_endpoint
            #Threads.@threads
            Threads.@threads for threadholder in 1:nthds
                c::Int64 = Threads.threadid()
                nowInd::Int64 = indexing[c]
                #new_triplet_corrections = nothing#::Array{Tuple, 1} = Tuple[]#Bookshelf[c]
                old_triplet_corrections::Vector{Tuple{Int64, Float64}} = Bookshelf[c]
                new_triplet_corrections::Vector{Tuple{Int64, Float64}} = newBookshelf[c]

                numthreads::Int = Threads.nthreads()
                threadnum::Int = c
                start::Int = m
                endpoint::Int = round(Int64,floor((m+r_endpoint)/2))
                interval::Int = endpoint-start+1

                threadrange::Int = ceil(Int64, interval/numthreads)
                s::Int = start + threadrange*(threadnum-1)
                e::Int = min(start + threadrange*threadnum -1, endpoint)

                for iIter::Int in s:e #m:round(Int64,floor((m+r_endpoint)/2)) #thread shplit
                    kIter::Int = r_endpoint-iIter+m
                    nowInd = iterate_X!(X, W, iIter, kIter, o, n, A, B, old_triplet_corrections, new_triplet_corrections, nowInd)
                end
                indexing[c] = nowInd
            end
        end
    end

end # function end



function iterate_X!(X::Matrix{Float64},W::Matrix{Float64}, iIter::Int, kIter::Int, o::Int, n::Int, A::Vector{Tuple{Int, Int}}, B::Vector{Tuple{Int, Int}}, now_corrections::Vector{Tuple{Int64, Float64}}, next_corrections::Vector{Tuple{Int64, Float64}}, nowInd::Int)
    epsi = 0
    iRange::Tuple{Int, Int} = A[iIter]
    kRange::Tuple{Int, Int} = B[kIter]

    i::Int = minimum(iRange)
    k::Int = maximum(kRange)

    C = Vector{Tuple{Int, Int}}()
    for j in k-1:-o:i+1
        a::Int = j-o+1
        b::Int = j
        if a < i+1
            a = i+1
        end
        push!(C, (a,b))
    end

    nextKey = now_corrections[nowInd]
    correctionsLength = length(now_corrections)

    for jIter::Int in 1:length(C)
        jRange::Tuple{Int, Int} = C[jIter]

        for i in iRange[1]:iRange[2]
            for j::Int in max(i+1, jRange[1]):jRange[2]


               Wij = W[j,i]

                for k in max(j+1,kRange[1]):kRange[2]
                Xik = X[k,i]
                  Xjk = X[k,j]
                  Xij = X[j,i]
                  Wik = W[k,i]
                  Wjk = W[k,j]

                  ### Check triangle i,j,k
                  ijkKey = (i-1)*n^2+(j-1)*n+k
                  # First see if this is the next triangle with a nonzero correction variable
                  if ijkKey == nextKey[1]

                          cor = nextKey[2]
                    # We need to scale the correction since we're projecting into the minimum
                    # weighted 2-norm direction

                    denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                    X[j,i] = Xij + cor*(Wik*Wjk/denom)
                    X[k,i] = Xik - cor*(Wij*Wjk/denom)
                    X[k,j] = Xjk - cor*(Wik*Wij/denom)
                    Xij = X[j,i]
                    Xik = X[k,i]
                    Xjk = X[k,j]
                    # Move along in the list of triplets with corrections
                    if nowInd < correctionsLength
                            nowInd +=1
                            nextKey = now_corrections[nowInd]
                    end
                  end

                  mu = (Xij - Xjk - Xik)

                  if mu > epsi
                    denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                    X[j,i] = Xij - mu*(Wik*Wjk/denom)
                    X[k,i] = Xik + mu*(Wij*Wjk/denom)
                    X[k,j] = Xjk + mu*(Wik*Wij/denom)
                    # Next time we see this triple we have to correct
                    push!(next_corrections,(ijkKey,mu))
                  end

                  ### Done checking triangle i,j,k

                  ### Check triangle i,k,j
                  Xij = X[j,i]
                  Xik = X[k,i]
                  Xjk = X[k,j]

                  ijkKey = (i-1)*n^2+(k-1)*n+j
                  # First see if this is the next triangle with a nonzero correction variable
                  if ijkKey == nextKey[1]

                          cor = nextKey[2]

                    denom = Wij*Wjk + Wik*Wij + Wjk*Wik
                    X[j,i] = X[j,i] - cor*(Wik*Wjk/denom)
                    X[k,i] = X[k,i] + cor*(Wij*Wjk/denom)
                    X[k,j] = X[k,j] - cor*(Wik*Wij/denom)
                    Xij = X[j,i]
                    Xik = X[k,i]
                    Xjk = X[k,j]
                    # Move along in the list of triplets with corrections
                    if nowInd < correctionsLength
                            nowInd +=1
                            nextKey = now_corrections[nowInd]
                    end
                  else
                    Xij = X[j,i]
                    Xik = X[k,i]
                    Xjk = X[k,j]
                  end
                  mu = (-Xij - Xjk + Xik)
                  if mu > epsi

                    denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                    X[j,i] = Xij + mu*(Wik*Wjk/denom)
                    X[k,i] = Xik - mu*(Wij*Wjk/denom)
                    X[k,j] = Xjk + mu*(Wik*Wij/denom)

                    # Next time we see this triple we have to correct
                    push!(next_corrections,(ijkKey,mu))
                  end
                  ### Done checking triangle i,k,j

                  ### Triangle j,k,i
                  Xij = X[j,i]
                  Xik = X[k,i]
                  Xjk = X[k,j]
                  ijkKey = (j-1)*n^2+(k-1)*n+i
                  # First see if this is the next triangle with a nonzero correction variable
                  if ijkKey == nextKey[1]

                          cor = nextKey[2]

                    denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                    X[j,i] = X[j,i] - cor*(Wik*Wjk/denom)
                    X[k,i] = X[k,i] - cor*(Wij*Wjk/denom)
                    X[k,j] = X[k,j] + cor*(Wik*Wij/denom)
                    Xij = X[j,i]
                    Xik = X[k,i]
                    Xjk = X[k,j]
                    # Move along in the list of triplets with corrections
                    if nowInd < correctionsLength
                            nowInd +=1
                            nextKey = now_corrections[nowInd]
                    end
                  else
                    Xij = X[j,i]
                    Xik = X[k,i]
                    Xjk = X[k,j]
                  end

                  mu = (-Xij + Xjk - Xik)

                  if mu > epsi
                    denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                    X[j,i] = Xij + mu*(Wik*Wjk/denom)
                    X[k,i] = Xik + mu*(Wij*Wjk/denom)
                    X[k,j] = Xjk - mu*(Wik*Wij/denom)

                    # Next time we see this triple we have to correct
                    push!(next_corrections,(ijkKey,mu))
                  end
                  ### Done checking triangle j,k,i

                end
            end
        end
    end
    nowInd
end
