

function cyclic_triple_loop!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},Bookshelf::Vector{Vector{Tuple{Int64, Float64}}}, newBookshelf::Vector{Vector{Tuple{Int64, Float64}}})
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
        unshift!(B, (a,b))
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
                    nowInd, BtDual = iterate(D, E, W, BtDual, iIter, kIter, o, n, A, B, old_triplet_corrections, new_triplet_corrections, nowInd)
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
                    nowInd, BtDual = iterate(D, E, W, BtDual, iIter, kIter, o, n, A, B, old_triplet_corrections, new_triplet_corrections, nowInd)
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


function iterate(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64}, BtDual::Float64, iIter::Int, kIter::Int, o::Int, n::Int, A::Vector{Tuple{Int, Int}}, B::Vector{Tuple{Int, Int}}, now_corrections::Vector{Tuple{Int64, Float64}}, next_corrections::Vector{Tuple{Int64, Float64}}, nowInd::Int)
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
                    #nowInd = ABCproj!(X,triplet_corrections,i,j,n,(i-1)*n^2+(j-1)*n + k,k,i,k,j, new_triplet_corrections, nowInd)

                    #  -X[j,i] + X[k,i] - X[k,j] <= 0
                    #nowInd = ABCproj!(X,triplet_corrections,i,k,n,(i-1)*n^2+(k-1)*n + j,j,i,k,j, new_triplet_corrections, nowInd)

                    #  -X[j,i] - X[k,i] + X[k,j] <= 0
                    #nowInd = ABCproj!(X,triplet_corrections,j,k,n,(j-1)*n^2+(k-1)*n + i,j,i,k,i, new_triplet_corrections, nowInd)
                    # @show nowInd
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




include("DykstraCC_Helper.jl")

function Cam_Dykstra_lamCC_TFA(A,GapTol::Float64=1e-3,ConTol::Float64=1e-3,
                lam::Float64=0.5,filename::String="DykstraLamCCoutput",gam::Float64=10.0,maxits::Int64=1000,statusFrequency::Int64=20,stagnationTol::Float64=1e-12)

        D,W = LamCC_DandW(A,lam)
        Dykstra_CC_TFA(A,W,D,GapTol,ConTol,filename,gam,maxits,statusFrequency,stagnationTol)

end



function Dykstra_CC_TFA(A,W::Matrix{Float64},D::Matrix{Float64},GapTol::Float64=1e-5,ConTol::Float64=1e-5,
                        filename::String="DykstraCCoutput",gam::Float64=10.0,maxits::Int64=1000,statusFrequency::Int64=20,stagnationTol::Float64=1e-12)
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
  #=
        println("\nPre iteration")
        @show length(current_corrections)
        @show length(current_corrections[1])

        @show length(next_corrections)
        @show length(next_corrections[1])
  =#
        # First time through constraints
        cyclic_triple_loop!(D,E,W,current_corrections,next_corrections)
        cyclic_double_loop!(E,F,P,Q,W)
  #=
        println("\nPost iteration")
        @show length(current_corrections)
        @show length(current_corrections[1])

        @show length(next_corrections)
        @show length(next_corrections[1])
        println("")
  =#

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
      #death
                next_corrections = Vector{Vector{Tuple{Int64, Float64}}}()

                for i in 1:nthds
                    Ac = Vector{Tuple{Int64, Float64}}()
                    push!(next_corrections, Ac)
                end

                for i in 1:nthds
                    Ac = current_corrections[i]
                    push!(Ac, (0, 0.0))
                end

                # Visiting the O(n^3) metric constraints is the bottleneck
        #=        println("\nPre iteration")
                @show length(current_corrections)
                @show length(current_corrections[1])

                @show length(next_corrections)
                @show length(next_corrections[1])
                =#

                tic()
                BtDual = cyclic_triple_loop!(D,E,W,current_corrections,next_corrections)
                TriTime = toq()
  #=
                println("\nPost iteration")
                @show length(current_corrections)
                @show length(current_corrections[1])

                @show length(next_corrections)
                @show length(next_corrections[1])
                println("")
  =#
                cyclic_double_loop!(E,F,P,Q,W)

                num_corrections = 0
                for i in 1:length(current_corrections)
                    num_corrections += length(current_corrections[i])
                end
                # Check convergence and give a progress update
                ConCheck, lastConCheck, lastTriCheck, objective, LPobjs, primals, duals, gaps, ConViolation,
                gap, stagnated, roundconverge, roundR, roundedGap, R, Bty =
                report_progress_CCDykstra(A,D,E,F,W,num_corrections,filename,iter,ConTol,GapTol,
                        lastConCheck,lastTriCheck,BtDual,TriTime,gam,LPobjs,
                        primals,duals,gaps,ConViolation,statusFrequency,stagnationTol)

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

                # If progress stagnates, return
                if stagnated
                        println("Progress stagnated, returned early")
                        X = E+D
                        Finalobj = LPobj(abs.(X-D),W)
                        FinalGap = gap
                        break
                end

                # If rounding the current iterate produces a solution, stop and
                # return that rounded solution
                if roundconverge
                        println("\t Converged to within desired tolerance by rounding the solution")
                        r = roundR
                        X = round.(E+D,r)
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




## This works for general (dense) correlation clustering problems


# Collect variables and report progress each iteration
function report_progress_CCDykstra(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},current_corrections::Int,
    filename::String,iter::Int64,ConTol::Float64,GapTol::Float64,lastConCheck::Float64,lastTriCheck::Float64,BtDual::Float64,triTime::Float64,
    gam::Float64,LPobjs::Vector{Float64},primals::Vector{Float64},duals::Vector{Float64},gaps::Vector{Float64},ConViolation::Vector{Float64},statusFrequency::Int64,stagnationTol::Float64)

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
        gapround = round(gap,5)
        PriRound = round(PrimalQP,5)
        DuRound = round(DualQP,5)
        ob = round(objective,3)
        time = round(triTime,3)

        if iter > 2
                stagnated = stagnationCheck(primals[end],PrimalQP,duals[end],DualQP,stagnationTol)
        else
                stagnated = false
        end
        # Check how dense the dual variable vector is--how many of (n choose 3)
        # triangle constraints have a nonzero correction term?
        nnz = current_corrections
        ysparsity = round(nnz*6/(n*(n-1)*(n-2)),5)

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

          #tic()
          if abs(gap) < GapRoundTol && lastConCheck < ConRoundTol
                # Round to a few decimals. If any one is good enough, stop
                for r = 5:-1:1
                        Xr = round.(E+D,r)
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
          #timeRound = toc()
        end

        # Save current progress
        push!(primals,PrimalQP)
        push!(duals,DualQP)
        push!(gaps,gap)
        push!(ConViolation,lastConCheck)


        println("$iter:\t Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio: $lastConCheck, TriVio = $lastTriCheck, 3Loop: $time, LPobj: $ob")
        if specialPrint
                PR = round(PrimalRound,5); rgap = round(roundedGap,5); rTri = round(roundTri,5)
                println("\t Rounded: \tNewPrimal = $PR, gap = $rgap, ConVio = $rTri")
        end

        open(filename, "a") do f
          write(f,"$iter:\t Dual = $DuRound, Primal = $PriRound, gap = $gapround, ConVio: $lastConCheck, TriVio = $lastTriCheck, 3Loop: $time, LPobj: $ob, nnz(y) = $nnz, sparsity = $ysparsity\n")
        end

        return ConCheck, lastConCheck, lastTriCheck, objective, LPobjs, primals, duals, gaps, ConViolation, gap, stagnated, roundconverge, roundR, roundedGap, R, Bty
end


# Visits double loop constraints for LambdaCC LP
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

# Visit each triangle constraint and perform necessary corrections and projections
# as a part of Dykstra's method



# This is the way that [Dax 2003} presents the Hildreth-like method,
# which he doesn't appropriately attribute to Hildreth.
function hildreth_triple_loop!(D::Matrix{Float64},E::Matrix{Float64},W::Matrix{Float64},
        now_corrections::Vector{Tuple{Int64,Float64}},next_corrections::Vector{Tuple{Int64,Float64}})
  n = size(D,1)
  correctionsLength = length(now_corrections)
  nowInd = 1
  # Grab next triplet in the list
  nextKey = now_corrections[nowInd]

  BtDual = 0.0
  @inbounds for i = 1:n-2
    for  j = i+1:n-1

     Dij = D[j,i]
     Wij = W[j,i]
    for  k = j+1:n

    Eij = E[j,i]
    Dik = D[k,i]
    Djk = D[k,j]
    Eik = E[k,i]
    Ejk = E[k,j]
    Wik = W[k,i]
    Wjk = W[k,j]
    denom = Wij*Wjk + Wik*Wij + Wjk*Wik
    aWinva = denom/(Wjk*Wik*Wij)
    ### Check triangle i,j,k

    ijkKey = (i-1)*n^2+(j-1)*n+k
    # First see if this is the next triangle with a nonzero correction variable

    if ijkKey == nextKey[1]
            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = Dik + Djk - Dij
    thetaI = (Eij - Ejk - Eik - b)/aWinva

    #delta = min(-thetaI,yI)
    delta = max(thetaI,-yI)

    if delta != 0
      E[j,i] = Eij - delta/Wij
      E[k,i] = Eik + delta/Wik
      E[k,j] = Ejk + delta/Wjk
    end

    if yI + delta != 0
        push!(next_corrections,(ijkKey,yI+delta))
        BtDual += b*(yI+delta)
    end

    ### Done checking triangle i,j,k

    Eij = E[j,i]
    Eik = E[k,i]
    Ejk = E[k,j]

    ### Check triangle i,k,j
    ijkKey = (i-1)*n^2+(k-1)*n+j
    if ijkKey == nextKey[1]

            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = -Dik + Djk + Dij
    thetaI = (-Eij - Ejk + Eik - b)/aWinva

    #delta = min(-thetaI,yI)
    delta = max(thetaI,-yI)

    if delta != 0
      E[j,i] = Eij + delta/Wij
      E[k,i] = Eik - delta/Wik
      E[k,j] = Ejk + delta/Wjk
    end

    # Update dual variable if it's not zero
    if yI + delta != 0
        push!(next_corrections,(ijkKey,yI+delta))
        BtDual += b*(yI+delta)
    end
    ### Done checking triangle i,k,j

    Eij = E[j,i]
    Eik = E[k,i]
    Ejk = E[k,j]
    ### Triangle j,k,i
    ijkKey = (j-1)*n^2+(k-1)*n+i
    if ijkKey == nextKey[1]

            yI = nextKey[2]
      # Move along in the list of triplets with corrections
      if nowInd < correctionsLength
              nowInd +=1
              nextKey = now_corrections[nowInd]
      end
    else
        yI = 0.0
    end

    b = Dik - Djk + Dij
    thetaI = (-Eij + Ejk - Eik - b)/aWinva

    #delta = min(-thetaI,yI)
    delta = max(thetaI,-yI)

    if delta != 0
      E[j,i] = Eij + delta/Wij
      E[k,i] = Eik + delta/Wik
      E[k,j] = Ejk - delta/Wjk
      # Next time we see this triple we have to correct
    end

    # Update dual variable if it's not zero
    if yI + delta != 0
        push!(next_corrections,(ijkKey,yI+delta))
        BtDual += b*(yI+delta)
    end
    ### Done checking triangle j,k,i

   end
  end
 end # end triple for loop

        return BtDual
end # end TripleLoop! function

function hildreth_double_loop!(E::Matrix{Float64},F::Matrix{Float64},
    P::Matrix{Float64},Q::Matrix{Float64},W::Matrix{Float64})

 n = size(P,1)
 @inbounds for i = 1:n-1
  for j = i+1:n

    #  -E - F <= 0
    thetaI = (-E[j,i] - F[j,i])*W[j,i]/2

    delta = min(-thetaI,P[j,i])
    E[j,i] = E[j,i] - delta/W[j,i]
    F[j,i] = F[j,i] - delta/W[j,i]
    P[j,i] = P[j,i] - delta

    #  E - F <= 0
    # corrections
   thetaI = (E[j,i] - F[j,i])*W[j,i]/2
   delta = min(-thetaI,Q[j,i])

   E[j,i] = E[j,i] + delta/W[j,i]
   F[j,i] = F[j,i] - delta/W[j,i]
   Q[j,i] = Q[j,i] - delta
  end
 end

end
