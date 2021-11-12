#cd("ADD Directory Here")

using Distributed

#Number of processors
@everywhere np = 1

if nprocs() < np
        addprocs(np-nprocs());
end

#Relaxation tightening level {0,1,2}
@everywhere RT=0 

#Number of pieces in piecewise relaxation {1,2,3,...}
@everywhere K=1

#Flag for GTEP co-optimization {0,1}
@everywhere G=1

#Optimality gap tolerance
@everywhere  gap_tol = 0.001

@everywhere begin
  using CSV, JuMP, CPLEX, NLsolve, LightGraphs, Ipopt, Statistics, SharedArrays, DataFrames
end

@everywhere include("load_data.jl");
@everywhere include("LBmodel.jl");
include("UBmodel.jl");
include("BT.jl");

t0 = time();
if G==1

        function Global(bounds, XL, XG)
                UB = 1e20
                LB = 0
                #Number of iterations
                niter = 1
                for iter = 1:niter
                        FBBT!(bounds)
                        t1 = time();
                        LB_trial,~, XL, XG, Qc, Ql = get_LB(bounds,RT,K,G)
                        println("LB time=  ", time()-t1,"\n")
                        LB = max(LB, LB_trial)
                        println("LB=    ",LB,"\n")
                        println("Qc="," ",Qc,"\n")
                        println("Ql="," ",Ql,"\n")
                        UB = min(UB, get_UB_AC_OPF_model(bounds, XL, G, XG, Qc, Ql))
                        gap=(UB-LB)/UB;
                        println("\n","UB=  ", UB,"\n")
                        println("gap= ",gap,"\n")
                        if gap < gap_tol
                                break;
                        end
                        if iter < niter
                                t2 = time();
                                OBBT!(bounds,UB,RT,K,G)
                                println("OBBT time  ", time()-t2)
                        end
                end
        end

        Global(bounds, XL, XG)

else

        function Global(bounds, XL)

                UB = 1e20
                LB = 0
                #Number of iterations
                niter = 1
                for iter = 1:niter
                        FBBT!(bounds)
                        t0 = time();
                        LB_trial,~, XL, Qc, Ql = get_LB(bounds,RT,K,G)
                        println("LB time=  ", time()-t0,"\n")
                        LB = max(LB, LB_trial)
                        println("LB=    ",LB,"\n")
                        println("Qc="," ",Qc,"\n")
                        println("Ql="," ",Ql,"\n")
                        UB = min(UB, get_UB_AC_OPF_model(bounds, XL, G, XG, Qc, Ql))
                        gap=(UB-LB)/UB;
                        println("\n","UB=  ", UB,"\n")
                        println("gap= ",gap,"\n")
                        if gap < gap_tol
                                break;
                        end
                        if iter < niter
                                t1 = time();
                                OBBT!(bounds,UB,RT,K,G)
                                println("OBBT time  ", time()-t1)
                        end
                end
        end

        Global(bounds, XL)

end

println("\n","Total runtime=  ", time()-t0);
