function get_UB_AC_OPF_model(bounds, XL, G, XG, Qc, Ql)


Theta_min = bounds.Theta_min
V_min = bounds.V_min
RE_min = bounds.RE_min
LE_min = bounds.LE_min
Theta_max = bounds.Theta_max
V_max = bounds.V_max
RE_max = bounds.RE_max
LE_max = bounds.LE_max



model=Model(Ipopt.Optimizer)
#set_optimizer_attribute(model, "print_level", 0)

@variable(model,-(ELcap*Branchdata[i,6]/baseMVA)<=PLe_S[i in 1:Num_Eline]<=(ELcap*Branchdata[i,6]/baseMVA))
@variable(model,-(ELcap*Branchdata[i,6]/baseMVA)<=QLe_S[i in 1:Num_Eline]<=(ELcap*Branchdata[i,6]/baseMVA))
@variable(model,-(ELcap*Branchdata[i,6]/baseMVA)<=PLe_R[i in 1:Num_Eline]<=(ELcap*Branchdata[i,6]/baseMVA))
@variable(model,-(ELcap*Branchdata[i,6]/baseMVA)<=QLe_R[i in 1:Num_Eline]<=(ELcap*Branchdata[i,6]/baseMVA))
@variable(model,-(Linedata_candidate[i,5]/baseMVA)<=PLc_S[i in 1:Num_Cline]<=(Linedata_candidate[i,5]/baseMVA))
@variable(model,-(Linedata_candidate[i,5]/baseMVA)<=QLc_S[i in 1:Num_Cline]<=(Linedata_candidate[i,5]/baseMVA))
@variable(model,-(Linedata_candidate[i,5]/baseMVA)<=PLc_R[i in 1:Num_Cline]<=(Linedata_candidate[i,5]/baseMVA))
@variable(model,-(Linedata_candidate[i,5]/baseMVA)<=QLc_R[i in 1:Num_Cline]<=(Linedata_candidate[i,5]/baseMVA))
@variable(model,0<=Pg[i in 1:Num_gen]<=Gendata[i,9], start=Gendata[i,9])
@variable(model,Gendata[i,5]<=Qg[i in 1:Num_gen]<=Gendata[i,4], start=Gendata[i,4])

if G==1
    @variable(model,Gendata_candidate[i,3]<=Pgc[i in 1:Num_Cgen]<=Gendata_candidate[i,2])
    @variable(model,Gendata_candidate[i,5]<=Qgc[i in 1:Num_Cgen]<=Gendata_candidate[i,4])
end

#original variables
@variable(model,V_min[i]<=Vi[i in 1:Num_bus]<=V_max[i], start=1)
@variable(model,-pi/2<=theta[i in 1:Num_bus]<=pi/2, start=0)
@variable(model,Theta_min[i]<=Theta[i in 1:Num_Eline+Num_Cline]<=Theta_max[i], start=0)

#auxiliary variables
@variable(model,RE_min[i]<=RE[i in 1:Num_Eline]<=RE_max[i], start=1)
@variable(model,LE_min[i]<=LE[i in 1:Num_Eline]<=LE_max[i], start=0)
@variable(model,-maximum(V_max)*maximum(V_max)<=RC[i in 1:Num_Cline]<=maximum(V_max)*maximum(V_max), start=1) 
@variable(model,-maximum(V_max)*maximum(V_max)<=LC[i in 1:Num_Cline]<=maximum(V_max)*maximum(V_max), start=0)  

@constraint(model,con_refbus,sum(theta[i] for i in 1:Num_bus if Busdata[i,2]==3)==0)
@constraint(model,con_Thetaexisitng[i=1:Num_Eline],Theta[i]-(theta[Branchdata[i,1]]-theta[Branchdata[i,2]])==0)
@constraint(model,con_Thetacandidate[i=1:Num_Cline],Theta[Num_Eline+i]-(theta[Linedata_candidate[i,1]]-theta[Linedata_candidate[i,2]])==0)

if G==1

    #Nodal power balance
    @constraint(model,con_nodalActivepowerbalance[i=1:Num_bus],-sum(Pg[j] for j in 1:Num_gen if Gendata[j,1]==i)-sum(Pgc[j] for j in 1:Num_Cgen if Gendata_candidate[j,1]==i)+
    baseMVA*sum(PLe_R[j] for j in 1:Num_Eline if Branchdata[j,2]==i)+baseMVA*sum(PLc_R[j] for j in 1:Num_Cline if Linedata_candidate[j,2]==i)+
    baseMVA*sum(PLe_S[j] for j in 1:Num_Eline if Branchdata[j,1]==i)+baseMVA*sum(PLc_S[j] for j in 1:Num_Cline if Linedata_candidate[j,1]==i)+
    sum(Loaddata[j,2] for j in 1:Num_load if Loaddata[j,1]==i)+Gsh[i]*Vi[i]^2==0)

    @constraint(model,con_nodalReactivepowerbalance[i=1:Num_bus],-sum(Qg[j] for j in 1:Num_gen if Gendata[j,1]==i)-sum(Qgc[j] for j in 1:Num_Cgen if Gendata_candidate[j,1]==i)+
    baseMVA*sum(QLe_R[j] for j in 1:Num_Eline if Branchdata[j,2]==i)+baseMVA*sum(QLc_R[j] for j in 1:Num_Cline if Linedata_candidate[j,2]==i)+
    baseMVA*sum(QLe_S[j] for j in 1:Num_Eline if Branchdata[j,1]==i)+baseMVA*sum(QLc_S[j] for j in 1:Num_Cline if Linedata_candidate[j,1]==i)+
    sum(Loaddata[j,3] for j in 1:Num_load if Loaddata[j,1]==i)-Bsh[i]*Vi[i]^2-Qc[i]+Ql[i]==0) 

    #Candidate generator limit
    @constraint(model,con_CgenActivelimit[i=1:Num_Cgen],Pgc[i]-Gendata_candidate[i,2]*XG[i]<=0)
    @constraint(model,con_CgenReactiveUPlimit[i=1:Num_Cgen],Qgc[i]-Gendata_candidate[i,4]*XG[i]<=0)
    @constraint(model,con_CgenReactiveDOWNlimit[i=1:Num_Cgen],Qgc[i]-Gendata_candidate[i,5]*XG[i]>=0)

else

    #Nodal power balance
    @constraint(model,con_nodalActivepowerbalance[i=1:Num_bus],-sum(Pg[j] for j in 1:Num_gen if Gendata[j,1]==i)+
    baseMVA*sum(PLe_R[j] for j in 1:Num_Eline if Branchdata[j,2]==i)+baseMVA*sum(PLc_R[j] for j in 1:Num_Cline if Linedata_candidate[j,2]==i)+
    baseMVA*sum(PLe_S[j] for j in 1:Num_Eline if Branchdata[j,1]==i)+baseMVA*sum(PLc_S[j] for j in 1:Num_Cline if Linedata_candidate[j,1]==i)+
    sum(Loaddata[j,2] for j in 1:Num_load if Loaddata[j,1]==i)+Gsh[i]*Vi[i]^2==0)
    
    @constraint(model,con_nodalReactivepowerbalance[i=1:Num_bus],-sum(Qg[j] for j in 1:Num_gen if Gendata[j,1]==i)+
    baseMVA*sum(QLe_R[j] for j in 1:Num_Eline if Branchdata[j,2]==i)+baseMVA*sum(QLc_R[j] for j in 1:Num_Cline if Linedata_candidate[j,2]==i)+
    baseMVA*sum(QLe_S[j] for j in 1:Num_Eline if Branchdata[j,1]==i)+baseMVA*sum(QLc_S[j] for j in 1:Num_Cline if Linedata_candidate[j,1]==i)+
    sum(Loaddata[j,3] for j in 1:Num_load if Loaddata[j,1]==i)-Bsh[i]*Vi[i]^2-Qc[i]+Ql[i]==0) 
end  


#Thermal constraint
@constraint(model,con_thermalExisitnglinelimitSend[i=1:Num_Eline],(PLe_S[i])^2+(QLe_S[i])^2<=(ELcap*Branchdata[i,6]/baseMVA)^2)
@constraint(model,con_thermalCandidatelinelimitSend[i=1:Num_Cline],(PLc_S[i])^2+(QLc_S[i])^2<=(Linedata_candidate[i,5]/baseMVA)^2)
@constraint(model,con_thermalExisitnglinelimitRec[i=1:Num_Eline],(PLe_R[i])^2+(QLe_R[i])^2<=(ELcap*Branchdata[i,6]/baseMVA)^2)
@constraint(model,con_thermalCandidatelinelimiRec[i=1:Num_Cline],(PLc_R[i])^2+(QLc_R[i])^2<=(Linedata_candidate[i,5]/baseMVA)^2)


#Exisitng line constraints
@constraint(model,con_ExisitngActivelineflowSend[i=1:Num_Eline],-PLe_S[i]+GE_ff[i]*Vi[Branchdata[i,1]]^2+GE_ft[i]*RE[i]+BE_ft[i]*LE[i]==0)
@constraint(model,con_ExisitngReactivelineflowSenc[i=1:Num_Eline],-QLe_S[i]-BE_ff[i]*Vi[Branchdata[i,1]]^2-BE_ft[i]*RE[i]+GE_ft[i]*LE[i]==0)
@constraint(model,con_ExisitngActivelineflowRec[i=1:Num_Eline],-PLe_R[i]+GE_tt[i]*Vi[Branchdata[i,2]]^2+GE_tf[i]*RE[i]-BE_tf[i]*LE[i]==0)
@constraint(model,con_ExisitngReactivelineflowRec[i=1:Num_Eline],-QLe_R[i]-BE_tt[i]*Vi[Branchdata[i,2]]^2-BE_tf[i]*RE[i]-GE_tf[i]*LE[i]==0)


#Candidate line constraints
for i=1:Num_Cline
    if XL[i]==0
    @constraint(model,[i],PLc_S[i]==0)
    @constraint(model,[i],PLc_R[i]==0)
   	@constraint(model,[i],QLc_S[i]==0)
   	@constraint(model,[i],QLc_R[i]==0)
    else
	@constraint(model,[i],-PLc_S[i]+GC_ff[i]*Vi[Linedata_candidate[i,1]]^2+GC_ft[i]*RC[i]+BC_ft[i]*LC[i]==0)
	@constraint(model,[i],-QLc_S[i]-BC_ff[i]*Vi[Linedata_candidate[i,1]]^2-BC_ft[i]*RC[i]+GC_ft[i]*LC[i]==0)
	@constraint(model,[i],-PLc_R[i]+GC_tt[i]*Vi[Linedata_candidate[i,2]]^2+GC_tf[i]*RC[i]-BC_tf[i]*LC[i]==0)
	@constraint(model,[i],-QLc_R[i]-BC_tt[i]*Vi[Linedata_candidate[i,2]]^2-BC_tf[i]*RC[i]-GC_tf[i]*LC[i]==0)
    end
end

#Definition for RE LE Ui
@NLconstraint(model,con_RE[i=1:Num_Eline],RE[i]-(Vi[Branchdata[i,1]]*Vi[Branchdata[i,2]]*cos(Theta[i]))==0)
@NLconstraint(model,con_LE[i=1:Num_Eline],LE[i]-(Vi[Branchdata[i,1]]*Vi[Branchdata[i,2]]*sin(Theta[i]))==0)
@NLconstraint(model,con_RC[i=1:Num_Cline],RC[i]-(Vi[Linedata_candidate[i,1]]*Vi[Linedata_candidate[i,2]]*cos(Theta[Num_Eline+i]))==0)
@NLconstraint(model,con_LC[i=1:Num_Cline],LC[i]-(Vi[Linedata_candidate[i,1]]*Vi[Linedata_candidate[i,2]]*sin(Theta[Num_Eline+i]))==0)

if G==1

    @objective(model,Min,sum(Pg[i] for i in 1:Num_gen)+sum(Pgc[i] for i in 1:Num_Cgen))

else

    @objective(model,Min,sum(Pg[i] for i in 1:Num_gen))

end

optimize!(model)

if termination_status(model) == MOI.OPTIMAL || termination_status(model) == MOI.LOCALLY_SOLVED
    if G==1
        return objective_value(model)+sum(Linedata_candidate[i,6]*XL[i] for i in 1:Num_Cline)+sum(Gendata_candidate[i,6]*XG[i] for i in 1:Num_Cgen)+sum(Qc[i] for i in 1:Num_bus)+sum(Ql[i] for i in 1:Num_bus)
    else
        return objective_value(model)+sum(Linedata_candidate[i,6]*XL[i] for i in 1:Num_Cline)+sum(Qc[i] for i in 1:Num_bus)+sum(Ql[i] for i in 1:Num_bus)
    end
else
    return 1e20
end

end