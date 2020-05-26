function [varout,delta,tau]=projvar(var,struc)
if any(isnan(var)|isinf(var))
    q=1:numel(var);
    disp(q(isnan(var)|isinf(var)));
    error('myApp:argChk','Check for NaN or Inf in observations listed above')
end

if isfield('esample',struc) && numel(var) == numel(struc.esample)
    var = var(struc.esample);
end

aux=sparse(struc.hhid,struc.tid,var.*struc.w,struc.N,struc.T,struc.obs);
Dy=full(sum(aux,2));
Ty=full(sum(aux,1))'; Ty=Ty((1:end-1)');
if struc.N<struc.T
    delta=struc.A*Dy+struc.B*Ty;
    tau=[struc.B'*(Dy-struc.C.invHHDH'*Ty)+struc.C.invHH*Ty;0];
else
    delta=struc.A.invDD*Dy+struc.B*(Ty-struc.A.invDDDH'*Dy);
    tau=[struc.B'*Dy+struc.C*Ty;0];
end
    
varout=(var-delta(struc.hhid)-tau(struc.tid')).*sqrt(struc.w);


