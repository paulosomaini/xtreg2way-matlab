function varout=projvar(var,struc)
if any(isnan(var)|isinf(var))
    q=1:numel(var);
    disp(q(isnan(var)|isinf(var)));
    error('myApp:argChk','Check for NaN or Inf in observations listed above')
end


aux=sparse(struc.hhid,struc.tid,var.*struc.w,struc.N,struc.T,struc.obs);
Dy=full(sum(aux,2));
Ty=full(sum(aux,1))'; Ty=Ty([1:end-1]');
if struc.N<struc.T;
    delta=struc.A*Dy+struc.B*Ty;
    tau=[struc.B'*Dy+struc.C.invHH*Ty+struc.C.invHHDH*struc.A*(struc.C.invHHDH'*Ty);0];
else
    delta=struc.A.invDD*Dy+struc.A.invDDDH*struc.C*(struc.A.invDDDH'*Dy)+struc.B*Ty;
    tau=[struc.B'*Dy+struc.C*Ty;0];
end
    
varout=(var-delta(struc.hhid)-tau(struc.tid')).*sqrt(struc.w);


