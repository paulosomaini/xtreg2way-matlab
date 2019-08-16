function struc=projdummies(hhid,tid,w)
if any(isnan(w)|isinf(w))
    q=1:numel(w);
    disp(q(isnan(w)|isinf(w)));
    error('myApp:argChk','Check for NaN or Inf in observations listed above')
end

struc.obs=numel(w); struc.w=w;
[~,~,hhid]=unique(hhid); [~,~,tid]=unique(tid);
struc.hhid=hhid; struc.tid=tid;
struc.N=max(hhid); struc.T=max(tid);
DH=sparse(struc.hhid,struc.tid,w,struc.N,struc.T,struc.obs);
DD=full(sum(DH,2));
HH=full(sum(DH,1)); HH=HH(1:end-1); DH=DH(:,1:end-1);
invHH=sparse(1:struc.T-1,1:struc.T-1,HH.^(-1),struc.T-1,struc.T-1,struc.T-1);
invDD=sparse(1:struc.N,1:struc.N,DD.^(-1),struc.N,struc.N,struc.N);
useinv=1;
if ~useinv;
    if struc.N<struc.T;
        struc.invA=diag(DD)-DH*invHH*DH';
        struc.B=-invA\DH*invHH;
        struc.C=invDD-invDD*DH/invA*DH'*invDD;
    else
        struc.invC=diag(HH)-DH'*invDD*DH;
        struc.A=invHH-invHH*DH'/invC*DH*invHH;
        struc.B=-A*DH*invHH;
    end
else
    if struc.N<struc.T;
        struc.A=inv(diag(DD)-DH*invHH*DH');
        struc.C.invHH=invHH;
        struc.C.invHHDH=invHH*DH';
        %struc.B=-struc.A*DH*invHH;
        struc.B=-struc.A*(struc.C.invHHDH)';
    else
        struc.C=inv(diag(HH)-DH'*invDD*DH);
        %invDD+invDD*DH*struc.C*DH'*invDD;
        struc.A.invDD=invDD;
        struc.A.invDDDH=invDD*DH;
        %Balt=-(struc.A.invDD*DH*invHH+struc.A.invDDDH*struc.C*(struc.A.invDDDH'*DH*invHH));
        %struc.B=-struc.A.invDD*DH*struc.C;
        struc.B=-struc.A.invDDDH*struc.C;
    end
        
end


    
