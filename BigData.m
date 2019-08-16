%%Unbalanced Panel, weights 1, heteroscedastic autocorrelation.
%Create data
clear
numgroups=4*10^6;
T=4*10^5;
obs=42*10^6;
probObs=obs/numgroups/T;

seed=14;
rand('seed',seed)
file='CrossCheckBalanced1000';

%Data to store
m=cell(numel(numgroups),1);
group=1;
%Each sample size
ng=numgroups; reps=1;
    e=rand(obs,1)*numgroups*T+1;
    e=sort(e);
    %isObs=rand(obs,1)>0.1;
    %e=1:obs;
    hhid=floor((e-1)/T+1);
    tid=floor(e-(hhid-1)*T);
    heff=randn(ng,1);
    teff=randn(T,1);
    w=rand(ng,1);w=w(hhid);
    x1=randn(obs,1)+0.5*heff(hhid)+0.25*teff(tid);
    x2=randn(obs,1)-0.25*heff(hhid)+0.5*teff(tid);
    autoc=rand(ng,1);
    initialv=randn(ng,1);
    
    betas=nan(reps,2);
    se_rob=nan(reps,2); % Robust to heteroskedasticity and serial correlation
    se_het=nan(reps,2); % Robust to heteroskedasticity
    se_cla=nan(reps,2); % Clasic assumptions
    tic;
    e=randn(obs,1); u=e;
    for o=1:obs
        if tid(o)>1, u_1=u(o-1); else u_1=initialv(hhid(o)); end
        u(o)=autoc(hhid(o))*u_1+e(o);
    end
    y=1+x1-x2+heff(hhid)+teff(tid)+u;
    data=[y x1 x2 hhid tid w];
    % Un balanced panel;
    data=data(isObs,:);
    fid = fopen([file '.csv'],'w');
    fprintf(fid, 'y,x1,x2,hhid,tid,w\n');
    fclose(fid);
    dlmwrite([file '.csv'], data,'-append', 'delimiter', ',', 'precision', 10);
   
%% Run the code;
%load the file;
data=dlmread([file '.csv'],',',1,0);
y=data(:,1);x1=data(:,2);x2=data(:,3);hhid=data(:,4);tid=data(:,5);w=data(:,6);
obs=numel(y); numgroups=numel(unique(hhid));
tic
    invSS=projdummies(hhid,tid,w);
    disp('Step 1'); toc;
    x1p=projvar(x1,invSS);
    x2p=projvar(x2,invSS);
    yp=projvar(y,invSS);
    disp('Step 2'); toc;
    reg=regress1(yp,[x1p x2p]);
    betas(1,:)=reg.beta';
    disp('Step 3'); toc;
    matCov=avar([x1p x2p],reg.res,invSS.hhid,reg.XX);
    stata_dof=((obs-1)/(obs-numel(reg.beta)-1))*(numgroups/(numgroups-1));
    se_rob(1,:)=sqrt(diag(matCov))*stata_dof;
    disp('Step 4a'); toc;
    matCov=avar([x1p x2p],reg.res,0,reg.XX);
    se_het(1,:)=sqrt(diag(matCov));
    disp('Step 4b'); toc;
    sig2hat=(reg.res'*reg.res)/(sum(invSS.w>0)-invSS.N-invSS.T+1-numel(reg.beta));
    se_cla(1,:)=sqrt(diag(sig2hat*inv(reg.XX)));
    disp('Step 4c'); toc;
    disp([betas(:)';se_cla(1,:);se_rob(1,:)])
toc;    

%% Use the wrapper instead of running each of steps separately:
disp('Use the wrapper!')
tic; [beta0,var0,yp,Xp,struc0] = xtreg2way(y,[x1 x2],hhid,tid,w,[],1); toc;
% now run a second regression taking advantage of the fact that the first
% step has been already performed (notice that the set of covariates in
% this regression is a subset of the previous set).
% This should return the same output as the previous regression.
tic; [beta1,var1] = xtreg2way(y,[x1 x2],hhid,tid,w,struc0,1); toc;
% Same here:
tic; [beta1,var1] = xtreg2way(yp,Xp,struc0); toc;
disp('*************************')
%Now use only a subset of covariates:
tic; [beta1,var1] = xtreg2way(y,[x1],hhid,tid,w,struc0,1); toc;
%Alternatively:
tic; [beta1,var1] = xtreg2way(yp,Xp(:,1),struc0); toc;



