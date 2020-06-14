disp('Use the wrapper!')
%se==0
diary table.txt
tic; xtreg2way(y,[x1 x2],hhid,tid,w,[],0,1); toc;
diary off
%se==1 (robust)
diary table.txt
tic; xtreg2way(y,[x1 x2],hhid,tid,w,[],1,1); toc;
diary off

%one control
%se==0
diary table.txt
tic; xtreg2way(y,x1,hhid,tid,w,[],0,1); toc;
diary off
%se==1
diary table.txt
tic; xtreg2way(y,x1,hhid,tid,w,[],1,1); toc;
diary off
