%% Copyright 2018, SueYeon Chung, All rights reserved.
%% Compute MFT manifold capacity of the general manifold.
%% Efficient calcaulation of MFT capacity, R_M, D_M. 
%% Ref: Chung, Lee,Sompolinsky. "Classification and geometry of general perceptual manifolds." PRX 8.3 (2018): 031003.

function [output] = manifold_analysis_outputanchor(XtotT, options)
% XtotT: cell file containing P manifolds. Each cell is M_i by N matrix, where M_i is number of
% samples in i_{th} manifold, and N is the feature dimension. 
% WARNING: N(the feature dimension) should be larger than P! The function will not run OTHERWISE!
% kappa : margin. default should be 0 
% n_t : number of gaussian vectors vec(T) sampled per manifold
% flag_NbyM : default should be set to false. 
% If input is given in a transposed form (N by M_i instead of M_i by N) set to true. 
% DEPENDENCIES: IBM ILOG CPLEX Optimization Studio is needed for the cplexqp function. 
% Installation of IBM ILOG CPLEX for Windows: Free educational version available on https://developer.ibm.com/academic/accessresources/#onthehub
% Installation of IBM ILOG CPLEX for Linux and MacOS: http://www-stat.wharton.upenn.edu/~josezubi/INSTALL
% If cplexqp is not available, matlab native function quadprog can be used instead. Comment out accordingly.

% Example Input. 
% P=20 manifolds with N=30 feature dimension sampled from the normal distribution.
% Each manifold is sampled M_i = 100 times. 
% Gaussian vectors vec(T) are sampled n_t = 200 times per manifold. 
% for i=1:20
%	XtotT{i}=randn(100,30)
% end
% options.kappa =0; options.n_t=200; options.flag_NbyM = false; 
% [output] = manifold_analysis_corr(XtotT,options)

kappa = options.kappa; 
n_t=options.n_t; % Default: n_t=500. If number of samples per manifold M>500, then set at M. 
flag_NbyM = options.flag_NbyM;
P=length(XtotT); 
for ii=1:P
    if ~flag_NbyM 
        Xtot{ii}=XtotT{ii}';    
    else
        Xtot{ii}=XtotT{ii}; 
    end
    N=size(Xtot{ii},1); 
end
clear XtotT; 
Xori=[]; 
for pp=1:P
    M_vec(pp)=size(Xtot{pp},2); %% #(samples) vec. 
    Xori=[Xori Xtot{pp}];  %% All data.
end
M0=sum(M_vec); % Total number of samples. 
X0=mean(Xori,2); %Global Mean 
clear Xori; 
centers=nan(N,P);
for pp=1:P
    Xtot0{pp}=Xtot{pp}-repmat(X0, [1 M_vec(pp)]);   % Previously X 
    centers(:,pp)=mean(Xtot0{pp},2); 
end
clear Xtot;

%% Center correlation analysis 
for ii=1:P
    M=M_vec(ii); 
    Xr=squeeze(Xtot0{ii}); % Data for each manifold, global mean subtracted (center not normalized yet) 
    Xr0=mean(Xr,2);
    %Xrr=(Xr-repmat(Xr0,[1 M])); % Don't Normalize center 
    Xrr= (Xr-repmat(Xr0,[1 M]))/norm(Xr0); % Normalize center 
    XtotInput{ii}=Xrr; 
end
clear Xtot0; 
%% First make (D+1) dimensional data 
for ii=1:P
    S_r=XtotInput{ii}; 
    [D,m]= size(S_r);   % D-dimensional data, 

    %% Project data into the smaller space (not necessary for small dataset) 
    if D>m
        [Q,R]=qr(S_r,0); % Q is [D, m]
        S_rNew=Q'*S_r;  
        S_rOld=S_r; D_old=D; 
        S_r=S_rNew; 
        [D,m]= size(S_r);
        fprintf('Reduced: D=%d, m=%d.\n', D,m)
    end 
    sD=zeros(D,m); 
    for kk=1:(D+1)
        if kk<D+1
            sD(kk,:)=S_r(kk,:);
        else
            sc = 1; 
        end
    end
    %% Make data D+1 dimensional, adding center dimension 
    sD1_0 = [sD; repmat(sc, [1 m])]; % (D+1) by m 
    sD1 = sD1_0/sc; 
    [a_Mfull, a_M, R_M,D_M,s_all] = each_manifold_analysis_D1(sD1, kappa, n_t); 
    R_M_vec(ii)=R_M; 
    D_M_vec(ii)=D_M; 
%     a_M_vec(ii)=a_M; 
    a_Mfull_vec(ii)=a_Mfull; 
     anchor_vec{ii}=s_all; 

    fprintf('%d th manifold: D=%d, m=%d, D_M=%.2f, R_M=%.2f, a_M_full=%.2f.\n',ii, D,m, D_M, R_M, a_Mfull)
end
fprintf('Average of %d manifold: <D_M>=%.2f, <R_M>=%.2f, 1/<1/a_M>=%.2f.\n',P, ...
    mean(D_M_vec), mean(R_M_vec), 1./mean(1./a_Mfull_vec))
fprintf('STD of %d manifold: std(D_M)=%.2f, std(R_M)=%.2f, std(a_M)=%.2f.\n',P, std(D_M_vec), std(R_M_vec), std(a_Mfull_vec))
output.a_Mfull_vec=a_Mfull_vec; 
output.R_M_vec=R_M_vec; 
output.D_M_vec=D_M_vec;
output.anchorpoints = anchor_vec;
end

function [a_Mfull, a_M, R_M,D_M,s_all] = each_manifold_analysis_D1(sD1, kappa, n_t)
[D1, m] = size(sD1);   % D+1-dimensional data
D = D1 - 1; 
sc = 1; 
c_hat=zeros(1,D+1)'; c_hat(end,1)=1; 
t_vec = randn(D+1,n_t);
[ss, gg] = maxproj(t_vec, sD1, sc);
ss_max=ss; 

s_all=zeros(D+1,n_t);
v_f_all=zeros(D+1,n_t); 
for jj=1:n_t
    if gg(jj)+kappa < 0 % Interior Points 
        v_f=t_vec(:,jj); 
        s_f=ss(:,jj); 
    else
        eps0=1e-8;
          [v_f, alpha, vminustsq, exitflag]=compute_v_allpt(t_vec(:,jj),ss(:,jj),eps0, kappa, sD1); 

         t0=t_vec(:,jj)'*c_hat; 
         v0=v_f'*c_hat; 
         if norm(t_vec(:,jj)-v_f)<1e-8 % Interior. 
            v_f=t_vec(:,jj); 
            s_f=ss(:,jj);
         else
            lambda=sum(alpha); 
            l_vec(jj)=lambda; 
            s_f=(t_vec(:,jj)-v_f)/lambda; 

            vD=v_f(1:D,1);tD=t_vec(1:D,jj);sD=s_f(1:D,1);
            ov=tD'*sD/(norm(tD)*norm(sD));
            norm(s_f(1:D,1));
            maxvs=max(v_f'*sD1);
         end
    end
    s_all(:,jj)=s_f; %% anchor point coordinates, [N,n_t]
    v_f_all(:,jj)=v_f; 
    a_all(jj)=1./(norm(v_f-t_vec(:,jj)).^2); 
end

% s0 = mean(ss_max,2);
s0 = mean(s_all,2);
ds0 = s_all-repmat(s0,[1 n_t]); 
ds = ds0(1:end-1,:)./repmat(s_all(end,:),[D 1]); 
R_M = sqrt(mean(sum(ds(:, :).^2,1)));
MW_M = mean(abs(sum(t_vec(1:end-1,:)'.*ds(:, :)',2)));  
tD_vec=t_vec(1:end-1,:); 
sD_vec= s_all(1:end-1,:); 
t_hat_vec = tD_vec./repmat(sqrt(sum(tD_vec.^2,1)),[D 1]);
s_hat_vec = sD_vec./repmat(sqrt(sum(sD_vec.^2,1)),[D 1]);
D_M = D*(mean(sum(t_hat_vec.*s_hat_vec,1)))^2;
k_M = R_M.*sqrt(D_M);
a_M=alphaB(0, R_M, D_M); 
a_Mfull = 1/mean((max(sum(t_vec.*s_all,1)+kappa,0)).^2 ./ sum(s_all.^2,1)); 
end

function [ output] = alphaB(kappa, radius, d)
k=kappa; R=radius; 
fun = @(r, t) A(k, R, r, t).*exp(-t.^2./2)./sqrt(2*pi).*P_d(r, d); 

L=50; % L is approximating infinity, which is practically impossible; large L can cause instability, small L can be too far from infinity. 
alphainv = integral2(fun, 0, L, -L, L);
output= 1./alphainv;     
    function out_p = P_d(r, d)
        %% Original p_d
        out_p = 2.^(1-d/2).*r.^(d-1).*exp(-0.5*r.^2)./gamma(d/2);      
        %% Approximate p_d
        %out_p = exp(-0.5.*d.*log(2)+d.*log(r)-0.5*r.^2-0.5*d*log(d./2)+d./2); 
    end 
    function out_r = A(k, R, r, t)         
        out_r = (R.*r-t+k).^2./(R.^2+1).*((t-(k-r./R))>0).*(((k+R.*r)-t)>0) ...
                + ((t-k).^2+r.^2).*(((k-r./R)-t)>0); 
    end 
end

function [v_f, alpha,vminustsq, exitflag]= compute_v_allpt(tt, sDi, eps0, kappa, sD1)
% tt=[D+1,1]; S=[D+1, m]; eps=small. 

[D1, m] = size(sDi); D=D1-1; 
Tk=sD1; 
sc= sDi(end,1); 
flag_con=1;  
k=1; vminustsqk=10000;  
Fk_old=vminustsqk; 
[v_k, vt_k , exitflag, alphak, vminustsqk]= minimize_vt_sq(tt, Tk, kappa); 
    
v_f=v_k; 
vminustsq=vminustsqk;
alpha=alphak;
end

function [s0, gt] = maxproj(t_vec, sD1, sc)
    n_t=size(t_vec,2);
    D1= size(t_vec,1); D=D1-1; 
    m= size(sD1,2); 
    for i=1:n_t
        tt= t_vec(:,i);
        [gt0, imax]= max(tt(1:D,1)'*sD1(1:D,:));
        sr = sD1(1:D,imax); 
        s0(:,i)=[sr; sc];
        gt(i)=tt'*s0(:,i); 
    end
end

function [v_f, vt_f , exitflag, alphar, normvt2] = minimize_vt_sq(t, sD1, kappa)
% normvt is the objective function. 
    [D1, m]=size(sD1); D=D1-1; 
    H=eye(D1,D1);
    f=-t; 
    Aineq= sD1'; %[m by D1];
    bineq= -kappa*ones(m,1); 
    
    %% cplexqp options 
%     tolerance=1e-10; 
%     options = cplexoptimset('cplex');
%     options.simplex.tolerances.feasibility = tolerance;
%     options.simplex.tolerances.optimality = tolerance;
%     options.display = 'off';
%     options.MaxIter = 1e25;
%     options.qpmethod = 2; 
    %%
%    [v_f, vt_f,exitflag, output, alpha0] = cplexqp(H,f,Aineq,bineq,[],[],[],[],[],options);
    [v_f, vt_f,exitflag, output, alpha0] = quadprog(H,f,Aineq,bineq);
    normvt=(vt_f+0.5*sum(t.^2))*2;
    normvt2= sum((v_f-t).^2);
    alphar=alpha0.ineqlin; % Vector of lagrange coefficients for each data 
end 

