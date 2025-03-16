% Basic equation of estimation
% z=Sig_ab*(Sig_bb^-1*z_obs)+noise

clear
close all
rng(10)

%% Generate synthetic manifold
syms u v w
% Manifold equation
% f(u,v,w)=3000*exp(-0.5*((u-50).^2/25^2+(v-50).^2/40^2))-5*u-20*w;
% f(u,v,w)=500*cos(2*pi*u/100)+250*sin(2*pi*v/50)-5*v-25*w;
% f(u,v,w)=11*u-5*v-25*w+sin(u.^2+v.^2+w.^2)*1e-9;
% f(u,v,w)=(-0.05*v+5).*u+(-0.05*u+5).*w+sin(u.^2+v.^2+w.^2)*1e-9;
f(u,v,w)=(500*cos(2*pi*(u+2*v-150)/160)+250*cos(4*pi*(u+v+w)/240)-5*v-20*w);
f_man=matlabFunction(f);

% Pole-to-plane equation (df/du, df/dv, df/dw)
L_u=matlabFunction(diff(f,u));
L_v=matlabFunction(diff(f,v));
L_w=matlabFunction(diff(f,w));

% Manifold gradient equation (dL/du, dL/dv, dL/dw)
f_u=@(u,v,w,bet) -bet*L_u(u,v,w)./...
    sqrt(L_u(u,v,w).^2+L_v(u,v,w).^2+L_w(u,v,w).^2);
f_v=@(u,v,w,bet) -bet*L_v(u,v,w)./...
    sqrt(L_u(u,v,w).^2+L_v(u,v,w).^2+L_w(u,v,w).^2);
f_w=@(u,v,w,bet) -bet*L_w(u,v,w)./...
    sqrt(L_u(u,v,w).^2+L_v(u,v,w).^2+L_w(u,v,w).^2);

%% Parameters
% Simulation type
con_unc=0;
if con_unc==0
    nd=50;
    MLengths=[50 30 100 5 20 40]';
end

% Domain
nx=100; % number of grids along x
ny=100; % number of grids along y
nz=100; % number of grids along z
nxyz=nx*ny*nz; % total number of grids
nl=6; % number of lithology
nc=4; % number of conditioning data
seed=1234;

% Estimation parameters
nsim=500;
bet=5; % ptp to manifold gradient coefficient

nabsc=5; % number of abscissae for numerical integration

%% Random path generation & reading conditioning data
rpath=randperm(nxyz)'; % random path of estimation

% Conditional or unconditional?
if con_unc==1
    % Reading conditioning data
    cdat=load('cdat.txt');
    nd=size(cdat,1);

    % Prioritized visiting of the location with conditioning information
    idd=cdat(:,1);
    jdd=cdat(:,2);
    kdd=cdat(:,3);
    vd=cdat(:,4);
    cond_loc=sub2ind([nx ny nz],idd,jdd,kdd);
    for ii=1:nd
        tmp1=find(rpath==cond_loc(ii));
        tmp2=rpath(ii);
        rpath(ii)=cond_loc(ii);
        rpath(tmp1)=tmp2;
    end
elseif con_unc==0
    vd=randi([1 nl],nd,1);
end

%
%
[id,jd,kd]=ind2sub([nx ny nz],rpath(1:nd+nsim));

% Initialize estimation vector
z_est=zeros(nxyz,1);

% Assigning conditional information to the estimation vector
d_loc=rpath(1:nd); % synthetic data location
%vd=randi([1 nl],nd,1);
z_est(d_loc)=vd(1:nd);

%% Synthetic TPM, Q, and p_mar generations
if con_unc==1
    TPM=load('TPM.txt');
    TPM=TPM;
    [V,~]=eigs(TPM');
    p_mar=V(:,1)'/sum(V(:,1));
    Q=logm(TPM);
elseif con_unc==0
    % [TPM,Q,p_mar]=TPM_gen1(MLengths);
    p_tgt=[0.2115 0.1304 0.3599 0.0330 0.0835 0.1817]';
    % [TPM,Q,p_mar]=TPM_gen1(MLengths);
    [TPM,Q,p_mar]=TPM_gen2(MLengths,p_tgt);
end

%% Hybrid sequential stochastic simulation
tic
rng(seed)
z_tmp=zeros(nxyz-nd,1);
% Sequential simulation
[ii0,jj0,kk0]=ind2sub([nx ny nz],rpath(nd+1:nd+nsim));
[ii1,jj1,kk1]=ind2sub([nx ny nz],rpath(1:nd));
x0=[ii0 jj0 kk0 f_man(ii0,jj0,kk0)];
x1=[ii1 jj1 kk1 f_man(ii1,jj1,kk1)];
dist_apr=pdist2(x0,x1);
dist_vec=dist_apr(1,:);
for ii=nd+1:nd+nsim
    % Locate current cell
    [i0,j0,k0]=ind2sub([nx ny nz],rpath(ii));

    [~,ord_apr]=sort(dist_vec');
    % Compute geodesic distances from the current to the conditioning cells
    dist=geodist([i0 j0 k0],[id(ord_apr(1:10)) jd(ord_apr(1:10)) kd(ord_apr(1:10))],f_man,nabsc,bet);

    % Find the nearer conditioning cells
    [dist,ord]=sort(dist);

    % Calculate the proportion of each lithology at the current computation
    q=sum(vd==1:nl)/nxyz;
    wght=p_mar./q;

    % Determine Lithology
    lith=lith_est(Q,vd(ord(1:nc)),dist(1:nc),0,wght); % check lith_est
    z_tmp(ii-nd)=lith;

    % Save the estimated lithology for the sequential estimation
    vd=[vd;lith];% vd(ord(1));
    [ii1,jj1,kk1]=ind2sub([nx ny nz],rpath(1:nd));
    x0=[ii0 jj0 kk0 f_man(ii0,jj0,kk0)];
    x1=[id(ii) jd(ii) kd(ii) f_man(id(ii),jd(ii),kd(ii))];
    dist_apr(:,ii)=pdist2(x0,x1);
    dist_vec=dist_apr(ii-nd,:);
end
toc

% Non-sequential estimation
parfor ii=nd+nsim+1:nxyz
    % Locate current cell
    [i0,j0,k0]=ind2sub([nx ny nz],rpath(ii));

    % Compute geodesic distances from the current to the conditioning cells
    dist=geodist([i0 j0 k0],[id jd kd],f_man,nabsc,bet);

    % Find the nearer conditioning cells
    [dist,ord]=sort(dist);

    % Determine Lithology
    lith=lith_est(Q,vd(ord(1)),dist(1));
    % lith=vd(ord(1));
    z_tmp(ii-nd)=lith;
end
% Assign stochastic estimates to the total estimation vector
for ii=nd+1:nxyz
    z_est(rpath(ii))=z_tmp(ii-nd);
end
toc

%% Figure
% Reshape 'dom' into 3D matrix for visualization
dom_3D=reshape(z_est,[ny, nx, nz]);

% Visualize slices of the 3D domain
figure('color','w','Position',[1000 700 1600 650])

subplot(1,2,1)
[X,Y,Z]=meshgrid(1:nx,ny/2+1:ny,1:nz);
slice(X,Y,Z,dom_3D(ny/2+1:ny,1:nx,1:nz),[1 nx],[ny/2+1 ny],[1 nz],'cubic');
hold on
[X,Y,Z]=meshgrid(nx/2+1:nx,1:ny,1:nz);
slice(X,Y,Z,dom_3D(1:ny,nx/2+1:nx,1:nz),[nx/2+1 nx],[1 ny],[1 nz],'cubic');
[X,Y,Z]=meshgrid(1:nx,1:ny,1:nz/2);
slice(X,Y,Z,dom_3D(1:ny,1:nx,1:nz/2),[1 nx],[1 ny],[1 nz/2],'cubic');

shading flat
axis equal
box on
set(gca,'fontsize',24,'fontname','times new roman','linewidth',2)
xlabel('u','fontsize',32,'fontweight','bold');
ylabel('v','fontsize',32,'fontweight','bold');
zlabel('w','fontsize',32,'fontweight','bold');
cb=colorbar;
set(cb,'fontsize',24,'fontname','times new roman','linewidth',2,'location','southoutside')
set(cb,'ytick',[1.4:(nl-1.8)/(nl-1):nl-0.4],'yticklabel',[1:nl])
camlight(-20,20)
set(gca,'cameraposition',[-479.7990 -312.8211  793.7013])
colormap(pink)
set(gca,'dataaspectratio',[1 1 2])

subplot(1,2,2)
[X,Y,Z]=meshgrid(1:nx,1:ny,1:nz/2);
slice(X,Y,Z,dom_3D(1:ny,1:nx,1:nz/2),[1 nx],[1 ny],[1 nz/2],'cubic');

shading flat
axis equal
box on
set(gca,'fontsize',24,'fontname','times new roman','linewidth',2)
xlabel('u','fontsize',32,'fontweight','bold');
ylabel('v','fontsize',32,'fontweight','bold');
zlabel('w','fontsize',32,'fontweight','bold');
cb=colorbar;
set(cb,'fontsize',24,'fontname','times new roman','linewidth',2,'location','southoutside')
set(cb,'ytick',[1.4:(nl-1.8)/(nl-1):nl-0.4],'yticklabel',[1:nl])
camlight(-20,20)
set(gca,'cameraposition',[-479.7990 -312.8211  793.7013])
colormap(pink)
[xx,yy,zz]=meshgrid(10:10:nx-10,10:10:ny-10,10:10:nz-10);
dLdu=f_u(xx,yy,zz,1);
dLdv=f_v(xx,yy,zz,1);
dLdw=f_w(xx,yy,zz,1);
hold on
% Pole-to-plane drawing
quiver3(xx(:,:,5),yy(:,:,5),ones(size(zz(:,:,5)))*nz/2,...
    dLdu(:,:,5),dLdv(:,:,5),dLdw(:,:,5),...
    'color',[0.75 0.75 0.75],'LineWidth', 1, 'AutoScaleFactor', 1);
axis([1 nx 1 ny 1 nz])
set(gca,'dataaspectratio',[1 1 2])
colormap(lines(nl))
%% Function for geodesic distance
function dist=geodist(uvw0,uvw_cond,f_man,n_absc,beta)
deluvw=1e-6;

u0=uvw0(1);
v0=uvw0(2);
w0=uvw0(3);
u=uvw_cond(:,1);
v=uvw_cond(:,2);
w=uvw_cond(:,3);
if beta~=0
n_cond=size(u,1);

uvw_lam=@(uvw1,uvw0,lambda) uvw0+(uvw1-uvw0)*lambda';

% n_absc % number of segments (abscisae)
% Legendre-Gauss Quadrature Weights and Nodes
[absc,wght]=lgwt(n_absc,-1,1);
lam=0.5*(absc+1);
ult=uvw_lam(u0,u,lam);% u(lambda)
vlt=uvw_lam(v0,v,lam);% v(lambda)
wlt=uvw_lam(w0,w,lam);% v(lambda)

% interpolation u_lambda and v_lambda
% Strike direction (=level-set curve direction)
df_du=-(f_man(ult(:)+deluvw,vlt(:),wlt(:))-f_man(ult(:),vlt(:),wlt(:)))/deluvw;
df_dv=-(f_man(ult(:),vlt(:)+deluvw,wlt(:))-f_man(ult(:),vlt(:),wlt(:)))/deluvw;
df_dw=-(f_man(ult(:),vlt(:),wlt(:)+deluvw)-f_man(ult(:),vlt(:),wlt(:)))/deluvw;
df_norm=sqrt(df_du.^2+df_dv.^2+df_dw.^2);
df_du=beta*df_du./df_norm;
df_dv=beta*df_dv./df_norm;
df_dw=beta*df_dw./df_norm;
% Computing distance of each segment and sum
g11=1+reshape(df_du,n_cond,n_absc).^2;
g12=reshape(df_du,n_cond,n_absc).*reshape(df_dv,n_cond,n_absc);
g13=reshape(df_du,n_cond,n_absc).*reshape(df_dw,n_cond,n_absc);
g22=1+reshape(df_dv,n_cond,n_absc).^2;
g23=reshape(df_dv,n_cond,n_absc).*reshape(df_dw,n_cond,n_absc);
g33=1+reshape(df_dw,n_cond,n_absc).^2;
% g=[g11 g12 g13;g12 g22 g23;g13 g23 g33];
% Integration
dist=sqrt((ult-u0).^2.*g11+(vlt-v0).^2.*g22+(wlt-w0).^2.*g33+...
    2*(ult-u0).*(vlt-v0).*g12+...
    2*(ult-u0).*(wlt-w0).*g13+...
    2*(vlt-v0).*(wlt-w0).*g23)*wght;
else
    dist=sqrt((u-u0).^2+(v-v0).^2+(w-w0).^2);
end
end

%% Function for Legendre-Gauss abscissae and weights
function [x,w]=lgwt(N,a,b)
N=N-1;
N1=N+1; N2=N+2;
xu=linspace(-1,1,N1)';
% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);
% Derivative of LGVM
Lp=zeros(N1,N2);
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=2;
% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    L(:,1)=1;
    Lp(:,1)=0;
    L(:,2)=y;
    Lp(:,2)=1;
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
    y0=y;
    y=y0-L(:,N2)./Lp;
end
% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;
% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end

%% Compute lithology (stochastic)
function lith_sim=lith_est(Q,lith,dist,flg,w)
if nargin==3 || flg==1
    flg=1;
    w=1;
end
n_lith=size(Q,1);
n_cnd=size(lith,1);
Lkh=ones(1,n_lith);
for ii=1:n_cnd
    T_ev=expm(Q*dist(ii));
    Lkh=Lkh.*(T_ev(lith(ii),:).*T_ev(:,lith(ii))');% Park, WRR, 2010
end
Lkh=Lkh.^(1/n_cnd);%%%%%%%%%%%%%%%
if flg==0
    Lkh=Lkh.*w; % Park, WRR, 20??
    p=cumsum(Lkh);
    p=p/p(end);
    pick=rand;
    lith_sim=find(p>=pick,1);
elseif flg==1
    [~,tmp]=sort(Lkh);
    lith_sim=tmp(end);
    % lith_sim=lith(1);
end
end


%% Generate transition probabilities

function [T,Q,p_mar]=TPM_gen1(L)

nl=size(L,1); % number of lithologies

% Off-diagonals of the transition count matrix
Tmp=randi([1 100],nl,nl);
TCM_off=full(spdiags(zeros(nl,1),0,Tmp));
% Embedded transition probability matrix (E)
E=TCM_off./sum(TCM_off,2);

% The pre-rate matrix (R)
R=E-eye(nl); % diagonals are set to -1

% Transition rate matrix (Q)
Q=R./L;

% Estimating the Transition Probability Matrix (T_est) using Q_est
T=expm(Q);

% Marginal probability vector (p_mar)
[V,~]=eigs(T');
p_mar=V(:,1)'/sum(V(:,1));
end


%% Generate transition probabilities constrained by marginal probabilities

function [T,Q,p_mar]=TPM_gen2(L,p_tgt)
% Validate input dimensions and normalize the target probabilities
nl=length(L);
p_tgt=p_tgt(:)'; % ensure row vector
p_tgt=p_tgt/sum(p_tgt); % normalize to sum to 1

% Prepare the initial guess for Q with proportional diagonal elements
Q0=zeros(nl); % rand(nl) 
Q0(~~eye(nl))=-1./L;

% Optimization options
options = optimoptions('fmincon',...
    'Algorithm','interior-point',...
    'Display','none',...
    'MaxIterations',10000,...
    'MaxFunctionEvaluations',10000,...
    'StepTolerance',1e-9, ...
    'OptimalityTolerance',1e-9, ...
    'ConstraintTolerance',1e-9);

% Flatten Q0 while excluding diagonal elements (to keep them fixed)
q0=Q0(~eye(nl));

% Define optimization problem
prob.options=options;
prob.solver='fmincon';
prob.objective=@(q) f_cost(q,p_tgt,Q0);
prob.x0=q0;
prob.lb=zeros(size(q0)); % Non-negative off-diagonal elements
prob.ub=[]; % No upper bounds
prob.nonlcon=@(q) T_con(q,Q0);

% Solve optimization
[q_opt,~]=fmincon(prob);

% Insert optimized off-diagonal elements back into Q
Q=Q0; % start with initial Q to retain diagonal values
Q(~eye(nl))=q_opt; % Update off-diagonal elements

% Generate T from Q and compute p_mar
T=expm(Q);
[V,~]=eigs(T');
p_mar=V(:,1)'/sum(V(:,1));


function cost=f_cost(q,p_tgt,Q0)
Q=Q0; % start with diagonal elements set
nl=size(Q0,1);
Q(~eye(nl))=q; % Update off-diagonal elements
T_tmp=expm(Q);
[V_tmp,~]=eigs(T_tmp');
p_tmp=V_tmp(:,1)'/sum(V_tmp(:,1));
cost=sum((p_tmp-p_tgt).^2); % cost function based on squared error
end

function [c,ceq]=T_con(q,Q0)
Q=Q0; % Start with diagonal elements set
nl=size(Q0,1);
Q(~eye(nl))=q; % Update off-diagonal elements
T_temp=expm(Q);
ceq=sum(T_temp, 2)-1; % Ensure row sums of T equal 1
c = []; % No inequality constraints
end
end