
% this file reads in the EGT best fit parameters
% computes the marginal statistics and correlation functions
% computes the spearmans correlation
% computes the underlying Gaussian correlation (using Lebrun Dutfoy)

clear ; clc; close all;

z = dlmread('results/EGT_BF_Best_Fit3D.txt',',');
XEG = z(:,2:7);
dtmp = dlmread('data/Small_IN100_20140317.csv',',',[2 0 2354 27]);
gorientation = dtmp(:,19:21);
Ngrain = length(z(:,1));
Centroid = z(:,2:4);
N1 = 189; N2 = 201; N3 = 117;
dx = .25;
v = z(:,5:7);
rotmat = z(:,8:16); 
% the rotation matrix is written as x1map, x2map, x3map where
% x1_rot = x1map(1)*x1 + x1map(2)*x2 + x1map(3)*x3

%I will convert the rotation matrix to axis angle, and then i will convert
%the axis to spherical coordinate system 
ellipse_orientation = zeros(Ngrain,3);

for j = 1:Ngrain
    % note that rotation matrix is given from data as transpose but
    % I will write it out in the same format as a transpose
    z = [rotmat(j,1) rotmat(j,2) rotmat(j,3); ...
        rotmat(j,4) rotmat(j,5) rotmat(j,6); ...
        rotmat(j,7) rotmat(j,8) rotmat(j,9)];
    z = real(vrrotmat2vec(z));
    ellipse_orientation(j,3) = z(4);
    [ellipse_orientation(j,1),ellipse_orientation(j,2),~]= ...
        cart2sph(z(1),z(2),z(3));

end

% ---------------------------------------------------------------------%
% get conditional PDF of grain velocities
    % this is conditioning on average distance of closest Np points
    Np = 5;
        Nbin = 5;
nn = zeros(Ngrain,1);
for j = 1:Ngrain
   z = setdiff(1:Ngrain,j);
   dsq =  (XEG(z,1)-XEG(j,1)).^2+ (XEG(z,2)-XEG(j,2)).^2+ (XEG(z,3)-XEG(j,3)).^2;
   z3 = sort(dsq);
   %z2 = z((dsq<= z3(10)));
   nn(j) = mean(z3(1:Np));
end
velocity_bins = cell(Nbin,1);
z = sort(nn);
z1 = floor(Ngrain/Nbin);
for j = 1:Nbin
   %velocity_bins{j} = nn(nn(:,2)<z(j*z1) & nn(:,2)>z((j-1)*z1+1),1 ); 
   velocity_bins{j} =   find(nn<z(j*z1) & nn>z((j-1)*z1+1)); 
end
grain_bin = zeros(Nbin,1);
for j = 1:Nbin
    grain_bin(j) = z(j*z1);
end

vcondY = cell(Nbin,1);
for j = 1:Nbin
    vcondY{j} = XEG(velocity_bins{j},4:6);
end


mean_and_var = zeros(Nbin,6);
for j = 1:Nbin
    for j1 = 1:3
    mean_and_var(j,2*j1-1:2*j1) = [mean(vcondY{j}(:,j1)) var(vcondY{j}(:,j1))];
    end
end
mean_and_var(:,[1 3 5])
mean_and_var(:,[2 4 6])


% ---------------------------------------------------------------------%
% look at covariance structure as a function of bins
param_mat = [v ellipse_orientation gorientation];
% get ranks for spearmans correlation
Nparam = length(param_mat(1,:));
mean_mat = zeros(Nbin,Nparam);
var_mat = zeros(Nbin,Nparam);
cov_mat = cell(Nbin,1);
cov_norm_mat = cell(Nbin,1);
for j = 1:Nbin
    for j1 = 1:Nparam
    mean_mat(j,j1) = mean(param_mat(velocity_bins{j},j1));
    var_mat(j,j1) = sum( (param_mat(velocity_bins{j},j1)-mean_mat(j,j1)).^2) ...
        /length(velocity_bins{j});
    end    
end
for j = 1:Nbin
cov_mat{j} = zeros(Nparam);
cov_norm_mat{j} = zeros(Nparam);
    for j1 = 1:Nparam
        for j2 = 1:Nparam
cov_mat{j}(j1,j2) = sum((param_mat(velocity_bins{j}(:),j1)-mean_mat(j,j1)).* ...
    (param_mat(velocity_bins{j}(:),j2)-mean_mat(j,j2)))/length(velocity_bins{j});      
cov_norm_mat{j}(j1,j2) = sum((param_mat(velocity_bins{j}(:),j1)-mean_mat(j,j1)).* ...
    (param_mat(velocity_bins{j}(:),j2)-mean_mat(j,j2)))/ ...
    sqrt(var_mat(j,j1)*var_mat(j,j2))/length(velocity_bins{j});      

        end
    end
    
end

% ---------------------------------------------------------------------%
% get Spearman's rho of each binned velocity distribution
param_mat_ranks = cell(Nbin,1);
for j = 1:Nbin
    Nv = length(velocity_bins{j});
param_mat_ranks{j} = zeros(Nv,Nparam);
for j1 = 1:Nparam
    param_mat_ranks{j}(:,j1) = tiedrank(param_mat(velocity_bins{j},j1));
end
end

ranks_mean = zeros(Nbin,Nparam);
ranks_var = zeros(Nbin,Nparam);
for j = 1:Nbin
for j1 = 1:Nparam
    ranks_mean(j,j1) = mean(param_mat_ranks{j}(:,j1));
end
for j1 = 1:Nparam
    ranks_var(j,j1) = sum( (param_mat_ranks{j}(:,j1)-ranks_mean(j,j1)).^2 );
end
end
ranks_rho = cell(Nbin,1);
for j = 1:Nbin
ranks_rho{j} = ones(Nparam);
for j1 = 2:Nparam
    for j2 = 1:(j1-1)
   ranks_rho{j}(j1,j2) = sum( (param_mat_ranks{j}(:,j1)-ranks_mean(j,j1)).* ...
       (param_mat_ranks{j}(:,j2)-ranks_mean(j,j2))) / ...
       sqrt(ranks_var(j,j1)*ranks_var(j,j2));
   ranks_rho{j}(j2,j1) = ranks_rho{j}(j1,j2);
    end
end
end

% compute underlying Gaussian correlation
Nataf_rho = cell(Nbin,1);
for j = 1:Nbin
Nataf_rho{j} = ones(Nparam);
for j1 = 2:Nparam
    for j2 = 1:(j1-1)
   Nataf_rho{j}(j1,j2) = 2*sin(pi/6*ranks_rho{j}(j1,j2));     
   Nataf_rho{j}(j2,j1) = Nataf_rho{j}(j1,j2);     
    end
end
end
R0inv = cell(Nbin,1); Gamma = cell(Nbin,1); Gammainv = cell(Nbin,1);
for j = 1:Nbin
R0inv{j} = Nataf_rho{j} \ eye(Nparam);
Gamma{j} = chol(R0inv{j});
Gammainv{j} = Gamma{j} \ eye(Nparam);
end

% ---------------------------------------------------------------------%
% get conditional densities of target

param_mat_CDF = cell(Nbin,1);
for j = 1:Nbin
    Nv = length(velocity_bins{j});
    param_mat_CDF{j} = zeros(1024,2*Nparam);
for j1 = 1:3
    dmin = 1;
dmax = max(param_mat(velocity_bins{j},j1))*1.2;
[~,~,param_mat_CDF{j}(:,2*j1-1),param_mat_CDF{j}(:,2*j1)]= ...
    kde(param_mat(velocity_bins{j},j1),1024,dmin,dmax);
end
for j1 = 4:Nparam
    dmin = min(param_mat(velocity_bins{j},j1));
    dmax = max(param_mat(velocity_bins{j},j1));
[~,~,param_mat_CDF{j}(:,2*j1-1),param_mat_CDF{j}(:,2*j1)]= ...
    kde(param_mat(velocity_bins{j},j1),1024,dmin,dmax);
end
end

% write out the conditional CDFs and sample space
fil = 'results/Target_CDFs.txt';
fid = fopen(fil,'w');
fprintf(fid,'5 conditional CDFs; [sample space CDF values] for 9 marks ');
fprintf(fid,'(v1,v2,v3,varphi1p,Phip,varphi2p,varphi1,Phi,varphi2p) ');
fprintf(fid,'first value is minimum of d_5 (average distance of closest 5 points)\n');
d5 = [0 z(z1) z(z1*2) z(z1*3) z(z1*4)];
for j = 1:5
    fprintf(fid,'%f \n',d5(j));
    for j1 = 1:1024
        for j2 = 1:18
            fprintf(fid,'%f,',param_mat_CDF{j}(j1,j2));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);

% write out the Gammainv for each conditional joint distribution
% it is written out as an upper matrix only
fil = 'results/Target_Gammainv.txt';
fid = fopen(fil,'w');
fprintf(fid,'Gamma inverse written out as an matrix of 9x9 for each ');
fprintf(fid,'conditional PDF written out successively. So its 81x5 values\n');
for j = 1:5
    for j1 = 1:9
        for j2 = 1:9
            fprintf(fid,'%f\n',Gammainv{j}(j2,j1));
        end
    end
end
fclose(fid);
% ---------------------------------------------------------------------%
% check to see if capturing the correlation and marginal CDFs

testq = 1; % 1 if skip this check, 0 if wanna check

if (testq == 0 )
% (1) simulate a lot of standard normals U, 
% (2) transfer to Y=Gammainv U, and then z=\Phi(Y)
% (3) then transfer to X (as F^-1 (z)
Nmc = 100000;
U = randn(Nmc,Nparam);
X = cell(Nbin,1);
for j = 1:Nbin
Y = U*(Gammainv{j}');
%Z = mvncdf(Y(1,:),zeros(1,Nparam),Nataf_rho);
%Z = mvncdf(Y,zeros(1,Nparam),Nataf_rho);
Z = normcdf(Y);
X{j} = zeros(Nmc,Nparam);
for j1 = 1:Nparam
    for j2 = 1:Nmc
        z2 = find(param_mat_CDF{j}(:,2*j1)<Z(j2,j1),1,'last');
        if (isempty(z2)==1) z2 = 1; end
    X{j}(j2,j1) = param_mat_CDF{j}(z2,2*j1-1);
    end
end
end

param_mat_ranks_sim = cell(Nbin,1);
ranks_mean_sim = cell(Nbin,1);
ranks_var_sim = cell(Nbin,1);
ranks_rho_sim = cell(Nbin,1);
mean_sim = cell(Nbin,1);
var_sim = cell(Nbin,1);
covvar_sim = cell(Nbin,1);
covvar_norm_sim = cell(Nbin,1);
for j = 1:Nbin    
% compute spearmans rho of simulated results
param_mat_ranks_sim{j} = zeros(Nmc,Nparam);
for j1 = 1:Nparam
    param_mat_ranks_sim{j}(:,j1) = tiedrank(X{j}(:,j1));
end

ranks_mean_sim{j} = zeros(Nparam,1);
ranks_var_sim{j} = zeros(Nparam,1);
mean_sim{j} = zeros(Nparam,1);
var_sim{j} = zeros(Nparam,1);
for j1 = 1:Nparam
    ranks_mean_sim{j}(j1) = mean(param_mat_ranks_sim{j}(:,j1));
    mean_sim{j}(j1) = mean(X{j}(:,j1));
    ranks_var_sim{j}(j1) = sum( (param_mat_ranks_sim{j}(:,j1)-ranks_mean_sim{j}(j1)).^2 )/Nmc;
    var_sim{j}(j1) = sum( (X{j}(:,j1)-mean_sim{j}(j1)).^2 )/Nmc;
end

ranks_rho_sim{j} = ones(Nparam);
covvar_norm_sim{j} = ones(Nparam);
covvar_sim{j} = zeros(Nparam);
covvar_sim{j}(1,1) = var_sim{j}(1);
for j1 = 2:Nparam
    covvar_sim{j}(j1,j1) = var_sim{j}(j1);
    for j2 = 1:(j1-1)
   ranks_rho_sim{j}(j1,j2) = sum( (param_mat_ranks_sim{j}(:,j1)-ranks_mean_sim{j}(j1)).* ...
       (param_mat_ranks_sim{j}(:,j2)-ranks_mean_sim{j}(j2))) / ...
       sqrt(ranks_var_sim{j}(j1)*ranks_var_sim{j}(j2))/Nmc;
   ranks_rho_sim{j}(j2,j1) = ranks_rho_sim{j}(j1,j2);
   covvar_norm_sim{j}(j1,j2) = sum( (X{j}(:,j1)-mean_sim{j}(j1)).* ...
       (X{j}(:,j2)-mean_sim{j}(j2))) / ...
       sqrt(var_sim{j}(j1)*var_sim{j}(j2))/Nmc;
   covvar_norm_sim{j}(j2,j1) = covvar_norm_sim{j}(j1,j2);
   covvar_sim{j}(j1,j2) = sum( (X{j}(:,j1)-mean_sim{j}(j1)).* ...
       (X{j}(:,j2)-mean_sim{j}(j2)))/Nmc;
   covvar_sim{j}(j2,j1) = covvar_sim{j}(j1,j2);
    end
    
end
end

cov_err = zeros(Nbin,1);
cov_norm_err = zeros(Nbin,1);
for j = 1:Nbin
    cov_err(j) = norm( cov_mat{j}-covvar_sim{j} )/norm(cov_mat{j});
    cov_norm_err(j) = norm( cov_norm_mat{j}-covvar_norm_sim{j} )/norm(cov_norm_mat{j});    
end


% compare error of marginal CDFs of conditional joint distributions
param_mat_CDF_sim = cell(Nbin,1);
for j = 1:Nbin
    param_mat_CDF_sim{j} = zeros(1024,2*Nparam);
for j1 = 1:3
    dmin = 1;
dmax = max(X{j}(:,j1))*1.2;
[~,~,param_mat_CDF_sim{j}(:,2*j1-1),param_mat_CDF_sim{j}(:,2*j1)]= ...
    kde(X{j}(:,j1),1024,dmin,dmax);
end
for j1 = 4:Nparam
    dmin = min(X{j}(:,j1));
    dmax = max(X{j}(:,j1));
[~,~,param_mat_CDF_sim{j}(:,2*j1-1),param_mat_CDF_sim{j}(:,2*j1)]= ...
    kde(X{j}(:,j1),1024,dmin,dmax);
end
end

% interpolate values over range
MCDF_error = zeros(Nbin,Nparam);
for j = 1:Nbin
    for j1 = 1:Nparam
x1 = param_mat_CDF_sim{j}(:,2*j1-1);
x2 = param_mat_CDF{j}(:,2*j1-1);
y1 = param_mat_CDF_sim{j}(:,2*j1);
y2 = param_mat_CDF{j}(:,2*j1);
xint1 = max([min(x1);min(x2)]);
xint2 = min([max(x1);max(x2)]);
xint = linspace(xint1,xint2,1024);
y1int = interp1(x1,y1,xint);
y2int = interp1(x2,y2,xint);
%MCDF_error(j,j1) = sum( (y1int-y2int).^2)/sum( y1int.^2);
MCDF_error(j,j1) = sqrt(sum((y1int-y2int).^2 / sum(y1int.^2)));

    end
end

end % if (testq ==0)

% ---------------------------------------------------------------------%
% simulate EGT parameters

% (1) simulate the centroids: Use Matern hard core proc 
N1 = 189; N2 = 201; N3 = 117;
Ngrain = 4500;
xmax = [N1*.25;N2*.25;N3*.25];
Centroid_simTmp = [xmax(2)*rand(Ngrain,1) xmax(2)*rand(Ngrain,1) xmax(2)*rand(Ngrain,1)];
z = find(Centroid_simTmp(:,1) < xmax(1) & Centroid_simTmp(:,3) < xmax(3));
Ngrain = length(z);
Centroid_sim = Centroid_simTmp(z,:);
%Centroid_sim = [xmax(1)*rand(Ngrain,1) xmax(2)*rand(Ngrain,1) xmax(3)*rand(Ngrain,1)];
markstmp = rand(Ngrain,1);
d = zeros(Ngrain,1);
% for j = 1:Ngrain
% z = setdiff(1:Ngrain,j);
%     d(j) = min( sqrt( (Centroid(j,1)-Centroid(z,1)).^2 + (Centroid(j,2)-Centroid(z,2)).^2 + ...
%         (Centroid(j,3)-Centroid(z,3)).^2));
% end
% R = min(d);
R = .226;

% Thinning loop: go through each point and find if distance to another 
% is closer. Remove one with smaller mark
grain_remain1 = [];
for j = 1:Ngrain
    d = ( sqrt( (Centroid_sim(j,1)-Centroid_sim(:,1)).^2 + (Centroid_sim(j,2)-Centroid_sim(:,2)).^2 + ...
        (Centroid_sim(j,3)-Centroid_sim(:,3)).^2));
    for j1 = 1:length(d)
        if (d(j1) < R)
            if (markstmp(j1) > markstmp(j)) grain_remain1 = [j grain_remain1]; end
            if (markstmp(j1) < markstmp(j)) grain_remain1 = [j1 grain_remain1]; end
        end % if
    end %j1

end %j
grain_remain2 = setdiff(1:Ngrain,grain_remain1);
Ngrain = length(grain_remain2);

% (2) simulate the grain parameters for every Poisson point
Nmc = Ngrain;
U = randn(Nmc,Nparam);
%Z = mvncdf(Y(1,:),zeros(1,Nparam),Nataf_rho);
%Z = mvncdf(Y,zeros(1,Nparam),Nataf_rho);
%Z = normcdf(Y);
X = zeros(Nmc,Nparam);
for j1 = 1:Nmc
% find 5 closest grains to bin grain parameters
z = setdiff(1:Ngrain,j1);
d = ( sqrt( (Centroid_sim(j,1)-Centroid_sim(grain_remain2(z),1)).^2 + ...
    (Centroid_sim(j,2)-Centroid_sim(grain_remain2(z),2)).^2 + ...
        (Centroid_sim(j,3)-Centroid_sim(grain_remain2(z),3)).^2));
d = sort(d);
b = mean(d(1:5));
jn = find(grain_bin > b,1,'first');
jn = 1;
Y = U*(Gammainv{jn}');
Y = normcdf(Y);    
    for j = 1:Nparam
        z2 = find(param_mat_CDF{jn}(:,2*j)<Y(j1,j),1,'last');
    if (isempty(z2)==1) z2 = 1; end    
    X(j1,j) = param_mat_CDF{jn}(z2,2*j-1);
    end
end

% (3) convert ellipse_orientation to rotation matrix
z = zeros(Ngrain,4);
rotmat = cell(Ngrain,1);
for j = 1:Ngrain
    rotmat{j} = zeros(3,3);
    z(j,4) = X(j,6);
    [z(j,1),z(j,2),z(j:3)] = sph2cart(X(j,4),X(j,5),1.0);
    rotmat{j} = vrrotvec2mat(z(j,:));
end

% write out simulated parameters
fil = 'results/EGT_Parameters_Simulation.txt';
fid = fopen(fil,'w');
XEG = [Centroid_sim(grain_remain2,:) X(:,1:3) X(:,7:9)];
for j = 1:Ngrain
fprintf(fid,'%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n',...
    j, XEG(j,1),XEG(j,2),XEG(j,3),XEG(j,4),XEG(j,5),XEG(j,6),rotmat{j}(1,1),...
    rotmat{j}(1,2),rotmat{j}(1,3),rotmat{j}(2,1),rotmat{j}(2,2),rotmat{j}(2,3),...
    rotmat{j}(3,1),rotmat{j}(3,2),rotmat{j}(3,3),XEG(j,7),XEG(j,8),XEG(j,9));
end
fclose(fid);

test1 = 0;
if (test1 ==0 )
% reconstruct microstructure
z = dlmread(fil);
XEG = z(:,2:7);
gorientation = z(:,17:19);
x1map = z(:,8:10);
x2map = z(:,11:13);
x3map = z(:,14:16);
N1 = 189; N2 = 201; N3 = 117;
dx = .25; xmin = [0;0;0]; xmax = [N1*dx;N2*dx;N3*dx];
[cellidEG, gneighborEG,grain_remainEG] = ...
    Tessellation_Simulation_EGT_3D_2(XEG,x1map,x2map,x3map,dx,xmin,xmax); 

% write cellidEG to file
fil = 'results/EGT_cellid_Simulation3D.txt';
fid = fopen(fil,'w');
fprintf(fid,'cellid of EGT simulation (N1=189,N2=201,N3=117). See files EGT_Simulation... for statistics');
for j = 1:(N1*N2*N3)
    fprintf(fid,'%d \n',cellidEG(j));
end
fclose(fid);

% write vtk file
fil = 'results/EGT_Sim_VTK.vtk';
Tessellation_VTK2(cellidEG,fil,xmin,dx);



elseif (test1==1)
% read cellidEG
N1 = 189; N2 = 201; N3 = 117;
cellidEG = zeros(N1,N2,N3);
fil = 'results/EGT_cellid_Simulation3D.txt';
fid = fopen(fil,'r');
tmp = fgets(fid);
for j = 1:(N1*N2*N3)
    cellidEG(j) = str2num(fgets(fid));
end
% compute various statistics (similar to tessellation_statistics2.m)    
grain_remainEG = unique(cellidEG);
NgrainEG = length(grain_remainEG);
EG_Centroid = zeros(NgrainEG,3);
EG_neigh = zeros(NgrainEG,1);
EG_vol = zeros(NgrainEG,1);
EG_m2 = zeros(NgrainEG,6); % I_xx, I_yy, I_zz, I_xy, I_xz, I_yz
for j = 1:NgrainEG
    z = find(cellidEG==grain_remainEG(j));
    EG_vol(j) = length(z)*dx^3;
    j3 = floor((z-1) / (N1*N2)  )+1;
j2 = floor( (z-1- N1*N2*(j3-1))/ N1) + 1;
j1 = z - N1*N2*(j3-1) - N1*(j2-1);
EG_Centroid(j,1) = sum((j1-.5)*dx+xmin(1))*dx^3/EG_vol(j);
EG_Centroid(j,2) = sum((j2-.5)*dx+xmin(2))*dx^3/EG_vol(j);
EG_Centroid(j,3) = sum((j3-.5)*dx+xmin(3))*dx^3/EG_vol(j); 
xdist = ((j1-.5)*dx+xmin(1) - EG_Centroid(j,1));
ydist = ((j2-.5)*dx+xmin(2) - EG_Centroid(j,2));
zdist = ((j3-.5)*dx+xmin(3) - EG_Centroid(j,3));
EG_m2(j,1) = sum( ydist.^2 + zdist.^2 )*dx^3;
EG_m2(j,2) = sum( xdist.^2 + zdist.^2 )*dx^3;
EG_m2(j,3) = sum( xdist.^2 + ydist.^2 )*dx^3;
EG_m2(j,4) = sum( xdist.*ydist )*dx^3;
EG_m2(j,5) = sum( xdist.*zdist )*dx^3;
EG_m2(j,6) = sum( ydist.*zdist )*dx^3;
end

% contiguous neighbors
for j = 1:NgrainEG
    EG_neigh(j) = length(gneighborEG{grain_remainEG(j)});
end

% grain shape distribution (fitting ellipsoid to each grain and outputting
% the three principal axes of ellipsoid). 
[EG_Shape,EG_Axis] = Get_Grain_Shape2(EG_Centroid,EG_m2,EG_vol);


fil = '/Users/kteferra/Documents/research/tessellation/results/EGT_Simulation_Vol_20140611.txt';        
fid1 = fopen(fil,'w');
fil = '/Users/kteferra/Documents/research/tessellation/results/EGT_Simulation_Neigh_20140611.txt';        
fid2 = fopen(fil,'w');
fil = '/Users/kteferra/Documents/research/tessellation/results/EGT_Simulation_Centroid_20140611.txt';        
fid3 = fopen(fil,'w');
fil = '/Users/kteferra/Documents/research/tessellation/results/EGT_Simulation_2ndMOI_20140611.txt';        
fid4 = fopen(fil,'w');
fil = '/Users/kteferra/Documents/research/tessellation/results/EGT_Simulation_Grain_Shape_20140611.txt';        
fid5 = fopen(fil,'w');
fprintf(fid1,'Grain ID (consistent with data grain id), Grain Volume in mm^3. (6/07/2014)\n');
fprintf(fid2,'Grain ID (consistent with data grain id),Number of neighbors per grain. (6/07/2014)\n');
fprintf(fid3,'Grain ID (consistent with data grain id), Centroid per grain (X,Y,Z). (6/07/2014)\n');
fprintf(fid4,'Grain ID (consistent with data grain id), 2nd Moment of Inertia per grain .');
fprintf(fid4,' (Ixx,Iyy,Izz,Ixy,Ixz,Iyz) (6/07/2014)\n');
fprintf(fid5,'Grain ID (consistent with data grain id), Best fit ellipsoid '); ...
fprintf(fid5,'(a,b,c). followed by directions for abc ');
fprintf(fid5, '(x_a,y_a,z_a,x_b_y_b,z_b,x_c,y_c,z_c) (6/07/2014) \n');
for j = 1:NgrainEG
fprintf(fid1,'%d, %f\n',grain_remainEG(j),EG_vol(j));
fprintf(fid2,'%d, %d\n',grain_remainEG(j),EG_neigh(j));
fprintf(fid3,'%d, %f, %f, %f\n',grain_remainEG(j),EG_Centroid(j,1),EG_Centroid(j,2),EG_Centroid(j,3));
fprintf(fid4,'%d, %f, %f, %f, %f, %f, %f\n',grain_remainEG(j),EG_m2(j,1),EG_m2(j,2),EG_m2(j,3), ...
    EG_m2(j,4),EG_m2(j,5),EG_m2(j,6));
fprintf(fid5,'%d, %f, %f, %f, ',grain_remainEG(j),EG_Shape(j,1),EG_Shape(j,2),EG_Shape(j,3));    
fprintf(fid5,'%f, %f, %f, %f, %f, %f, %f, %f, %f\n', EG_Axis{j}(1,1), ...
    EG_Axis{j}(2,1),EG_Axis{j}(3,1),EG_Axis{j}(1,2), EG_Axis{j}(2,2), ...
    EG_Axis{j}(3,2), EG_Axis{j}(1,3),EG_Axis{j}(2,3),EG_Axis{j}(3,3));
end %j
fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); fclose(fid5);

    
end %(if (test1 ==0 )

test2 = 0;
% compare simulated tessellation statistics with best fit
if (test2==0)
    
Ngrain = 2353;
bound = [1 0 Ngrain 1];
Dat_Vol = dlmread('results/Data_Vol_20140317.txt',',',bound);  
bound = [1 0 2334 1];
EG_Vol = dlmread('results/EGT_BF_Vol_20140317.txt',',',bound);
bound = [1 0 2403 1];
EG_Sim_Vol = dlmread('results/EGT_Simulation_Vol_20140611.txt',',',bound);

Nbinhist=80;
[dn dx] = hist(Dat_Vol(:,2),Nbinhist);
[en ex]  = hist(EG_Vol(:,2),Nbinhist);
[enSim exSim]  = hist(EG_Sim_Vol(:,2),Nbinhist);
dx = dx - (max(Dat_Vol(:,2))-min(Dat_Vol(:,2)))/Nbinhist/2;
ex = ex - (max(EG_Vol(:,2))-min(EG_Vol(:,2)))/Nbinhist/2;
exSim = exSim - (max(EG_Sim_Vol(:,2))-min(EG_Sim_Vol(:,2)))/Nbinhist/2;
for j = 1:Nbinhist-1
dn(j) = dn(j)/(dx(j+1)-dx(j));
en(j) = en(j)/(ex(j+1)-ex(j));
enSim(j) = enSim(j)/(exSim(j+1)-exSim(j));
end
dn(Nbinhist) = dn(Nbinhist)/(dx(Nbinhist)-dx(Nbinhist-1));
en(Nbinhist) = en(Nbinhist)/(ex(Nbinhist)-ex(Nbinhist-1));
enSim(Nbinhist) = enSim(Nbinhist)/(exSim(Nbinhist)-exSim(Nbinhist-1));

figure
stairs(dx, dn/length(Dat_Vol(:,1)), 'k','LineWidth',3);
hold on
stairs(ex, en/length(EG_Vol(:,1)), 'b','LineWidth',3);
stairs(exSim, enSim/length(EG_Sim_Vol(:,1)), 'r--','LineWidth',3);
xlabel('grain volume (\mu m^3)','FontSize',20)
axis([0.0 500.0 0.0 .1])
set(gca,'FontSize',20)
print -dpng results/Grain_Volume_BF_vs_Sim_EGT_20140611.png
figure
stairs(dx, dn/length(Dat_Vol(:,1)), 'k','LineWidth',3);
hold on
stairs(ex, en/length(EG_Vol(:,1)), 'b','LineWidth',3);
stairs(exSim, enSim/length(EG_Sim_Vol(:,1)), 'r--','LineWidth',3);
axis([0 100 0 .1])
set(gca,'FontSize',20)
hlegend=legend('Data','EGT BF','EGT Sim');
%set(hlegend,'FontSize',12)
print -dpng results/Grain_Volume_BF_vs_Sim_EGT_2_20140611.png

Ngrain = 2353;
bound = [1 0 Ngrain 3];
Dat_Shape = dlmread('results/Data_Grain_Shape_20140317.txt',',',bound);
bound = [1 0 2334 3];
EG_Shape = dlmread('results/EGT_BF_Grain_Shape_20140317.txt',',',bound);
bound = [1 0 2403 3];
EG_Sim_Shape = dlmread('results/EGT_Simulation_Grain_Shape_20140611.txt',',',bound);
EG_Shape2 = [EG_Shape(:,1) EG_Shape(:,4) EG_Shape(:,3) EG_Shape(:,2)];
EG_Shape = EG_Shape2;
EG_Shape2 = [EG_Sim_Shape(:,1) EG_Sim_Shape(:,4) EG_Sim_Shape(:,3) EG_Sim_Shape(:,2)];
EG_Sim_Shape = EG_Shape2;

Nbinhist=80;

Ndat = find(Dat_Shape(:,2) >.01 & Dat_Shape(:,3)>.01 & Dat_Shape(:,4)>.01 & Dat_Vol(:,2) > 2);
NEG = find(EG_Shape(:,2) >.01 & EG_Shape(:,3)>.01 & EG_Shape(:,4)>.01 & EG_Vol(:,2) > 2);
NEG_Sim = find(EG_Sim_Shape(:,2) >.01 & EG_Sim_Shape(:,3)>.01 & EG_Sim_Shape(:,4)>.01 & EG_Sim_Vol(:,2) > 2);

% I am scaling each term by its volume ratio
Dat_ShapeA = Dat_Shape(Ndat,4)./Dat_Shape(Ndat,2);
Dat_ShapeB = Dat_Shape(Ndat,4)./Dat_Shape(Ndat,3);
EG_ShapeA = EG_Shape(NEG,4)./EG_Shape(NEG,2);
EG_ShapeB = EG_Shape(NEG,4)./EG_Shape(NEG,3);
EG_Sim_ShapeA = EG_Sim_Shape(NEG_Sim,4)./EG_Sim_Shape(NEG_Sim,2);
EG_Sim_ShapeB = EG_Sim_Shape(NEG_Sim,4)./EG_Sim_Shape(NEG_Sim,3);

[dnA dxA] = hist(Dat_ShapeA,Nbinhist);
[dnB dxB] = hist(Dat_ShapeB,Nbinhist);
[enA exA] = hist(EG_ShapeA,Nbinhist);
[enB exB] = hist(EG_ShapeB,Nbinhist);
[enA_Sim exA_Sim] = hist(EG_Sim_ShapeA,Nbinhist);
[enB_Sim exB_Sim] = hist(EG_Sim_ShapeB,Nbinhist);

dnA = dnA/sum(dnA);
dnB = dnB/sum(dnB);
enA = enA/sum(enA);
enB = enB/sum(enB);
enA_Sim = enA_Sim/sum(enA_Sim);
enB_Sim = enB_Sim/sum(enB_Sim);
dxA = dxA - (max(Dat_ShapeA)-min(Dat_ShapeA))/Nbinhist/2;
dxB = dxB - (max(Dat_ShapeB)-min(Dat_ShapeB))/Nbinhist/2;
exA = exA - (max(EG_ShapeA)-min(EG_ShapeA))/Nbinhist/2;
exB = exB - (max(EG_ShapeB)-min(EG_ShapeB))/Nbinhist/2;
exA_Sim = exA_Sim - (max(EG_Sim_ShapeA)-min(EG_Sim_ShapeA))/Nbinhist/2;
exB_Sim = exB_Sim - (max(EG_Sim_ShapeB)-min(EG_Sim_ShapeB))/Nbinhist/2;

for j = 1:Nbinhist    
dnA(j) = dnA(j)/(dxA(2)-dxA(1));
dnB(j) = dnB(j)/(dxB(2)-dxB(1));
enA(j) = enA(j)/(exA(2)-exA(1));
enB(j) = enB(j)/(exB(2)-exB(1));
enA_Sim(j) = enA_Sim(j)/(exA_Sim(2)-exA_Sim(1));
enB_Sim(j) = enB_Sim(j)/(exB_Sim(2)-exB_Sim(1));
end

figure
stairs(dxA, dnA, 'k','LineWidth',3);
hold on
stairs(exA, enA, 'b','LineWidth',3);
stairs(exA_Sim, enA_Sim, 'r--','LineWidth',3);
xlabel('Best-fit ellipsoidal ratio (A/C)','FontSize',20)
ylabel('Probability','FontSize',20)
axis([1 10 0 1])
hlegend=legend('Data','EGT BF','EGT Sim');
set(gca,'FontSize',20)
print -dpdf results/Grain_Shape_BF_vs_Sim_EGTA_20140611.pdf

figure
stairs(dxB, dnB, 'k','LineWidth',3);
hold on
stairs(exB, enB, 'b','LineWidth',3);
stairs(exB_Sim, enB_Sim, 'r--','LineWidth',3);
xlabel('Best-fit ellipsoidal ratio (A/B)','FontSize',20)
ylabel('Probability','FontSize',20)
hlegend=legend('Data','EGT BF','EGT Sim');
axis([1 5 0 3])
set(gca,'FontSize',20)
print -dpdf results/Grain_Shape_BF_vs_Sim_EGTB_20140611.pdf

    
% number of neighbors

end % if (test2==0)




% ---------------------------------------------------------------------%
% analyze correlation strucures of velocities

% (1) see correlation as a function of distance

% find distances
dc = 1;
d = zeros(2353^2,1);
for j = 1:Ngrain
    z2 = setdiff(1:Ngrain,1:j);
    d(dc:(length(z2)+dc-1)) = sqrt( (Centroid(j,1)- Centroid(z2,1)).^2 + ...
        (Centroid(j,2)- Centroid(z2,2)).^2 + (Centroid(j,3)- Centroid(z2,3)).^2 );
    dc = dc + length(z2);

end
d10 = zeros(2353*10,1); % distribution of closest 10
d1loc = zeros(2353,2);
for j = 1:Ngrain
   z2 = setdiff(1:Ngrain,j);
   d =  sqrt( (Centroid(z2,1)-Centroid(j,1)).^2+ ...
       (Centroid(z2,2)-Centroid(j,2)).^2+ (Centroid(z2,3)-Centroid(j,3)).^2 );
   z3 = sort(d);
   %z2 = z((dsq<= z3(10)));
   d10(10*(j-1)+1:10*j) = (z3(1:10));
   d1loc(j,:) = [find(d==min(d)) j];
end
d1 = zeros(2353,1);
for j = 1:2353
    d1(j) = min(d10( 10*(j-1)+1:10*j) );
end

% find covariance of closest
v1 = 0;
z2 = velocity_bins{1};
% find distances between these
dv1 = zeros(length(z2)-1);
for j = 1:length(z2)
   z3 = setdiff(z2,z2(j));
   dv1(:,j) =  sqrt( (Centroid(z3,1)-Centroid(z2(j),1)).^2+ ...
       (Centroid(z3,2)-Centroid(z2(j),2)).^2+ (Centroid(z3,3)-Centroid(z2(j),3)).^2 );
   
end

% ---------------------------------------------------------------------%
% compute correlation of grain volume as a function of closest neighbors
% from data. My aim is to show that the data is uncorrelated spatially
bound = [1 0 Ngrain 1];
Dat_Vol = dlmread('results/Data_Vol_20140317.txt',',',bound);

% find grain id of closest 10 grains per grain
dId = zeros(2353,100);
for j = 1:Ngrain
    d =  sqrt( (Centroid(:,1)-Centroid(j,1)).^2+ ...
       (Centroid(:,2)-Centroid(j,2)).^2+ (Centroid(:,3)-Centroid(j,3)).^2 );
   [~,z3] = sort(d);
   dId(j,:) = z3(1:100);
end
Nclose = 11;
Vol_Corr = zeros(Nclose,1);
Euler_Corr = zeros(Nclose,9);
for j1 = 1:Nclose
    for j = 1:Ngrain
        Vol_Corr(j1) = (Dat_Vol(j,2)-mean(Dat_Vol(:,2)))*(Dat_Vol(dId(j,j1),2)-mean(Dat_Vol(:,2)))/var(Dat_Vol(:,2)) + Vol_Corr(j1);
        for j2 = 1:3
            for j3 = 1:3
        Euler_Corr(j1, 3*(j2-1)+j3) = gorientation(j,j2)*...
            gorientation(dId(j,j1),j3) + Euler_Corr(j1,3*(j2-1)+j3);
            end
        end
    end
    Vol_Corr(j1) = Vol_Corr(j1)/Ngrain;
    Euler_Corr(j1,:)=Euler_Corr(j1,:)/Ngrain;
end

sqmean = mean(Dat_Vol(:,2))^2*ones(Nclose,1);
plot(0:(Nclose-1),Vol_Corr,'x','MarkerSize',8,'LineWidth',5)
hold on
plot(0:(Nclose-1),sqmean,'r','LineWidth',2)
xlabel('Closest 10 grains','FontSize',20)
ylabel('Correlation of Volume','FontSize',20)
hlegend=legend('data', 'Square of mean');
%axis([0 5 0 3])
set(gca,'FontSize',20)
print -dpdf results/Volume_Correlation.pdf


j2 = 1; j3 = 1;
sqmean = mean(gorientation(:,j2))*mean(gorientation(:,j3))*ones(Nclose,1);
plot(0:(Nclose-1),Euler_Corr(:,3*(j2-1)+j3),'x','MarkerSize',8,'LineWidth',5)
hold on
plot(0:(Nclose-1),sqmean,'r','LineWidth',2)
xlabel('Closest 10 grains','FontSize',20)
ylabel('Correlation of 1st Euler angle','FontSize',20)
hlegend=legend('data', 'Square of mean');
%axis([0 5 0 3])
set(gca,'FontSize',20)
print -dpdf results/Euler1_Correlation.pdf



    

