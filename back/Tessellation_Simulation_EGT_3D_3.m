%%% functions

% this file generates a microstructure with high curvature GBs
% by randomly perturbing the ellipsoid velocity in any direction 
% this is somewhat arbitrary

function cellid= ...
    Tessellation_Simulation_EGT_3D_3(X,x1map,x2map,x3map,dx,xmin,xmax)



% xmin and xmax are the bounds of the domain
% Eorientation is the ellipsoidal orientation
% for each grain it gives the Euler angles wrt a reference
N1 = round((xmax(1)-xmin(1))/dx);
N2 = round((xmax(2)-xmin(2))/dx);
N3 = round((xmax(3)-xmin(3))/dx);
Ngrain = length(X(:,1));
gcompete = 1:Ngrain;
v1 = X(gcompete,4);
v2 = X(gcompete,5);
v3 = X(gcompete,6);
cellid = zeros(N1,N2,N3);
voxelid = (1:(N1*N2*N3))';
badgrain = zeros(length(voxelid),15);
% find grain that arrives to each voxel first
for j1 = 1:N1
    for j2 = 1:N2
        for j3 =1:N3
node = [(j1-.5)*dx+xmin(1); (j2-.5)*dx+xmin(2);(j3-.5)*dx+xmin(3)];
% must establish angle of each grains nucleation site to node of interest
x1 = node(1)-X(gcompete,1);
x2 = node(2)-X(gcompete,2);
x3 = node(3)-X(gcompete,3);
d = sqrt(x1.^2 + x2.^2 + x3.^2);

xp1 = x1map(gcompete,1).*x1 + x1map(gcompete,2).*x2 + x1map(gcompete,3).*x3;
xp2 = x2map(gcompete,1).*x1 + x2map(gcompete,2).*x2 + x2map(gcompete,3).*x3;
xp3 = x3map(gcompete,1).*x1 + x3map(gcompete,2).*x2 + x3map(gcompete,3).*x3;

th = atan( abs(xp2)./abs(xp1)    );
ph = atan(  sqrt((xp1).^2 + (xp2).^2) ./ abs(xp3));
vang = sqrt(v1.^2.*v2.^2.*v3.^2./( v2.^2.*v3.^2.*(cos(th)).^2.*(sin(ph)).^2 + ...
v1.^2.*v3.^2.*(sin(th)).^2.*(sin(ph)).^2 + v1.^2.*v2.^2.*(cos(ph)).^2 ) );

t1 = d./vang;
cellid(j1,j2,j3) = gcompete(find(t1==min(t1),1));  
        end %j3
    end %j2
end %j1
j3
% c = 0;
% % 
% % loop through following until voxelid is empty
% while (isempty(voxelid)==0)
% % find all voxels disconnected from grain
% voxelid = [];
% c = c + 1;
% for j = 1:Ngrain
% gc = find(cellid==j);
% j1 = round((X(j,1)-xmin(1))/dx);
% j2 = round((X(j,2)-xmin(2))/dx);
% j3 = round((X(j,3)-xmin(3))/dx);
% baccounted = N1*N2*(j3-1)+N1*(j2-1)+j1;
% z = length(baccounted)+1;
% while (z ~= length(baccounted))
% z = length(baccounted);
% potneigh = [baccounted+1; baccounted+N1; baccounted-1;...
% baccounted-N1; baccounted+N1*N2; baccounted-N1*N2];
% baccounted = unique([intersect(gc,potneigh); baccounted]);
% end %while
% z = setdiff(gc,baccounted);
% badgrain(z,c) = j;  
% voxelid = [z; voxelid];
%   
% end % j
% for j = 1:length(voxelid)
% j3 = floor((voxelid(j)-1) / (N1*N2)  )+1;
% j2 = floor( (voxelid(j)-1- N1*N2*(j3-1))/ N1) + 1;
% j1 = voxelid(j) - N1*N2*(j3-1) - N1*(j2-1);
% gcompete = setdiff(1:Ngrain,badgrain(voxelid(j),:))';    
% node = [(j1-.5)*dx+xmin(1); (j2-.5)*dx+xmin(2);(j3-.5)*dx+xmin(3)];
% v1 = X(gcompete,4);
% v2 = X(gcompete,5);
% v3 = X(gcompete,6);
% x1 = node(1)-X(gcompete,1);
% x2 = node(2)-X(gcompete,2);
% x3 = node(3)-X(gcompete,3);
% d = sqrt(x1.^2 + x2.^2 + x3.^2);
% xp1 = x1map(gcompete,1).*x1 + x1map(gcompete,2).*x2 + x1map(gcompete,3).*x3;
% xp2 = x2map(gcompete,1).*x1 + x2map(gcompete,2).*x2 + x2map(gcompete,3).*x3;
% xp3 = x3map(gcompete,1).*x1 + x3map(gcompete,2).*x2 + x3map(gcompete,3).*x3;
% th = atan( abs(xp2)./abs(xp1)    );
% ph = atan(  sqrt((xp1).^2 + (xp2).^2) ./ abs(xp3));
% vang = sqrt(v1.^2.*v2.^2.*v3.^2./( v2.^2.*v3.^2.*(cos(th)).^2.*(sin(ph)).^2 + ...
% v1.^2.*v3.^2.*(sin(th)).^2.*(sin(ph)).^2 + v1.^2.*v2.^2.*(cos(ph)).^2 ) );
% t1 = d./vang;
% cellid(j1,j2,j3) = gcompete(find(t1==min(t1),1));  
% end %j
% 
% end % while
% 
% % find all grains that have only 1 neighbor (means enclosed)
% % this grain is eliminated and its enclosing grain takes its voxels
% gneighbor = cell(Ngrain,1);
% for j1 = 1:N1-1
%     for j2 = 1:N2-1
%         for j3 = 1:N3-1
%  if (cellid(j1,j2,j3) ~= cellid(j1+1,j2,j3))
% gneighbor{cellid(j1,j2,j3)} = [cellid(j1+1,j2,j3);gneighbor{cellid(j1,j2,j3)}];
% gneighbor{cellid(j1+1,j2,j3)} = [cellid(j1,j2,j3);gneighbor{cellid(j1+1,j2,j3)}];
% end % if
%  if (cellid(j1,j2,j3) ~= cellid(j1,j2+1,j3))
% gneighbor{cellid(j1,j2,j3)} = [cellid(j1,j2+1,j3);gneighbor{cellid(j1,j2,j3)}];
% gneighbor{cellid(j1,j2+1,j3)} = [cellid(j1,j2,j3);gneighbor{cellid(j1,j2+1,j3)}];
%  end % if        
%   if (cellid(j1,j2,j3) ~= cellid(j1,j2,j3+1))
% gneighbor{cellid(j1,j2,j3)} = [cellid(j1,j2,j3+1);gneighbor{cellid(j1,j2,j3)}];
% gneighbor{cellid(j1,j2,j3+1)} = [cellid(j1,j2,j3);gneighbor{cellid(j1,j2,j3+1)}];
%  end % if        
%         end %j3
%     end %j2
% end %j1
% grain_remain = 1:Ngrain;
% for j = 1:Ngrain
%     gneighbor{j} = unique(gneighbor{j});
%     if (length(gneighbor{j})<=1)
%         gc = (cellid==j);
%         cellid(gc) = gneighbor{j};        
%         grain_remain = setdiff(grain_remain,j);
%     end %if
% end %j
% % find boundary surface faces on voxels
% for j = 1:Ngrain
%     gneighbor{j}
% end
end % function
