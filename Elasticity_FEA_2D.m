% Elasticity_FEA_2D.m
% Solves for the stress & strain on a 2D surface due to various loads,
%   boundary conditions, and material properties.
% This file uses the companion file Input.dat to define the problem.
% by Sean Ehle
% 2011-12-18

clear
%% acquire data from the file Input.dat to define the problem
fid = fopen('Input.dat'); % open a file path to the input data
Stress = ~cell2mat(textscan(fid,'%d','Headerlines',4)); %stress/strain flag
NodeInfo = cell2mat(textscan(fid,'%f%f%f','Headerlines',2)); % node locs
ElemInfo = cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f%f%f','Headerlines',2));
PntLoads = cell2mat(textscan(fid,'%f%f%f%f','Headerlines',2)); % point load
SrfLoads = cell2mat(textscan(fid,'%f%f%f%f%f%f','Headerlines',2)); %surface
BdyLoads = cell2mat(textscan(fid,'%f%f%f%f','Headerlines',2)); % body force
BConInfo = cell2mat(textscan(fid,'%f%f%f%f','Headerlines',3)); % B.C. info
fclose(fid); % close the file path, once all info is acquired
clear fid % clear unneeded variable from memory

%% initialize global system variables
m = 2*size(NodeInfo,1); % multi-use integer variable
K = sparse(zeros(m)); % stiffness matrix
u = NaN(m,1); % nodal displacement vector
f = zeros(m,1); % nodal force vector

%% initialize elemental variables to prevent growth inside loops
m = size(ElemInfo,1); % multi-use integer variable
D = zeros(3,3,m); % D matrix
B = zeros(3,6,m); % gradient matrix
GlAd = zeros(m,6); % global address of element nodes
epT = zeros(3,m); % element thermal strain (eps_xx eps_yy eps_xy)
eps = epT; % element strain (eps_xx eps_yy gam_xy)
sig = epT; % element stress (sig_xx sig_yy tau_xy)
sigp = epT; % principle stresses (sig_1 sig_2 sig_3)
x = epT; % x-matrix for plotting
y = epT; % y-matrix for plotting
sigz = zeros(1,m); % out of plane stress for strain analysis (sigma_zz)

%% assemble the global system
for m = 1:size(PntLoads,1) % loop through each point force component
  n = 2*PntLoads(m,2) - PntLoads(m,3); % global address index
  f(n) = f(n) + PntLoads(m,4); % add force to node / create {p}
end % for point force

for m = 1:size(BConInfo,1) % loop through each boundary condition
  n = 2*BConInfo(m,2) - BConInfo(m,3); % global address index
  u(n) = BConInfo(m,4); % define the displacement / {u}
end % for BC

for m = 1:size(ElemInfo,1) % loop through each element
  i = ElemInfo(m,2);  j = ElemInfo(m,3);  k = ElemInfo(m,4); % find i,j,k
  t = ElemInfo(m,5); % thickness in z-dir
  E = ElemInfo(m,6); % modulus of Elasticity
  nu = ElemInfo(m,7); % Poisson's ratio
  alpha = ElemInfo(m,8); % alpha*delta_T
  if Stress % build D matrix for stress
    D(:,:,m) = E/(1 - nu^2)*[ 1 nu      0     ;
                             nu  1      0     ;
                              0  0 (1 - nu)/2];
  else % build D matrix for strain
    D(:,:,m) = E/(1 + nu)*[(1-nu)/(1-2*nu)   nu/(1-2*nu)   0   ; 
                             nu/(1-2*nu)   (1-nu)/(1-2*nu) 0   ;
                                 0                 0       0.5];
  end % if D matrix
  A = 0.5*det([1 NodeInfo(i,2:3); % Area of element
               1 NodeInfo(j,2:3);
               1 NodeInfo(k,2:3)]);
  % Calculate b,c coefficients
  bi = diff(NodeInfo([k j],3));  ci = diff(NodeInfo([j k],2));
  bj = diff(NodeInfo([i k],3));  cj = diff(NodeInfo([k i],2));
  bk = diff(NodeInfo([j i],3));  ck = diff(NodeInfo([i j],2));
  % create gradient matrix
  B(:,:,m) = 1/(2*A)*[bi  0 bj  0 bk  0 ;
                       0 ci  0 cj  0 ck ;
                      ci bi cj bj ck bk];
  Ke = B(:,:,m)'*D(:,:,m)*B(:,:,m)*t*A; % element stiffness matrix
  epT(:,m) = [alpha;alpha;0]; % thermal strain
  fe = B(:,:,m)'*D(:,:,m)*epT(:,m)*t*A; % thermal load
  
  for n =1:size(SrfLoads,1)  % search all surface loads
    if SrfLoads(n,2)==m  % if load is on this element
      L = norm(diff(NodeInfo(SrfLoads(n,3:4),2:3))); % length between nodes
      if SrfLoads(n,5) % if x-direction
        temp = [SrfLoads(n,3)==ElemInfo(2:4)|SrfLoads(n,4)==ElemInfo(2:4);
                0 0 0]; % find ijk for nodes, add zeros after for y-vals
        temp = temp(:); % unravel the 2x3 into a 6x1
        fe = fe + t*L/2*SrfLoads(n,6)*temp; % add half load to each node
      else % if y-direction
        temp = [0 0 0; % find ijk for nodes, add zeros before for x-vals
                SrfLoads(n,3)==ElemInfo(2:4)|SrfLoads(n,4)==ElemInfo(2:4)];
        temp = temp(:); % unravel the 2x3 into a 6x1
        fe = fe + t*L/2*SrfLoads(n,6)*temp; % add half load to each node
      end % if direction
    end % if element
  end % surface load loop
  
  for n =1:size(BdyLoads,1)  % search all body loads
    if BdyLoads(n,2)==m || BdyLoads(n,2)==inf % if on this element
      if BdyLoads(n,3) % if x-direction
        fe = fe + t*A/3*BdyLoads(n,4)*[1;0;1;0;1;0];
      else % if y-direction
        fe = fe + t*A/3*BdyLoads(n,4)*[0;1;0;1;0;1];
      end % if direction
    end % if element
  end % body load loop
  
  GlAd(m,:) = [2*i-1 2*i 2*j-1 2*j 2*k-1 2*k]; % global addresses
  K(GlAd(m,:),GlAd(m,:)) = K(GlAd(m,:),GlAd(m,:)) + Ke; % element stiffness
  f(GlAd(m,:)) = f(GlAd(m,:)) + fe; % element forces
end % for primary element loop
clear PntLoads SrfLoads BdyLoads BConInfo Ke fe % clear unneeded variables
clear i j k t E nu alpha A bi bj bk ci cj ck L % clear unneeded variables

%% solve for unknown nodal displacements and element strains and stresses
temp = find(isnan(u)); % find address indices of unknown displacements
u(temp) = K(temp,temp)\f(temp); % trim, solve, and replace into global
for m = 1:size(ElemInfo,1) % loop through each element
  eps(:,m) = B(:,:,m)*u(GlAd(m,:)); % solve for total strain
  sig(:,m) = D(:,:,m)*(eps(:,m) - epT(:,m)); % solve for elastic stress
  if ~Stress % if strain analysis
    sigz(m) = ElemInfo(m,7)*(sig(1,m) + sig(2,m)) - ...
      ElemInfo(m,6)*ElemInfo(m,8);
    temp = eig([sig(1,m) sig(3,m)     0    ; % principle stresses
                sig(3,m) sig(2,m)     0    ;
                   0        0     sigz(m)]);
    sigp(1,m) = max(temp); % sigma_1
    sigp(3,m) = min(temp); % sigma_3
    sigp(2,m) = sum(temp) - sigp(1,m) - sigp(3,m); % sigma_2
  end % if strain
end % for element loop
if Stress
  fprintf('\n-----------------------------\nS T R E S S   A N A L Y S I S')
  fprintf('\n-----------------------------\n')
  C = (sig(1,:) + sig(2,:))/2; % centers of Mohr's circles
  R = sqrt(((sig(1,:) - sig(2,:))/2).^2 + sig(3,:).^2); % radii of circles
  sigp(1,:) = C + R; % principle stresses sigma_1
  sigp(2,:) = C - R; % principle stresses sigma_2  ;  sigma_3 = 0
  Tresca = (max(sigp(1,:),0) - min(sigp(2,:),0))/2; % Tresca equivalent
else
  fprintf('\n-----------------------------\nS T R A I N   A N A L Y S I S')
  fprintf('\n-----------------------------\n')
  Tresca = (sigp(1,:) - sigp(3,:))/2; % Tresca equivalent stress
end
vMises = sqrt(((sigp(1,:)-sigp(2,:)).^2 + (sigp(2,:)-sigp(3,:)).^2 + ...
  (sigp(1,:)-sigp(3,:)).^2)/2); % von Mises equivalent stress
clear B D GlAd n temp Stress % clear unneeded variables from memory

%% output analysis to the screen
fprintf('\n\tvon Mises analysis\n')
fprintf('The maximum stress   (%9.3e Pa) is in element(s) %d.\n',...
  max(vMises),find(vMises==max(vMises)))
FSv = ElemInfo(:,9)'./vMises; % calculate factor of safety in each element
fprintf('The minimum factor of safety (%4.3f) is in element(s) %d.\n',...
  min(FSv),find(FSv==min(FSv)))
fprintf('\n\tTresca analysis\n')
fprintf('The maximum shear    (%9.3e Pa) is in element(s) %d.\n',...
  max(Tresca),find(Tresca==max(Tresca)))
FSt = ElemInfo(:,10)'./Tresca; % calculate factor of safety in each element
fprintf('The minimum factor of safety (%4.3f) is in element(s) %d.\n',...
  min(FSt),find(FSt==min(FSt)))

%% plot the solution
for m = 1:3 % loop through each polygon vertex (3 for triangles)
  x(m,:) = NodeInfo(ElemInfo(:,m+1),2)'; % element x-coords for i;j;k nodes
  y(m,:) = NodeInfo(ElemInfo(:,m+1),3)'; % element x-coords for i;j;k nodes
end % for polygon vertices
figure(1) % open figure 1
set(gcf,'Color','w') % set background color to white
subplot(1,2,1)
fill(x,y,vMises) % plot VonMises stress in each element
axis equal % retain geometric shape (a square is a square)
title('von Mises') % title the graph
colorbar % add legend
subplot(1,2,2)
fill(x,y,Tresca) % plot Tresca stress in each element
axis equal % retain geometric shape (a square is a square)
title('Tresca') % title the graph
colorbar % add legend
clear NodeInfo ElemInfo m x y % clear unneeded variables from memory
