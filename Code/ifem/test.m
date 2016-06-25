function  [N,equilError, resError, mixedError, errH1] = test(method, afem, mixed)
  savedir = 'figures';
  %clear all;
  %close all; 
  %profile on;
  %% Parameters
  theta = 0.3;    uh =0;

  global refinemethod;
  refinemethod = @uniformrefine;

  maxN =  2e5;
  maxIt = 12;
  
  switch method
  case 'square_sin'
    [node, elem, bdFlag, pde] = squaresin();
  case 'square_ana'
    [node, elem, bdFlag, pde] = squareana(10);
  case 'square_peak'
    [node, elem, bdFlag, pde] = squarepeak(10, 0.51, 0.117);
  case 'square_one'
    [node, elem, bdFlag, pde] = squareone();
  case 'lshape_corner'
    [node, elem, bdFlag, pde] = lshapecorner();
  case 'lshape_one'
    [node, elem, bdFlag, pde] = lshapeone();
  case 'crack_one'
    [node, elem, bdFlag, pde] = crackone();
  end

  % AFEM suffix
  if afem
    method =sprintf('%s_afem', method);
  end
  if mixed
    method = sprintf('%s_mixed', method);
  end
  mkdir(sprintf('%s/%s', savedir, method)); clf;

  showmesh(node, elem);
  saveas(gca, sprintf('%s/%s/mesh_initial.png',savedir, method));

  plotapproxh1(method,node, elem, pde, bdFlag,pde.Du, 10)
  return
  % one uniform refinement
  [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
  if afem
    CompareAfem(method, node, elem, pde, bdFlag, theta,maxN)
  else
    [N,equilError, resError, mixedError, errH1] = CompareUniform(method, node, elem, pde, bdFlag, maxN)
  end
  %profile viewer

  return;

  %figure(1);  showresult(node,elem,uh,[-50,12]);    
  fluxpatch(node, elem, uh, Duh, pde.f, 2);
  sig = flux(node,elem, uh, Duh, pde.f);
  divfluxelem(node,elem, sig);
  sigelem = fluxelem(node, elem, sig);
  errorestimate(node, elem, Duh, sig, pde.f)
  %errorestimate(node, elem, Duh, sig, pde.f)
  fluxequil(node, elem, sig, pde.f);
  [y,i] = max(fluxequil(node, elem, sig,pde.f))
  return;
  %figure(1);  showresult(node,elem,uh,[-50,12]);    
  %showmesh(node,elem);
  %findelem(node, elem);

  NE = size(T.edge,1)
  N = size(node,1)
  NT = size(elem, 1)

  % Print
  return






  
  

  %showmesh(node, elem)
  return;
  [bdNode, bdEdge, isBdNode] = findboundary(elem);

  % Generate mapping from elements to vertices
  t2v = sparse([1:NT, 1:NT, 1:NT], elem, 1, NT, N)
  % Find all elements using vertex 9
  nodeStar = find(t2v(:,9));
  % Find create list of the elements in our patch
  elemStar = elem([nodeStar],:)
  % Find boundary of this patch
  [~, bdEdgeStar, ~] = findboundary(elemStar)
  % Create index matrix belonging to this patch
  % isBdEdgeStar = sparse(bdEdgeStar(:,1), bdEdgeStar(:,2), 1);
  % Create vector of element rows
  bla = [sort(elemStar(1,[2,3]));sort(elemStar(1,[1,3]));sort(elemStar(1,[1,2]))];
  % Check which of the edges is on the boundary
  isbdEdgeElem = ismember(bla, bdEdgeStar, 'rows')
  % Get coordinates of this element
  coord = node(elem(1,:),:)
  % Create local mass matrix
  massB = stimaB(coord')
  % Remove boundary edge
  massB(isbdEdgeElem, :) = 0;
  massB(:, isbdEdgeElem) = 0;
  massB
  % Remove boundary coordinates
  % massB(
  return
  [B_K, ~, B_K_det] = affine_transformations(node, elemStar);
  [elems2edges, edges2nodes] = get_edges(elemStar);
  elems2edges
  edges2nodes
  signs_e = signs_edges(elemStar)
  K_RT0 = stiffness_matrix_RT0(elems2edges,B_K_det,signs_e);
  M_RT0 = mass_matrix_RT0(elems2edges,B_K,B_K_det,signs_e);
  return


  %% Set up PDE data
  pde.f = 0;
  pde.g_D = @exactu;
  pde.Du=[];% used for recoverFlux;
  pde.d=[];
  %%  Adaptive Finite Element Method
  % *SOLVE* -> *ESTIMATE* -> *MARK* -> *REFINE*
  for k = 1:2
      % Step 1: SOLVE
      [uh,Duh,eqn] = Poisson(node,elem,pde,bdFlag);
      % Plot mesh and solution
      figure(1);  showresult(node,elem,uh,[-50,12]);    
      % Step 2: ESTIMATE
      eta = estimateresidual(node,elem,uh,pde);    % residual type
      sig = flux(node,elem, uh)
      % Step 3: MARK
      markedElem = mark(elem,eta,theta);
      % Step 4: REFINE
      [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
  end
end

function err = exactH1error(node, elem, bdFlag, pde, uh)
  if isempty(pde.Du)
    err = approxH1error(node, elem, bdFlag, pde, uh,3);
  else
    err = getH1error(node,elem,pde.Du,uh);
  end
end


function CompareAfem(method, nodeOri, elemOri, pde, bdFlagOri, theta,maxN) 
  function plot
    f2 = figure(2);clf;
    loglog(Nu,  erru, 'b-o',Nr, errr, 'r-x'); hold on;
    loglog(Ne,  erre, '-s', 'color', [0 0.5 0]); 
    legend({'uniform($U_k$)','residual($U_k$)', 'equilibrated($U_k, \zeta$)'}, 'interpreter', 'latex');
    title('Comparison AFEM performance')
    xlabel('Number of vertices');
    ylabel('Exact error');
    saveas(f2, sprintf('%s/%s/norm_%d_%g.png',savedir, method, maxN, theta));
    saveas(f2, sprintf('%s/%s/norm_%d_%g.fig',savedir, method, maxN, theta));
    loglog(Nm, errm, 'm-d');
    legend({'uniform($U_k$)','residual($U_k$)', 'equilibrated($U_k, \zeta$)', 'mixed($U_k, \sigma$)'}, 'interpreter', 'latex');
    saveas(f2, sprintf('%s/%s/norm_mixed_%d_%g.png',savedir, method, maxN, theta));
    saveas(f2, sprintf('%s/%s/norm_mixed_%d_%g.fig',savedir, method, maxN, theta));
  end

  savedir = 'figures/';
  % Uniform refinements
  Nu = [];
  erru = [];
  % AFEM residual estimator
  Nr = [];
  errr = [];
  % AFEM mixed estimator
  Nm = [];
  errm = [];
  % AFEM equil estimator
  Ne = [];
  erre = [];


  % Uniform refinement
  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
  while (size(node,1) < maxN)
    uh = Poisson(node, elem, pde, bdFlag);
    % Calculate real error
    Nu(end+1) = size(node, 1);
    erru(end+1) = exactH1error(node, elem, bdFlag, pde, uh);

    [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
  end
  plot
  erru

  % Residual refinements
  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
  while (size(node,1) < maxN)
    uh = Poisson(node, elem, pde, bdFlag);

    % Calculate the real error
    Nr(end+1) = size(node, 1);
    errr(end+1) = exactH1error(node, elem, bdFlag, pde,uh);

    % Perform afem step
    eta = estimateresidual(node, elem, uh, pde);
    markedElem = mark(elem,eta,theta);
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
  end
  plot
  errr

  % Mixed refinements
  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
  while (size(node,1) < maxN)
    uh = Poisson(node, elem, pde, bdFlag);
    [Duh,~] = gradu(node, elem, uh);

    % Calculate the real error
    Nm(end+1) = size(node, 1);
    errm(end+1) = exactH1error(node, elem, bdFlag, pde,uh);

    % Perform afem step
    [~,sigma] = PoissonRT0(node, elem, pde, bdFlag);
    [~, eta]  = getL2errorRT0(node, elem,  Duh, sigma);
    
    markedElem = mark(elem,sqrt(eta),theta);
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
  end
  plot
  errm


  % Equilibrated refinements
  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
  while (size(node, 1) < maxN)
    plot
    uh = Poisson(node, elem, pde, bdFlag);
    [Duh,~] = gradu(node, elem, uh);

    % Calculate the real error
    Ne(end+1) = size(node, 1);
    erre(end+1) = exactH1error(node, elem, bdFlag, pde,uh);

    % Perform afem step
    sig = flux(node,elem,  Duh, pde.f);
    [~, eta, ~] =  equilresidualestimate(node, elem, Duh, sig, pde.f);
    markedElem = mark(elem,eta,theta);
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
  end
  erre
  plot
end

function [N,equilError, resError, mixedError, errH1] = CompareUniform(method, node, elem, pde, bdFlag, maxN)
  errH1 = [];
  mixedError = [];
  equilError = [];
  resError = [];
  N = [];
  % Uniform refinement
  t = 1;
  while  (size(node, 1) < maxN)
    t = t+1;
    N(end+1) = size(node,1);

    % Solve the poisson problem
    uh = Poisson(node, elem, pde, bdFlag);
    [Duh,~] = gradu(node, elem, uh);

    % Calculate real error
    if isempty(pde.Du)
      errH1(end+1) = approxH1error(node, elem, bdFlag, pde, uh,3)
    else
      errH1(end+1) = getH1error(node,elem,pde.Du,uh)
    end
    [~,sigma] = PoissonRT0(node, elem, pde, bdFlag);
    mixedError(end+1) = getL2errorRT0(node, elem,  Duh, sigma);

    size(elem)

    % Calculate the flux
    sig = flux(node,elem,  Duh, pde.f);

    %Calculate the error esimate
    [equilError(end+1), ~, ~] =  equilresidualestimate(node, elem, Duh, sig, pde.f);

    %Calculate the residual eror
    eta = estimateresidual(node, elem, uh, pde);
    resError(end+1) = sqrt(sum(eta.^2));
    %ploterror(method,node, elem, uh, N,equilError, resError,mixedError, errH1, pde.theorate);

    [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
  end
end

function plotuh
  for j =1:9
    uh = Poisson(node, elem, pde, bdFlag);
    %f1 = figure(1); showmesh(node,elem);
    f2 = figure(2); %showsolution(node, elem, uh,2); colorbar;
    colormap('jet')
    trisurf(elem, node(:,1), node(:,2), uh', 'FaceColor', 'interp', 'EdgeColor', 'interp');
    caxis([-1,1])
    h = colorbar;
    view(2)
    saveas(f2, sprintf('%s/result%d.png', savedir, j));
    [node,elem, bdFlag] = uniformbisect(node, elem);
  end
  f2 = figure(2);
  showsolution(node,elem,uh,3)
  saveas(f2, 'figures/result.png');
end

function plotapproxh1(method,node, elem, pde, bdFlag,Du, maxIt)
  savedir = 'figures';
  mkdir(sprintf('%s/%s', savedir, method));
  errH1 = [];
  approx1H1 = [];
  approx2H1 = [];
  approx3H1 = [];
  N = [];
  t = 1;
  while  (t <= maxIt) 
    [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
    N(end+1) = size(node,1);
    uh = Poisson(node, elem, pde, bdFlag);
    %figure(1);  showresult(node,elem,uh,[-50,12]);    
    [Duh,~] = gradu(node, elem, uh);

    % Calculate real error
    errH1(end+1) = getH1error(node,elem,Du,uh)
      approx1H1(end+1) = approxH1error(node, elem, bdFlag, pde, uh,1)
    if (t <= maxIt - 1)
      approx2H1(end+1) = approxH1error(node, elem, bdFlag, pde, uh,2)
    end
    if (t <= maxIt - 2)
      approx3H1(end+1) = approxH1error(node, elem, bdFlag, pde, uh,3)
    end
    % Draw figuress!
    f2 = figure(2); clf;
    len1 = size(approx1H1,2);
    len2 = size(approx2H1,2);
    len3 = size(approx3H1,2);

    loglog(N(1:len1),  approx1H1 ./ errH1(1:len1), 'r-x',N(1:len2),  approx2H1 ./ errH1(1:len2), 'k-s', N(1:len3), approx3H1 ./ errH1(1:len3), 'm-d');
    xlabel('Number of vertices');
    ylabel('Relative error: $\|\nabla{U_{k+i}} - \nabla{ U_k}\| / \|\nabla{u} - \nabla{U_k}\|$', 'interpreter', 'latex');
    legend({'$i=1$',' $i=2$', '$i=3$'}, 'interpreter', 'latex');

    f3 = figure(3);clf;
    loglog(N,  errH1, 'b-o',N(1:len1), approx1H1, 'r-x', N(1:len2), approx2H1, 'k-s', N(1:len3), approx3H1, 'm-d'); 
    legend({'$\|\nabla{u} - \nabla {U_k}\|$','$\|\nabla{U_{k+1}} - \nabla{U_k}\|$',' $\|\nabla{U_{k+2}} - \nabla{U_k}\|$', '$\|\nabla{U_{k+3}} - \nabla{U_k}\|$'}, 'interpreter', 'latex');
    xlabel('Number of vertices');
    ylabel('Error');

    saveas(f3, sprintf('%s/%s/approx_H1_%d.fig',savedir, method, t) );
    saveas(f2, sprintf('%s/%s/approx_H1_rel_%d.fig',savedir, method, t));

    saveas(f3, sprintf('%s/%s/approx_H1_%d.eps',savedir, method, t), 'epsc');
    saveas(f2, sprintf('%s/%s/approx_H1_rel_%d.eps',savedir, method, t), 'epsc');

    t = t+1;
  end
end

% Approximate H1 error by two times bisecting
function err = approxH1error(node, elem, bdFlag, pde, uh, it);
  % Refine the grid it times; interpolate uh to this new grid
  for i=1:it
    [node, elem, bdFlag, HB] = uniformbisect(node, elem, bdFlag);
    uh = nodeinterpolate(uh, HB);
  end
    
  % Calculate the `real' solution for this grid
  u = Poisson(node, elem, pde, bdFlag);

  % uh - u lives on the same triangulation,  and Duh is constant on each triangle
  [Derr,area] = gradu(node, elem, uh - u);

  % error elementwise
  err = sqrt(sum( area.*sum(Derr.^2,2)));

  % Calculate the error
  %err = getH1error(node, elem, @(p) 0*p, uh - u)
end

%  Returns global information about the edges of a triange
%    [edges, signs, patchedge] = elem2edge(t, patchvert, T)
%      Given element t, auxstructure T and the global index
%      of this patch vertex, it returns the global indices of 
%      the edges, the global signs of these edges (-1 or +1), and
%      indicates which edge is connected to the patchvertex
function [edges,signs,patchedge] = elem2edge(t, patchvert, T)
  % Indices of the edges used by triang t
  edges = T.elem2edge(t,:);
  % Returns the list of `main' triangles belonging to these edges
  tmp = T.edge2elem(edges,1);
  % If triangle t is the first index, we get positive sign, else -1
  signs(1:3) = -1;
  signs(find(tmp == t)) = 1;

  % Find all the vertices belonging to the edges of t
  tmp = T.edge(edges,:);
  % any returns per column, transpose for rowwise
  patchedge = any(tmp' == patchvert);
end

function int = quadmid(nodes, area, f) % apply midpoint quadrature on triangle with nodes
  midpoints = 0.5 * [nodes(3,:) + nodes(2,:); nodes(1,:) + nodes(3,:); nodes(1,:) + nodes(2,:)];
  int = area / 3.0 * (f(midpoints(1,:)) + f(midpoints(2,:)) + f(midpoints(3,:)));
end

% Calculates the integral of u_h over the entire domain
function int = intdomain(node, elem, uh)
  NT = size(elem, 1);
  area = simplexvolume(node, elem);
  int = 0;
  for t=1:NT
    % Values of uh at the nodes
    val = uh(elem(t,:));

    % Values of uh at the midpoints of the edges
    valmid = 0.5 * [val(3) + val(2); val(1) + val(3); val(1) + val(2)];

    % Quadrature the integral
    intelem = area(t) / 3.0 * (valmid(1) + valmid(2) + valmid(3));
    int = int + intelem;
  end
end

function lengths = edgelengths(coords)
  edges = [coords(3,:) - coords(2,:); coords(3,:) - coords(1,:); coords(2,:) - coords(1,:)];
  lengths = [norm(edges(1,:)), norm(edges(2,:)), norm(edges(3,:))];
end

% Returns the three RT basis functions on a given triangle
function basis = basisRT(coords, area, signs)
  lengths = edgelengths(coords);
  % Construct the three basis functions
  basis = {@(x) signs(1) * lengths(1) / ( 2 * area) * (x - coords(1,:)),
          @(x) signs(2) * lengths(2) / ( 2 * area) * (x - coords(2,:)),
          @(x) signs(3) * lengths(3) / ( 2 * area) * (x - coords(3,:))};
end

% Calculates the inner products of the three local raviart
% thomas basis functions with signs on a given triangle
function M = localMassRT(coords, area, signs)
  basis = basisRT(coords, area, signs);
  % eveluate each basis function in the midpoints
  midpoints = 0.5 * [coords(3,:) + coords(2,:); coords(1,:) + coords(3,:); coords(1,:) + coords(2,:)];
  M = zeros(3,6);
  for i = 1:3
     M(i, :) = [basis{i}(midpoints(1,:)), basis{i}(midpoints(2,:)), basis{i}(midpoints(3,:))];
  end
  % calculate the final interactions
  M = area/ 3.0 * M * M';
%  M = zeros(3,3);
%  for i = 1:3
%    for j = 1:3
%      M(i,j) = quadmid(coords, @(x) dot( basis{i}(x) , basis{j}(x)));
%    end
%  end

end

% Calculates the inner products between Qh basis functions
% and divergence of RT basis functions on a given triangle
%
% Basis for Qh is given by the functions with |T_n| on x \in T_i
% and -|T_i| on x \in T_N
function B = localQhVh(volScale, coords, area, signs) 
  lengths = edgelengths(coords);
  B = volScale  * [signs(1) * lengths(1); signs(2) * lengths(2); signs(3) * lengths(3)];
end

% Calculates the RHS for the first set of equations,
% for a single triangle
% 
% Coordinates, signs of the edge, and patchedge which
% indicates which (LOCAL!) edges of the triangle arae
% connected to the patch vertex
function B = localRHS1(Duh, coords, area, signs, patchedge) 
  midpoints = 0.5 * [coords(3,:) + coords(2,:); coords(3,:) + coords(1,:); coords(2,:) + coords(1,:)];

  % Construct the three RT basis functions
  basis = basisRT(coords,area, signs);

  % Apply midpoint quadrature
  % hat_a dot(grad uh, \phi_j)
  % hat_a(midpoint) = 1/2 or 0. Depending on whether the edge is connected
  % to the patch vertex
  B = zeros(3,1);
  for j=1:3
    %f = @(x) dot(Duh, basis{j}(x));
    DuhT = Duh';
    int = area / 3.0 * (...
           patchedge(1) * 0.5 * basis{j}(midpoints(1,:)) * DuhT + ...
           patchedge(2) * 0.5 * basis{j}(midpoints(2,:)) * DuhT+ ...
           patchedge(3) * 0.5 * basis{j}(midpoints(3,:)) * DuhT);
    B(j) = int;
  end
end

% Calculates the RHS for the second set of equations,
% for a single triangle
% 
% On a single triangle 1 .. NT-1 only 1 basis function lives
%
% Duh is the gradient  of Uh
% Dphi_a is the gradient of phi_a
% f is the RHS
% patchedge indicates whether the edge is connected to the patchver


function int = localRHS2(Duh, Dphi_a, f, coords, area, volScale, patchedge) 
  edges = [coords(3,:) - coords(2,:); coords(3,:) - coords(1,:); coords(2,:) - coords(1,:)];
  midpoints = 0.5 * [coords(3,:) + coords(2,:); coords(1,:) + coords(3,:); coords(2,:) + coords(1,:)];
  Dprod = dot(Duh, Dphi_a);
  % Apply midpoint quadrature
  % hat_a(midpoint) = 1/2 or 0. Depending on whether the edge is connected
  % to the patch vertex
  int = volScale * area / 3.0 * (...
        patchedge(1) * 0.5 * f(midpoints(1,:)) - Dprod + ...
        patchedge(2) * 0.5 * f(midpoints(2,:)) - Dprod + ...
        patchedge(3) * 0.5 * f(midpoints(3,:)) - Dprod);
end

function A = localStiffRT(coords, signs)
  % coords are sorted counter clockwise
  % edge_i is the edge oppositve i
  edges = [coords(3,:) - coords(2,:); coords(1,:) - coords(3,:); coords(2,:) - coords(1,:)]
  % rotate 90 degrees clockwise
  normals = [edges(:,2), -edges(:,1)]
  % normalize
  normals(1,:) = 1.0 /  norm(normals(1,:)) * normals(1,:)
  normals(2,:) = 1.0 /  norm(normals(2,:)) * normals(2,:)
  normals(3,:) = 1.0 /  norm(normals(3,:)) * normals(3,:)
end

function sig = fluxpatch(T, node, elem,area, Dlambda,  Duh, f, patchvert, elemPatch, isBdPatch)
  % Number of edges (global DOFs)
  NE = size(T.edge,1);

  % Substructures of this patch
  NTP = size(elemPatch,1);

  dofQh = NTP - ~isBdPatch;

  % Assemble massmatrix M11 for zeta vs vh
  M11 = sparse(NE, NE);
  for k=1:NTP
    t = elemPatch(k);
    % Find edges and global signs belonging to this edge
    [edges,signs, ~] = elem2edge(t, patchvert, T);
    %Calculate localmass matrix
    Mt = localMassRT(node(elem(t,:),:),area(t), signs);
    %Add result to big matrix
    M11(edges, edges) = M11(edges, edges) + Mt;
  end

  % For calculation of mixed matrix M12 for rh vs dih vh
  M12 = sparse(NE, dofQh);


  % In case of a boundary patch, we don't have to rescale
  if (isBdPatch)
    areaSlave = 1;
  else
    % Calculate interactions on triangles 0..NTP-1 inside patch
    % Each triangle also induces an interaction on NTP, called the slave
    slaveElem = elemPatch(end);
    [edgesSlave, signsSlave, patchedgeSlave] = elem2edge(slaveElem, patchvert, T);
    coordsSlave = node(elem(slaveElem,:),:);
    areaSlave = area(slaveElem);
  end

  for k=1:dofQh
    t = elemPatch(k);

    % Calculate interactions on triangle t
    [edges,signs,~] = elem2edge(t, patchvert, T);
    Mt = localQhVh(areaSlave, node(elem(t,:),:), area(t), signs);

    % Add in the big matrix
    M12(edges, k) = M12(edges, k) + Mt;

    % In case of a boundary patch we don't have a slave
    if (isBdPatch)
      continue;
    end

    % Calculate interactions on slave with basis function t (has value -area(t))
    Mt = localQhVh(-area(t), coordsSlave, areaSlave, signsSlave);
    % Add to the big matrix
    M12(edgesSlave, k) = M12(edgesSlave,k) + Mt;
  end

  % Second set of equations, is just the tranpose of the latter
  M21 = M12';


  % Calculate the right hand sides, first set of eq
  %   b1 = (psi_a nabla uh, phi_j)
  b1 = zeros(NE,1);

  % First set of equations
  %   Again we assemble per triangle in patch
  for k=1:NTP
    t = elemPatch(k);

    % On each triangle there live three basis functions of RT
    [edges,signs, patchedges] = elem2edge(t, patchvert, T);

    % Calculate local in products on this triangle
    bt = localRHS1(Duh(t,:), node(elem(t,:),:), area(t), signs, patchedges);

    % Add to value in the big vector
    b1(edges) = b1(edges) + bt;
  end

  % Calculate the right hand sides, second set of eq
  %  b2 = (psi_a f - nabla psi_a nabla uh)
  b2 = zeros(dofQh, 1);

  if (~isBdPatch) 
    % Local coordinate of the patch vertex on the slave
    locpatchvertSlave = find(elem(slaveElem,:) == patchvert);
  end

  % Assemble second set of equations
  for k= 1:dofQh
    t = elemPatch(k);
    % Find out which edges of this patch are connected to the patch vert
    [~,~, patchedges] = elem2edge(t, patchvert, T);

    % local coordinate of the patchvert 
    locpatchvert = find(elem(t,:) == patchvert);

    % Calculate (g^a, x_t)
    bt = localRHS2(Duh(t,:), Dlambda(t,:, locpatchvert), f, node(elem(t,:),:),  area(t), areaSlave, patchedges);

    % Save
    b2(k) = bt;

    if (isBdPatch)
      continue;
    end

    % Calculate contribution on slavetriangle
    bt = localRHS2(Duh(slaveElem,:), Dlambda(slaveElem,:, locpatchvertSlave), f, coordsSlave, areaSlave,  -area(t), patchedgeSlave);

    b2(k) = b2(k) + bt;
  end

  % Construct a list of edge INDICES used in this patch
  indPatch = reshape(T.elem2edge(elemPatch, :), 3*NTP,1);
  edgePatch = accumarray(indPatch, ones(3*NTP,1));
  % edgePatch now holds the edge indices of our patch, with val == 2 if it is an interior edge


  % In case of an interior vertex, we have dof's at every interior edge
  % Otherwise we have dof's at every interior edge and at all edges that are also exterior edges
  if (~isBdPatch)
    dofEdge = find(edgePatch == 2);
  else
    intEdges = find(edgePatch == 2);
    extEdges = find(edgePatch == 1);
    
    isBdEdge = T.edge2elem(extEdges, 1) == T.edge2elem(extEdges,2);

    dofEdge = sort([intEdges; extEdges(isBdEdge)]);
  end



  % Create total matrix of the lhs M
  M = [M11(dofEdge, dofEdge), - M12(dofEdge,:);
       M21(:, dofEdge),        zeros(dofQh,dofQh)];      

  % Calculate big right hand side
  b = [-b1(dofEdge); b2];

  % Solve
  sol = M\b;

  %The first dofEdge coefficients define zeta
  sig = sparse(NE,1);
  sig(dofEdge) = sol(1:length(dofEdge));

  return
  % Print
  patchvert
  f
  printmat(node, 'nodes')
  printmat(elem, 'elements')
  printmat(uh, 'uh')
  printmat(Duh, 'Duh')
  printmat(T.edge, 'edges')
  printmat(elemPatch, 'elements in patch')
  printmat(elemPatch(1:dofQh), 'elements with dof for  Qh')
  printmat(dofEdge, 'edges with dof for Vh')
  printmat(full(M11), 'zeta vs vh')
  printmat(full(M12), 'rh vs div vh')
  printmat(full(M21), 'div zeta vs qh')
  printmat(b1, 'rhs 1')
  printmat(b2, 'rhs 2')

  printmat(full(M), 'Full matrix')
  printmat(b, 'Full RHS')
  printmat(sig, 'Solution')
end

% Calculate flux over entire domain
function sig = flux(node, elem, Duh, f)
  disp('Calculating flux')
  tic
  % Calculate auxiliary datastructure
  T = auxstructure(elem);

  % Number of nodes, number of edges, number of triangles
  N = size(node,1);
  NE = size(T.edge,1);
  NT = size(elem,1);

  % Fix the order of the triangulation
  [elem, ~, area] = fixorder(node, elem);

  % Calculate the boundary vertices
  [~, ~, isBdNode] = findboundary(elem);

  % Calculate gradient of basis
  [Dlambda, ~, ~] = gradbasis(node,elem);

  % NT \times N matrix
  t2v = sparse([1:NT, 1:NT, 1:NT], elem, 1, NT, N);

  % Crate solution
  sig = zeros(NE,1);
  for v=1:N
    % The indices of triangles inside this patch
    elemPatch = find(t2v(:,v));

    % Find solution for this patch
    isBdPatch = isBdNode(v);

    sigPatch = fluxpatch(T, node, elem, area, Dlambda,  Duh, f, v,elemPatch, isBdPatch);
    sig = sig + sigPatch;
  end
  toc
end

% Calculate the divergence of the flux on each element 
%   Note that the divergence is constant on each element
function divsig = divfluxelem(node, elem, flux) 
  T = auxstructure(elem);
  NT = size(elem,1);
  divsig = zeros(NT,1);
  area = simplexvolume(node, elem);

  for t=1:NT
    coords = node(elem(t,:),:);
    [edges,signs,~] = elem2edge(t, -1, T);

    lengths = edgelengths(coords);

    basis = 1.0 / area(t) * lengths .* signs;
    divsig(t) = dot(basis, flux(edges));
  end
end

% Returns an arrary of functions representing
%  the flux on each element
function sigelem =  fluxelem(node, elem, flux)
  T = auxstructure(elem);
  NT = size(elem,1);
  sigelem = cell(NT,1);
  area = simplexvolume(node, elem);

  for t=1:NT
    coords = node(elem(t,:),:);
    [edges,signs,~] = elem2edge(t, -1, T);
    basis = basisRT(coords, area(t), signs);

    % construct the function
    sigelem{t} = @(x) basis{1}(x) * flux(edges(1)) +  ...
                      basis{2}(x) * flux(edges(2)) +  ...
                      basis{3}(x) * flux(edges(3));
  end
end

% Approximate the norm \|zeta + Duh\|^2 using midpoint quadrature
function e = fluxerror(coords, area, signs, flux, duh)
    % Gather RT basis for this element
    basis = basisRT(coords, area, signs);
    midpoints = 0.5 * [coords(3,:) + coords(2,:); coords(1,:) + coords(3,:); coords(1,:) + coords(2,:)];

    % Evaluate flux at midpoints given as the linear combination of flux at edges
    M = zeros(1,6);
    for i = 1:3
      M = M + flux(i) * [basis{i}(midpoints(1,:)), basis{i}(midpoints(2,:)), basis{i}(midpoints(3,:))];
    end

    % Add the Duh
    M = M + [duh, duh, duh];

    % Result, edge-midpoint quadrature
    e = sqrt(area / 3.0 * M * M');
%    % Flux on this triangle is the linear combination, given by value of flux at edges
%
%    % construct the function
%    sigelem = @(x) basis{1}(x) * flux(1) +  ...
%                   basis{2}(x) * flux(2) +  ...
%                   basis{3}(x) * flux(3);
%
%    % Calculate (Duh + Sigma_h) inner product itself
%    error_func = @(x) dot(duh + sigelem(x), duh + sigelem(x));
%
%    % create norm
%    e = sqrt(quadmid(coords, area, error_func));
end

% Returns an upper bound for \|nabla(u - u_h)\|^2
function [toterr, errelem, oscelem] = equilresidualestimate(node, elem, duh, flux, f)
  tic
  disp('Equilibrated residual estimate')
  NT = size(elem,1);
  T = auxstructure(elem);
  area = simplexvolume(node, elem);
  %sigelem = fluxelem(node, elem, flux);
  divsigelem =  divfluxelem(node, elem, flux);
  errelem = zeros(NT,1);
  oscelem = zeros(NT,1);
  toterr = 0;
  for t=1:NT
    % Gather triangle specific information
    coords = node(elem(t,:),:);
    [edges,signs,~] = elem2edge(t, -1, T);

    % Calculate (Duh + Sigma_h) inner product itself
    %error_func = @(x) dot(duh(t,:) + sigelem{t}(x), duh(t,:) + sigelem{t}(x));

    % create norm
    %error_1 = sqrt(quadmid(coords, area(t), error_func))
    error_1 = fluxerror(coords, area(t), signs, flux(edges), duh(t,:));
    % error_t holds \|nabla uh + sigma_h\|_K

    % calculate h_K/ Pi \|f - div sigma\|_K
    h_K = max(edgelengths(coords));

    error_2 = sqrt(quadmid(coords,area(t), @(x) (f(x)-divsigelem(t))^2));

    % add contribution of this triangle to the total
    errelem(t) = error_1;
    oscelem(t) = error_2;

    toterr = toterr + (error_1 + h_K/pi * error_2)^2;
  end
  toterr = sqrt(toterr);
  toc
end

% The (squared) error indicator for every triangle (eta)
function etaellem =etaestimate(node, elem, duh, f)
  T = auxstructure(elem);
  NT = size(elem,1);
  etaellem = zeros(NT, 1);
  area = simplexvolume(node,elem);
  for t=1:NT
    % Calculate the jumps over edges

    % Find indices of edges, edges, and their normals
    edgeInd = T.elem2edge(t,:);
    edges = node(T.edge(edgeInd,2),:) - node(T.edge(edgeInd,1),:);
    lengths = [norm(edges(1,:));
               norm(edges(2,:));
               norm(edges(3,:))];
    normals = [edges(:,2), -edges(:,1)];

    % Find triangles on both sides belonging to these edges
    t1 = T.edge2elem(edgeInd,1);
    t2 = T.edge2elem(edgeInd,2);

    % Calculate jump^2 on each edge (note that t1=t2 in case of boundary!)
    jump = zeros(3,1);
    for e=1:3
      % Notice that norm(normals(e,:)) = vol(e)
      jump(e) = 1.0/lengths(e) * dot(duh(t1(e),:) -duh(t2(e),:), normals(e,:))^2;
    end

    % Calculate erros
    %h_t = diameter = max length
    h = max(lengths);
    err_a = h^2 * quadmid(node(elem(t,:),:), area(t), @(x) f(x)^2);
    err_j = h* sum(jump);

    etaellem(t) = err_a + err_j;
  end
end

% The (squared) oscillation for every triangle
function oscelem = oscestimate(node, elem, f)
  NT = size(elem, 1);
  oscelem = zeros(NT, 1);
  area = simplexvolume(node, elem);

  for t =1:NT
    coords = node(elem(t,:),:);
    % Calculate the projection of f onto the constants (is just average)
    int_f = quadmid(coords, area(t), f);
    c = 1.0  / area(t) * int_f;

    % Calculate the oscilation term
    oscelem(t) = max(edgelengths(coords)) * quadmid(coords, area(t), @(x) (f(x) - c)^2);
  end
end

% Classical residual error estimate (the total error): eta + osc
function [toterr, errelem] = residualesimate(node, elem, duh, f)
  errelem = etaestimate(node, elem, duh, f) + oscestimate(node, elem, f);
  toterr = sum(errelem);
end

% Results vector with (div flux, 1)_K- (f,1)_k for all K \in Th
function result = fluxequil(node, elem, flux, f)
  divsig = divfluxelem(node, elem, flux);
  area = simplexvolume(node,elem);
  NT = size(elem,1);

  result = zeros(NT,1);
  for t=1:NT
    int_f = quadmid(node(elem(t,:),:), f);
    int_flux = area(t) * divsig(t);
    result(t) = int_flux - int_f;

  end
end

% Results vector with [[nabla u_h]] - [[sigma]] over all the interior edges
function result = fluxjump(node, elem, Duh, flux)
  T = auxstructure(elem);
  NE = size(T.edge,1);
  result = zeros(NE, 1);
  bdEdge = T.edge2elem(:,1) == T.edge2elem(:,2);
  for e=1:NE
    if (bdEdge(e))
      continue;
    end

    t1 = T.edge2elem(e,1);
    t2 = T.edge2elem(e,2);
    edge = node(T.edge(e,1),:) - node(T.edge(e,2),:);
    normal = [edge(:,2), -edge(:,1)];
    % normalize
    normal = normal * 1.0/ norm(normal);

    % Calculate first jump
    j1 = dot(Duh(t1,:), normal) - dot(Duh(t2,:), normal)
    % Calculate second jump
    flux(e)
    result(e) = j1 - flux(e);
  end
end

function vismesh(node, elem)
  T = auxstructure(elem);
  figure(1);
  subplot(1,3,1);
  showmesh(node, elem);
  findedge(node, T.edge);
  subplot(1,3,2);
  showmesh(node, elem);
  findnode(node);
  subplot(1,3,3);
  showmesh(node, elem);
  findelem(node, elem);
end

function [node, elem, bdFlag, pde, Du, theorate] = lshapeone() 
  global refinemethod;
  %%  Generate an initial mesh
  [node,elem] = squaremesh([-1,1,-1,1],1);
  [node,elem] = delmesh(node,elem,'x<0 & y<0');
  bdFlag = setboundary(node,elem,'Dirichlet');
  %node = [-1,0;-1,1;0,1; 1,1;1,0;1,-1;0,-1; 0,0];
  %elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
  %bdFlag = setboundary(node,elem,'Dirichlet');
  elem = fixorientation(node,elem);   % counter-clockwise oritentation
  %[node, elem, bdFlag ] = refinemethod(node, elem, bdFlag);
  %% Set up PDE data
  pde.f = @(p)  ones(size(p,1),1);
  pde.g_D = 0;
  pde.theorate = 1.0/3.0;
  pde.Du = [];
end

function [node, elem, bdFlag, pde] = lshapecorner() 
  [node, elem] = squaremesh([-1,1,-1,1],1);
  [node, elem] = delmesh(node, elem, 'x>0& y<0');
  bdFlag = setboundary(node,elem,'Dirichlet');
  elem = fixorientation(node,elem);   % counter-clockwise oritentation
  % Symbolic symbls
  syms r t x y;
  p = sym('p', 1:2);

  r = sqrt(x^2 + y^2);
  t = atan(y/(x+0.0000001));

  % Symbolic functions
  u = (x+1)*x*(x-1)*(y+1)*y*(y-1) * r^(2/3) * sin(2/3*t);
  lap(x,y) =  (diff(u,x,2) + diff(u,y,2));
  grad(x,y) = [diff(u,x,1), diff(u,y,1)];

  % Convert to matlab shizzle
  pde.f = matlabFunction(-lap(p(1), p(2)), 'vars', {p});
  pde.Du = matlabFunction(grad(p(1), p(2)), 'vars' , {p});
  pde.g_D = 0;
end

% u = 2^{4a} x^a(1-x)^a y^a (1-y)^a.
function [node, elem, bdFlag, pde] = squareana(a)
  [node, elem] = squaremesh([0,1,0,1],0.5);
  [node, elem] = squaremesh([0,1,0,1],1);
  bdFlag = setboundary(node, elem, 'Dirichlet');
  elem = fixorientation(node,elem);   % counter-clockwise oritentation

  [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
  % Symbolic symbls
  syms x y;
  p = sym('p', 1:2);

  % Symbolic functions
  u(x,y) = 2^(4*a)*x^a*(1-x)^a* y^a* (1-y)^a;
  lap(x,y) =  (diff(u,x,2) + diff(u,y,2));
  grad(x,y) = [diff(u,x,1), diff(u,y,1)];

  % Convert to matlab shizzle
  pde.f = matlabFunction(-lap(p(1), p(2)), 'vars', {p});
  pde.Du = matlabFunction(grad(p(1), p(2)), 'vars' , {p});
  pde.g_D = 0;
  pde.theorate = 0.5;
end

% u = x(x-1)y(y-1)exp(-a((x-x_c)^2 + (y-y_c)^2))
%   Peak centered at (x_c, y_c) with intensity a
function [node, elem, bdFlag, pde] = squarepeak(a,x_c, y_c)
  [node, elem] = squaremesh([0,1,0,1],0.5);
  [node, elem] = squaremesh([0,1,0,1],1);
  bdFlag = setboundary(node, elem, 'Dirichlet');
  elem = fixorientation(node,elem);   % counter-clockwise oritentation
  % Symbolic symbls
  syms x y;
  p = sym('p', 1:2);

  % Symbolic functions
  u(x,y) = x*(x-1)*y*(y-1)*exp(-a*((x-x_c)^2 + (y-y_c)^2));
  lap(x,y) =  (diff(u,x,2) + diff(u,y,2));
  grad(x,y) = [diff(u,x,1), diff(u,y,1)];

  % Convert to matlab shizzle
  pde.f = matlabFunction(-lap(p(1), p(2)), 'vars', {p});
  pde.Du = matlabFunction(grad(p(1), p(2)), 'vars' , {p});

  pde.g_D = 0;
end

function [node, elem, bdFlag, pde] = squaresin() 
  % Exact derivative on squaremesh
  function z = DuSquare(p)
    x = p(:,1); y = p(:,2);
    z(:,1) = 2*pi*cos(2*pi*x).*sin(2*pi*y);
    z(:,2) = 2*pi*sin(2*pi*x).*cos(2*pi*y);
  end
  %%  Generate an initial mesh
  node = [[0,0]; [0.5, 0.5]; [1,0];[1,1];[0,1]];
  elem = [[1,2,3]; [3,2,4]; [4,2,5]; [2,1,5]];
  [node, elem] = squaremesh([0,1,0,1],1);
  [elem, ~, ~] = fixorder(node, elem);
  bdFlag = setboundary(node,elem,'Dirichlet');

  [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
  %pde.f = @(p) 2*p(:,1) + p(:,2);
  pde.f = @(p) 8*pi^2*sin(2*pi*p(:,1)).*sin(2*pi*p(:,2));
  %pde.f = @(p) 2*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2));
  pde.g_D = 0;
  pde.Du = @DuSquare;
  pde.theorate = 0.5;
end

function [node, elem, bdFlag, pde] = squareone() 
  % Exact derivative on squaremesh
  function z = DuSquare(p)
    x = p(:,1); y = p(:,2);
    z(:,1) = 2*pi*cos(2*pi*x).*sin(2*pi*y);
    z(:,2) = 2*pi*sin(2*pi*x).*cos(2*pi*y);
  end
  %%  Generate an initial mesh
  node = [[0,0]; [0.5, 0.5]; [1,0];[1,1];[0,1]];
  elem = [[1,2,3]; [3,2,4]; [4,2,5]; [2,1,5]];
  [node, elem] = squaremesh([0,1,0,1],1);
  [elem, ~, ~] = fixorder(node, elem);
  bdFlag = setboundary(node,elem,'Dirichlet');
  %pde.f = @(p) 2*p(:,1) + p(:,2);
  pde.f = @(p) 8*pi^2*sin(2*pi*p(:,1)).*sin(2*pi*p(:,2));
  %pde.f = @(p) 2*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2));
  pde.g_D = 0;
  pde.Du = @DuSquare;
end

function [node, elem, bdFlag, pde] = crackone() 
  %%  Generate an initial mesh
  node = [1,0; 0,1; -1,0; 0,-1; 0,0; 1,0];        % nodes
  elem = [5,1,2; 5,2,3; 5,3,4; 5,4,6];            % elements
  % Node 6 == Node 5, to ensure the crack property
  elem = fixorientation(node, elem);
  bdFlag = setboundary(node,elem,'Dirichlet');    % Dirichlet boundary condition
  %[node, elem, bdFlag ] = refinemethod(node, elem, bdFlag);
  %% Set up PDE data
  pde.f = @(p)  ones(size(p,1),1);
  pde.g_D = 0;
  pde.theorate = 1.0/4.0;
  pde.Du = [];
end
