function test
  close all; 
  profile on;
  %% Parameters
  theta = 0.5;    uh =0;

  global refinemethod;
  refinemethod = @uniformrefine;

  maxN =  10000;
  maxIt = 20;
  %%  Generate an initial mesh
  %[node, elem, bdFlag, pde, Du] = squaresin();
  %[node, elem, bdFlag, pde, Du] = lshapeone();
  %[node, elem, bdFlag, pde, Du] = lshapecorner();
  [node, elem, bdFlag, pde, Du] = squareana(10);
  %[node, elem, bdFlag, pde, Du] = squarepeak(10, 0.51, 0.117);
  errH1 = zeros(maxIt,1);
  equilError = zeros(maxIt,1);
  resError = zeros(maxIt,1);
  N = zeros(maxIt, 1);
  N(1) = size(node, 1);
  % Uniform refinement
  t = 1;
  while  (t <= maxIt) && (N(t) < maxN) 
    N(t) = size(node,1);
    uh = Poisson(node, elem, pde, bdFlag);
    figure(1);  showresult(node,elem,uh,[-50,12]);    
    [Duh,~] = gradu(node, elem, uh);
    if isreal(Du)
      if Du == 0
        errH1(t) = 0;
      else
        % Du should hold \|Nabla u\|^2_{Omega}
        errH1(t) = Du - intdomain(node, elem, uh);
      end
    else
      errH1(t) = getH1error(node,elem,Du,uh)^2;
    end
    sig = flux(node,elem, uh, Duh, pde.f);
    [equilError(t), eta] =  equilresidualestimate(node, elem, Duh, sig, pde.f);
    [resError(t), ~] = residualesimate(node, elem, Duh, pde.f);
    markedElem = mark(elem,eta,theta);
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
    %[node,elem, bdFlag] = refinemethod(node, elem, bdFlag);
    t = t+1;
  end
  t = t-1;
  equilError
  resError
  errH1
  figure(2);
  N = N(1:t); errH1 = errH1(1:t); equilError = equilError(1:t); resError = resError(1:t);
  x = (1:t) - 1;
  loglog(N, errH1, 'r-',N, resError, 'g-', N, equilError, 'b-')
  legend({'$\|\nabla{(u - u_h)}\|^2$','res($uh$)', 'equil($u_h, \sigma_h$)'}, 'interpreter', 'latex');
  title('Comparison equilibrated error estimator');
  xlabel('Number of degrees of freedom');
  ylabel('Error');
  figure(3);
  loglog(N, sqrt(equilError) ./ sqrt(errH1));
  title('Efficiency index');
  legend({'$C_{eff}$'}, 'interpreter', 'latex');
  profile viewer
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

function int = quadmid(nodes, f) % apply midpoint quadrature on triangle with nodes
  area = polyarea(nodes(:, 1), nodes(:,2));
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
function basis = basisRT(coords, signs)
  lengths = edgelengths(coords);
  area = polyarea(coords(:,1), coords(:,2));
  % Construct the three basis functions
  basis = {@(x) signs(1) * lengths(1) / ( 2 * area) * (x - coords(1,:)),
          @(x) signs(2) * lengths(2) / ( 2 * area) * (x - coords(2,:)),
          @(x) signs(3) * lengths(3) / ( 2 * area) * (x - coords(3,:))};
end

% Calculates the inner products of the three local raviart
% thomas basis functions with signs on a given triangle
function M = localMassRT(coords, signs)
  basis = basisRT(coords, signs);
  M = zeros(3,3);
  for i = 1:3
    for j = 1:3
      M(i,j) = quadmid(coords, @(x) dot( basis{i}(x) , basis{j}(x)));
    end
  end
end

% Calculates the inner products between Qh basis functions
% and divergence of RT basis functions on a given triangle
%
% Basis for Qh is given by the functions with |T_n| on x \in T_i
% and -|T_i| on x \in T_N
function B = localQhVh(volScale, coords, signs) 
  lengths = edgelengths(coords);
  area = polyarea(coords(:,1), coords(:,2));
  B = volScale  * [signs(1) * lengths(1); signs(2) * lengths(2); signs(3) * lengths(3)];
end

% Calculates the RHS for the first set of equations,
% for a single triangle
% 
% Coordinates, signs of the edge, and patchedge which
% indicates which (LOCAL!) edges of the triangle arae
% connected to the patch vertex
function B = localRHS1(Duh, coords, signs, patchedge) 
  area = polyarea(coords(:,1), coords(:,2));
  midpoints = 0.5 * [coords(3,:) + coords(2,:); coords(3,:) + coords(1,:); coords(2,:) + coords(1,:)];

  % Construct the three RT basis functions
  basis = basisRT(coords, signs);

  % Apply midpoint quadrature
  % hat_a dot(grad uh, \phi_j)
  % hat_a(midpoint) = 1/2 or 0. Depending on whether the edge is connected
  % to the patch vertex
  B = zeros(3,1);
  for j=1:3
    f = @(x) dot(Duh, basis{j}(x));
    int = area / 3.0 * (...
           patchedge(1) * 0.5 * f(midpoints(1,:)) + ...
           patchedge(2) * 0.5 * f(midpoints(2,:)) + ...
           patchedge(3) * 0.5 * f(midpoints(3,:)));
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


function int = localRHS2(Duh, Dphi_a, f, coords, volScale, patchedge) 
  edges = [coords(3,:) - coords(2,:); coords(3,:) - coords(1,:); coords(2,:) - coords(1,:)];
  midpoints = 0.5 * [coords(3,:) + coords(2,:); coords(1,:) + coords(3,:); coords(2,:) + coords(1,:)];
  area = polyarea(coords(:,1), coords(:,2));
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

function sig = fluxpatch(node, elem,uh, Duh, f, patchvert)
  [~, ~, isBdNode] = findboundary(elem);
  isBdPatch = isBdNode(patchvert);
  T = auxstructure(elem);
  N = size(node,1);
  NT = size(elem, 1);
  NE = size(T.edge,1);

  [elem, ~, area] = fixorder(node, elem);
  % Calculate gradient of the barycentric basis
  [Dlambda, ~, ~] = gradbasis(node,elem);


  % NT \times N matrix
  t2v = sparse([1:NT, 1:NT, 1:NT], elem, 1, NT, N);
  % The indices of triangles inside this patch
  elemPatch = find(t2v(:,patchvert));

  % Substructures of this patch
  TPatch = auxstructure(elem(elemPatch,:));
  NTP = size(elemPatch,1);
  NEP = size(TPatch.edge,1);
  dofQh = NTP - ~isBdPatch;

  % Assemble massmatrix M11 for zeta vs vh
  M11 = sparse(NE, NE);
  for k=1:NTP
    t = elemPatch(k);
    % Find edges and global signs belonging to this edge
    [edges,signs, ~] = elem2edge(t, patchvert, T);
    %Calculate localmass matrix
    Mt = localMassRT(node(elem(t,:),:),signs);
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
    Mt = localQhVh(areaSlave, node(elem(t,:),:), signs);

    % Add in the big matrix
    M12(edges, k) = M12(edges, k) + Mt;

    % In case of a boundary patch we don't have a slave
    if (isBdPatch)
      continue;
    end

    % Calculate interactions on slave with basis function t (has value -area(t))
    Mt = localQhVh(-area(t), coordsSlave, signsSlave);
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
    bt = localRHS1(Duh(t,:), node(elem(t,:),:), signs, patchedges);

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
    bt = localRHS2(Duh(t,:), Dlambda(t,:, locpatchvert), f, node(elem(t,:),:), areaSlave, patchedges);

    % Save
    b2(k) = bt;

    if (isBdPatch)
      continue;
    end

    % Calculate contribution on slavetriangle
    bt = localRHS2(Duh(slaveElem,:), Dlambda(slaveElem,:, locpatchvertSlave), f, coordsSlave,  -area(t), patchedgeSlave);

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
function sig = flux(node, elem,uh, Duh, f)
  T = auxstructure(elem);
  N = size(node,1);
  NE = size(T.edge,1);
  % Crate solution
  sig = zeros(NE,1);
  for v=1:N
    % Find solution for this patch
    sigPatch = fluxpatch(node, elem, uh, Duh, f, v);
    sig = sig + sigPatch;
  end
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
  for t=1:NT
    coords = node(elem(t,:),:);
    [edges,signs,~] = elem2edge(t, -1, T);
    basis = basisRT(coords, signs);

    % construct the function
    sigelem{t} = @(x) basis{1}(x) * flux(edges(1)) +  ...
                      basis{2}(x) * flux(edges(2)) +  ...
                      basis{3}(x) * flux(edges(3));
  end
end


% Returns an upper bound for \|nabla(u - u_h)\|^2
function [toterr, errelem] = equilresidualestimate(node, elem, duh, flux, f)
  NT = size(elem,1);
  sigelem = fluxelem(node, elem, flux);
  divsigelem =  divfluxelem(node, elem, flux);
  errelem = zeros(NT,1);
  for t=1:NT
    coords = node(elem(t,:),:);
    % Calculate (Duh + Sigma_h) inner product itself
    error_func = @(x) dot(duh(t,:) + sigelem{t}(x), duh(t,:) + sigelem{t}(x));

    % create norm
    error_1 = sqrt(quadmid(coords, error_func));
    % error_t holds \|nabla uh + sigma_h\|_K

    % calculate h_K/ Pi \|f - div sigma\|_K
    h_K = max(edgelengths(coords));

    error_2 = h_K/pi*sqrt(quadmid(coords,@(x) (f(x)-divsigelem(t))^2));

    % add contribution of this triangle to the total
    errelem(t) = (error_1 + error_2)^2;
  end
  toterr = sum(errelem);
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
    err_a = h^2 * quadmid(node(elem(t,:),:), @(x) f(x)^2);
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
    int_f = quadmid(coords, f);
    c = 1.0  / area(t) * int_f;

    % Calculate the oscilation term
    oscelem(t) = max(edgelengths(coords)) * quadmid(coords, @(x) (f(x) - c)^2);
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

function [node, elem, bdFlag, pde, Du] = lshapeone() 
  global refinemethod;
  %%  Generate an initial mesh
  node = [-1,0;-1,1;0,1; 1,1;1,0;1,-1;0,-1; 0,0];
  elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];    % elements
  bdFlag = setboundary(node,elem,'Dirichlet');
  elem = fixorientation(node,elem);   % counter-clockwise oritentation
  %[node, elem, bdFlag ] = refinemethod(node, elem, bdFlag);
  %% Set up PDE data
  pde.f = @(p) 1;
  pde.g_D = 0;
  Du = 0.2140758036240825;
end

function [node, elem, bdFlag, pde, Du] = lshapecorner() 
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
  Du = matlabFunction(grad(p(1), p(2)), 'vars' , {p});

  pde.g_D = 0;
end

% u = 2^{4a} x^a(1-x)^a y^a (1-y)^a.
function [node, elem, bdFlag, pde, Du] = squareana(a)
  [node, elem] = squaremesh([0,1,0,1],0.5);
  bdFlag = setboundary(node, elem, 'Dirichlet');
  elem = fixorientation(node,elem);   % counter-clockwise oritentation
  % Symbolic symbls
  syms x y;
  p = sym('p', 1:2);

  % Symbolic functions
  u(x,y) = 2^(4*a)*x^a*(1-x)^a* y^a* (1-y)^a;
  lap(x,y) =  (diff(u,x,2) + diff(u,y,2));
  grad(x,y) = [diff(u,x,1), diff(u,y,1)];

  % Convert to matlab shizzle
  pde.f = matlabFunction(-lap(p(1), p(2)), 'vars', {p});
  Du = matlabFunction(grad(p(1), p(2)), 'vars' , {p});

  pde.g_D = 0;
end

% u = x(x-1)y(y-1)exp(-a((x-x_c)^2 + (y-y_c)^2))
%   Peak centered at (x_c, y_c) with intensity a
function [node, elem, bdFlag, pde, Du] = squarepeak(a,x_c, y_c)
  [node, elem] = squaremesh([0,1,0,1],0.5);
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
  Du = matlabFunction(grad(p(1), p(2)), 'vars' , {p});

  pde.g_D = 0;
end

function [node, elem, bdFlag, pde, Du] = squaresin() 
  % Exact derivative on squaremesh
  function z = DuSquare(p)
    x = p(:,1); y = p(:,2);
    z(:,1) = pi*cos(pi*x).*sin(pi*y);
    z(:,2) = pi*sin(pi*x).*cos(pi*y);
  end
  %%  Generate an initial mesh
  node = [[0,0]; [0.5, 0.5]; [1,0];[1,1];[0,1]];
  elem = [[1,2,3]; [3,2,4]; [4,2,5]; [2,1,5]];
  [elem, ~, ~] = fixorder(node, elem);
  bdFlag = setboundary(node,elem,'Dirichlet');
  %pde.f = @(p) 2*p(:,1) + p(:,2);
  pde.f = @(p) 2*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2));
  pde.g_D = 0;
  Du = @DuSquare;
end