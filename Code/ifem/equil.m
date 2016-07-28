% Returns an upper bound for \|nabla(u - u_h)\|^2
function [toterr, errelem, oscelem] = equil(node, elem, duh,  f)
  % Calculate the flux
  flux = flux(node,elem,  duh, f);
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
