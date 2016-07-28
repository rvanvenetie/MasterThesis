function   test(method, estimators, afem)
  savedir = 'figures';
  %clear all;
  %close all; 
  %profile on;
  %% Parameters
  theta = 0.3;    uh =0;

  global refinemethod;
  refinemethod = @uniformrefine;

  maxN =  1e4;
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
  %plotuh(node, elem, pde, bdFlag);
  %return;

  % AFEM suffix
  if afem
    return
    method =sprintf('%s_afem', method);
  end

  % Estimators suffix
  method = sprintf('%s_',method);
  for estimator = estimators
    % append the first letter to the method
    method = sprintf('%s%s',method, estimator{1}(1));
  end

  disp(method)

  % Make dir
  mkdir(sprintf('%s/%s', savedir, method)); clf;

  showmesh(node, elem);
  export_fig(gca, sprintf('%s/%s/mesh_initial.pdf',savedir, method), '-painters');

  %plotapproxh1(method,node, elem, pde, bdFlag,pde.Du, 10)
  %return
  % one uniform refinement
  [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
  if afem
    CompareAfem(method, node, elem, pde, bdFlag, theta,maxN)
  else
     compareuniform(method, estimators, node, elem, pde, bdFlag, maxN);
  end
  %profile viewer
end



%function [N,equilError, resError, mixedError, errH1] = CompareUniform(method, node, elem, pde, bdFlag, maxN)
%  errH1 = [];
%  mixedError = [];
%  equilError = [];
%  resError = [];
%  N = [];
%  % Uniform refinement
%  t = 1;
%  while  (size(node, 1) < maxN)
%    t = t+1;
%    N(end+1) = size(node,1);
%
%    % Solve the poisson problem
%    uh = Poisson(node, elem, pde, bdFlag);
%    [Duh,~] = gradu(node, elem, uh);
%
%    % Calculate real error
%    if isempty(pde.Du)
%      errH1(end+1) = approxH1error(node, elem, bdFlag, pde, uh,3)
%    else
%      errH1(end+1) = getH1error(node,elem,pde.Du,uh)
%    end
%    [~,sigma] = PoissonRT0(node, elem, pde, bdFlag);
%    mixedError(end+1) = getL2errorRT0(node, elem,  Duh, sigma);
%
%    size(elem)
%
%    % Calculate the flux
%    sig = flux(node,elem,  Duh, pde.f);
%
%    %Calculate the error esimate
%    [equilError(end+1), ~, ~] =  equil(node, elem, Duh, sig, pde.f);
%
%    %Calculate the residual eror
%    eta = estimateresidual(node, elem, uh, pde);
%    resError(end+1) = sqrt(sum(eta.^2));
%    %ploterror(method,node, elem, uh, N,equilError, resError,mixedError, errH1, pde.theorate);
%
%    [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
%  end
%end

function plotuh(node, elem, pde, bdFlag)
  savedir = 'figures/uh';
  for j =1:8
    [node,elem, bdFlag] = uniformbisect(node, elem);
    uh = Poisson(node, elem, pde, bdFlag);
    %f1 = figure(1); showmesh(node,elem);
    if 0
      f2 = figure(2); %showsolution(node, elem, uh,2); colorbar;
      clf;
      colormap('jet')
      trisurf(elem, node(:,1), node(:,2), uh', 'FaceColor', 'interp', 'EdgeColor', 'k');
      caxis([-1,1])
      h = colorbar;
      view(2)
      export_fig(f2, sprintf('%s/result%d.pdf', savedir, j), '-painters');
    end
  end
  f2 = figure(3);
  clf;
  colormap('jet')
  trisurf(elem, node(:,1), node(:,2), uh', 'FaceColor', 'interp', 'EdgeColor', 'interp');
  caxis([-1,1])
  h = colorbar;
  view(3)
  set(f2, 'color', 'none');
  export_fig(f2, sprintf('%s/result.png', savedir));
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
