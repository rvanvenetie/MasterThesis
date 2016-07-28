function [node, elem, bdFlag,pde] = example(method)
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
