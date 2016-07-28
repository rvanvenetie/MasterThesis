function err = exactH1error(node, elem, bdFlag, pde, uh)
  if isempty(pde.Du)
    err = approxH1error(node, elem, bdFlag, pde, uh,3);
  else
    err = getH1error(node,elem,pde.Du,uh);
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
