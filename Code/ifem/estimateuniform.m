function [N, err] = estimateuniform(estimator, node, elem, pde, bdFlag, maxN)
  err = [];
  N = [];
  % Uniform refinement
  t = 1;
  while  (size(node, 1) < maxN)
    t = t+1;
    N(end+1) = size(node,1);

    % Solve the poisson problem
    uh = Poisson(node, elem, pde, bdFlag);
    [Duh,area] = gradu(node, elem, uh);

    % Calculate the error according to the estimator type
    switch estimator
    case 'exact'
      % Calculate real error
      err(end+1) =  exactH1error(node, elem, bdFlag, pde, uh)
    case 'mixed'
      % Compute the mixed FEM estimator
      [~,sigma] = PoissonRT0(node, elem, pde, bdFlag)
      err(end+1) = getL2errorRT0(node, elem,  Duh, sigma)
    case 'equilibrated'
      % Calculate the equilibrated flux estimator
      %sig = flux(node,elem,  Duh, pde.f);

      %Calculate the error esimate
      [erre,eta, ~] =  equil(node, elem, Duh, pde.f)
      err(end+1) = erre % sqrt(sum(eta.^2))
    case 'residual'
      %Calculate the residual eror
      eta = estimateresidual(node, elem, uh, pde);
      err(end+1) = sqrt(sum(eta.^2))
    case 'zz'
      % Calculate the zienkewicz-zhu estimator
      err(end+1) = zzestimate(node, elem, Duh, area, pde.f)
    otherwise
      disp('You dumb boy, select a valid estimator')
      return
    end
    %ploterror(method,node, elem, uh, N,equilError, resError,mixedError, errH1, pde.theorate);
    [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
  end
end
