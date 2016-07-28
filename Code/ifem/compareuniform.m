% Calculates various estimators for different estimators and plots
function [N,equilError, resError, mixedError, errH1] = compareuniform(method, estimators, node, elem, pde, bdFlag, maxN)
  savedir = 'figures';

  % Create new figures 
  f2 = figure(2);clf; hold on;% absolute errors
  f3 = figure(3);clf; hold on;% relative errors

  % Iteratively create legends
  legendAbs = {};
  legendRel = {};

  % First, calculate the exact error and plot
  [N, errH1] = EstimateUniform('exact', node, elem, pde, bdFlag, maxN)
  figure(f2); loglog(N,  errH1, 'b-o'); 
  legendAbs{end + 1} = '$\|\nabla{\tilde U_k} - \nabla {U_k}\|$';

  % Calculate & plot the result for the various estimators
  for estimator = estimators

    % estimate this thingy
    disp(estimator{1})
    [~, err] = EstimateUniform(estimator{1}, node, elem, pde, bdFlag, maxN)

    % Plot it to a figure, depending on the estimator used
    switch estimator{1}
    case 'residual'
      legendAbs{end+1} = 'residual($U_k$)';
      figure(f2); loglog(N, err, 'r-o'); % abs
      figure(f3); semilogx(N, errH1 ./ err, 'r-o'); % rel
    case 'mixed'
      legendAbs{end+1} = 'mixed($U_k, \sigma$)';
      figure(f2); loglog(N, err, 'm-o'); % abs
      figure(f3); semilogx(N, errH1 ./ err, 'm-o'); % rel
    case 'equilibrated'
      legendAbs{end+1} = 'equilibrated($U_k, \zeta$)';
      figure(f2); loglog(N, err, '-o', 'color', [0 0.5 0]); % abs
      figure(f3); semilogx(N, errH1 ./ err, '-o', 'color', [0 0.5 0]); % rel
    case 'zz'
      legendAbs{end+1} = 'ZZ($U_k$)';
      figure(f2); loglog(N, err, 'k-o'); %abs
      figure(f3); semilogx(N, errH1 ./ err, 'k-o'); %rel
    end
    % Same for every estimator
    legendRel{end+1} = legendAbs{end};
  end

  % add theoretical convergence triangle
  len = max(size(N))
  if (len > 2)
    (log(errH1(end)) - log(errH1(end -1)))/ (log(N(end)) - log(N(end-1)))
    x1 = N(end-2); x2 = N(end-1);
    y12 = errH1(end-1); y11 = exp(log(y12) - pde.theorate);

    figure(f2); loglog([x1 x1 x2 x1], [y12 y11 y11 y12], 'color', [0,0,0]);
  end

  % add legends and shit for absolute figure
  figure(f2);
  title('Comparison error estimators');
  xlabel('Number of vertices');
  ylabel('Error');
  legend(legendAbs, 'interpreter', 'latex')

  % Save, todo
  %saveas(f2, sprintf('%s/%s/norm_slope_%d.png',savedir, method, len));


  % add legends and shit for relative figure
  f3 = figure(f3); 
  title('Efficiency index');
  xlabel('Number of vertices');
  legend(legendRel, 'interpreter', 'latex', 'location', 'best');
  
  % save, todo
  %saveas(f3, sprintf('%s/%s/efficiency_%d.fig',savedir, method, len));
end
function [N, err] = EstimateUniform(estimator, node, elem, pde, bdFlag, maxN)
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
