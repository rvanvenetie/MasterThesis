% Calculates various estimators for different estimators and plots
function compareuniform(method, estimators, node, elem, pde, bdFlag, maxN)
  savedir = 'figures';

  % Create new figures 
  f2 = figure(2);clf; hold on;% absolute errors
  f3 = figure(3);clf; hold on;% relative errors

  % Iteratively create legends
  legendAbs = {};
  legendRel = {};

  % First, calculate the exact error and plot
  [N, errH1] = estimateuniform('exact', node, elem, pde, bdFlag, maxN)
  figure(f2); loglog(N,  errH1, 'b-o'); 
  legendAbs{end + 1} = '$\|\nabla{U_{\star,k}} - \nabla {U_k}\|$';

  % Calculate & plot the result for the various estimators
  for estimator = estimators
    if strcmp(estimator{1}, 'exact')
      continue
    end

    % estimate this thingy
    disp(estimator{1})
    [~, err] = estimateuniform(estimator{1}, node, elem, pde, bdFlag, maxN)

    % Plot it to a figure, depending on the estimator used
    switch estimator{1}
    case 'residual'
      legendAbs{end+1} = 'residual($U_k$)';
      figure(f2); loglog(N, err, 'r-o'); % abs
      figure(f3); semilogx(N, err ./ errH1, 'r-o'); % rel
    case 'mixed'
      legendAbs{end+1} = 'mixed($U_k, \sigma$)';
      figure(f2); loglog(N, err, 'm-o'); % abs
      figure(f3); semilogx(N, err ./errH1, 'm-o'); % rel
    case 'equilibrated'
      legendAbs{end+1} = 'equilibrated($U_k, \zeta$)';
      figure(f2); loglog(N, err, '-o', 'color', [0 0.5 0]); % abs
      figure(f3); semilogx(N, err ./errH1, '-o', 'color', [0 0.5 0]); % rel
    case 'zz'
      legendAbs{end+1} = 'ZZ($U_k$)';
      figure(f2); loglog(N, err, 'k-o'); %abs
      figure(f3); semilogx(N, err ./ errH1, 'k-o'); %rel
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
  export_fig(f2, sprintf('%s/%s/norm_slope_%d',savedir, method, len),'-pdf', '-painters'); % Export to PDF using export_fig


  % add legends and shit for relative figure
  f3 = figure(f3); 
  title('Efficiency index');
  xlabel('Number of vertices');
  legend(legendRel, 'interpreter', 'latex', 'location', 'best');
  
  % save, todo
  %saveas(f3, sprintf('%s/%s/efficiency_%d.fig',savedir, method, len));
  export_fig(f3, sprintf('%s/%s/efficiency_%d',savedir, method, len),'-pdf', '-painters'); % Export to PDF using export_fig
end
