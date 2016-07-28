function efficiencyZZ
  maxN = 2e3
  meths = {'square_sin', 'square_ana', 'lshape_one' , 'crack_one'};
  colors = {'r', 'b', [0 0.5 0], 'm'};
  markers = ['o','x','s','d'];
  handles = [0,0,0,0];
  fig = figure(6); clf;
  for i=1:4
    % Retreive the example data
    [node, elem, bdFlag, pde] = examples(meths{i});

    % Calculate the exact estimator
    [N, errH1] = estimateuniform('exact', node, elem, pde, bdFlag, maxN);
    % Calculate the ZZ estimator
    [N, errZZ] = estimateuniform('zz', node, elem, pde, bdFlag, maxN);
    
    figure(fig);
    loglog(N, errH1 ./ errZZ, '-', 'color', colors{i}, 'marker', markers(i)); hold on;
  end
  figure(fig); hold on;
  title('Efficiency index comparison');
  xlabel('Number of vertices');
  legend({'Square (sin)', 'Square (poly)', 'L-shape (one)', 'Crack (one)'}, 'location', 'best');

  saveas(fig, 'figures/efficiency_ZZ.fig');
  set(fig, 'Color', 'none'); % Set background transparent
  export_fig(fig, 'figures/efficiency_ZZ','-pdf', '-painters'); % Export to PDF using export_fig
end
