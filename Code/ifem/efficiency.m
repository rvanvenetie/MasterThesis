function efficiency
  meths = {'square_sin', 'square_ana', 'lshape_one' , 'crack_one'};
  colors = {'r', 'b', [0 0.5 0], 'm'};
  markers = ['o','x','s','d'];
  handles = [0,0,0,0];
  figure(6); clf;
  title('Efficiency index comparisons');
  xlabel('Number of vertices');
  for i=1:4
    figure(1);
    [N,equilError, resError, mixedError, errH1] = test(meths{i}, false, true)

    figure(6); hold on;
    handles(i) = semilogx(N, errH1 ./ mixedError, '-', 'color', colors{i}, 'marker', markers(i)); hold on;
    semilogx(N, errH1 ./ equilError, ':', 'color', colors{i}, 'marker', markers(i)); hold on;
  end
  legend(handles, {'Square (sin)', 'Square (poly)', 'L-shape (one)', 'Crack (one)'}, 'location', 'best');
  saveas(gca, 'figures/efficiency_compare.png');
  saveas(gca, 'figures/efficiency_compare.fig');
end
