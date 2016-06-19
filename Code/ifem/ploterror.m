function ploterror(method,node, elem, uh,N,equilError, resError, errH1)    
  savedir = 'figures';
  mkdir(sprintf('%s/%s', savedir, method));
  len = max(size(N));

  N
  N.^(-0.5)
  2.^(- (1:len))
  errH1
  resError
  equilError

  f2 = figure(2);clf;
  loglog(N,  errH1, 'b-o',N, resError, 'r-o'); hold on;
  loglog(N, equilError, '-o', 'color', [0 0.5 0]); 
  legend({'$\|\nabla{u} - \nabla {U_k}\|$','residual($U_k$)', 'equilibrated($U_k, \zeta$)'}, 'interpreter', 'latex');
  title('Comparison error estimators');
  xlabel('Number of vertices');
  ylabel('Error');
  saveas(f2, sprintf('%s/%s/norm_%d.png',savedir, method, len));
  saveas(f2, sprintf('%s/%s/norm_%d.fig',savedir, method, len));
  if (len > 2)
    x1 = N(end-2)
    x2 = N(end-1)

    y12 = errH1(end-1);
    y11 = exp(log(y12) - 1.0/2.0);

    loglog([x1 x1 x2 x1], [y12 y11 y11 y12], 'color', [0,0,0]);
  end
  saveas(f2, sprintf('%s/%s/norm_slope_%d.png',savedir, method, len));
  saveas(f2, sprintf('%s/%s/norm_slope_%d.fig',savedir, method, len));

  f3 = figure(3);
  semilogx(N, errH1 ./ resError, 'r-o'); hold on;
  semilogx(N,  errH1 ./ equilError, '-o', 'color', [0 0.5 0]);
  title('Efficiency index');
  xlabel('Number of vertices');
  legend({'residual($U_k$)', 'equilibrated($U_k, \zeta$)'}, 'interpreter', 'latex', 'location', 'northwest');


  f4 = figure(4); colormap('jet');
  showsolution(node,elem,uh,3)
  f5 = figure(5); colormap('jet'); 
  showsolution(node,elem,uh,2)

  saveas(f3, sprintf('%s/%s/efficiency_%d.png',savedir, method, len));
  saveas(f3, sprintf('%s/%s/efficiency_%d.fig',savedir, method, len));

  saveas(f4, sprintf('%s/%s/uh_3d_%d.png',savedir, method, len));
  saveas(f5, sprintf('%s/%s/uh_2d_%d.png',savedir, method, len));
end
