% Performs AFEM driven by the estimator presented, and plots the results
function compareafem(method, estimators, node, elem, pde, bdFlag, theta, maxN)
  savedir = 'figures';

  % Create new plot figures
  f2 = figure(2); clf; hold on;

  % Iteratively build legend
  leg = {}

  % Calculate & plot the result for the various estimators
  for estimator = estimators
    % display the current estimator
    disp(estimator{1})

    % calculate the shit for this estimator
    [N,err] =  estimateafem(estimator{1}, node, elem, pde, bdFlag, theta, maxN)

    legend({'uniform($U_k$)','residual($U_k$)', 'equilibrated($U_k, \zeta$)', 'mixed($U_k, \sigma$)'}, 'interpreter', 'latex');
    % plot the results
    switch estimator{1}
    case 'uniform'
      leg{end+1} = 'uniform($U_k$)'
      figure(f2); loglog(N, err, 'b-o')
    case 'residual'
      leg{end+1} = 'residual($U_k$)';
      figure(f2); loglog(N, err, 'r-x'); % abs
    case 'mixed'
      leg{end+1} = 'mixed($U_k, \sigma$)';
      figure(f2); loglog(N, err, 'm-d'); % abs
    case 'equilibrated'
      leg{end+1} = 'equilibrated($U_k, \zeta$)';
      figure(f2); loglog(N, err, '-s', 'color', [0 0.5 0]); % abs
    case 'zz'
      leg{end+1} = 'ZZ($U_k$)';
      figure(f2); loglog(N, err, 'k-*'); %abs
    end
  end
  figure(f2);
  legend(leg, 'interpreter', 'latex');
  title('Comparison AFEM performance')
  xlabel('Number of vertices');
  ylabel('Exact error');
end

function [N,err] = estimateafem(estimator, node, elem, pde, bdFlag, theta, maxN)
  if strcmp(estimator, 'uniform')
    [N,err] = estimateuniform('exact',node,elem,pde,bdFlag,maxN)
    return
  end
  N = []
  err = []
  while (size(node,1) < maxN)
    % Calculate the discrete solution and its derivative
    uh = Poisson(node, elem, pde, bdFlag);
    [Duh,area] = gradu(node, elem, uh);

    % Store the number of nodes and the `real' error
    N(end+1) = size(node, 1);
    err(end+1) = exactH1error(node, elem, bdFlag, pde,uh)

    % Perform afem step
    switch estimator
    case 'mixed'
      % Compute the mixed FEM estimator
      [~,sigma] = PoissonRT0(node, elem, pde, bdFlag);
      [~, eta]  = getL2errorRT0(node, elem,  Duh, sigma);
      eta = sqrt(eta)
    case 'equilibrated'
      % Compute the equilibrated flux estimator
      [~, eta, ~] =  equil(node, elem, Duh, sig, pde.f);
    case 'residual'
      % Compute the classical residual estimator
      eta = estimateresidual(node, elem, uh, pde);
    case 'zz'
      % Calculate the zienkewicz-zhu estimator
      [~,eta] = zzestimate(node, elem, Duh, area, pde.f);
    otherwise
      disp('You dumb boy, select a valid estimator')
      return
    end

    % Apply dorfler marking
    markedElem = mark(elem,eta,theta);

    % Refine the marked elements
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
  end
end

%function CompareAfem(method, estimators, nodeOri, elemOri, pde, bdFlagOri, theta,maxN) 
%  function plot
%    f2 = figure(2);clf;
%    loglog(Nu,  erru, 'b-o',Nr, errr, 'r-x'); hold on;
%    loglog(Ne,  erre, '-s', 'color', [0 0.5 0]); 
%    legend({'uniform($U_k$)','residual($U_k$)', 'equilibrated($U_k, \zeta$)'}, 'interpreter', 'latex');
%    title('Comparison AFEM performance')
%    xlabel('Number of vertices');
%    ylabel('Exact error');
%    saveas(f2, sprintf('%s/%s/norm_%d_%g.png',savedir, method, maxN, theta));
%    saveas(f2, sprintf('%s/%s/norm_%d_%g.fig',savedir, method, maxN, theta));
%    loglog(Nm, errm, 'm-d');
%    legend({'uniform($U_k$)','residual($U_k$)', 'equilibrated($U_k, \zeta$)', 'mixed($U_k, \sigma$)'}, 'interpreter', 'latex');
%    saveas(f2, sprintf('%s/%s/norm_mixed_%d_%g.png',savedir, method, maxN, theta));
%    saveas(f2, sprintf('%s/%s/norm_mixed_%d_%g.fig',savedir, method, maxN, theta));
%  end
%
%  savedir = 'figures/';
%  % Uniform refinements
%  Nu = [];
%  erru = [];
%  % AFEM residual estimator
%  Nr = [];
%  errr = [];
%  % AFEM mixed estimator
%  Nm = [];
%  errm = [];
%  % AFEM equil estimator
%  Ne = [];
%  erre = [];
%
%
%  % Uniform refinement
%  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
%  while (size(node,1) < maxN)
%    uh = Poisson(node, elem, pde, bdFlag);
%    % Calculate real error
%    Nu(end+1) = size(node, 1);
%    erru(end+1) = exactH1error(node, elem, bdFlag, pde, uh);
%
%    [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
%  end
%  plot
%  erru
%
%  % Residual refinements
%  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
%  while (size(node,1) < maxN)
%    uh = Poisson(node, elem, pde, bdFlag);
%
%    % Calculate the real error
%    Nr(end+1) = size(node, 1);
%    errr(end+1) = exactH1error(node, elem, bdFlag, pde,uh);
%
%    % Perform afem step
%    eta = estimateresidual(node, elem, uh, pde);
%    markedElem = mark(elem,eta,theta);
%    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
%  end
%  plot
%  errr
%
%  % Mixed refinements
%  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
%  while (size(node,1) < maxN)
%    uh = Poisson(node, elem, pde, bdFlag);
%    [Duh,~] = gradu(node, elem, uh);
%
%    % Calculate the real error
%    Nm(end+1) = size(node, 1);
%    errm(end+1) = exactH1error(node, elem, bdFlag, pde,uh);
%
%    % Perform afem step
%    [~,sigma] = PoissonRT0(node, elem, pde, bdFlag);
%    [~, eta]  = getL2errorRT0(node, elem,  Duh, sigma);
%    
%    markedElem = mark(elem,sqrt(eta),theta);
%    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
%  end
%  plot
%  errm
%
%
%  % Equilibrated refinements
%  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
%  while (size(node, 1) < maxN)
%    plot
%    uh = Poisson(node, elem, pde, bdFlag);
%    [Duh,~] = gradu(node, elem, uh);
%
%    % Calculate the real error
%    Ne(end+1) = size(node, 1);
%    erre(end+1) = exactH1error(node, elem, bdFlag, pde,uh);
%
%    % Perform afem step
%    sig = flux(node,elem,  Duh, pde.f);
%    [~, eta, ~] =  equil(node, elem, Duh, sig, pde.f);
%    markedElem = mark(elem,eta,theta);
%    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
%  end
%  erre
%  plot
%end
%
