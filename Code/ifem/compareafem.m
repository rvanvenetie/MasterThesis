function CompareAfem(method, estimators, nodeOri, elemOri, pde, bdFlagOri, theta,maxN) 
  function plot
    f2 = figure(2);clf;
    loglog(Nu,  erru, 'b-o',Nr, errr, 'r-x'); hold on;
    loglog(Ne,  erre, '-s', 'color', [0 0.5 0]); 
    legend({'uniform($U_k$)','residual($U_k$)', 'equilibrated($U_k, \zeta$)'}, 'interpreter', 'latex');
    title('Comparison AFEM performance')
    xlabel('Number of vertices');
    ylabel('Exact error');
    saveas(f2, sprintf('%s/%s/norm_%d_%g.png',savedir, method, maxN, theta));
    saveas(f2, sprintf('%s/%s/norm_%d_%g.fig',savedir, method, maxN, theta));
    loglog(Nm, errm, 'm-d');
    legend({'uniform($U_k$)','residual($U_k$)', 'equilibrated($U_k, \zeta$)', 'mixed($U_k, \sigma$)'}, 'interpreter', 'latex');
    saveas(f2, sprintf('%s/%s/norm_mixed_%d_%g.png',savedir, method, maxN, theta));
    saveas(f2, sprintf('%s/%s/norm_mixed_%d_%g.fig',savedir, method, maxN, theta));
  end

  savedir = 'figures/';
  % Uniform refinements
  Nu = [];
  erru = [];
  % AFEM residual estimator
  Nr = [];
  errr = [];
  % AFEM mixed estimator
  Nm = [];
  errm = [];
  % AFEM equil estimator
  Ne = [];
  erre = [];


  % Uniform refinement
  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
  while (size(node,1) < maxN)
    uh = Poisson(node, elem, pde, bdFlag);
    % Calculate real error
    Nu(end+1) = size(node, 1);
    erru(end+1) = exactH1error(node, elem, bdFlag, pde, uh);

    [node, elem, bdFlag] = uniformbisect(node, elem, bdFlag);
  end
  plot
  erru

  % Residual refinements
  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
  while (size(node,1) < maxN)
    uh = Poisson(node, elem, pde, bdFlag);

    % Calculate the real error
    Nr(end+1) = size(node, 1);
    errr(end+1) = exactH1error(node, elem, bdFlag, pde,uh);

    % Perform afem step
    eta = estimateresidual(node, elem, uh, pde);
    markedElem = mark(elem,eta,theta);
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
  end
  plot
  errr

  % Mixed refinements
  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
  while (size(node,1) < maxN)
    uh = Poisson(node, elem, pde, bdFlag);
    [Duh,~] = gradu(node, elem, uh);

    % Calculate the real error
    Nm(end+1) = size(node, 1);
    errm(end+1) = exactH1error(node, elem, bdFlag, pde,uh);

    % Perform afem step
    [~,sigma] = PoissonRT0(node, elem, pde, bdFlag);
    [~, eta]  = getL2errorRT0(node, elem,  Duh, sigma);
    
    markedElem = mark(elem,sqrt(eta),theta);
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
  end
  plot
  errm


  % Equilibrated refinements
  node = nodeOri; elem = elemOri; bdFlag = bdFlagOri;
  while (size(node, 1) < maxN)
    plot
    uh = Poisson(node, elem, pde, bdFlag);
    [Duh,~] = gradu(node, elem, uh);

    % Calculate the real error
    Ne(end+1) = size(node, 1);
    erre(end+1) = exactH1error(node, elem, bdFlag, pde,uh);

    % Perform afem step
    sig = flux(node,elem,  Duh, pde.f);
    [~, eta, ~] =  equil(node, elem, Duh, sig, pde.f);
    markedElem = mark(elem,eta,theta);
    [node,elem,bdFlag] = bisect(node,elem,markedElem,bdFlag);
  end
  erre
  plot
end

function [N, err] = EstimateAfem(estimator, node, elem, pde, bdFlag, maxN)
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
