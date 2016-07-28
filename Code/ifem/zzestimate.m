function [toterr, errelem] = zzestimate(node, elem, Duh, area, f)
  % Recover the gradient, i.e. smoothen it
  Ru = recovery(node, elem, Duh, area);

  % Calculate L2 error for both components
  [~, errElem1] = getL2error(node,elem,Duh(:,1),Ru(:,1));
  [~, errElem2] = getL2error(node,elem,Duh(:,2),Ru(:,2));

  % Calculate the Zienkewicz-Zhu estimator for each element
  errelem = sqrt(errElem1 + errElem2);

  % total error
  toterr = sqrt(sum(errElem1 + errElem2));
end
