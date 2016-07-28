[node ,elem] = circlemesh(0, 0,1, 0.65);
node(:,2) = .7*node(:,2);
groen = [200, 238, 200] ./255;
rood = [236, 124, 111] ./255;

T = auxstructure(elem);
fig = figure(1); clf;
showmesh(node, elem, 'facecolor' ,groen);
findnode(node, 'all', 'noindex');
set(fig, 'Color', 'none'); % Set background transparent
export_fig(fig, 'figures/mesh/ex_mesh.pdf', '-painters');

% Highlight some stars, 3 = bdr, 7 = int
bdrVtx = 3;
intVtx = 7;
NT = size(elem, 1);
N = size(node,1);
t2v = sparse([1:NT, 1:NT, 1:NT], elem, 1, NT, N)
bdrPatch = [find(t2v(:,bdrVtx))];
intPatch = [find(t2v(:,intVtx))];
%bdrpatch
fig = figure(2); clf;
showmesh(node, elem, 'facecolor' ,groen);
findelem(node, elem, bdrPatch, 'noindex', 'facecolor' , rood);
findnode(node, 'all', 'noindex', 'markersize', 15);
findnode(node, bdrVtx, 'noindex','color' , 'k', 'markersize', 30);
text(node(bdrVtx,1)+0.07, node(bdrVtx,2) + 0.065, 'a', 'Fontsize', 16, 'fontweight', 'bold')
set(fig, 'Color', 'none'); % Set background transparent
export_fig(fig, 'figures/mesh/ex_bdpatch.pdf', '-painters');

%intpatch
fig = figure(3); clf;
showmesh(node, elem, 'facecolor' , groen);
findelem(node, elem, intPatch, 'noindex', 'facecolor' ,rood);
findnode(node, 'all', 'noindex', 'markersize', 15);
findnode(node, intVtx, 'noindex','color' , 'k', 'markersize', 30);
text(node(intVtx,1)+ 0.05, node(intVtx,2) + 0.05, 'a', 'Fontsize', 16, 'fontweight', 'bold')
set(fig, 'Color', 'none'); % Set background transparent
export_fig(fig, 'figures/mesh/ex_intpatch.pdf', '-painters');
