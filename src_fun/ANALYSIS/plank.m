%choose n

graph=tril(Models.models{n}.structure);
edge_pot = tril(Models.models{n}.theta.edge_potentials);
G = Models.models{n}.theta.G;

[G_all] = getEdgePotAll(graph,G);
[S_all,T_all,Edge_wt_all] = vectorize_model(graph,G_all);
Edge_wt_all = transpose(Edge_wt_all);

%CRF for G_11 (that is, Phi 11)
[G_on] = getEdgePot(graph,G,4);
[S_11,T_11,Edge_wt_11] = vectorize_model(graph,G_on);
Edge_wt_11 = transpose(Edge_wt_11);

%Edge_wt_all = normalize(Edge_wt_all,'range');
%Edge_wt_11 = normalize(Edge_wt_11,'range');

MODEL_ALL = digraph(S_all,T_all,Edge_wt_all);
MODEL_11  = digraph(S_11,T_11,Edge_wt_11);

Lwidths = 5;
colorweight_all = MODEL_ALL.Edges.Weight/max(MODEL_ALL.Edges.Weight);
colorweight_11 = MODEL_11.Edges.Weight/max(MODEL_11.Edges.Weight);

figure
MDL_all = plot(MODEL_ALL, 'LineWidth', Lwidths, 'Layout', 'subspace3');
MDL_all.ShowArrows = 'off';
MDL_all.Marker = 's';
colormap winter
MDL_all.EdgeCData = colorweight_all;
colorbar;
title('All');

figure
MDL_11 = plot(MODEL_11, 'LineWidth', Lwidths, 'Layout', 'subspace3');
MDL_11.ShowArrows = 'off';
MDL_11.Marker = 's';
colormap winter
MDL_11.EdgeCData = colorweight_11;
colorbar;
title('Phi 11');