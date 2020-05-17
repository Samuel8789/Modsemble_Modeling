function [max_idx,g_max] = find_graph_cut_max(fittedmodel,XData)
yfitted=feval(fittedmodel,XData);
[g_max,idx]=max(yfitted);
max_idx=XData(idx);
 end