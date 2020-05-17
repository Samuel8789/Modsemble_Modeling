function [indices] = all_but_me(low,high)
N = high-low+1;
tmp = repmat((low:high)', 1, N);
tmp = tmp(~eye(size(tmp)));
tmp = reshape(tmp,N-1,N)';
indices = num2cell(tmp,2)';
end
