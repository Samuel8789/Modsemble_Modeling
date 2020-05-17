function [s,t] = ster(graph)
N = length(graph);
K = N;
a=1;
b=1;
s=[];
t=[];

for a = 1:N
    for b = 1:N
        if graph(a,b) == 1
            s = [s a];
            t = [t b];
            b=b+1;
        end
    end
    a=a+1;
end
end

