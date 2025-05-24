function [vec] = TurnCell2Vector(cells)

n = length(cells);
vec = zeros(n,1);
for i = 1 : n
   vec(i) = cells{i};
end


end