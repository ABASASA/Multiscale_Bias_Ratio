function [SNR] = ComputeSNRNumericllySPD(sigma)

Ntrail = 2500;
func =  @(x,y) [sin(2*pi*y) .* cos(2*pi*x), y.^2, x.*y;  y.^2, 1, 0;  x.*y, 0,  cos(pi * x)];

res = zeros(Ntrail, 2); %[clean, noised]

for i = 1 : Ntrail
    point = rand(1,2);
    clean = func(point(1), point(2));
    res(i,:) = [norm(clean)^2, ...
                norm(AddNoiseSPDInner(clean, sigma))^2];
end

SNR = mean(res(:,1) ./ res(:,2));
end


function [value] = AddNoiseSPDInner(X, sigma)
     noise = randn(3,3);
     noise = triu(noise) + triu(noise).' - diag(diag(noise));
     value = X(:) * sigma .* noise(:).';
end