function y = spline_kernel(t,n)
%   y = spline_kernel(t,n)
%   Evaluate the B-spline function of order n on the points t

y = zeros(size(t));

for k = 0:(n+1);
    y = y + (n+1)*(-1)^k/(factorial(k)*factorial(n+1-k))*(max(0,(n+1)/2+t-k)).^n;
end

end