
function v = evaluate_spline_interpolate_from_coef(u,X,Y,order,bdc)
% v = evaluate_spline_interpolate_from_coef(u,X,Y,order,bdc)
% Evaluate the spline interpolate with interpolation coefficient on the
% grid [X,Y].
% 
% Inputs:
%   - u: interpolation coefficients
%   - X: first coordinates of the evaluation grid
%   - Y: second coordinates of the evaluation grid
%   - order: interpolation order
%   - bdc: boundary condition
% 
% Available interpolations are:
%    order = 0           nearest neighbor (default)
%    order = 1           bilinear
%    order = 3,5,7,9,11  spline with specified order
% Available boundary conditions are:
%    bdc = 0     zero at the border
%    bdc = 1     symmetric boundary condition (default)
%    bdc = 2     periodic boundary condition

if nargin < 5
    bdc = 1; % symmetric boundary condition by default
end
if nargin < 4
    order = 0; % nearest neighbor interpolation by default
end
[ny,nx] = size(u);
[my,mx] = size(X);
if order==0 % nearest neighbor interpolation
    iX = floor(X+0.5);
    iY = floor(Y+0.5);
    n1 = 0;
    n2 = 0;
else % other interpolations:
    iX = floor(X);
    iY = floor(Y);
    X = X-iX;
    Y = Y-iY;
    if sum([1,3,5,7,9,11]==order)==1 % spline interpolation
        n2 = (order+1)/2;
        c = mat_coeff_splinen(order);
    else
        error('Unrecognized interpolation order.');
    end
    n1 = 1-n2;
end
% boundary conditions:
supx1 = max(0,1-min(iX(:)))-n1;
supx2 = max(0,max(iX(:))-nx)+n2;
iX = iX + supx1;
supy1 = max(0,1-min(iY(:)))-n1;
supy2 = max(0,max(iY(:))-ny)+n2;
iY = iY + supy1;
switch bdc
    case 0 % zero at border
        u = [zeros(size(u,1), supx1),u,zeros(size(u,1), supx2)];
        u = [zeros(supy1, size(u,2));u;zeros(supy2, size(u,2))];
    case 1 % symmetric boundary conditions
        u = [u(:,1+supx1:-1:2),u,u(:,end-1:-1:end-supx2)];
        u = [u(1+supy1:-1:2,:);u;u(end-1:-1:end-supy2,:)];
    case 2 % periodic boundary conditions
        u = [u(:,(end-supx1+1):end),u,u(:,1:supx2)];
        u = [u((end-supy1+1):end,:);u;u(1:supy2,:)];
end
nny = size(u,1);  % that is, nny = ny + supy1 + supy2
if order==0 % nearest neighbor interpolation
    v = reshape(u(iY+nny*(iX-1)),my,mx);
else % other interpolations
    % compute interpolation coefficients
    cx = c*( (ones(order+1,1)*(X(:)')).^((0:order)'*ones(1,mx*my)) );
    cy = c*( (ones(order+1,1)*(Y(:)')).^((0:order)'*ones(1,mx*my)) );
    v = zeros(size(mx*my));
    for dx = n1:n2
        for dy = n1:n2
            v = v + cy(n2+1-dy,:).*cx(n2+1-dx,:).*reshape(u(iY+dy+nny*(iX+dx-1)),1,[]);
        end
    end
    v = reshape(v,my,mx);
end
end



% coefficients of piecewise polynomial spline function with order n
function c = mat_coeff_splinen(n)
c = zeros(n+1,n+1);
a = zeros(1,n+2);
a(1) = 1/prod(1:n);
for k=1:n+1
    a(k+1) = -a(k)*(n+2-k)/k;
end
for k=0:n+1
    for p=0:n
        xp = prod((n-p+1:n)./(1:p))*k^(n-p);
        for i=k:n
            c(i+1,p+1) = c(i+1,p+1)+a(i-k+1)*xp;
        end
    end
end
end




