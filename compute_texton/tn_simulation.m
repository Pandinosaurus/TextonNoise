function F = tn_simulation(alpha, order, mni, X, Y, seed)
% F = tn_simulation(alpha, order, mni, X, Y, seed)
% samples a texton noise f whose kernel function
% is the convolution of the interpolation coefficients alpha 
%	and a B-spline function.
%
% INPUT:
%  - alpha: interpolation coefficients (graylevel or color)
%  - order: order for the spline interpolation
%  - mni: mean number of impacts per point (gives the intensity lambda of
%  the Poisson process)
%  - [X,Y] grid of the domain so simulate
%  - seed (optional): seed for the pseudo-random number generator.
% OUTPUT:
%  - F a matrix such that F(i,j,:) = f(X(i,j), Y(i,j), :)

% PNRG seed:
if nargin==6
    RandStream.setDefaultStream(RandStream('mt19937ar','Seed',seed));
end

% Support of impulse function for the grid size:
ma = size(alpha,1);
na = size(alpha,2);
C = size(alpha,3);
m = floor((order+1)/2);
mh = ma + 2*m;
nh = na + 2*m;
lambda = mni/(mh*nh); % Poisson process intensity

% Generation of the Poisson process P=[Px,Py]
% Computation of the domain:
xmin = min(X(:)) - mh;
xmax = max(X(:)) + mh;
ymin = min(Y(:)) - mh;
ymax = max(Y(:)) + mh;
nbpoints = poissrnd(lambda*(xmax-xmin)*(ymax-ymin));
Px = xmin + (xmax-xmin)*rand(nbpoints,1);
Py = ymin + (ymax-ymin)*rand(nbpoints,1);

F = zeros([size(X) C]);
disp(['Number of points is ',num2str(nbpoints)]);
t = 0.1;
dt = 0.1;
for k=1:nbpoints
    xk = Px(k);
    yk = Py(k);
    Ik = (X>=(xk-m)) & (X<=(xk+m+ma)) & (Y>=(yk-m)) & (Y<=(yk+m+ma));
    for c=1:C
        Fc = F(:,:,c);
        Fc(Ik) = Fc(Ik) + evaluate_spline_interpolate_from_coef(alpha(:,:,c),X(Ik)-Px(k),Y(Ik)-Py(k),order,0);
        F(:,:,c) = Fc;
    end
    if(k/nbpoints >=t)
        disp([num2str(100*t),' % of points done']);
        t = t +dt;
    end
end
% CLT normalization:
% Need to remove lambda x integral of impulse function h 
for c=1:C
    F(:,:,c) = (F(:,:,c)-lambda*sum(sum(alpha(:,:,c))))*1/sqrt(lambda);
end

end










