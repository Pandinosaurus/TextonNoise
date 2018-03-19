function beta = tn_color_correction(alpha, order, u)
%   beta = tn_color_correction(alpha, order, u)
%   Apply a color transformation matrix to the
%   texton coefficients alpha in order to force the
%   color covariance of the texton noise to be the
%   empirical color covariance of u.   
%
% Warning: alpha and u are supposed to have the same size 
%   (apply correction before cropping!)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Computation of the cross-correlation matrix of input image u:
% Image size
M = size(u,1);
N = size(u,2);
nc = size(u,3);
u = reshape(u, [M*N, nc]);
mu = mean(u);
u = u - repmat(mu,[M*N,1]);
% Cross-correlation matrix of input image u:
B = 1/(M*N)*u'*u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Computation of the covariance matrix of the noise associated with the
% coefficients alpha
% alpha is surrounded by zero coefficient so a circshift corresponds to a
% shift with zero boundary condition...

% Need to compute the sampling of the kernel of order 2*order+1 that has
% support ]-2*order, 2*order[
A = zeros(nc,nc);
for tx = (-2*order):(2*order)
    for ty = (-2*order):(2*order)
        S = reshape(alpha, [M*N, nc])' * reshape(circshift(alpha, [-tx, -ty]), [M*N, nc]);
        A = A + spline_kernel(tx,2*order+1)*spline_kernel(ty,2*order+1)*S;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Computation of squareroots matrices and color correction

SA = sqrtm(A);
invSA = inv(SA);
SB = sqrtm(B);
C = SB*invSA;
beta = reshape( reshape(alpha,[M*N, nc]) * C' ,[M,N,nc]);


end



