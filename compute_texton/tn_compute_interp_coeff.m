function alpha = tn_compute_interp_coeff(u,support,order, nit)
%   alpha = tn_compute_interp_coeff(u,support,order, nit)
%   Compute the interpolation coefficients for the texton
%   noise associated to the exemplar texture u.
%
%   Inputs:
%   - u: exemplar texture (graylevel or color image)
%   - order: interpolation order
%   - support: binary image giving the desired support of alpha
%   - nit: number of iterations 
%		Output:
%		- alpha: interpolation coefficients

% Image sizes
[M,N,C] = size(u);

% Extend the binary image support into an image of size support
tmp = zeros(M,N);
tmp(1:size(support,1),1:size(support,2)) = support;
support = repmat(tmp,[1,1,C]);

% Compute fft for spectral projection:
tu = zeros(M,N,C);
fproj = zeros(M,N,C);
modb = dft_of_sampled_spline(2*order+1, M, N);
for c=1:C
    tu(:,:,c) = 1/sqrt(M*N)*(u(:,:,c)-mean(mean(u(:,:,c))));
    fproj(:,:,c) = fft2(tu(:,:,c))./sqrt(modb);
end

% Initialization:
alpha = randn(M,N,C);

% Iterative alternating projections:
for it=1:nit
    % spectral projection:
    falpha = fft2(alpha);
    sp = sum(falpha.*conj(fproj),3);
    sp = repmat(exp(1i*angle(sp)),[1,1,3]);
    falpha = sp.*fproj;
    alpha = real(ifft2(falpha));
    
    % support projection:
    alpha = support.*alpha;
end

end


function fftb = dft_of_sampled_spline(n, M, N)

% M and N are supposed to be larger then 2*n+1
b = zeros(M,N);
for t1 = -(n-1)/2:(n-1)/2
    if t1<0
        indt1 = 1+t1+M;
    else
        indt1 = 1+t1;
    end
    
    for t2 = -(n-1)/2:(n-1)/2
        if t2<0
            indt2 = 1+t2+N;
        else
            indt2 = 1+t2;
        end
        b(indt1, indt2) = spline_kernel(t1,n)*spline_kernel(t2,n);
    end
end

fftb = real(fft2(b));

end



