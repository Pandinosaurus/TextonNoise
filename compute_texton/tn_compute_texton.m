function tn_compute_texton(input_image, output_folder, s)
% tn_compute_texton(input_image, output_folder)
% computes the texton associated to an exemplar texture.
% The texton is written in a .texton file on the hard disk.
%   Optional input: texton size s

if nargin<3
    s = 64;
end

% Read input image
u = double(imread(input_image));
% get filename without root folder and extension:
ind = strfind(input_image, '/');
tmp = input_image(ind(end)+1:end);
imgname = strrep(tmp, '.png', '');

% mean value
m = mean(mean(u,2));

% replace the input by its periodic component
%   (in order to reduce the border-to-border discontinuity in the exemplar)
u = perdecomp(u); 

% Parameters
order = 1; % interpolation order
nit = 50; % number of alternating projections

% Compute the interpolation coefficients
%   and appply the color correction
support = ones(s);
tic;
alpha = tn_compute_interp_coeff(u, support, order, nit);
beta = tn_color_correction(alpha, order, u);
beta = beta(1:s,1:s,:);
toc

filename = [imgname '_s' num2str(s) '.texton'];
write_texton_file([output_folder '/' filename],1,order,m/255,beta/255,'name',size(u));

% quit matlab
exit
end

