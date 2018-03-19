%% EXAMPLE OF TEXTON NOISE

% Exemplar texture
u = double(imread('input_textures/0070_4.png')); name = '0070_4';
% u = double(imread('input_textures/MarbleGranite0007_5_thumbhuge_128_a.png')); name = 'MarbleGranite0007_5_thumbhuge_128_a';
m = mean(mean(u,2));

% replace the input by its periodic component
%   (in order to reduce the border-to-border discontinuity in the exemplar)
u = perdecomp(u); 

figure(1);
imshow(uint8(u));
title('Original texture');

s = 32; % support size
order = 1; % interpolation order
nit = 50; % number of alternating projections
mni = 30; % mean number of impacts per pixel
seed = 2021; % random seed for synthesis

% Compute the interpolation coefficients
%   and appply the color correction
support = ones(s);
tic;
alpha = tn_compute_interp_coeff(u, support, order, nit);
beta = tn_color_correction(alpha, order, u);
alpha = alpha(1:s,1:s,:);
beta = beta(1:s,1:s,:);
toc

% Domain of evaluation of the procedural noise:
h = 0.25;
MF = 150;
NF = 200;
x = 1:h:NF;
y = 1:h:MF;
Y = y'*ones(size(x));
X = ones(size(y'))*x;

% %A pixel grid of size 300 x 400
% MF = 300;
% NF = 400;
% Y = (1:MF)'*ones(1,NF);
% X = ones(MF,1)*(1:NF);

% Simulation for corrected coefficients
F = tn_simulation(beta, order, mni, X, Y, seed);
figure(2);
clf;
imshow(uint8(repmat(m,[size(F,1) size(F,2) 1])+F));
title('Spline noise, corrected coeff (true contrast)');

%% Write the texton in a file

filename = ['output/' name '_s' num2str(s) '.texton'];
write_texton_file(filename,1,1,m,beta,'name',size(u));
visut = visualize_texton(filename,s+4,s+4);
imshow(uint8(visut))
