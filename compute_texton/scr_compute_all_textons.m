%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to compute the textons of all images in the folder input_textures/
% For each image of size MxN, we compute the textons of sizes 
%   16x16, 32x32, 64x64, ..., 2^pmax x 2^pmax 
% where pmax = ceil(log2(min(M,N)))-1 saved as a .texton file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To run the script in a terminal:
% nohup matlab -nojvm < scr_compute_all_textons.m > log_scr_compute_all_textons.txt 2>&1 &

inputFolder = 'input_textures/';
outputFolder = 'output_textons/';


% Input images
imgFiles = dir([inputFolder,'*.png']); % run for all images
% imgFiles = dir([inputFolder,'wall1021.png']); % test with one image
nFiles = length(imgFiles);

% Parameters
order = 1; % interpolation order
nit = 50; % number of iterations in the texton computation algorithm

if nFiles>1
    disp(['Computation of all textons: ' num2str(nFiles) ' images']);
else
    disp(['Computation of all textons: ' num2str(nFiles) ' image']);
end

for k = 1:nFiles
    tic;
    % image filename
    imgname = imgFiles(k).name;
    disp(['Texton computation for image k = ',num2str(k),' over ',num2str(nFiles) ]);
    disp(imgname)
    u = double(imread([inputFolder,imgname]));
    imgname = strrep(imgname, '.png', '');
    
    % Compute periodic component
    u = perdecomp(u);

    % Compute max size for texton
    M = size(u,1);
    N = size(u,2);
    disp(['Image size: ' num2str(N) 'x' num2str(M)])
    pmax = ceil(log2(min(M,N)))-1;

    % Compute the mean
    m = mean(mean(u,2));
    
    for pwoftwo = (pmax-2):pmax
        
        s = 2^pwoftwo - 2*order; % so that final size is 2^pwoftwo 
                        % after added the 0 coefficients around the values
        disp(['Compute texton of size ' num2str(s) 'x' num2str(s)]);
        
        % Texton computation
        support = ones(s,s);
        alpha = tn_compute_interp_coeff(u, support, order, nit);
        beta = tn_color_correction(alpha, order, u);
        beta = beta(1:s,1:s,:);

        % Write texton files
        filetexton = [outputFolder,imgname,'_s',num2str(2^pwoftwo),'.texton'];
        filevisutexton = [outputFolder,imgname,'_s',num2str(2^pwoftwo),'_visutexton.png'];
        write_texton_file(filetexton, ...
            1, order, m/255, beta/255, imgFiles(k).name, size(u));
        imwrite(uint8(255*visualize_texton(filetexton, M, N)), filevisutexton);
        
        
    
    end
    t = toc;
    disp(['Total time to compute the ', num2str(pmax-3),' textons: ', num2str(t) 's']);
    fprintf('\n\n');
end





