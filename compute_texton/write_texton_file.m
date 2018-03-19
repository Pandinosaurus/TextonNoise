function write_texton_file(filename, version, order, mu, texton, name_ori, size_ori)
% write_texton_file(filename, version, order, mu, texton, name_ori, size_ori)

if(nargin < 7)
    error('There is not enough input arguments');
end

switch version
    case 1
        write_texton_file_v1(filename, order, mu, texton, name_ori, size_ori);
    otherwise
        error('version number not supported: only supported versions are 1');
end



end


function write_texton_file_v1(filename, order, mu, texton, name_ori, size_ori)

% Open file
fileid = fopen(filename, 'w');

% First line: filetype
fprintf(fileid, 'TEXTON');

% version number
fprintf(fileid, '\n1');

% interpolation order
fprintf(fileid, ['\n',num2str(order)]);

% RGB mean
fprintf(fileid, ['\n',num2str(mu(1), '%e'),' ',num2str(mu(2), '%e'),' ',num2str(mu(3), '%e')]);


% Linear TEXTON must have a one pixel frame of zero values! 
%   This must be added when writing the file
texton = add_black_frame_if_not_present(texton);

% Texton size
fprintf(fileid, ['\n',num2str(size(texton,2)),' ',num2str(size(texton,1))]);

% Coefficients
texton = list_texton_coefficients(texton);
fprintf(fileid, '\n%e %e %e', texton);
% for m=size(texton,1):-1:1
%     for n=1:size(texton,2)
%         fprintf(fileid, ['\n',num2str(texton(m,n,1), '%e'),' ',num2str(texton(m,n,2), '%e'),' ',num2str(texton(m,n,3), '%e')]);
%     end
% end


% Name and size of original image
fprintf(fileid, ['\nOriginal image: ',name_ori,' (',num2str(size_ori(2)),' x ', num2str(size_ori(1)),')']);


end

function B = add_black_frame_if_not_present(A)
% add a zero frame if the texton is not already zero at the border
%if(any(A(:,1))||any(A(:,end))||any(A(1,:))||any(A(end,:)))
    M = size(A,1);
    N = size(A,2);
    B = zeros(M+2,N+2,3);
    B(2:(M+1),2:(N+1),:) = A;
%end
end

function B = list_texton_coefficients(A)
% List the texton coefficient in the good order for opengl textures (line
% by line and starting at bottom left.
M = size(A,1);
N = size(A,2);
B = zeros(size(A,3),M*N);
for c=1:size(A,3)
    B(c,:) = reshape(A(end:-1:1,:,c)', [1,M*N]);
end
end







