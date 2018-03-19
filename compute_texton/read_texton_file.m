function [order, mu, texton] = read_texton_file(filename)
% [order, mu, texton] = read_texton_file(filename)

% Open file
fileid = fopen(filename, 'r');

% First line: Should be TEXTON
tline = fgetl(fileid);
if(strcmp(tline, 'TEXTON')==0)
    error('Not a TEXTON file');
end
% Second line: Version
version = str2num(fgetl(fileid));
switch version
    case 1
        [order, mu, texton] = read_texton_file_v1(fileid);
    otherwise
        error('version number not supported: only supported versions are 1');
end

end


function [order, mu, texton] = read_texton_file_v1(fileid)

% pursue scan file:
% 3rd line: order
order = str2num(fgetl(fileid));
% 4th line: mean
mu = str2num(fgetl(fileid));
mu = reshape(mu, [1,1,3]);
% 5th line: size
s = str2num(fgetl(fileid));
M = s(2);
N = s(1);
% M*N subsequent lines: texton coefficients
texton = fscanf(fileid, '\n%e %e %e',[3, M*N]);
texton = from_coefficients_list_to_texton(texton, M, N);

% close file
fclose(fileid);

end

function B = from_coefficients_list_to_texton(A, M, N)
B = zeros(M,N,3);
for c=1:3
    C = reshape(A(c,:),[N,M]);
    C = C';
    B(:,:,c) = C(end:-1:1,:);
end
end



