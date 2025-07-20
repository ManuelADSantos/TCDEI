% =========================================================================
%                       Manuel Santos   2019231352
% =========================================================================

clear; close all; clc;

load("Yale.mat");

% Select train faces
trainFaces = allFaces_Yale(:,1:numPhotosPerSubject_Yale*(numSubjects_Yale-1));

% Calculate average face
avgFace = mean(trainFaces,2);

% Plot average face
figure(1)
imagesc(reshape(avgFace(:),h_Yale,w_Yale));
colormap("gray"); axis off;

% Calculate SVD
X = double(trainFaces)-avgFace.*ones(size(trainFaces));
[U,S,V] = svd(X,'econ');

% Show Eigenfaces
sizeTile = 8;
EigenFaces = zeros(h_Yale*sizeTile,w_Yale*sizeTile);
count = 1;
for i=1:sizeTile
    for j=1:sizeTile
        EigenFaces(1+(i-1)*h_Yale:i*h_Yale,1+(j-1)*w_Yale:j*w_Yale) = reshape(U(:,count),h_Yale,w_Yale);
        count = count + 1;
    end
end
figure(2)
imagesc(EigenFaces.*-1),colormap('gray');
axis off

% Show singular values
figure(3)
semilogy(diag(S))
grid on