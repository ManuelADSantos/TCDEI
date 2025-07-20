% =========================================================================
%                       Manuel Santos   2019231352
% =========================================================================

clear; close all; clc;

load("AT&T.mat");

% Select train faces
trainFaces = allFaces_ATeT(:,1:numPhotosPerSubject_ATeT*(numSubjects_ATeT-1));

% Calculate average face
avgFace = mean(trainFaces,2);

% Plot average face
figure(1), axis off;
imagesc(reshape(avgFace(:),h_ATeT,w_ATeT));
colormap("gray");

% Calculate SVD
X = double(trainFaces)-avgFace.*ones(size(trainFaces));
[U,S,V] = svd(X,'econ');

% Show Eigenfaces
sizeTile = 8;
EigenFaces = zeros(h_ATeT*sizeTile,w_ATeT*sizeTile);
count = 1;
for i=1:sizeTile
    for j=1:sizeTile
        EigenFaces(1+(i-1)*h_ATeT:i*h_ATeT,1+(j-1)*w_ATeT:j*w_ATeT) = reshape(U(:,count),h_ATeT,w_ATeT);
        count = count + 1;
    end
end
figure(2),axis off
imagesc(EigenFaces.*-1),colormap('gray');

% Show singular values
figure(3),grid on
semilogy(diag(S))