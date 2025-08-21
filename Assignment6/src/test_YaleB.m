% =========================================================================
%                       Manuel Santos   2019231352
% =========================================================================

clear; close all; clc;

load("YaleB.mat");

% Show Faces of First Individual
Faces_One = zeros(h_YaleB*8,w_YaleB*8);
count = 1;
for i=1:8
    for j=1:8
        Faces_One(1+(i-1)*h_YaleB:i*h_YaleB,1+(j-1)*w_YaleB:j*w_YaleB) = reshape(allFaces_YaleB(:,count),h_YaleB,w_YaleB);
        count = count + 1;
    end
end
figure
imagesc(Faces_One),colormap('gray');
axis off; 

% Show First facef of All Individorizeuals
Face_One_All = zeros(h_YaleB*5,w_YaleB*8);
count = 1;
for i=1:5
    for j=1:8
        Face_One_All(1+(i-1)*h_YaleB:i*h_YaleB,1+(j-1)*w_YaleB:j*w_YaleB) = reshape(allFaces_YaleB(:,(count-1)*64 + 1),h_YaleB,w_YaleB);
        count = count + 1;
        if count == 40
            break;
        end
    end
end
figure
imagesc(Face_One_All),colormap('gray');
axis off;

% Select train faces
trainFaces = allFaces_YaleB(:,1:numPhotosPerSubject_YaleB*(numSubjects_YaleB-8));

% Calculate average face
avgFace = mean(trainFaces,2);

% Plot average face
figure
imagesc(reshape(avgFace(:),h_YaleB,w_YaleB));
colormap("gray"); axis off; title("Average Face");

% Calculate SVD
X = double(trainFaces)-avgFace.*ones(size(trainFaces));
[U,S,V] = svd(X,'econ');

% Show Eigenfaces
sizeTile = 8;
EigenFaces = zeros(h_YaleB*sizeTile,w_YaleB*sizeTile);
count = 1;
for i=1:sizeTile
    for j=1:sizeTile
        EigenFaces(1+(i-1)*h_YaleB:i*h_YaleB,1+(j-1)*w_YaleB:j*w_YaleB) = reshape(U(:,count),h_YaleB,w_YaleB);
        count = count + 1;
    end
end

% Plot Eigenfaces
figure
imagesc(EigenFaces.*-1),colormap('gray');
axis off; title("Faces Space - Eigenfaces");

% Show singular values
figure
semilogy(diag(S))
grid on; title("Singular Values");


