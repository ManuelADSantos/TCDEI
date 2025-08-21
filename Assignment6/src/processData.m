% =========================================================================
%                       Manuel Santos   2019231352
% =========================================================================

%% ========== AT&T dataset preparation ==========
clear, close all; clc;

dataset_ATeT = uigetdir(pwd, 'Select the folder containing the AT&T dataset');

w_ATeT = 92; h_ATeT = 112; 
numSubjects_ATeT = 40;
numPhotosPerSubject_ATeT = 10;

allFaces_ATeT = zeros(w_ATeT*h_ATeT,numSubjects_ATeT*numPhotosPerSubject_ATeT,'uint8');

for i = 1:numSubjects_ATeT
    for j = 1:numPhotosPerSubject_ATeT         
        image = im2gray(imread([dataset_ATeT '/s' num2str(i) '/' num2str(j) '.pgm']));
        allFaces_ATeT(:,(i-1)*numPhotosPerSubject_ATeT + j) = reshape(image,w_ATeT*h_ATeT,1);
        disp(['AT&T Subject ' num2str(i) ' photo ' num2str(j) ' ==> Done saving']);
    end
end

figure
imshow(image)
hold on

save("AT&T.mat","allFaces_ATeT","numSubjects_ATeT","numPhotosPerSubject_ATeT","h_ATeT","w_ATeT");


%% ========== Yale dataset preparation ==========
clear, close all; clc;

dataset_Yale = uigetdir(pwd, 'Select the folder containing the Yale dataset');

w_Yale = 320; h_Yale = 243; 
numSubjects_Yale = 15;
numPhotosPerSubject_Yale = 11;

allFaces_Yale = zeros(w_Yale*h_Yale,numSubjects_Yale*numPhotosPerSubject_Yale,'uint8');

images = dir([dataset_Yale '/archive/subject*']);
for i = 1:size(images)
    image = im2gray(imread([dataset_Yale '/archive/' images(i).name]));
    allFaces_Yale(:,i) = reshape(image,w_Yale*h_Yale,1);
    disp(['Yale ' images(i).name ' ==> Done saving']);
end

figure
imshow(image)

save("Yale.mat","allFaces_Yale","numSubjects_Yale","numPhotosPerSubject_Yale","h_Yale","w_Yale");


%% ========== YaleB dataset preparation ==========
clear, close all; clc;

dataset_YaleB = uigetdir(pwd, 'Select the folder containing the Yale B dataset');

w_YaleB = 168; h_YaleB = 192; 
numSubjects_YaleB = 38;
numPhotosPerSubject_YaleB = 64;

allFaces_YaleB = zeros(w_YaleB*h_YaleB,numSubjects_YaleB*numPhotosPerSubject_YaleB,'uint8');

for i = 1:numSubjects_YaleB+1
    if(i==14)
        continue;
    end
    i_form = sprintf('%02d',i);
    images = dir([dataset_YaleB '/yaleB' i_form '/*.pgm']);
    for j = 1:size(images)     
        image = im2gray(imread([dataset_YaleB '/yaleB' i_form '/' images(j).name]));
        allFaces_YaleB(:,(i-1)*numPhotosPerSubject_YaleB + j) = reshape(image,w_YaleB*h_YaleB,1);
        disp(['YaleB Subject ' num2str(i) ' photo ' num2str(j) ' ==> Done saving']);
    end
end

figure
imshow(image)
hold on

save("YaleB.mat","allFaces_YaleB","numSubjects_YaleB","numPhotosPerSubject_YaleB","h_YaleB","w_YaleB");

% cols_with_all_zeros = find(all(allFaces_YaleB==0))
