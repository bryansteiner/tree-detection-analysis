% Image - insert the name
I = imread('img/trees1.png');
[x,y,z] = size(I);


% RGB
RGB = I(:,:,:);
red = RGB(:,:,1);
green = RGB(:,:,2);
blue = RGB(:,:,3);

% LAB
LAB = rgb2lab(RGB);
L = LAB(:,:,1);
a = LAB(:,:,2);
b = LAB(:,:,3);

% % Entropy
% neighborhoods1 = cell(x,y);
% entropy1Vals = cell(x,y);
% neighborhoods2 = cell(x,y);
% entropy2Vals = cell(x,y);
% neighborhoods3 = cell(x,y);
% entropy3Vals = cell(x,y);
%
% for i = 9:(x-9)
%     for j = 9:(y-9)
%         neighborhoods1{i,j} = L((i-2):(i+2), (j-2):(j+2));
%         entropy1Vals{i,j} = entropy(neighborhoods1{i,j});
%         neighborhoods2{i,j} = L((i-4):(i+4), (j-4):(j+4));
%         entropy2Vals{i,j} = entropy(neighborhoods2{i,j});
%         neighborhoods3{i,j} = L((i-8):(i+8), (j-8):(j+8));
%         entropy3Vals{i,j} = entropy(neighborhoods3{i,j});
%     end
% end

% Texture
texture_sigma1_theta1 = cell(x,y);
texture_sigma1_theta2 = cell(x,y);
texture_sigma1_theta3 = cell(x,y);
texture_sigma1_theta4 = cell(x,y);
texture_sigma1_theta5 = cell(x,y);
texture_sigma1_theta6 = cell(x,y);
texture_sigma2_theta1 = cell(x,y);
texture_sigma2_theta2 = cell(x,y);
texture_sigma2_theta3 = cell(x,y);
texture_sigma2_theta4 = cell(x,y);
texture_sigma2_theta5 = cell(x,y);
texture_sigma2_theta6 = cell(x,y);
texture_sigma3_theta1 = cell(x,y);
texture_sigma3_theta2 = cell(x,y);
texture_sigma3_theta3 = cell(x,y);
texture_sigma3_theta4 = cell(x,y);
texture_sigma3_theta5 = cell(x,y);
texture_sigma3_theta6 = cell(x,y);

% 6 thetas (0, 30, 60, 90, 120, 150)
theta = zeros(6, 1);
for i = 1:6
    theta(i,1) = (i-1)*30;
end

% Filter Banks
filter_sigma1_theta1 = zeros(19, 19);
filter_sigma1_theta2 = zeros(19, 19);
filter_sigma1_theta3 = zeros(19, 19);
filter_sigma1_theta4 = zeros(19, 19);
filter_sigma1_theta5 = zeros(19, 19);
filter_sigma1_theta6 = zeros(19, 19);
filter_sigma2_theta1 = zeros(19, 19);
filter_sigma2_theta2 = zeros(19, 19);
filter_sigma2_theta3 = zeros(19, 19);
filter_sigma2_theta4 = zeros(19, 19);
filter_sigma2_theta5 = zeros(19, 19);
filter_sigma2_theta6 = zeros(19, 19);
filter_sigma3_theta1 = zeros(19, 19);
filter_sigma3_theta2 = zeros(19, 19);
filter_sigma3_theta3 = zeros(19, 19);
filter_sigma3_theta4 = zeros(19, 19);
filter_sigma3_theta5 = zeros(19, 19);
filter_sigma3_theta6 = zeros(19, 19);

% Filter Banks: sigma = 1 and 6 thetas
sigma = 1;
for i = -9:9
    for j = -9:9
        G0_0 = (1/(2*pi*sigma^2)).*exp(-(i.^2+j.^2)/(2*sigma^2));
        G2_90 = ((j.^2-sigma^2)/(2*pi*sigma^6)).*exp(-(i.^2+j.^2)/(2*sigma^2));

        filter_sigma1_theta1(i+10, j+10) = G0_0*cos(theta(1)) + G2_90*sin(theta(1));
        filter_sigma1_theta2(i+10, j+10) = G0_0*cos(theta(2)) + G2_90*sin(theta(2));
        filter_sigma1_theta3(i+10, j+10) = G0_0*cos(theta(3)) + G2_90*sin(theta(3));
        filter_sigma1_theta4(i+10, j+10) = G0_0*cos(theta(4)) + G2_90*sin(theta(4));
        filter_sigma1_theta5(i+10, j+10) = G0_0*cos(theta(5)) + G2_90*sin(theta(5));
        filter_sigma1_theta6(i+10, j+10) = G0_0*cos(theta(6)) + G2_90*sin(theta(6));
    end
end

% Filter Banks: sigma = sqrt(2) and 6 thetas
sigma = sqrt(2);
for i = -9:9
    for j = -9:9
        G0_0 = (1/(2*pi*sigma^2)).*exp(-(i.^2+j.^2)/(2*sigma^2));
        G2_90 = ((j.^2-sigma^2)/(2*pi*sigma^6)).*exp(-(i.^2+j.^2)/(2*sigma^2));

        filter_sigma2_theta1(i+10, j+10) = G0_0*cos(theta(1)) + G2_90*sin(theta(1));
        filter_sigma2_theta2(i+10, j+10) = G0_0*cos(theta(2)) + G2_90*sin(theta(2));
        filter_sigma2_theta3(i+10, j+10) = G0_0*cos(theta(3)) + G2_90*sin(theta(3));
        filter_sigma2_theta4(i+10, j+10) = G0_0*cos(theta(4)) + G2_90*sin(theta(4));
        filter_sigma2_theta5(i+10, j+10) = G0_0*cos(theta(5)) + G2_90*sin(theta(5));
        filter_sigma2_theta6(i+10, j+10) = G0_0*cos(theta(6)) + G2_90*sin(theta(6));
    end
end

% Filter Banks: sigma = sqrt(2) and 6 thetas
sigma = 2;
for i = -9:9
    for j = -9:9
        G0_0 = (1/(2*pi*sigma^2)).*exp(-(i.^2+j.^2)/(2*sigma^2));
        G2_90 = ((j.^2-sigma^2)/(2*pi*sigma^6)).*exp(-(i.^2+j.^2)/(2*sigma^2));

        filter_sigma3_theta1(i+10, j+10) = G0_0*cos(theta(1)) + G2_90*sin(theta(1));
        filter_sigma3_theta2(i+10, j+10) = G0_0*cos(theta(2)) + G2_90*sin(theta(2));
        filter_sigma3_theta3(i+10, j+10) = G0_0*cos(theta(3)) + G2_90*sin(theta(3));
        filter_sigma3_theta4(i+10, j+10) = G0_0*cos(theta(4)) + G2_90*sin(theta(4));
        filter_sigma3_theta5(i+10, j+10) = G0_0*cos(theta(5)) + G2_90*sin(theta(5));
        filter_sigma3_theta6(i+10, j+10) = G0_0*cos(theta(6)) + G2_90*sin(theta(6));
    end
end

% Convultion of Filter Bank and LChannels
conv_sigma1_theta1 = conv2(filter_sigma1_theta1, L);
conv_sigma1_theta2 = conv2(filter_sigma1_theta2, L);
conv_sigma1_theta3 = conv2(filter_sigma1_theta3, L);
conv_sigma1_theta4 = conv2(filter_sigma1_theta4, L);
conv_sigma1_theta5 = conv2(filter_sigma1_theta5, L);
conv_sigma1_theta6 = conv2(filter_sigma1_theta6, L);
conv_sigma2_theta1 = conv2(filter_sigma2_theta1, L);
conv_sigma2_theta2 = conv2(filter_sigma2_theta2, L);
conv_sigma2_theta3 = conv2(filter_sigma2_theta3, L);
conv_sigma2_theta4 = conv2(filter_sigma2_theta4, L);
conv_sigma2_theta5 = conv2(filter_sigma2_theta5, L);
conv_sigma2_theta6 = conv2(filter_sigma2_theta6, L);
conv_sigma3_theta1 = conv2(filter_sigma3_theta1, L);
conv_sigma3_theta2 = conv2(filter_sigma3_theta2, L);
conv_sigma3_theta3 = conv2(filter_sigma3_theta3, L);
conv_sigma3_theta4 = conv2(filter_sigma3_theta4, L);
conv_sigma3_theta5 = conv2(filter_sigma3_theta5, L);
conv_sigma3_theta6 = conv2(filter_sigma3_theta6, L);

for i = 1:x
    for j= 1:y
        texture_sigma1_theta1{i,j} = conv_sigma1_theta1(i+9, j+9);
        texture_sigma1_theta2{i,j} = conv_sigma1_theta2(i+9, j+9);
        texture_sigma1_theta3{i,j} = conv_sigma1_theta3(i+9, j+9);
        texture_sigma1_theta4{i,j} = conv_sigma1_theta4(i+9, j+9);
        texture_sigma1_theta5{i,j} = conv_sigma1_theta5(i+9, j+9);
        texture_sigma1_theta6{i,j} = conv_sigma1_theta6(i+9, j+9);
        texture_sigma2_theta1{i,j} = conv_sigma2_theta1(i+9, j+9);
        texture_sigma2_theta2{i,j} = conv_sigma2_theta2(i+9, j+9);
        texture_sigma2_theta3{i,j} = conv_sigma2_theta3(i+9, j+9);
        texture_sigma2_theta4{i,j} = conv_sigma2_theta4(i+9, j+9);
        texture_sigma2_theta5{i,j} = conv_sigma2_theta5(i+9, j+9);
        texture_sigma2_theta6{i,j} = conv_sigma2_theta6(i+9, j+9);
        texture_sigma3_theta1{i,j} = conv_sigma3_theta1(i+9, j+9);
        texture_sigma3_theta2{i,j} = conv_sigma3_theta2(i+9, j+9);
        texture_sigma3_theta3{i,j} = conv_sigma3_theta3(i+9, j+9);
        texture_sigma3_theta4{i,j} = conv_sigma3_theta4(i+9, j+9);
        texture_sigma3_theta5{i,j} = conv_sigma3_theta5(i+9, j+9);
        texture_sigma3_theta6{i,j} = conv_sigma3_theta6(i+9, j+9);
    end
end

allInfo = [red(:)';
    (green(:)');
    (blue(:)');
    (L(:)');
    (a(:)');
    (b(:)');
    (cell2mat(texture_sigma1_theta1(:)'));
    (cell2mat(texture_sigma1_theta2(:)'));
    (cell2mat(texture_sigma1_theta3(:)'));
    (cell2mat(texture_sigma1_theta4(:)'));
    (cell2mat(texture_sigma1_theta5(:)'));
    (cell2mat(texture_sigma1_theta6(:)'));
    (cell2mat(texture_sigma2_theta1(:)'));
    (cell2mat(texture_sigma2_theta2(:)'));
    (cell2mat(texture_sigma2_theta3(:)'));
    (cell2mat(texture_sigma2_theta4(:)'));
    (cell2mat(texture_sigma2_theta5(:)'));
    (cell2mat(texture_sigma2_theta6(:)'));
    (cell2mat(texture_sigma3_theta1(:)'));
    (cell2mat(texture_sigma3_theta2(:)'));
    (cell2mat(texture_sigma3_theta3(:)'));
    (cell2mat(texture_sigma3_theta4(:)'));
    (cell2mat(texture_sigma3_theta5(:)'));
    (cell2mat(texture_sigma3_theta6(:)'))];

% Load in Boosted Trees Classifier
load('classifier.mat');
newStruct = struct('predictFcn', predictFcn, 'ClassificationEnsemble', ClassificationEnsemble, 'About', About, 'HowToPredict', HowToPredict);

values = newStruct.predictFcn(allInfo);
valueArray = reshape(values, x, y);

valueArray=(valueArray-2)*-255;
image(valueArray);
percent = sum(sum(valueArray))/255/(x*y);
strcat(num2str(percent*100), '% is foliage')
