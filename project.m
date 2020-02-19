% Image and Sample Size
I = imread('img/trees1.png');
samples = 75;

% Locations
randX = xlsread('data/treeplaces.xlsx', 'A62:BW62');
randY = xlsread('data/treeplaces.xlsx', 'A63:BW63');
pixel = [randX; randY];
% RGB
pixelRGB = cell(1, samples);
red = cell(1, samples);
green = cell(1, samples);
blue = cell(1, samples);
% L* A* B*
pixelLAB = cell(1, samples);
L = cell(1, samples);
a = cell(1, samples);
b = cell(1, samples);
% Entropy
neighborhoods1 = cell(1, samples);
entropy1Vals = cell(1, samples);
neighborhoods2 = cell(1, samples);
entropy2Vals = cell(1, samples);
neighborhoods3 = cell(1, samples);
entropy3Vals = cell(1, samples);
% Texture
texture_sigma1_theta1 = cell(1, samples);
texture_sigma1_theta2 = cell(1, samples);
texture_sigma1_theta3 = cell(1, samples);
texture_sigma1_theta4 = cell(1, samples);
texture_sigma1_theta5 = cell(1, samples);
texture_sigma1_theta6 = cell(1, samples);
texture_sigma2_theta1 = cell(1, samples);
texture_sigma2_theta2 = cell(1, samples);
texture_sigma2_theta3 = cell(1, samples);
texture_sigma2_theta4 = cell(1, samples);
texture_sigma2_theta5 = cell(1, samples);
texture_sigma2_theta6 = cell(1, samples);
texture_sigma3_theta1 = cell(1, samples);
texture_sigma3_theta2 = cell(1, samples);
texture_sigma3_theta3 = cell(1, samples);
texture_sigma3_theta4 = cell(1, samples);
texture_sigma3_theta5 = cell(1, samples);
texture_sigma3_theta6 = cell(1, samples);

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

% LChannels
LabChannel = rgb2lab(I);
LChannel = LabChannel(:,:,1);

% Convultion of Filter Bank and LChannels
conv_sigma1_theta1 = conv2(filter_sigma1_theta1, LChannel);
conv_sigma1_theta2 = conv2(filter_sigma1_theta2, LChannel);
conv_sigma1_theta3 = conv2(filter_sigma1_theta3, LChannel);
conv_sigma1_theta4 = conv2(filter_sigma1_theta4, LChannel);
conv_sigma1_theta5 = conv2(filter_sigma1_theta5, LChannel);
conv_sigma1_theta6 = conv2(filter_sigma1_theta6, LChannel);
conv_sigma2_theta1 = conv2(filter_sigma2_theta1, LChannel);
conv_sigma2_theta2 = conv2(filter_sigma2_theta2, LChannel);
conv_sigma2_theta3 = conv2(filter_sigma2_theta3, LChannel);
conv_sigma2_theta4 = conv2(filter_sigma2_theta4, LChannel);
conv_sigma2_theta5 = conv2(filter_sigma2_theta5, LChannel);
conv_sigma2_theta6 = conv2(filter_sigma2_theta6, LChannel);
conv_sigma3_theta1 = conv2(filter_sigma3_theta1, LChannel);
conv_sigma3_theta2 = conv2(filter_sigma3_theta2, LChannel);
conv_sigma3_theta3 = conv2(filter_sigma3_theta3, LChannel);
conv_sigma3_theta4 = conv2(filter_sigma3_theta4, LChannel);
conv_sigma3_theta5 = conv2(filter_sigma3_theta5, LChannel);
conv_sigma3_theta6 = conv2(filter_sigma3_theta6, LChannel);

for i = 1:length(randX)
   % RGB
   pixelRGB{i} = impixel(I, randX(i), randY(i));
   red{i} = (pixelRGB{i}(1));
   green{i} = (pixelRGB{i}(2));
   blue{i} = (pixelRGB{i}(3));
   % L* A* B*
   pixelLAB{i} = rgb2lab(pixelRGB{i}/255);
   L{i} = (pixelLAB{i}(1));
   a{i} = (pixelLAB{i}(2));
   b{i} = (pixelLAB{i}(3));
   % Entropy
%    neighborhoods1{i} = imcrop(LChannel, [randX(i)-2 randY(i)-2 5 5]);
   neighborhoods1{i} = LChannel((randX(i)-2):(randX(i)+2), (randY(i)-2):(randY(i)+2));
   entropy1Vals{i} = entropy(neighborhoods1{i});
%    neighborhoods2{i} = imcrop(LChannel, [randX(i)-4 randY(i)-4 9 9]);
   neighborhoods2{i} = LChannel((randX(i)-4):(randX(i)+4), (randY(i)-4):(randY(i)+4));
   entropy2Vals{i} = entropy(neighborhoods2{i});
%    neighborhoods3{i} = imcrop(LChannel, [randX(i)-8 randY(i)-8 17 17]);
   neighborhoods3{i} = LChannel((randX(i)-8):(randX(i)+8), (randY(i)-8):(randY(i)+8));
   entropy3Vals{i} = entropy(neighborhoods3{i});
   % Texture
   texture_sigma1_theta1{i} = conv_sigma1_theta1((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma1_theta2{i} = conv_sigma1_theta2((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma1_theta3{i} = conv_sigma1_theta3((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma1_theta4{i} = conv_sigma1_theta4((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma1_theta5{i} = conv_sigma1_theta5((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma1_theta6{i} = conv_sigma1_theta6((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma2_theta1{i} = conv_sigma2_theta1((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma2_theta2{i} = conv_sigma2_theta2((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma2_theta3{i} = conv_sigma2_theta3((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma2_theta4{i} = conv_sigma2_theta4((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma2_theta5{i} = conv_sigma2_theta5((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma2_theta6{i} = conv_sigma2_theta6((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma3_theta1{i} = conv_sigma3_theta1((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma3_theta2{i} = conv_sigma3_theta2((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma3_theta3{i} = conv_sigma3_theta3((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma3_theta4{i} = conv_sigma3_theta4((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma3_theta5{i} = conv_sigma3_theta5((pixel(2, i)+9), (pixel(1, i)+9));
   texture_sigma3_theta6{i} = conv_sigma3_theta6((pixel(2, i)+9), (pixel(1, i)+9));
end

% 27 dimensional feature vector for color, texture, and entropy values
%       Color features (1-6):       R,G,B and L*,A*,B*
%       Entropy features (7-9):     5�5, 9�9, and 17�17 search windows
%       Texture features (10-27):   3 sigmas (sigma= 1, sqrt(2), 2) and 6 orientations (0, 30, 60, 90, 120, 150)
features = [red; green; blue; L; a; b;
            entropy1Vals; entropy2Vals; entropy3Vals;
            texture_sigma1_theta1; texture_sigma1_theta2; texture_sigma1_theta3; texture_sigma1_theta4; texture_sigma1_theta5; texture_sigma1_theta6;
            texture_sigma2_theta1; texture_sigma2_theta2; texture_sigma2_theta3; texture_sigma2_theta4; texture_sigma2_theta5; texture_sigma2_theta6;
            texture_sigma3_theta1; texture_sigma3_theta2; texture_sigma3_theta3; texture_sigma3_theta4; texture_sigma3_theta5; texture_sigma3_theta6];
