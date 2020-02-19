I = imread('img/reference.png');
[x,y,z] = size(I);

red = I(:,:,1);
green = I(:,:,2);
bblue = I(:,:,3);

lab = rgb2lab(I);
L = lab(:,:,1);
A = lab(:,:,2);
B = lab(:,:,3);

neighborhoods1 = cell(x,y);
entropy1Vals = cell(x,y);
neighborhoods2 = cell(x,y);
entropy2Vals = cell(x,y);
neighborhoods3 = cell(x,y);
entropy3Vals = cell(x,y);

for i = 1:x
    for j = 1:y
        neighbors(1) = img(r-1,c-1);    % Upper left.  r = row, c = column.
        neighbors(2) = img(r-1,c);      % Upper middle.  r = row, c = column.
        neighbors(3) = img(r-1,c+1);    % Upper right.  r = row, c = column.
        neighbors(4) = img(r,c-1);      % left.  r = row, c = column.
        neighbors(5) = img(r,c+1);      % right. r = row, c = column.
        neighbors(6) = img(r+1,c+1);    % Lowerleft.  r = row, c = column.
        neighbors(7) = img(r+1,c);      % lower middle.  r = row, c = column.
        neighbors(8) = img(r+1,c-1);    % Lower left.  r = row, c = column.
    end
end
