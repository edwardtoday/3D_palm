function extractShapeF()

sizeH = 200; %
sizeW = 200; %
lnum = 8;

pixDis = 0.28;
% load facesMatrix64.mat;
temp = pixDis * ((sizeH-2)/2+0.5);
[X,Y] = meshgrid(-temp : pixDis : temp);

fid = fopen('F:\myPalm\BigSub3D_correct\Sub3D_I_3_0.dat', 'r');
% fid = fopen('Sub3D_II_100_0.dat', 'r');
Z = fread(fid, [sizeH,sizeW], 'double');
fclose(fid);
[fx, fy] = gradient(Z);
fxy = fx.^2 + fy.^2;
noiseP = find(fxy>0.1); % 0.1 used for corrected and smoothed Sub3D, 1 used for original Sub3D
flag = ones(sizeH, sizeW); %1 for valid point, 0 for invalid point
flag(noiseP) = 0;
flag = 1 - flag;
se = strel('disk',5);  
flag = imdilate(flag, se);
flag = 1 - flag;
noiseP = find(flag==0);
meanZ = sum(sum(Z .* flag)) / (sizeH*sizeW - length(noiseP));
Z(noiseP) = 5;
%neend't smooth
% Z = smooth(Z, 7);

%%for show the mask which get rid of bad quality region
% flag(1,:) = 0; flag(end,:) = 0;
% flag(:,1) = 0; flag(:,end) = 0;
% imshow(flag');

%calculate the reference 0 plane
% use Z(6:35, 65:136) for calculate the mean value as reference 0
refRect = Z(6:35, 65:136);
refFlag = flag(6:35, 65:136);
refVal = sum(sum(refRect .* refFlag)) / sum(sum(refFlag));
Z = Z - refVal;

%search the min point in Z(65:190, 41:160)
rectforMin = Z(65:190, 41:160);
rectforMin(1:35, 1:25) = 5;
rectforMin(1:35, end-25:end) = 5;
% flagforMin = Z(65:190, 41:160);
palmH = min(min(rectforMin));
ind = find(Z == palmH);
% Z(ind) = 5;

% %find the level regions from 0 to deepest point 
% %find the region > 0
% palmH = -palmH;
% Z = -Z;
% Ln = 8;
% step = palmH / Ln;
% levelH = [0:step:palmH];
% for i = 1:Ln
%     L0 = zeros(sizeH, sizeW);
%     L0( find( Z>=levelH(Ln-i+1) ) ) = 1;
%     %the 1st level
%     if i==1
%         [L,num] = bwlabel(L0);
%         if num > 1
%             for j = 1:num
%                 indL = find(L==j);
%                 if length(find(indL==ind)) > 0
%                     L0(:,:) = 0;
%                     L0(indL) = 1;
%                 end
%             end
%         end
%         Lp = L0;
%         L0 = logical(L0');  
%         [x, y] = find(L0);
%         temp = find(L0);
%         mx = mean(x);
%         my = mean(y);
% 
% %         saveIm(L0, i);
%         figure;
%         imshow(L0);
%     else
%         %dilate the Lp (previous level) and than & with L0
%         se = strel('disk',35 - 3*i);  
%         L1 = imdilate(Lp, se);
%         L0 = L0 & L1;
%         Lp = L0;
%         L0 = L0';   
% %         saveIm(L0, i);
%         figure;
%         imshow(L0);
%     end    
% end

%for show
im_level = zeros(sizeH, sizeW); %for show the levels
%find the level regions from 0 to deepest point 
%find the region > 0
palmH = -palmH;
Z = -Z';
Ln = 8;
step = palmH / Ln;
levelH = [0:step:palmH];
for i = 1:Ln
    L0 = zeros(sizeH, sizeW);
    L0( find( Z>=levelH(Ln-i+1) ) ) = 1;
    %the 1st level
    if i==1
        [L,num] = bwlabel(L0);
        if num > 1
            for j = 1:num
                indL = find(L==j);
                if length(find(indL==ind)) > 0
                    L0(:,:) = 0;
                    L0(indL) = 1;
                end
            end
        end
        Lp = L0;
        [x, y] = find(L0);
        temp = find(L0);
        mx = round(mean(x));
        my = round(mean(y));
        im_level(temp) = 155;

%         saveIm(L0, i);
%         figure;
%         imshow(L0);
    else
        %dilate the Lp (previous level) and than & with L0
        se = strel('disk',35 - 3*i);  
        L1 = imdilate(Lp, se);
%         figure; imshow(Lp);
%         figure; imshow(L1);
        L0 = L0 & L1;
%         figure; imshow(L0);
        Ltemp = L0 - Lp;
%         figure; imshow(Ltemp);
        Lp = L0; 
        temp = find(Ltemp);
        im_level(temp) = 255-(i-1)*30;
%         saveIm(L0, i);
%         figure;
%         imshow(L0);
    end    
end

im_level = uint8(im_level);
% imshow(im_level);
% imwrite(im_level, 'F:\report\2010-1-1\pic\im_levels.bmp');

line_401 = imread('F:\report\2010-1-1\pic\line_4\line_all.bmp');
line_200 = line_401(201-mx+1:201-mx+1+199, 201-my+1:201-my+1+199);
% imshow(line_200);
temp = find(line_200);
im_level(temp) = 20;
im_level(mx, my) = 255;
figure;
imshow(im_level);
imwrite(im_level, 'F:\report\2010-1-1\pic\im_levels_line_4.bmp');
a = 0;

% %%% Z(6:35, 65:136) = 0;
% Z(6:35, 65) = 1;
% Z(6:35, 136) = 1;
% Z(6, 65:136) = 1;
% Z(35, 65:136) = 1;

% figure;
% mesh(X,Y,Z);
% daspect([1 1 1]);
% view([90 90]);

a = 0;

%-----------------------------------------------------------------
function saveIm(data, nlevel)
filename = ['G:\Area\Sub3D_I_4_9', '_' num2str(nlevel), '.bmp'];
imwrite(data, filename);


function [z] = smooth(z, fs)
% z -- 128*128
% n -- filter size
t = floor(fs/2)-1;
for i = 0:t;
    [m,n] = size(z);
    z = [z(:, 2*i+1) z(:,:) z(:, n-2*i)];
end

for i = 0:t;
    [m,n] = size(z);
    z = [z(2*i+1, :); z(:,:); z(m-2*i, :)];
end

h = ones(fs,fs)/(fs*fs);
z = filter2(h, z, 'valid');