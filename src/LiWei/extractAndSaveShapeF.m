function extractAndSaveShapeF()

sizeH = 200; %
sizeW = 200; %
Ln = 4; %level num: 4  8  12  16
lineN = 8; %line num of different angles: 4  8  16  32  which means the ray: 8  16  32  64
lineMask = gen_line_mask(lineN); %lineMask is a [n, 401, 401] matrix

% readPath = 'F:\liwei\DB\BigROI\BigSub3D_356_10_correct_I\';
readPath = 'F:\liwei\DB\BigROI\BigSub3D_356_10_correct_I_Sel\';
savePath = 'F:\liwei\DB\BigROI\shapeF\F_Level_4_line_8\';

FL = Ln+1; %feature length 1+8
fileList = dir([readPath '*.dat']);
fileNum = size(fileList, 1);

for i = 1:fileNum    
    level_A = zeros(1, Ln); %store the area of every level
    ray_len = zeros(Ln, lineN*2); %store the ray length of every level
    ray_flag = ones(Ln, lineN*2); %if the ray line touch the boundary, the flag = 0, else flag = 1.
    
    fn = [readPath fileList(i).name];
    fid = fopen(fn, 'r');
    Z = fread(fid, [sizeH,sizeW], 'double');
    fclose(fid);
    
    %%%get the ray and area feature%%%%%%%%%%%%%%%%%%%
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
    palmH = min(min(rectforMin));
    
    %find the level regions from 0 to deepest point 
    %find the region > 0
    palmH = -palmH;
    Z = -Z';
    ind = find(Z == palmH);

    step = palmH / Ln;
    levelH = [0:step:palmH];
    for k = 1:Ln
        L0 = zeros(sizeH, sizeW);
        L0( find( Z>=levelH(Ln-k+1) ) ) = 1;
        %the 1st level
        if k==1
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
            mx = round(mean(x)); %mx will be used in the following
            my = round(mean(y)); %my will be used in the following
            level_A(k) = length(find(L0));

            %find the max length ray 
            tempL = zeros(1, lineN);
            tempRay = zeros(1, lineN*2);
            tempflag = ones(1, lineN*2);
            for j = 1:lineN
                line_200 = zeros(sizeH, sizeW);
                line_200(:,:) = lineMask(j, 201-mx+1:201-mx+1+199, 201-my+1:201-my+1+199);
                tempIm = L0 & line_200;
                [x, y] = find(tempIm);
                tempL(j) = length(x);
                %get every ray length
                min_x = min(x); min_y = min(y); max_x = max(x); max_y = max(y);
                if j<lineN/4
                    tempRay(j) = max_y - my;
                    if max_y == 200 || min_x == 1
                        tempflag(j) = 0;
                        tempRay(j) = 0;
                    end
                    tempRay(j+lineN) = my - min_y;
                    if min_y == 1 || max_x == 200
                        tempflag(j+lineN) = 0;
                        tempRay(j+lineN) = 0;
                    end
                elseif j<2*lineN/4
                    tempRay(j) = mx - min_x;
                    if max_y == 200 || min_x == 1
                        tempflag(j) = 0;
                        tempRay(j) = 0;
                    end
                    tempRay(j+lineN) = max_x - mx;
                    if min_y == 1 || max_x == 200
                        tempflag(j+lineN) = 0;
                        tempRay(j+lineN) = 0;
                    end
                elseif j<3*lineN/4
                    tempRay(j) = mx - min_x;
                    if min_x == 1 || min_y == 1
                        tempflag(j) = 0;
                        tempRay(j) = 0;
                    end
                    tempRay(j+lineN) = max_x - mx;
                    if max_x == 200 || max_y == 200
                        tempflag(j+lineN) = 0;
                        tempRay(j+lineN) = 0;
                    end
                else
                    tempRay(j) = my - min_y;
                    if min_y == 1 || min_x == 1
                        tempflag(j) = 0;
                        tempRay(j) = 0;
                    end
                    tempRay(j+lineN) = max_y - my;
                    if max_y == 200 || max_x == 200
                        tempflag(j+lineN) = 0;
                        tempRay(j+lineN) = 0;
                    end
                end
            end
            [maxL, maxI] = max(tempL); %maxI will be used in the following

            %store the ray length from the start direction (max length ray)
            %start direction is between -pi/2 and pi/2
            if maxI>lineN/2
                temp = tempRay(maxI+lineN : end);
                temp = [temp tempRay(1:maxI+lineN-1)];
                ray_len(k, :) = temp(:);

                temp = tempflag(maxI+lineN : end);
                temp = [temp tempflag(1:maxI+lineN-1)];
                ray_flag(k, :) = temp(:);
            else
                temp = tempRay(maxI : end);
                temp = [temp tempRay(1:maxI-1)];
                ray_len(k, :) = temp(:);

                temp = tempflag(maxI : end);
                temp = [temp tempflag(1:maxI-1)];
                ray_flag(k, :) = temp(:);
            end
        else
            %dilate the Lp (previous level) and than & with L0
            se = strel('disk',round(35 - 3*k*(8.0/Ln)));  
            L1 = imdilate(Lp, se);
            L0 = L0 & L1;
            Lp = L0;   
            level_A(k) = length(find(L0));

            %find the max length ray 
            tempL = zeros(1, lineN);
            tempRay = zeros(1, lineN*2);
            tempflag = ones(1, lineN*2);
            for j = 1:lineN
                line_200 = zeros(sizeH, sizeW);
                line_200(:,:) = lineMask(j, 201-mx+1:201-mx+1+199, 201-my+1:201-my+1+199);
                tempIm = L0 & line_200;
                [x, y] = find(tempIm);

                %get every ray length
                min_x = min(x); min_y = min(y); max_x = max(x); max_y = max(y);
                if j<lineN/4
                    tempRay(j) = max_y - my;
                    if max_y == 200 || min_x == 1
                        tempflag(j) = 0;
                        tempRay(j) = 0;
                    end
                    tempRay(j+lineN) = my - min_y;
                    if min_y == 1 || max_x == 200
                        tempflag(j+lineN) = 0;
                        tempRay(j+lineN) = 0;
                    end
                elseif j<2*lineN/4
                    tempRay(j) = mx - min_x;
                    if max_y == 200 || min_x == 1
                        tempflag(j) = 0;
                        tempRay(j) = 0;
                    end
                    tempRay(j+lineN) = max_x - mx;
                    if min_y == 1 || max_x == 200
                        tempflag(j+lineN) = 0;
                        tempRay(j+lineN) = 0;
                    end
                elseif j<3*lineN/4
                    tempRay(j) = mx - min_x;
                    if min_x == 1 || min_y == 1
                        tempflag(j) = 0;
                        tempRay(j) = 0;
                    end
                    tempRay(j+lineN) = max_x - mx;
                    if max_x == 200 || max_y == 200
                        tempflag(j+lineN) = 0;
                        tempRay(j+lineN) = 0;
                    end
                else
                    tempRay(j) = my - min_y;
                    if min_y == 1 || min_x == 1
                        tempflag(j) = 0;
                        tempRay(j) = 0;
                    end
                    tempRay(j+lineN) = max_y - my;
                    if max_y == 200 || max_x == 200
                        tempflag(j+lineN) = 0;
                        tempRay(j+lineN) = 0;
                    end
                end
            end

            %store the ray length from the start direction (max length ray)
            %start direction is between -pi/2 and pi/2
            if maxI>lineN/2
                temp = tempRay(maxI+lineN : end);
                temp = [temp tempRay(1:maxI+lineN-1)];
                ray_len(k, :) = temp(:);

                temp = tempflag(maxI+lineN : end);
                temp = [temp tempflag(1:maxI+lineN-1)];
                ray_flag(k, :) = temp(:);
            else
                temp = tempRay(maxI : end);
                temp = [temp tempRay(1:maxI-1)];
                ray_len(k, :) = temp(:);

                temp = tempflag(maxI : end);
                temp = [temp tempflag(1:maxI-1)];
                ray_flag(k, :) = temp(:);
            end
        end    
    end
    H_A = [palmH level_A];
    ray_flag = int8(ray_flag);
    %%%%%%%%%end of get the ray and area feature%%%%%%%%%%%%%%%%%%%%%%%
    
    %save features
    savename = [savePath fileList(i).name];
    fid = fopen(savename, 'wb');
    fwrite(fid, H_A, 'double');  %H_A is 1*9
    fwrite(fid, ray_len', 'double'); %ray_len is 8*32
    fwrite(fid, ray_flag', 'int8'); %ray_flag is 8*32
    fclose(fid);
    disp(i);
end