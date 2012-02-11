function palmShapeF_H_A()

% readpath = 'F:\liwei\DB\BigROI\BigSub3D_356_20_correct\';
% writepath = 'F:\liwei\DB\BigROI\BigSub3DNew\';
% savepath = 'F:\liwei\DB\BigROI\results\';

% fidtxt = fopen([savepath 'deepInf.txt'], 'wt');

sizeH = 200; %
sizeW = 200; %
Ln = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

readPath = 'F:\myPalm\BigSub3D_correct\';

FL = Ln+1; %feature length 1+8
fileList = dir([readPath '*.dat']);
fileNum = size(fileList, 1);
fAll = ones(fileNum, FL);

nameIndex = zeros(fileNum, 3);
for i = 1:fileNum
    if(fileList(i).name(8) == 'I')
        nameIndex(i, 1) = 2;
    else
        nameIndex(i, 1) = 1;
    end
    tpos = findstr(fileList(i).name, '_');
    nameIndex(i, 2) = str2num( fileList(i).name(tpos(2)+1 : tpos(3)-1) );
    nameIndex(i, 3) = str2num( fileList(i).name(tpos(3)+1) );
    
    
    fn = [readPath fileList(i).name];
    fid = fopen(fn, 'r');
    Z = fread(fid, [sizeH,sizeW], 'double');
    fAll(i, :) = getPalmH_A(Z);
    fclose(fid);
end

ct1 = clock;

ClientScore = [];
ImpostorScore = [];
ClIndex = [];
ImIndex = [];
for i = 1:fileNum-1
    for j = i+1:fileNum
        f1 = zeros(1, Ln);
        f2 = zeros(1, Ln);
        f1 = fAll(i, :);
        f2 = fAll(j, :);
        
        currDis = abs(f1 - f2);

        if(nameIndex(i,2) == nameIndex(j,2))
            ClientScore = [ClientScore; currDis];
            ClIndex = [ClIndex; nameIndex(i,:) nameIndex(j,:) currDis];
        else
            ImpostorScore = [ImpostorScore; currDis];
            ImIndex = [ImIndex; nameIndex(i,:) nameIndex(j,:) currDis];
        end
    end
end
save('ClientScore', 'ClientScore');
save('ImpostorScore', 'ImpostorScore');
save('ClIndex', 'ClIndex');
save('ImIndex', 'ImIndex');
for i = 1:FL
    plotROC(ClientScore(:,i),ImpostorScore(:,i));
end

ct2 = clock;
ct = ct2-ct1

MyEnd = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H_A] = getPalmH_A(Z)

Ln = 8;
level_A = zeros(1, Ln);
sizeH = 200; %
sizeW = 200; %
[fx, fy] = gradient(Z);
fxy = fx.^2 + fy.^2;
noiseP = find(fxy>0.1);
flag = ones(sizeH, sizeW); %1 for valid point, 0 for invalid point
flag(noiseP) = 0;
meanZ = sum(sum(Z .* flag)) / (sizeH*sizeW - length(noiseP));
Z(noiseP) = 5;
%neend't smooth
% Z = smooth(Z, 7);

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

%find the level regions from 0 to deepest point 
%find the region > 0
palmH = -palmH;
Z = -Z;

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
        L0 = L0';  
        level_A(i) = length(find(L0));
%         saveIm(L0, i);
    else
        %dilate the Lp (previous level) and than & with L0
        se = strel('disk',35 - 3*i);  
        L1 = imdilate(Lp, se);
        L0 = L0 & L1;
        Lp = L0;
        L0 = L0';   
        level_A(i) = length(find(L0));
%         saveIm(L0, i);
    end    
end

H_A = [palmH level_A];

%------------------------------------------------------------------
function plotROC(ClientScore, ImpostorScore)
% ClientScore = ClientScore;
% ImpostorScore = ImpostorScore;

ClientScore = sort(ClientScore);
ImpostorScore = sort(ImpostorScore);
TotalClient = length(ClientScore);
TotalImpostor = length(ImpostorScore);

count = 0;

count = count+1;
Index = find(ImpostorScore<=ImpostorScore(1));
FAR(count) = length(Index)*100/TotalImpostor;        
Index = find(ClientScore<ImpostorScore(1));        
GAR(count) = length(Index)*100/TotalClient;

for i=-8:0.05:2    
    ThIndex = round(TotalImpostor*(10^i)/100);
    if(ThIndex>1)                
        count = count+1;
        Index = find(ImpostorScore<=ImpostorScore(ThIndex));
        FAR(count) = length(Index)*100/TotalImpostor;
        Index = find(ClientScore<ImpostorScore(ThIndex));
        GAR(count) = length(Index)*100/TotalClient;
    end
end

% compute EER
Diff = abs((100-GAR)-FAR);
[minDiff,idx] = min(Diff);             
[(FAR(idx)+100-GAR(idx))/2, FAR(idx), 100-GAR(idx)];
if(GAR(idx) ~= 100)
    EER = (FAR(idx)+100-GAR(idx))/2  %EER
else
    EER = 0
end

EERcount = 1;
for i=-8:0.05:2    
    ThIndex = round(TotalImpostor*(10^i)/100);
    if(ThIndex>1)                
        EERcount = EERcount+1;
        if(EERcount == idx)
            EERThreshold = ImpostorScore(ThIndex);
        end
    end
end
% d_prime = abs(mean(ImpostorScore)-mean(ClientScore))/sqrt(var(ImpostorScore)/2+var(ClientScore)/2)

% figure, semilogx(FAR,GAR);
% xlabel('False Acceptance Rate')
% ylabel('Genuine Acceptance Rate')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot Distribution
% Step = 0.005;
% ImpostorScore(find(ImpostorScore>3000)) = 3000;
% ClientScore(find(ClientScore>3000)) = 3000;
minScore = min([min(ClientScore),min(ImpostorScore)]);
maxScore = max([max(ClientScore),max(ImpostorScore)]);
Step = (maxScore-minScore)/1000;
X = minScore:Step:maxScore+Step;
n = hist(ClientScore,X);
nfreq = n/sum(n)*100;
figure,plot(X,nfreq,'r-', 'LineWidth',2);
n = hist(ImpostorScore,X);
nfreq = n/sum(n)*100;
hold on
plot(X,nfreq,'g-.', 'LineWidth',2);
xlabel('Distance');
ylabel('Frequence %')
legend('Genuine', 'Impostor');
% text(0.25,4,'Genuine')
% text(0.8,6,'Impostor')
title('Genuine Red Curve; Impostor Green Curve')
aa=0;
