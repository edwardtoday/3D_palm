function [MD, HCA] = extract_global_feature(zonly_file, N)    
    %% Read binary data
    datfile = fopen(zonly_file);
    matdata = fread(datfile, [768,576], 'single');
    fclose(datfile);

    %% Crop ROI, reference plane and max region
    roi_size = 400;
    roi = matdata(235:(234+roi_size),69:(68+roi_size));

    refplane = roi(1:(roi_size*0.25),(roi_size*0.25):(roi_size*0.9));
    maxregion = roi((roi_size*0.4):(roi_size*0.95),(roi_size*0.3):(roi_size*0.7));

    % mesh(roi)
    % contourf(roi)

    %% Calculate max depth (MD)
    d_ref = max(max(refplane));
    d_max = min(min(maxregion));

    MD = d_ref - d_max;

    %% Discretize ROI to different levels
    num_levels = N;
    step = MD/(num_levels+1);
    G_k = zeros(roi_size,roi_size,'int32');

    for row = 1:roi_size
        for col = 1:roi_size
            d_pixel = roi(row,col);
            if(d_pixel > d_ref || d_pixel < d_max)
                G_k(row,col) = -1;
            elseif(d_pixel == d_max) 
                if(row >= roi_size * 0.4 && row <= roi_size * 0.95 && col >= roi_size * 0.3 && col <= roi_size * 0.7)
                    G_k(row,col) = 0;
                else
                    G_k(row,col) = -1;
                end
            else
                G_k(row,col) = 1 + (d_pixel-d_max)/step;
            end
        end
    end
    
    %% Calculate horizontal cross-sectional area (HCA)
    L = zeros(roi_size,roi_size,num_levels+1,'int32');
    L_k = zeros(roi_size,roi_size,'int32');
    for k = 0:num_levels
        result = (G_k==k);
        if(k > 0)
            se = strel('disk',70-5*k,4);
            result = result & (imdilate(L(:,:,k),se));
        end
        L(:,:,k+1) = int32(result);
        L_k = L_k + int32(result)*k;
    end
    % contourf(L_k);

    HCA = zeros(num_levels,1,'int32');
    parfor k = 1:num_levels
        HCA(k) = sum(sum(L_k==k));
    end
    % HCA

    %% Calculate radial line length (RLL)

end