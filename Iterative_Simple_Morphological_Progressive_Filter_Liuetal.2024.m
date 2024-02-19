
%% -2024/02/19
% Author contact, Dr Yinxue Liu, yinxue.liu@outlook.com
% This code is the algorithm presented in the paper, Liu et al., 2024. 
% A bare-earth version of GoogleDEM to simulate flooding in New Delhi, India 

% The DEM is separated into tiles to avoid memory issue
% set overlapping areas between tiles to avoid inconsistent elevation
% buffer_pix = 100



% input include
% path for result file
% building mask 
% smrf based xyz file


% other parameters include smrf, smrf_its1, smrf_its2 parameters: s,w,et,es.
% the current setting smrf_its1_s = smrf_s/2, smrf_its1_w = smrf_w*2, other
% parameters are the same.
% smrf_its2_w = smrf_its_w*2, other parameters are the same.

% change path to the path of the code
cd 'Delhi/smrf_bdmask/'
%path of the output data
path = 'Delhi/Yamuna/smrf_withbd/';
% mkdir(path); 

function [dem_f,obj] = tile_processing(
    x_2d, y_2d, z_2d, obj_ori, ...
    pixel_r, pixel_c, tile_r,tile_c,path, ...
    s, w, et, es, iteration, ...
    info1.GeoTIFFTags.GeoKeyDirectoryTag)

    % prelocate the output array
    arr = spalloc(row,col,row*col);
    arr_obj = spalloc(row,col,row*col);

    for num_r = 1:tile_r
        for num_c = 1:tile_c
            disp(['The row number of the tile is', string(num_r),'and the col number of the tile is', string(num_c)]);
            
            rr_range = pixel_r*(num_r-1)+1-(logical(num_r>1))*buffer_pix:pixel_r*num_r;%rows allowing overlapping 100 pixels
            cc_range = pixel_c*(num_c-1)+1-(logical(num_c>1))*buffer_pix:pixel_c*num_c;%record the end location of column of tile
            
            % if the range beyonds or less than the size of the matrix, in the case
            % residual~=0
            if num_r == tile_r
                if pixel_r*num_r~=row %check the edge of row equals to the length of matrix
                    rr_range = pixel_r*(num_r-1)+1-(logical(num_r>1))*buffer_pix:row;
                end
            end
            if num_c == tile_c
                if pixel_c*num_c~=col %check the edge of column equals to the length of matrix
                    cc_range = pixel_c*(num_c-1)+1-(logical(num_c>1))*buffer_pix:col;
                end
            end

            if (num_r> 1)&&(num_c==1)
                rr_overlap = pixel_r*(num_r-1)+1-buffer_pix:pixel_r*(num_r-1);
                cc_overlap = cc_range;
                
                ii_overlap = repelem(rr_overlap,length(cc_overlap),1)';
                ii_overlap = ii_overlap(:);
                jj_overlap = repelem(cc_overlap,1,length(rr_overlap));
                tile_arr_overlap = sparse(ii_overlap,jj_overlap,1,row,col);
                
            elseif (num_r==1)&&(num_c >1)


                rr_overlap = rr_range;
                cc_overlap = pixel_c*(num_c-1)+1-buffer_pix:pixel_c*(num_c-1);

                ii_overlap = repelem(rr_overlap,length(cc_overlap),1)';
                ii_overlap = ii_overlap(:);
                jj_overlap = repelem(cc_overlap,1,length(rr_overlap));
                tile_arr_overlap = sparse(ii_overlap,jj_overlap,1,row,col);
                
            elseif (num_r>1)&&(num_c>1)


                %large sparse array
                ii = repelem(rr_range,length(cc_range),1)';
                ii = ii(:);
                jj = repelem(cc_range,1,length(rr_range));
                L_arr = sparse(ii,jj,1,row,col);
                
                %small sparse array

                rr_non_overlap = rr_range(buffer_pix+1:end);
                cc_non_overlap = cc_range(buffer_pix+1:end);
                ii_non_overlap = repelem(rr_non_overlap,length(cc_non_overlap),1)';
                ii_non_overlap = ii_non_overlap(:);
                jj_non_overlap = repelem(cc_non_overlap,1,length(rr_non_overlap));
                tile_arr_non_overlap = sparse(ii_non_overlap,jj_non_overlap,1,row,col);
                %difference between two array represent the overlap array 
                tile_arr_overlap = L_arr - tile_arr_non_overlap;
                        
            else

                rr_overlap = 0;
                cc_overlap = 0;
                tile_arr_overlap = [];
            end
            
            %% conduct fetch specific row/col and smrf
            x = reshape((x_2d(rr_range,cc_range)).',[],1);
            y = reshape((y_2d(rr_range,cc_range)).',[],1);
            z = reshape((z_2d(rr_range,cc_range)).',[],1);
            z = double(z);
            
            % bd of overlapping pixels between tiles might need to be separated
            % from the non-ovelapping areas
            bd_t = bd(rr_range,cc_range);
            

            % rd_t = rd(rr_range,cc_range);
            % % get xi yi as the corrdinate boundary
            xi = min(x):0.5:max(x);
            yi = max(y):-0.5:min(y);
            
            % smrf
            [ZI R isobject ZIpro ZImin isObjectCell] = smrf_road(x,y,z,'c',c,'s',s,'w',w,'et',et,'es',es,'cutnet',w,'xi',xi,'yi',yi,...
                'inpaintmethod',4,'objectmask',bd_t)
            geotiffwrite(fullfile(path,strcat('Tile_rc',num2str(num_r),num2str(num_c)...
                ,'s',num2str(s),'_w',num2str(w),'_et',num2str(et),'_es',num2str(es),...
                '_ZIpro.tif')),ZIpro,R,'GeoKeyDirectoryTag', info1.GeoTIFFTags.GeoKeyDirectoryTag);

            geotiffwrite(fullfile(path,strcat('Tile_rc',num2str(num_r),num2str(num_c)...
                ,'s',num2str(s),'_w',num2str(w),'_et',num2str(et),'_es',num2str(es),...
                '_obj.tif')),isObjectCell,R,'GeoKeyDirectoryTag', info1.GeoTIFFTags.GeoKeyDirectoryTag);
            
            % check output value of ZIpro
            if ~isempty(find(ZIpro>300))
                error('There is something wrong with the overlapping pixels!')
            end

            %% save the output tiles and combine back to the full image
            z_smrf = reshape(ZIpro,[],1);
            ii = repelem(rr_range,length(cc_range),1)';
            ii = ii(:);
            jj = repelem(cc_range,1,length(rr_range));
            tile_arr = sparse(ii,jj,z_smrf,row,col);
            
            % create sparse array at overlapping pixels

            % find the index at array and current tile_array overlap 
            if ~isempty(tile_arr_overlap)
                idx = find(tile_arr_overlap==1);
                arr(idx) = 0;
                arr_obj(idx) = 0;
            end

            % add the tile_arr to arr
            arr = arr+tile_arr;
            
            %isObjectCell
            obj = reshape(isObjectCell,[],1);

            arr_obj = arr_obj+sparse(ii(:),jj,obj,row,col);
            
            %% R
            if (num_r==1) &&(num_c==1)
                x_corner = R(3,1);
                y_corner = R(3,2);
            end
            
            % update the value of pixels with smrf-applied
            arr_final = full(arr);    
            z_2d(rr_range,cc_range) = arr_final(rr_range,cc_range);
            clear arr_final;
            
        end
    end

    %% save ZIpro, isObjectCell from above
    ZIpro = full(arr);
    isObjectCell = full(arr_obj);
    R(3,1) = x_corner;
    R(3,2) = y_corner;

    % write the filtered array to tif
    dem_f = fullfile(path,strcat('s',num2str(s),'_w',num2str(w),'_et',num2str(et),'_es',num2str(es),...
        '_ZIpro.tif'))
    geotiffwrite(dem_f,ZIpro,R,'GeoKeyDirectoryTag', info1.GeoTIFFTags.GeoKeyDirectoryTag)
    % write the object array to tif
    obj = fullfile(path,strcat('s',num2str(s),'_w',num2str(w),'_et',num2str(et),'_es',num2str(es),...
        '_obj.tif'))
    geotiffwrite(obj,isObjectCell,R,'GeoKeyDirectoryTag', info1.GeoTIFFTags.GeoKeyDirectoryTag)
end

%define the base parameter for filter
c = 0.5;
s = 0.05; % a slightly larger start and decrease with iterations
w = 15;  % a slightly smaller start and increase with iterations
et = 0.1; % the same value for all the steps of the filter
es = 1.25; % the same value for all the steps of the filter

% to be processed GoogleDEM in Geotiff format
referfilename = 'Delhi/Yamuna/GoogleDEM.tif';
info1 = geotiffinfo(referfilename);% % get the georeference from tif file
z_2d = imread(referfilename)
n = ndims(z_2d);
z_2d = double(z_2d);

if n==3
    [row col band]=size(z_2d);
elseif n==2
    [row col] = size(z_2d);
end

% to be processed GoogleDEM in xyz format
filename = 'Delhi/Yamuna/Yamuna_mosaic_masked.xyz';
M = dlmread(filename);
x_2d = reshape(M(:,1),[col,row]).';
y_2d = reshape(M(:,2),[col,row]).';
clear M;

%% chunk the data into tiles
pixel_r = 5000;% 5000 by 5000 pixels of each run
pixel_c = 5000;
% buffer pixels between tiles
buffer_pix = 100; %temporary value
thres = 1000; %exceeding this number will define a new tile of the residual pixels
% decide numbers of tile
r_residual = rem(row,pixel_r);
c_residual = rem(col,pixel_c);
tile_r = fix(row/pixel_r);
tile_c = fix(col/pixel_c);
if (r_residual>thres)
    tile_r = tile_r+1;
end
if (c_residual>thres)
    tile_c = tile_c +1;
end

%% smrf iteration-1
% z_2d is the original GoogleDEM
z_2d = imread(dem);
% obj is the initally defined object, in this case building footprint from OpenStreetMap
% building layer
build = 'Delhi/Yamuna/Yamuna_rec_bd_final.tif'; % here the raster was generated from the openstreetmap building footprint
bd = imread(build);
obj_ori = logical(mod(bd,-128)); %make it logical

function [dem_f2,obj_f2] = tile_processing(
    x_2d, y_2d, z_2d, obj_ori, ...
    pixel_r, pixel_c, tile_r,tile_c, path, ...
    s, w, et, es, iteration, ...
    info1.GeoTIFFTags.GeoKeyDirectoryTag)

%% smrf iteration-2
s = 0.025;
w = 60;
iteration = 2;
% z_2d is the fitered dem from the above smrf
z_2d = imread(dem_f1);
% obj1 is the detected object dem from the above smrf
obj_ori = imread(obj1)

function [dem_f2,obj_f2] = tile_processing(
    x_2d, y_2d, z_2d, obj_ori, ...
    pixel_r, pixel_c, tile_r,tile_c, path, ...
    s, w, et, es, iteration, ...
    info1.GeoTIFFTags.GeoKeyDirectoryTag)

%% smrf iteration-3
w = w*2;
iteration = 3;

% read the GoogleDEM-F2, Object_f2
% z_2d is the fitered dem from the above smrf
z_2d = imread(dem_f2);
% obj_f2 is the detected object dem from the above smrf
obj_ori = imread(obj_f2)

function [dem_f3,obj_f3] = tile_processing(
    x_2d, y_2d, z_2d, obj_ori, ...
    pixel_r, pixel_c, tile_r,tile_c, path, ...
    s, w, et, es, iteration, ...
    info1.GeoTIFFTags.GeoKeyDirectoryTag)
