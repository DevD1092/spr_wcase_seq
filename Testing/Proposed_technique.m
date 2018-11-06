% This script is to test any proposed star recognition algorithm


clear all;

clc;

close all;

clear classes;


% Initialize camera parameter

% SI STS parameters

FOV = 15;

img_height = 1024;

img_width = 1024;

pixel_size = 13.3/1024;

f = (img_height)*pixel_size /2/ tand(FOV/2);

% g=4;
% 
% pr=ceil(384/g);

img_dimension=[img_width img_height];


%Load SPD

% pattern_database_file='pattern_catalogue_proposed_v4_v384.txt';


%Load SC

catalog_path='..\simulate\SKY2000_Magnitude6_doublestars_0.12.txt';


% Initialize camera attitude

int=2; %2.5: 10368, 3.75: 4608, 5:2592, 7.5: 1152

RAini=1;%from [0, 360];

RAend=359;

% Generate RA in a sequential manner.
RA=[RAini:int:RAend]';

% Generate RA randomly

% a = 50;
% b = 100;
% r = (b-a).*rand(1000,1) + a;

% RA = (RAend - RAini).*rand(1000,1) + RAini;

DECini=-89; %from(-90,90);

DECend=89;

% Generate DEC in a sequntial manner.

DEC=[DECini:int:DECend]; %Column-wise

% Generate DEC randomly.

% DEC = (DECend - DECini).*rand(1000,1) + DECini;


angle=45;

var_test = 0;

%feature_vector tolerance

% tol=1;

% fn_threshold=0.7;


%Simulator parameters

cent_variance=0;

no_ran_star=0;

SNR=10; % SNR = 10 : ideal 

background_noise=0.0; % Background noise = 0.0 : ideal

PSF_set=0;          % Change this setting to (1,2 or 3) if you want to calculate the centroid. PSF_set = 0 means exact co-ordinate of the star given with no magnitude consideration.


%centroider parameters

thresh=0.3;


% Read star coordinates in Earth reference frame from star catalog

[SKYMAP_No,star_RA,star_DEC,star_MAG]= textread(catalog_path,'%d %f %f %f');

% Si is coordiante of star in Earth reference frame,

% the 3 column are X, Y,and Z

Si = [cosd(star_DEC).*cosd(star_RA) cosd(star_DEC).*sind(star_RA)  sind(star_DEC)];

% for i=1: length(star_RA)

%     % Convert RA, and DEC into ECI unit vector

%     ECI_vector =[cosd(star_DEC(i))* cosd(star_RA(i)) cosd(star_DEC(i))* sind(star_RA(i))  sind(star_DEC(i))];

%     Si(i,:)= ECI_vector;

% end


catalog=struct('SKYMAP_No',SKYMAP_No,'star_RA',star_RA,'star_DEC',star_DEC,'star_MAG',star_MAG,'Si',Si);

 

% Testing   variables
acc = 0;
match_acc=0;
count_acc = 0;
time_elapsed = [];
Star_ref_empty_cases = 0;
Star_ref_empty_cases_list = [];
not_successful_RA_empty = [];
not_successful_DEC_empty = [];
multiple_star_id = 0;
multiple_star_id_list = [];
failed_images_number_of_stars = [];
false_match_RA = [];
false_match_DEC = [];
avg_number_star_return = [];
PSID_avg_length = [];
id_shift_list = [];
id_shift_list_gre_1 = [];
% starID=zeros(size(RA,1),size(DEC,2));

% star_q_fov=zeros(size(RA,1),size(DEC,2));

% star_q_ad=zeros(size(RA,1),size(DEC,2));

% star_num_returned_match=zeros(size(RA,1),size(DEC,2));

% starnum_ad_frm_sra=zeros(size(RA,1),size(DEC,2));

% starnum_fov_frm_ad=zeros(size(RA,1),size(DEC,2));

%% THE LOOP FOR GENERATING IMAGES


for ii = 1:size(RA,1)
    
    for jj = 1:size(DEC,2)

        id=[];

        fn=[];

        temp=[];
        
%              RA = 1;
%              DEC = -59;       
       % Generate sky image at predetermined attitude

        Reci2body= Convert_Axis_2_AttitudeMatrix(RA(ii),DEC(jj),angle);

        [star_matrix, I]= Plot_sky_images(Reci2body, FOV, img_height, img_width, pixel_size,cent_variance, no_ran_star, SNR, background_noise, PSF_set,catalog);
    
        [row_star_matrix,col_star_matrix] = size(star_matrix);
        
        for x = 1:row_star_matrix
            id = [id star_matrix(x,1)];
        end
        
%         id
        
        len_id = length(id);
        
      if (len_id > 3)  
 %% Got the image I simulated ; Start the algo from here.

    %% Calculate star centroid
    
    %% ------- Uncomment this for sub - pixel centroiding process ---------
    tic
%     [Sc centroid magnitude] = centroider_mine(I, FOV, img_height, img_width, pixel_size); % Change the PSF_set to (1,2 or 3) if you are using this.

 %% ------- Uncomment this for sub - pixel centroiding process ---------
 
%     centroider_time(RA,DEC+91)=toc;
%     centroider_star(RA,DEC+91)=size(Sc,1);
   
% --- This is for exact co-ordinates of the star with PSF_setting = 0 ----

centroid = [];
Sc = [];


A = (1 : len_id)' ;
for i = 1 : len_id
    centroid(i,1) = star_matrix(i,11);
    centroid(i,2) = star_matrix(i,12);
    Sc(i,1) = star_matrix(i,8);
    Sc(i,2) = star_matrix(i,9);
    Sc(i,3) = star_matrix(i,10);
end

% --- This is for exact co-ordinates of the star with PSF_setting = 0 ----
    
no_star= size(Sc,1);
    
% choose a reference star that is nearest to the center
d=zeros(no_star,1);
for i=1:no_star
    temp= centroid(i,:)-[img_height/2 img_width/2];
    d(i)= sum(temp.*temp);    
end
k= find(d==min(d));
k=k(1);

S1= Sc(k,:);
distance=zeros(no_star,1);

% Co-ordinates of the star nearest to the center of the image.
% Shifting the co-oridnates
Star_ref = [];
PSID_list = [];
centroid_new = [];
centroid_newer = [];
dist_new = [];
center_star_x_cord = (centroid(k,1));
center_star_y_cord = (centroid(k,2));
shift_x_cord = abs((img_height/2) - center_star_x_cord);
shift_y_cord = abs((img_height/2) - center_star_y_cord);
j = 1;
for i = 1 : length(centroid)
    if(i ~= k)
    if(center_star_x_cord < 512)
            centroid_new(i,1) = centroid(i,1) + shift_x_cord;
    end
    if (center_star_x_cord > 512)
            centroid_new(i,1) = centroid(i,1) - shift_x_cord;
    end
    if(center_star_y_cord < 512)
            centroid_new(i,2) = centroid(i,2) + shift_y_cord;
    end
    if(center_star_y_cord > 512)
           centroid_new(i,2) = centroid(i,2) - shift_y_cord;
    end
    if(center_star_x_cord == 512)
           centroid_new(i,1) = centroid(i,1);
    end
     if(center_star_y_cord == 512)
          centroid_new(i,2) = centroid(i,2);
     end
    end
    
    if(i == k)
        centroid_new(i,1) = 512;
        centroid_new(i,2) = 512;
    end
   
    dist_new(i) = sqrt((centroid_new(i,1) - 512) ^ 2 + (centroid_new(i,2) - 512) ^ 2); 
    
    if(centroid_new(i,1) >= 0 && centroid_new(i,1) <= 1024)
        if(centroid_new(i,2) >= 0 && centroid_new(i,2) <= 1024)
            if(dist_new(i) < 510)     % Just to be on the safe side.
            centroid_newer(j,1) = centroid_new(i,1);
            centroid_newer(j,2) = centroid_new(i,2);
            j = j + 1;
        end
        end
    end    
end

if(length(centroid_newer) > 2)
% Listing the angles of the stars in the image the same anticlockwise order
% These angles are of the stars in the image.
angle_hor_list = [];
sorted_angle_hor_list = [];
angle_final_list = [];
star_ids = [];
star_ids_sorted = [];

for i=1: length(centroid_newer)    
  
   if(centroid_newer(i,1) ~= (img_height/2) || centroid_newer(i,2) ~= (img_height/2))
       
        x_cord = centroid_newer(i,1);
        y_cord = centroid_newer(i,2);
       
        % x_cord corresponds to the row number of the star in the image.
        % y_cord corresponds to the column number of the star in the image.
        
        dim_x = abs(x_cord - (img_height/2));
        dim_y = abs(y_cord - (img_height/2));
        angle_hor = atand(dim_x/dim_y);
        
        % See the documentation for quadrant information.
        % Star lies in the first quadrant. 
        if( x_cord > 0 && x_cord < (img_height/2))
            if( y_cord > (img_height/2) && y_cord < (img_height))
                angle_hor_list = [angle_hor_list angle_hor];            
            end
        end
        
        % Star lies in the second quadrant.
        if(x_cord > 0 && x_cord < (img_height/2))
            if(y_cord > 0 && y_cord <  (img_height/2))
                angle_hor = 180 - angle_hor;
                angle_hor_list = [angle_hor_list angle_hor];            
            end
        end
        
        % Star lies in the third quadrant.
        if( x_cord >  (img_height/2) && x_cord < img_height)
            if( y_cord > 0 && y_cord <  (img_height/2))
                angle_hor = 180 + angle_hor;
                angle_hor_list = [angle_hor_list angle_hor];            
            end
        end
        
        % Star lies in the fourth quadrant.
        if( x_cord >  (img_height/2) && x_cord < img_height)
            if( y_cord >  (img_height/2) && y_cord < (img_height))
                angle_hor = 360 - angle_hor;
                angle_hor_list = [angle_hor_list angle_hor];            
            end
        end    
        
        % Covering the extreme conditions: 0,90,180,270.
        
        if(x_cord ==  (img_height/2) && y_cord >  (img_height/2))
            angle_hor = 0;
            angle_hor_list = [angle_hor_list angle_hor];
        end
        
        if(x_cord ==  (img_height/2) && y_cord <  (img_height/2))
            angle_hor = 180;
            angle_hor_list = [angle_hor_list angle_hor];
        end
         
        if(y_cord ==  (img_height/2) && x_cord <  (img_height/2))
            angle_hor = 90;
            angle_hor_list = [angle_hor_list angle_hor];
        end
        
         if(y_cord ==  (img_height/2) && x_cord >  (img_height/2))
             angle_hor = 270;
            angle_hor_list = [angle_hor_list angle_hor];
         end
          star_ids = [star_ids; i angle_hor];
   end
end
   
    sorted_angle_hor_list = sort(angle_hor_list);    
    
    for l = 1 : length(sorted_angle_hor_list)
        
        if(l ~= length(sorted_angle_hor_list))
            angle_between = sorted_angle_hor_list(l+1) - sorted_angle_hor_list(l);
        end
        
        if(l == length(sorted_angle_hor_list))
            angle_between = 360 - sorted_angle_hor_list(l) + sorted_angle_hor_list(1);
        end        
        angle_final_list(l) =  angle_between;   
    end
    
    count_star_id = 1;
    for jjj = 1 : length(sorted_angle_hor_list)
        for kk = 1 : (length(centroid_newer) - 1)
            if(star_ids(kk,2) == sorted_angle_hor_list(jjj))
                star_ids_sorted(count_star_id) = star_ids(kk,1);
                 count_star_id =  count_star_id  + 1;
            end
        end        
    end
    
    
    % Calculating the radial distances in the same order.
    dist_im_shift_list = [];
    for i = 1 : length(star_ids_sorted)
        x_cord_shift_im = abs((img_height/2) - centroid_newer(star_ids_sorted(i),1));
        y_cord_shift_im = abs((img_height/2) - centroid_newer(star_ids_sorted(i),2));
        dist_im_shift = sqrt((x_cord_shift_im).^2 + (y_cord_shift_im).^2);
        dist_im_shift_list = [dist_im_shift_list dist_im_shift];
    end
%% Load the Star pattern databases that were generated

SPD_angle = dlmread('..\SPD\SPD_angles_Mv_6.txt');

SPD_dist = dlmread('..\SPD\SPD_distance.txt');

PSID_list = [];
voting_list = [];
for i = 1 : length(SPD_dist)
        PSID = 0;
        ind_SPD_dist = 0;
        dist_matched_image = 0;
        ind_dist_matched_1_1 = 0;
        ind_dist_matched_1_2 = 0;
        ind_dist_matched_2_1 = 0;
        ind_dist_matched_2_2 = 0;
        ind_dist_matched_1 = 0;
        ind_dist_matched_2 =0;
        ind_dist_matched_image = [];  
        diff_image_list_1 = [];
        diff_image_list_2 = [];
        diff_indices = [];
        diff_indices_1 = [];
        diff_indices_2 = [];
        
            % Getting the matches in the SPD
            for j = 2 : length(SPD_dist(i,:))
                if(SPD_dist(i,j) ~= 0)
                    diff_image_1 = abs(dist_im_shift_list(1) - SPD_dist(i,j));
                    diff_image_list_1  = [diff_image_list_1 diff_image_1];
                    diff_image_2 = abs(dist_im_shift_list(2) - SPD_dist(i,j));
                    diff_image_list_2  = [diff_image_list_2 diff_image_2];
                end
            end
           diff_image_list_1 = sort(diff_image_list_1);
           diff_image_list_2 = sort(diff_image_list_2);
           if(diff_image_list_1(1) < 2.5) % 2.15 - initial
               for j = 2 : length(SPD_dist(i,:))
                   diff_image_1 = abs(dist_im_shift_list(1) - SPD_dist(i,j));
                   if(diff_image_1 == diff_image_list_1(1))
                       ind_dist_matched_1_1 = j;
                   end
               end
                if(diff_image_list_1(2) < 2.5) % 2.15 - initial
                    for j = 2 : length(SPD_dist(i,:))
                        diff_image_1 = abs(dist_im_shift_list(1) - SPD_dist(i,j));
                        if(diff_image_1 == diff_image_list_1(2))
                            ind_dist_matched_1_2 = j;
                        end
                    end
                end
                PSID = i;
            end
            if(diff_image_list_2(1) < 2.5) % 2.15 - initial
                for j = 2 : length(SPD_dist(i,:))
                    diff_image_2 = abs(dist_im_shift_list(2) - SPD_dist(i,j));
                    if(diff_image_2 == diff_image_list_2(1))
                        ind_dist_matched_2_1 = j;
                    end
                end
                if(diff_image_list_2(2) < 2.4) % 2.0 - initial
                    for j = 2 : length(SPD_dist(i,:))
                        diff_image_2 = abs(dist_im_shift_list(2) - SPD_dist(i,j));
                        if(diff_image_2 == diff_image_list_2(2))
                            ind_dist_matched_2_2 = j;
                        end
                    end
                end
                PSID = i;
            end
            
             if(PSID == 0)
                continue;
             end
            
            if(ind_dist_matched_1_1 == 0 && ind_dist_matched_1_2 == 0)
                continue;
            end
            
            if(ind_dist_matched_2_1 == 0 && ind_dist_matched_2_2 == 0)
                continue;
            end
            
            if(ind_dist_matched_2_2 ~= 0 && ind_dist_matched_1_1 ~= 0 && ind_dist_matched_1_2 ~= 0 && ind_dist_matched_2_1 ~= 0)
                diff_indices = [abs(ind_dist_matched_1_1 - ind_dist_matched_2_1) abs(ind_dist_matched_1_1 -  ind_dist_matched_2_2) abs(ind_dist_matched_1_2 - ind_dist_matched_2_1) abs(ind_dist_matched_1_2 - ind_dist_matched_2_2)];
                temp_ind = [];
                diff_indices_up = [];
                for zz = 1 : length(diff_indices)
                    if(diff_indices(zz) ~= 0)
                        diff_indices_up = [diff_indices_up diff_indices(zz)];
                        temp_ind = [temp_ind zz];
                    end
                end
                
                min_diff_indices = find(diff_indices_up == min(diff_indices_up));
                
                if(length(min_diff_indices) == 0)
                    continue;
                end
                
                if(length(min_diff_indices) > 1)
                    if(min(diff_indices_up) ~= 1)
                        for zz = 1 : length(temp_ind)
                            if(temp_ind(zz) == 1)
                                ind_dist_matched_1 = ind_dist_matched_1_1;
                                ind_dist_matched_2 = ind_dist_matched_2_1;
                                break;
                            end
                            if(temp_ind(zz) == 2)
                                ind_dist_matched_1 = ind_dist_matched_1_1;
                                ind_dist_matched_2 = ind_dist_matched_2_2;
                                break;
                            end
                            if(temp_ind(zz) == 3)
                                ind_dist_matched_1 = ind_dist_matched_1_2;
                                ind_dist_matched_2 = ind_dist_matched_2_1;
                                break;
                            end
                            if(temp_ind(zz) == 4)
                                ind_dist_matched_1 = ind_dist_matched_1_2;
                                ind_dist_matched_2 = ind_dist_matched_2_2;
                                break;
                            end
                        end
                    end
                end
                
                if(length(min_diff_indices) > 1)
                    if(min(diff_indices_up) == 1)
                       if((SPD_angle(PSID,ind_dist_matched_1_1) - 0.7) < angle_final_list(1) && angle_final_list(1) < (SPD_angle(PSID,ind_dist_matched_1_1) + 0.7))
                            ind_dist_matched_1 = ind_dist_matched_1_1;
                       end
                       if((SPD_angle(PSID,ind_dist_matched_1_2) - 0.7) < angle_final_list(1) && angle_final_list(1) < (SPD_angle(PSID,ind_dist_matched_1_2) + 0.7))
                            ind_dist_matched_1 = ind_dist_matched_1_2;
                       end
                       if(diff_indices(1) == 1 && ind_dist_matched_1 == ind_dist_matched_1_1)
                            ind_dist_matched_2 = ind_dist_matched_2_1;
                       end
                       if(diff_indices(2) == 1 && ind_dist_matched_1 == ind_dist_matched_1_1)
                            ind_dist_matched_2 = ind_dist_matched_2_2;
                       end
                       if(diff_indices(3) == 1 && ind_dist_matched_1 == ind_dist_matched_1_2)
                            ind_dist_matched_2 = ind_dist_matched_2_1;
                       end
                       if(diff_indices(4) == 1 && ind_dist_matched_1 == ind_dist_matched_1_2)
                            ind_dist_matched_2 = ind_dist_matched_2_2;
                       end
                    end
                end
                
                if(length(min_diff_indices) == 1)
                    for zz = 1 : length(diff_indices)
                        if(diff_indices_up(min_diff_indices) == diff_indices(zz))
                            ind_matched = zz;
                            break;
                        end
                    end
                    if(ind_matched == 1)
                        ind_dist_matched_1 = ind_dist_matched_1_1;
                        ind_dist_matched_2 = ind_dist_matched_2_1;
                    end
                    if(ind_matched == 2)
                        ind_dist_matched_1 = ind_dist_matched_1_1;
                        ind_dist_matched_2 = ind_dist_matched_2_2;
                    end
                    if(ind_matched == 3)
                        ind_dist_matched_1 = ind_dist_matched_1_2;
                        ind_dist_matched_2 = ind_dist_matched_2_1;
                    end
                    if(ind_matched == 4)
                        ind_dist_matched_1 = ind_dist_matched_1_2;
                        ind_dist_matched_2 = ind_dist_matched_2_2;
                    end
                end
            end
            
            if(ind_dist_matched_1_2 == 0 && ind_dist_matched_2_2 == 0)
                ind_dist_matched_1 = ind_dist_matched_1_1;
                ind_dist_matched_2 = ind_dist_matched_2_1;
            end
            
            if(ind_dist_matched_2_2 == 0 && ind_dist_matched_1_2 ~= 0)
                ind_dist_matched_2 = ind_dist_matched_2_1;
                count_zeros = 0;
                for j = 1 : length(SPD_angle(PSID,:))
                    if(SPD_angle(PSID,j) == 0)
                        count_zeros = count_zeros + 1;
                    end
                end
                if(ind_dist_matched_2 == 2)
                    if((68 - count_zeros) == ind_dist_matched_1_1 || (68 - count_zeros) == ind_dist_matched_1_2)
                        ind_dist_matched_1 = 68 - count_zeros;
                    end
                end
                if(ind_dist_matched_2 ~= 2)
                    if(ind_dist_matched_1_1 < ind_dist_matched_2 && ind_dist_matched_1_2 < ind_dist_matched_2)
                        diff_indices_1 = [abs(ind_dist_matched_1_1 - ind_dist_matched_2) abs(ind_dist_matched_1_2 - ind_dist_matched_2)];
                        min_diff_indices_1 = find(diff_indices_1 == min(diff_indices_1));
                        if(min_diff_indices_1 == 1)
                            ind_dist_matched_1 = ind_dist_matched_1_1;
                        end
                        if(min_diff_indices_1 == 2)
                            ind_dist_matched_1 = ind_dist_matched_1_2;
                        end
                    else
                        diff_indices_1 = [(ind_dist_matched_1_1 - ind_dist_matched_2) (ind_dist_matched_1_2 - ind_dist_matched_2)];
                        if(diff_indices_1(1) == 0)
                            ind_dist_matched_1 = ind_dist_matched_1_2;
                            diff_indices_1 = [(ind_dist_matched_1_2 - ind_dist_matched_2)];
                        elseif(diff_indices_1(2) == 0)
                            ind_dist_matched_1 = ind_dist_matched_1_1;
                            diff_indices_1 = [(ind_dist_matched_1_1 - ind_dist_matched_2)];
                        end
                        if(length(diff_indices_1) == 2)
                            if(diff_indices_1(1) == diff_indices_1(2))
                                continue;
                            end
                            min_diff_indices_1 = find(diff_indices_1 == min(diff_indices_1));
                            if(min_diff_indices_1 == 1)
                                ind_dist_matched_1 = ind_dist_matched_1_1;
                            end
                            if(min_diff_indices_1 == 2)
                                ind_dist_matched_1 = ind_dist_matched_1_2;
                            end
                        end
                    end
                end
            end
            
           if(ind_dist_matched_1_2 == 0 && ind_dist_matched_2_2 ~= 0)
                ind_dist_matched_1 = ind_dist_matched_1_1;
                count_zeros = 0;                
                for j = 1 : length(SPD_angle(PSID,:))
                    if(SPD_angle(PSID,j) == 0)
                        count_zeros = count_zeros + 1;
                    end
                end
                if((68 - count_zeros) == ind_dist_matched_1)
                     ind_dist_matched_2 = min(ind_dist_matched_2_1,ind_dist_matched_2_2);
                end                
                if((68 - count_zeros) ~= ind_dist_matched_1)
                    diff_indices_2 = [(ind_dist_matched_2_1 - ind_dist_matched_1) (ind_dist_matched_2_2 - ind_dist_matched_1)];
                    if(diff_indices_2(1) == 0)
                        ind_dist_matched_2 = ind_dist_matched_2_2;
                        diff_indices_2 = [(ind_dist_matched_2_2 - ind_dist_matched_1)];
                    elseif(diff_indices_2(2) == 0)
                        ind_dist_matched_2 = ind_dist_matched_2_1;
                        diff_indices_2 = [(ind_dist_matched_2_2 - ind_dist_matched_1)];
                    end
                    if(length(diff_indices_2) == 2)
                        if(diff_indices_2(1) == diff_indices_2(2))
                            continue;
                        end
                        diff_indices_2 = [abs(ind_dist_matched_2_1 - ind_dist_matched_1) abs(ind_dist_matched_2_2 - ind_dist_matched_1)];
                            if(ind_dist_matched_2_1 > ind_dist_matched_1)
                                ind_dist_matched_2 = ind_dist_matched_2_1;
                            end
                            if(ind_dist_matched_2_2 > ind_dist_matched_1)
                                ind_dist_matched_2 = ind_dist_matched_2_2;
                            end
                            min_diff_indices_2 = find(diff_indices_2 == min(diff_indices_2));
%                         if(length(min_diff_indices_2) == 1)
%                             if(min_diff_indices_2 == 1)
%                                 ind_dist_matched_2 = ind_dist_matched_2_1;
%                             end
%                             if(min_diff_indices_2 == 2)
%                                 ind_dist_matched_2 = ind_dist_matched_2_2;
%                             end
%                         end
                        if(length(min_diff_indices_2) == 2)
                            if(ind_dist_matched_2_1 > ind_dist_matched_1)
                                ind_dist_matched_2 = ind_dist_matched_2_1;
                            end
                            if(ind_dist_matched_2_2 > ind_dist_matched_1)
                                ind_dist_matched_2 = ind_dist_matched_2_2;
                            end
                        end
                    end
                end
            end
            
            if(ind_dist_matched_1 == 0 || ind_dist_matched_2 == 0)
                continue;
            end
            
            if(abs(ind_dist_matched_1 - ind_dist_matched_2 ) == 1)
                if(ind_dist_matched_1 > ind_dist_matched_2)
                    temp = ind_dist_matched_1;
                    ind_dist_matched_1 = ind_dist_matched_2;
                    ind_dist_matched_2 = temp;                    
                end
            end
            
            % Re-arranging the SPD according to the matches got.
            % No need to re-arrange the angle_final_list.
            
            % Calculating the number of zeros in the SPD of that particular Star_ID
                    count_zeros = 0;
                    
                    for j = 1 : length(SPD_angle(PSID,:))
                        if(SPD_angle(PSID,j) == 0)
                            count_zeros = count_zeros + 1;
                        end
                    end
                    
                    SPD_new = [];
                    ind_incr_SPD = 0;
                    ind_incr_SPD_1 = 2;
                    limit = length(SPD_angle(PSID,:)) - count_zeros - 1;
                    for j = 1 : limit
                        if(ind_dist_matched_1 + ind_incr_SPD > 68)
                            SPD_new(j) = SPD_angle(PSID,(ind_incr_SPD_1));
                            ind_incr_SPD_1 = ind_incr_SPD_1 + 1;
                            continue;
                        end
                        if(SPD_angle(PSID,ind_dist_matched_1 + ind_incr_SPD) ~= 0)
                            SPD_new(j) = SPD_angle(PSID,(ind_dist_matched_1 + ind_incr_SPD));
                        end
                        if(SPD_angle(PSID,ind_dist_matched_1 + ind_incr_SPD) == 0)
                            SPD_new(j) = SPD_angle(PSID,(ind_incr_SPD_1));
                            ind_incr_SPD_1 = ind_incr_SPD_1 + 1;
                        end
                        ind_incr_SPD = ind_incr_SPD + 1;
                    end
                    
                    % Start the comparison of the SPD and the image angles
                    % from here. Voting takes place.
                    
                    ind_to_comp = 1;
                    add_SPD_angles = 0;
                    no_matches = 0;
                    
                    for m = 1 : length(angle_final_list)
                        if(ind_to_comp > length(SPD_new))
                            break;
                        end
                        if((SPD_new(ind_to_comp) - 0.85) < angle_final_list(m) && angle_final_list(m) < (SPD_new(ind_to_comp) + 0.85)) % 0.85 - best
                            ind_to_comp = ind_to_comp + 1;
                            no_matches = no_matches + 1;
                            continue;
                        else
                            add_SPD_angles = SPD_new(ind_to_comp);
                            for n = (ind_to_comp + 1) : length(SPD_new)
                                add_SPD_angles = add_SPD_angles + SPD_new(n);
                                if((add_SPD_angles - 1.1) < angle_final_list(m) && angle_final_list(m) < (add_SPD_angles + 1.1)) % 1.1 - best
                                    ind_to_comp = n + 1;
                                    no_matches = no_matches + 1;
                                    break;
                                end
                                if(add_SPD_angles >  (angle_final_list(m) + 3))
                                    ind_to_comp = ind_to_comp + 1;
                                    break;
                                end
                            end
                        end
                    end
                    
                    %------------------- Highest vote match ---------------------%
                    
                    voting_list = [voting_list no_matches];
                    PSID_list = [PSID_list PSID];
                    
                    %------------------- Highest vote match ---------------------%
end
end

if(length(PSID_list) ~= length(voting_list))
    voting_list = [];
    PSID_list = [];
end

PSID_avg_length = [PSID_avg_length length(PSID_list)];

max_vote = max(voting_list);
ind = find(voting_list == max_vote);
Star_ref = PSID_list(ind);
 

      
% Star_ref
toc
time_elapsed = [time_elapsed toc];


if(length(Star_ref) <= 10)
    for p = 1 : length(Star_ref)
        if(ismember(Star_ref(p),id) == 1)
            match_acc = match_acc + 1;
            break;
        end
    end
end

% ll = intersect(Star_ref,id);
% if(length(Star_ref) >= 1)
%     if(length(ll) == 0)
%         false_match_RA = [false_match_RA RA(ii)];
%         false_match_DEC = [false_match_DEC DEC(jj)];
%     end
% end
% 
% if(length(Star_ref) > 1)
%     multiple_star_id = multiple_star_id + 1;
%     multiple_star_id_list = [multiple_star_id_list length(Star_ref)];
%     id_shift_list_gre_1 = [id_shift_list_gre_1 length(centroid_newer)];
% end
% 
% if(length(Star_ref) == 1)
%     id_shift_list = [id_shift_list length(centroid_newer)];
% end
% 
% 
% if(length(Star_ref) == 0)
%     not_successful_RA_empty = [not_successful_RA_empty RA(ii)];
%     not_successful_DEC_empty = [not_successful_DEC_empty DEC(jj)];
% end


count_acc = count_acc + 1;
end

if(count_acc == 100)
    return;
end
      
% imshow(I);

pause(0.2)

        
    end
end
% % save('test.mat');
