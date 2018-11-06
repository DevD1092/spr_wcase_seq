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

f = (img_height)*pixel_size /2/ tand(FOV/2);        % Change this according to the FOV taken for SPD generation

% g=4;
% 
% pr=ceil(384/g);

img_dimension=[img_width img_height];


%Load SPD

% pattern_database_file='pattern_catalogue_proposed_v4_v384.txt';


%Load SC

catalog_path='..\simulate\SKY2000_Magnitude6_doublestars_0.12.txt';

% Initialize camera attitude

int=0.3; %2.5: 10368, 3.75: 4608, 5:2592, 7.5: 1152

RAini=1;%from [0, 360];

RAend=359;

RA=[RAini:int:RAend]'; %row-wise

DECini=-89; %from(-90,90);

DECend=89;

DEC=[DECini:int:DECend]; %Column-wise

angle=45;

var_test = 0;

%feature_vector tolerance

% tol=1;

% fn_threshold=0.7;


%Simulator parameters

cent_variance=0;

no_ran_star=0;

% SNR=7;

% background_noise=0.1;

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
failed_images_number_of_stars = [];
false_match_RA = [];
false_match_DEC = [];
avg_number_star_return = [];
failed_cases_RA = [];
failed_cases_DEC = [];
x = 1;

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
        
%          RA = 1;
%          DEC = -69.8;     % -69.8 (empty return) ; -72.8 (false match)  
       % Generate sky image at predetermined attitude

        Reci2body= Convert_Axis_2_AttitudeMatrix(RA(ii),DEC(jj),angle);

        [star_matrix, I, x_cord, y_cord]= Plot_sky_images_mine_random_stars(Reci2body, FOV, img_height, img_width, pixel_size,cent_variance, no_ran_star, PSF_set,catalog);
    
        [row_star_matrix,col_star_matrix] = size(star_matrix);
        
        for x = 1:row_star_matrix
            id = [id star_matrix(x,1)];
        end
        
        id
        
        len_id = length(id);
        
      if (len_id > 3)  
 %% Got the image I simulated ; Start the algo from here.

    %% Calculate star centroid
    
    %% ------- Uncomment this for sub - pixel centroiding process ---------
    tic
%     [Sc centroid magnitude] = centroider_mine(I, FOV, img_height, img_width, pixel_size);

 %% ------- Uncomment this for sub - pixel centroiding process ---------
 
%     centroider_time(RA,DEC+91)=toc;
%     centroider_star(RA,DEC+91)=size(Sc,1);
   
% --- This is for exact co-ordinates of the star with PSF_setting = 0 ----

centroid = [];
Sc = [];


A = (1 : len_id)' ;
for i = 1 : length(x_cord)
    centroid(i,1) = x_cord(i);
    centroid(i,2) = y_cord(i);
%     Sc(i,1) = star_matrix(i,8);
%     Sc(i,2) = star_matrix(i,9);
%     Sc(i,3) = star_matrix(i,10);
end
      
 d=zeros(length(centroid),1);
for i=1:length(centroid)
    temp= centroid(i,:)-[img_height/2 img_width/2];
    d(i)= sum(temp.*temp);    
end
k= find(d==min(d));
k=k(1);

if(k <= (size(star_matrix,1)))
    
% Co-ordinates of the star nearest to the center of the image.
            % Shifting the co-oridnates
            
            
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
            
            if(length(centroid_newer) <= 2)
                  false_match_RA = [false_match_RA RA(ii)];
                  false_match_DEC = [false_match_DEC DEC(jj)];
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
                            angle_hor_list = [angle_hor_list  angle_hor];
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
                
                % Calculate the radial distances in the same order as that of the angles generated.
                
                dist_im_shift_list = [];
                for i = 1 : length(star_ids_sorted)
                    x_cord_shift_im = abs((img_height/2) - centroid_newer(star_ids_sorted(i),1));
                    y_cord_shift_im = abs((img_height/2) - centroid_newer(star_ids_sorted(i),2));
                    dist_im_shift = sqrt((x_cord_shift_im).^2 + (y_cord_shift_im).^2);
                    dist_im_shift_list = [dist_im_shift_list dist_im_shift];
                end
                
                
                %% Load the Star pattern database that was generated
                
                SPD_angle = dlmread('..\SPD\SPD_angles_Mv_6.txt');
                SPD_dist = dlmread('..\SPD\SPD_dist.txt');
                
                %% Starting comparing angles in the image with that of the SPD.
                
                voting_list = []; % To store the votes of PSIDs.
                PSID_list = [];
                %Starting the comparison with the SPD.
                
                for i = 1 : length(SPD_angle)
                    
                    %Step 1:
                    % The below loop is to find the match that should surely exist in the
                    % SPD and then starting the comparison sequence from this index in SPD
                    PSID = 0;
                    ind_SPD = 0;
                    angle_matched_image = 0;
                    ind_angle_matched_image = 0;
                    
                    for q = 1 : length(angle_final_list)
                        diff_angle_list = [];
                        diff_radial_dist  = [];
                        for j = 2 : length(SPD_angle(i,:))
                            if(SPD_angle(i,j) ~= 0)
                                diff_angle = abs(angle_final_list(q) - SPD_angle(i,j));
                                diff_angle_list = [diff_angle_list diff_angle];
                            end
                        end
                        diff_angle_list = sort(diff_angle_list);
                        if(diff_angle_list(1) < 0.1)    
                            for j = 2 : length(SPD_angle(i,:))
                                diff_angle = abs(angle_final_list(q) - SPD_angle(i,j));
                                if(diff_angle == diff_angle_list(1))
                                    pros_ind_SPD_1 = j;
                                end
                                if(diff_angle == diff_angle_list(2))
                                    pros_ind_SPD_2 = j;
                                end
                                if(diff_angle == diff_angle_list(3))
                                    pros_ind_SPD_3 = j;
                                end
                            end
                            radial_dist_im =  dist_im_shift_list(q);
                            radial_dist_1_SPD = SPD_dist(i,pros_ind_SPD_1);
                            radial_dist_2_SPD = SPD_dist(i, pros_ind_SPD_2);
                            radial_dist_3_SPD = SPD_dist(i, pros_ind_SPD_3);
                            diff_radial_dist = [abs(radial_dist_im - radial_dist_1_SPD) abs(radial_dist_im - radial_dist_2_SPD) abs(radial_dist_im - radial_dist_3_SPD)];
                            min_diff_radial_dist = min(diff_radial_dist);
                            ind_matched = find(diff_radial_dist == min_diff_radial_dist);
                            if(ind_matched == 1)
                                angle_matched_image = angle_final_list(q);
                                ind_angle_matched_image = q;
                                ind_SPD = pros_ind_SPD_1;
                                PSID = i;
                                break;
                            end
                            if(ind_matched == 2)
                                angle_matched_image = angle_final_list(q);
                                ind_angle_matched_image = q;
                                ind_SPD = pros_ind_SPD_2;
                                PSID = i;
                                break;
                            end
                            if(ind_matched == 3)
                                angle_matched_image = angle_final_list(q);
                                ind_angle_matched_image = q;
                                ind_SPD = pros_ind_SPD_3;
                                PSID = i;
                                break;
                            end
                        end
                        if(PSID ~= 0)
                            break;
                        end
                    end
                    
                    if(PSID == 0)
                        continue;      % If no angle is matched then continue testing other stars.
                    end
                    
                    % Rearranging the angle list with the sequence starting from the angle
                    % in the image that got matched and now to be compared with the SPD sequence.
                    angle_final_list_new = [];
                    ind_incr = 0;
                    ind_incr_1 = 1;
                    for w = 1 : length(angle_final_list)
                        if(ind_angle_matched_image + ind_incr <= length(angle_final_list))
                            angle_final_list_new(w) = angle_final_list(ind_angle_matched_image + ind_incr);
                        end
                        if(ind_angle_matched_image + ind_incr > length(angle_final_list))
                            angle_final_list_new(w) = angle_final_list(ind_incr_1);
                            ind_incr_1 = ind_incr_1 + 1;
                        end
                        ind_incr = ind_incr + 1;
                    end
                    
                    % Calculating the number of zeros in the SPD of that particular Star_ID
                    count_zeros = 0;
                    
                    for j = 1 : length(SPD_angle(PSID,:))
                        if(SPD_angle(PSID,j) == 0)
                            count_zeros = count_zeros + 1;
                        end
                    end
                    
                    %Step 2:
                    % Rearranging the SPD sequence of that particular Star_ID
                    
                    SPD_new = [];
                    ind_incr_SPD = 0;
                    ind_incr_SPD_1 = 2;
                    limit = length(SPD_angle(PSID,:)) - count_zeros - 1;
                    for j = 1 : limit
                        if(ind_SPD + ind_incr_SPD > 68)
                            SPD_new(j) = SPD_angle(PSID,(ind_incr_SPD_1));
                            ind_incr_SPD_1 = ind_incr_SPD_1 + 1;
                            continue;
                        end
                        if(SPD_angle(PSID,ind_SPD + ind_incr_SPD) ~= 0)
                            SPD_new(j) = SPD_angle(PSID,(ind_SPD + ind_incr_SPD));
                        end
                        if(SPD_angle(PSID,ind_SPD + ind_incr_SPD) == 0)
                            SPD_new(j) = SPD_angle(PSID,(ind_incr_SPD_1));
                            ind_incr_SPD_1 = ind_incr_SPD_1 + 1;
                        end
                        ind_incr_SPD = ind_incr_SPD + 1;
                    end
                    
                    %Step 3:
                    % Start the sequence comparison from here
                    ind_to_comp = 2;
                    add_SPD_angles = 0;
                    no_matches = 0;
                    
                    for m = 2 : length(angle_final_list_new)
                        if(ind_to_comp > length(SPD_new))
                            break;
                        end
                        if((SPD_new(ind_to_comp) - 0.25) < angle_final_list_new(m) && angle_final_list_new(m) < (SPD_new(ind_to_comp) + 0.25))
                            ind_to_comp = ind_to_comp + 1;
                            no_matches = no_matches + 1;
                            continue;
                        else
                            add_SPD_angles = SPD_new(ind_to_comp);
                            for n = (ind_to_comp + 1) : length(SPD_new)
                                add_SPD_angles = add_SPD_angles + SPD_new(n);
                                if((add_SPD_angles - 0.3) < angle_final_list_new(m) && angle_final_list_new(m) < (add_SPD_angles + 0.3)) 
                                    ind_to_comp = n + 1;
                                    no_matches = no_matches + 1;
                                    break;
                                end
                                if(add_SPD_angles >  (angle_final_list_new(m) + 3))
                                    ind_to_comp = ind_to_comp + 1;
                                    break;
                                end
                            end
                        end
                    end
                    
                    %------------------- Perfect match ---------------------%
                    
                    % if (no_matches ==  (length(angle_final_list_new) - 1))
                    %     Star_ref = [Star_ref PSID];
                    
                    %------------------- Perfect match ---------------------%
                    
                    %------------------- Highest vote match ---------------------%
                    
                    voting_list = [voting_list no_matches];
                    PSID_list = [PSID_list PSID];
                    
                    %------------------- Highest vote match ---------------------%
                    
                end
            end
        
      
max_vote = max(voting_list);
ind = find(voting_list == max_vote);
Star_ref = PSID_list(ind);
      
      
Star_ref
toc
time_elapsed = [time_elapsed toc];



for p = 1 : length(Star_ref)
    if(ismember(Star_ref(p),id) == 1)
        match_acc = match_acc + 1;
        break;
    end
end

ll = intersect(Star_ref,id);
if(length(Star_ref) >= 1)
    if(length(ll) == 0)
        false_match_RA = [false_match_RA RA(ii)];
        false_match_DEC = [false_match_DEC DEC(jj)];
    end
end


if(length(Star_ref) > 1)
    multiple_star_id = multiple_star_id + 1;
end

if(length(Star_ref) == 0)
    not_successful_RA_empty = [not_successful_RA_empty RA(ii)];
    not_successful_DEC_empty = [not_successful_DEC_empty DEC(jj)];
end


count_acc = count_acc + 1;
end

if(count_acc == 250)
    return;
end
      
imshow(I);

pause(0.2)

        
    end
     end
end
% % save('test.mat');
    
    