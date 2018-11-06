% This script is to generate the SPD for the proposed idea to prevent


clear all;
clc;
%% Initialize camera parameter
img_height = 1024;
img_width = 1024;
pixelsize= 13.3/1024;
FOV =15;
angle = 45;
test = 0;
% Max_N = 15; 

file_path='..\simulate\SKY2000_Magnitude6_3_doublestars_0.1.txt';
[SKYMAP_No,star_RA,star_DEC,star_MAG]= textread(file_path,'%d %f %f %f');

% Si is coordiante of star in Earth reference frame, 
% the 3 column are X, Y,and Z
no_stars=length(star_RA);
Si = zeros(no_stars, 3);
catalog = zeros(no_stars, 4);
angle_final_list = [];

for i=1: no_stars
    % Convert RA, and DEC into ECI unit vector
%     i = 31;
    ECI_vector =[cosd(star_DEC(i))* cosd(star_RA(i)) cosd(star_DEC(i))* sind(star_RA(i))  sind(star_DEC(i))];
    Si(i,:)= ECI_vector;
end

for i=1: no_stars

   angle_hor_list = [];
   sorted_angle_hor_list = [];
   RA= star_RA(i);
   DEC= star_DEC(i);
   
   C= Convert_Axis_2_AttitudeMatrix(RA,DEC,angle);
   [R_camera_to_earth,star_matrix]= Find_neighbor_star_half_FOV(C,FOV, img_height, img_width, pixelsize);
   
  
    for j = 1 : size(star_matrix,1)
        
        
        if(star_matrix(j,1) ~= i)
        x_cord = star_matrix(j,11);
        y_cord = star_matrix(j,12);
       
        % x_cord corresponds to the row number of the star in the image.
        % y_cord corresponds to the column number of the star in the image.
        
        dim_x = abs(x_cord - (img_height/2));   
        dim_y = abs(y_cord - (img_width/2));
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
            if(y_cord > 0 && y_cord < (img_height/2))
                angle_hor = 180 - angle_hor;
                angle_hor_list = [angle_hor_list angle_hor];            
            end
        end
        
        % Star lies in the third quadrant.
        if( x_cord > (img_height/2) && x_cord < img_height)
            if( y_cord > 0 && y_cord < (img_height/2))
                angle_hor = 180 + angle_hor;
                angle_hor_list = [angle_hor_list angle_hor];            
            end
        end
        
        % Star lies in the fourth quadrant.
        if( x_cord > (img_height/2) && x_cord < img_height)
            if( y_cord > (img_height/2) && y_cord < (img_height))
                angle_hor = 360 - angle_hor;
                angle_hor_list = [angle_hor_list angle_hor];            
            end
        end    
        
        % Covering the extreme conditions: 0,90,180,270.
        
        if(x_cord == (img_height/2) && y_cord > (img_height/2))
            angle_hor_list = [angle_hor_list 0];
        end
        
        if(x_cord == (img_height/2) && y_cord < (img_height/2))
            angle_hor_list = [angle_hor_list 180];
        end
         
        if(y_cord == (img_height/2) && x_cord < (img_height/2))
            angle_hor_list = [angle_hor_list 90];
        end
        
         if(y_cord == (img_height/2) && x_cord > (img_height/2))
            angle_hor_list = [angle_hor_list 270];
         end  
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
        
        angle_final_list(i,l) =  angle_between;    
    
    end
end