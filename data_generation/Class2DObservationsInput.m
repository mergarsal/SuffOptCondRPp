classdef Class2DObservationsInput
   properties
        pt_number = 100
        noise = 0.1                   % in pixels 
        number_outliers = 0
        FoV_par = 120                 % in degrees
        max_parallax = 2.0            % in meters
        min_depth = 1.0               % in meters
        max_depth = 8.0               % in meters 
        focal_length = 500.0          % in pixels, if 0 => compute from image size 
        image_size = [1280, 720]      % [rows, cols] in pixels  
   end
end
