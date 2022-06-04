%getAPH(2001, 300, 900, false)

function[AP] = getAPH_func(year,dayofyear,UTseconds,update_files)
 %% APH DATA RETRIEVAL FUNCTION
 % -----------------------------------------------------------------------
 % This functions retrieves APH data from the ftp server of the National
 % Centers for Environmental Information (NCEI) to be implemented in the
 % MATLAB function atmosnrlmsise00. This function should go in conjunction
 % with getf107.m. FTP Server and corresponding directory can be found at
 % ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/
 % -----------------------------------------------------------------------
 % Input: - year : Year [yr]
 % - dayofyear : Day of Year [days]
 % - UTseconds : Seconds as given in Universal Time [s]
 % - update_files : Boolean for File Update [true/false]
 % -----------------------------------------------------------------------
 % Output: - AP : Contains APH info, see the MATLAB function
 % atmosnrlmsise00 for exact information.
 % -----------------------------------------------------------------------
 % Developed by David Ju 09-02-17 Copyright (c)

 %% Check Input Values
 if year < 1932 || (year == 1932 && dayofyear <= 2)
     disp('ERROR: No (complete) data exists before 1932/01/02')
     return
 end
 if dayofyear < 0
     disp('ERROR: Please enter a positive day of year (real integer)')
     return
 end
 if UTseconds < 0
     disp('ERROR: Please enter a positive UT seconds (real)')
     return
 end 
 %% Retrieve Magnetic Index Data
 % Update Files
 if update_files == true
     disp('Downloading takes a few minutes')
     ngdc = ftp('ftp.ngdc.noaa.gov');
     contents = '/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/*';
     mget(ngdc,contents);
     close(ngdc);
     disp('Downloading Complete')
 end
 %% Convert Time to Serial Date Number
 % Requested Times
 date_n = datenum(year,1,dayofyear) + UTseconds/(24*3600);

 %% Check Time 57H Before First Time
 time_57h = date_n - 57/24;
 [year_p,~,~,~,~] = datevec(time_57h);

 if year_p < year
     prev_year = true;
     year = [year_p; year];
 else
     prev_year = false;
 end

 %% Get Data File(s)
 % Convert Number to String
 year_str = cellstr(num2str(year));

 % Allocation Matrix
 data_set = [];

 % Start Data Selection
 for i = 1:length(year)
     file_name = fullfile(pwd,'STP','GEOMAGNETIC_DATA','INDICES','KP_AP',...
     year_str{i});
     % Open File
     fid = fopen(file_name,'r'); % Open File
     data = fscanf(fid,'%c'); % Read File
     fclose(fid); % Close File

     % Reshape to Matrix Format
     data = reshape(data,72,length(data)/72)';

     % Take Last 3 Days of Previous Year
     if (prev_year == true && i == 1)
         data_set = data(end-2:end,:);
     else
         data_set = [data_set; data];
     end
 end

 % Find Current 3H Interval
 UTinterval = floor(UTseconds/(3*3600));
 APinterval = 32:34; % Location in Folder

 % Find Date
 if prev_year == true
     n_d = dayofyear + 3;
 else
     n_d = dayofyear;
 end

 % Retrieve "Useful" Column Data
 AP0 = str2double(data_set(n_d,56:58));

 for i = 1:20
     APi(i) = str2double(data_set(n_d,APinterval+3*UTinterval));
     UTinterval = UTinterval -1;

     if UTinterval < 0
         UTinterval = 7;
         n_d = n_d -1;
     end
 end

 APmean1 = mean(APi(5:12));
 APmean2 = mean(APi(13:end));

 AP = [AP0,APi(1:4),APmean1,APmean2];
end
