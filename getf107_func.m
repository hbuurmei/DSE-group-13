%[a, b] = getf107(2020, 45, false)

function[f107a,f107d] = getf107_func(year,dayofyear,update_files)
 %% 10.7cm SOLAR RADIO FLUX DATA RETRIEVAL FUNCTION
 % -----------------------------------------------------------------------
 % This functions retrieves 10.7 solar radio data from the ftp server of
 % the National Centers for Environmental Information (NCEI) to be
 % implemented in the MATLAB function atmosnrlmsise00. This function
 % should go in conjunction with getAPH.m. FTP Server and corresponding
 % directory can be found at
 % ftp://ftp.ngdc.noaa.gov/STP/space-weather/solar-data/solar-features/
 % solar-radio/noontime-flux/penticton/penticton_observed/listings/
 % Non existent data values will be discarded for averaging and a previous
 % existing value will be taken for f107d.
 % -----------------------------------------------------------------------
 % Input: - year : Year [yr]
 % - dayofyear : Day of Year [days]
 % - update_files : Boolean for File Update [true/false]
 % -----------------------------------------------------------------------
 % Output: - f107a : Averaged Flux Data over 81 days
 % - f107d : Previous Day Fulx Data
 % ----------------------------------------------------------------------- % Developed by David Ju 09-02-17 Copyright (c)

 %% Check Input Values
 if year < 1947 || (year == 1947 && dayofyear <= 44)
     disp('ERROR: No (complete) data exists before 1947/02/14')
     return
 end
 if dayofyear < 0
     disp('ERROR: Please enter a positive day of year (real integer)')
     return
 end

 %% Retrieve Magnetic Index Data
 % Update Files
 if update_files == true
     disp('Download initiated...')
     ngdc = ftp('ftp.ngdc.noaa.gov');
     contents = ['/STP/space-weather/solar-data/solar-features/solar-radio/' ...
              'noontime-flux/penticton/penticton_observed/listings/' ...
              'listing_drao_noontime-flux-observed_daily.txt'];
     mget(ngdc,contents);
     close(ngdc);
     disp('Downloading Complete')
 end

 %% Retrieve Data File
 % Start Data Selection
 file_name = fullfile(pwd,'STP','space-weather','solar-data', ...
 'solar-features', 'solar-radio', 'noontime-flux' ,'penticton',...
 'penticton_observed','listings',...
 'listing_drao_noontime-flux-observed_daily.txt');

 % Open File
 fid = fopen(file_name,'r'); % Open File
 data = textscan(fid,'%4d %2d %2d %f','treatasempty','.'); % Read File
 fclose(fid); % Close File

 %% Find Data Location
 % Time Conversion
 date_n = datenum(year,1,dayofyear); % Convert to Serial Date Number
 [yy,mm,dd,~,~,~] = datevec(date_n); % Convert to Year, Month and Day.

 % Find Date
 idty = find(yy == data{1});
 idtm = find(mm == data{2});
 idtd = find(dd == data{3});

 idt1 = intersect(idty,idtm);
 idt = intersect(idt1,idtd);

 if idt == 0
     disp('ERROR: No matching dates were found.')
     return
 end

 %% Parse Values
 % Averaged 81 Days Value
 f107a = mean(data{4}((idt-40):(idt+40)));

 if isnan(f107a)
     disp('ERROR: No data available for 81 consecutive days.')
 end

 % f107 Previous Day Value
 f107d = data{4}(idt-1);
 % Take Closest Point if NaN
 if isnan(f107d)
     disp('Data incomplete, the closest previous value will be taken.')
 end

 while isnan(f107d)
     idt = idt -1;
     f107a = (idt-1);
 end
end
