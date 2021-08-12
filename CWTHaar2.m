function [exitcode] = CWTHaar2(instrument,location)
% exitcode = instrument;
% cloud1and3 ID in cl51 may need && max_cloud_index(cbj(1,1))<400 line?
% Instrument
% added gradient record 

rwpfiles = dir('*.bln');
files = dir('*.nc');
nof = length(files);

for K = 1:nof
contlim = 200;    
stdvlim = 200; 
if (strcmp(instrument,'15k'))
    instrument = '15k'; 
    bsc_variable_name = 'beta_raw';
    time_variable_name = 'time';
    alt_variable_name =  'range';
    precip_threshold1 = 230000;
    precip_threshold2 = 400000;
    cloud_threshold = 5000;
    min_cloud_threshold_height = 1;
    cloud2fraction = 2; 
    cloud3fraction = 2;
    vres = 15;
    amax= 100; %max a dilations
    ind = 14; %min vertical height
    sgol = 15; 
    sgolclouds = 5; 
    cla=200;
    
elseif (strcmp(instrument,'CL51'))
    instrument = 'CL51'; 
    bsc_variable_name = 'beta_att';
    time_variable_name = 'time';
    precip_threshold1 = 2000;
    precip_threshold2 = 8000;
    cloud_threshold = 10;
    min_cloud_threshold_height = 10;
    cloud2fraction = 5; 
    cloud3fraction = 1.3;
    vres = 10;
    amax= 150; %max a dilations
    ind = 11; %min vertical height
    sgol = 31; 
    sgolclouds = 7; 
    range = 15400;
    cla=400;
    
elseif (strcmp(instrument,'CL31'))
    instrument = 'CL31'; 
    bsc_variable_name = 'Daily_Profile';
    time_variable_name = 'Daily_Time';
    precip_threshold1 = 2000;
    precip_threshold2 = 8000;
    cloud_threshold = 10;
    min_cloud_threshold_height = 10;
    cloud2fraction = 5; 
    cloud3fraction = 1.3;
    vres = 10;
    amax= 200; %max a dilations
    ind = 11; %min vertical height
    sgol = 31; 
    sgolclouds = 7; 
    range = 7700;
    cla=400;
    
elseif (strcmp(instrument,'8k'))
    instrument = '8k'; 
    bsc_variable_name = 'beta_att';
    time_variable_name = 'time';
    alt_variable_name =  'range';
    precip_threshold1 = 230000;
    precip_threshold2 = 400000;
    cloud_threshold = 3e-08;
    min_cloud_threshold_height = 1;
    cloud2fraction = 1.05; 
    cloud3fraction = 1.5;
    vres = 15;
    amax= 150; %max a dilations
    ind = 8; %min vertical height
    sgol = 31; 
    sgolclouds = 7; 
    cla=150;
else
    return;
end

%% Read in Files - will need to divide into Lufft vs Vaisala 

%Lufft lidar data
if (strcmp(instrument,'15k')) || (strcmp(instrument,'8k'))
%     endname = [location 'rolling.nc']; %removed for historical 
    endname = '20210501_ESX_CHM178003_000.nc'; 
    Time = ncread(endname,'time');
    alt = ncread(endname,'range');
    Braw = double(ncread(endname,bsc_variable_name));
    dateoffile = endname(1,1:8);
    %remove last 20 points 
    % Time = Time(1:end-20,1);
    % Braw = Braw(:,1:end-20);

    %for weird daily archived files modification for testing 
%     Time = Time(41:end,1);
%     Braw = Braw(:,41:end);
    load('OverlapFuctions.mat');

    %% Correct with new Lufft Overlap fuction 
%     for b=1:length(Time)
%          Braw_raw(:,b) = Braw(:,b) .* Old_Overlap(:);
%          Profile(:,b) =  Braw_raw(:,b) ./ New_Overlap(:); 
%     end
%      clearvars b Braw_raw New_Overlap Old_Overlap  
%     Braw=Profile; 
    
elseif (strcmp(instrument,'CL51')) || (strcmp(instrument,'CL31'))
  InFile = files(K).name;
  Braw = double(ncread(InFile,'beta_att'));
  Time = double(ncread(InFile,'time'));
  dateoffile = InFile(1,10:17);
    %make alt arrays 
    altstart=10;
    dalt = 10;
    N = range/vres;
    alt = transpose(altstart + (0:N-1)*dalt);
end   

    %Read-in settings files 
    setfileid = [location '_Settings.txt'];
    Settings = readtable(setfileid);
    mpbl = Settings.Max_PBL;
    srvec = datevec(Settings.Sunrise);
    srdt = datetime(Settings.Sunrise);
    srhrs = round(srvec(:,4) + (srvec(:,5)/60) + (srvec(:,6)/3600));
    ssvec = datevec(Settings.Sunset);
    ssdt = datetime(Settings.Sunset);
    sshrs = round(ssvec(:,4) + (ssvec(:,5)/60) + (ssvec(:,6)/3600));
    
    %modied for historical processing 
    filedatetime = datetime(dateoffile,'InputFormat', 'MMddyyyy');
    today = datevec(filedatetime); 
%     today = datevec('Jul-13-2020 00:00');
    
    %this grabs settings for today
    for s=1:length(ssvec)
        if srvec(s,1:3) == today(1,1:3)
            sr = srdt(s,1);
            ss = ssdt(s,1);
            ssprev = ssdt(s-1,1);
            ssnext = ssdt(s+1,1);
            srprev = srdt(s-1,1);
            srnext = srdt(s+1,1);
        end 
    end 
     
 clearvars s            
%% Format timestamps 
%Lufft Instruments 
%% Format Time Stamps
if (strcmp(instrument,'15k')) || (strcmp(instrument,'8k'))
    timevecUTC = datevec(Time/(3600*24) + datenum(1904, 1, 1, 0, 0, 0)); %edited for historical daily files 
    timevecUTC = datevec(Time/(3600*24) + datenum(1904, 1, 1, 0, 0, 0));
    timewaveUTC = datetime(timevecUTC,'InputFormat','YYYY-MM-dd HH:mm:SS');
    TimehrsUTC = timevecUTC(:,4) + (timevecUTC(:,5)/60) + (timevecUTC(:,6)/3600);
    timenumUTC = datenum(timewaveUTC);
end 

%change this for vaisala formatting 
if (strcmp(instrument,'CL51')) || (strcmp(instrument,'CL31'))
    timevecUTC = datevec(Time);
    timewaveUTC = datetime(timevecUTC);
    TimehrsUTC = timevecUTC(:,4) + (timevecUTC(:,5)/60) + (timevecUTC(:,6)/3600);
    timenumUTC = datenum(timewaveUTC);
end

%time expand 
t0 = datetime([timevecUTC(1,1:3) 0 0 0]); %grab day and zero it to 00:00:00 time 
t1 = datetime([timevecUTC(1,1:3) 23 59 59]); %grab day and end it at 23:59:59 time 

tExnum = (floor(abs(seconds(t0-t1))/16)); %how many 15s intervals should be expected from t0 to t1
timeexUTC = transpose(t0 + seconds(0:16:tExnum*16)); %create expected 24 hour - 15s intervals 
TimeExUTCnum = datenum(timeexUTC); %make it into a numerical format for time expansion
timeEXvec = datevec(timeexUTC);

for te = 1:length(TimeExUTCnum) 
    tendifex = abs(timenumUTC - TimeExUTCnum(te)); %find original (timenumUTC) time closest to new time (TimeExUTCnum)
    [iex,jex] = min(tendifex); %index closest original time 
    startex(te) = jex;
    if min(tendifex) > 8.1017e-05 %flag if closest profile is more than 7 seconds apart 
        nanflag(te) = 1; 
    else 
        nanflag(te) = 0; 
    end 
end 

%create index for calculating mean
stopex = startex(2:end)-1;
stopex = [stopex length(timenumUTC)]; 
stopex(stopex==0) = 1 ;
%time expand array
for it=1:length(startex)
        if nanflag(it) == 0
            pex = nanmean(Braw(:,startex(it):stopex(it)),2);
            Profile_TEx(:,it) = pex;
        else 
            Profile_TEx(:,it) = nan(N,1);
        end 
end 
%% 10-min averages and time formats 

%10min time 
% t3 = datetime('04-Dec-2016 00:00:00');
t3 = datetime([timeEXvec(1,1:6)]);
% t4 = datetime('04-Dec-2016 23:59:59');
t4 = datetime([timeEXvec(end,1:6)]);
tf=length(timeexUTC);
tfile = (floor(abs(seconds(t4-t3)/600)))*10;
% tfile = (floor(abs(seconds(timewaveUTC(1)-timewaveUTC(tf)))/600))*10;
UTC10 = transpose(t3 + minutes(0:10:tfile));
UTC10num = datenum(UTC10);
timevecUTC10 = datevec(UTC10);
TimehrsUTC10 = timevecUTC10(:,4) + (timevecUTC10(:,5)/60) + (timevecUTC10(:,6)/3600);
 
% 10 min backscatter array 
for t = 1:length(UTC10num) 
    tendif = abs(timenumUTC - UTC10num(t)); 
    [i10,j10] = min(tendif); 
    start(t) = j10;
end 
stop = start(2:end)-1;
stop = [stop length(timenumUTC)]; 

for ii=1:length(start)
        p = nanmean(Braw(:,start(ii):stop(ii)),2);
        bsc10min(:,ii) = p;
end 

%%  Processing: Precipitation Flag
%Precip array 1=precip
n= length(UTC10num);

for r=1:n
    % Define sample data and a threshold value.
    A=transpose(bsc10min(:,r));
    B=transpose(bsc10min(:,r));
    % Find logical vector where bsc > threshold
    binaryVector = A > precip_threshold1;
    binaryVector2 = B > precip_threshold2;
    % Label each region with a label - an "ID" number.
    [labeledVector, numRegions] = bwlabel(binaryVector);
    [labeledVector2, numRegions2] = bwlabel(binaryVector2);
    % Measure lengths of each region and the indexes
    measurements = regionprops(labeledVector, A, 'Area', 'PixelValues','PixelIdxList');
    measurements2 = regionprops(labeledVector2, B, 'Area', 'PixelValues','PixelIdxList');
    % Find regions where the area (length) are 3 or greater and put the values into a cell of a cell array
    [c,e]=max([measurements.Area]);[cc,ee]=max([measurements2.Area]);
    Area1min = floor(300/vres); %300meters
    Area2min = floor(50/vres); %50meters
    if numRegions > 0 && measurements(e).Area >= Area1min && measurements(e).PixelIdxList(1,1)==1
        if r==1 
            precip(r)= 1; precip(r+1)= 1;
        elseif r==144 %Exclude last profile
            precip(r)= 1; precip(r-1)= 1; 
        else 
            precip(r)= 1; precip(r-1)= 1; precip(r+1)= 1;
        end 
    elseif numRegions2 >0 && measurements2(ee).Area >= Area2min && measurements2(ee).PixelIdxList(1,1)==1
        if r==1 
            precip(r)= 1; precip(r+1)= 1;
        elseif r==144 %Exclude last profile
            precip(r)= 1; precip(r-1)= 1; 
        else 
            precip(r)= 1; precip(r-1)= 1; precip(r+1)= 1;
        end
    elseif r==1
        precip(r) = 0; 
    elseif r>1 && length(precip)==r-1
        precip(r) = 0; 
    elseif r>1 && length(precip)==r
        precip(r)*1;
    end  
end 

%% APPLY HAAR for clouds and PBL 
b= length(alt); %define vertical extent of dilations 

%for ml times 

for ifile=1:n
    %Define time-dependant settings 
%     if  (TimehrsUTC10(ifile) >= sr+5) && (TimehrsUTC10(ifile) <= ss+1) 
    if (UTC10(ifile) >= sr+hours(5)) && (UTC10(ifile) <= ss+hours(1)) || (UTC10(ifile) <= ss+hours(1)) && (UTC10(ifile) <= ssprev+hours(1))
        a= amax;
        hl = floor(mpbl(1,1)/vres); 
        ma(ifile)=a; %to index transitionts 
        detectlim(ifile) = hl;
%     elseif (TimehrsUTC10(ifile) >= sr+3)  && (TimehrsUTC10(ifile) < sr+5)
    elseif (UTC10(ifile) >= sr+hours(3))  && (UTC10(ifile) < sr+hours(5))
        a = floor(amax/2);
        hl = floor((mpbl(1,1)/1.5)/vres); 
        ma(ifile)=a; %to index transitionts
        detectlim(ifile) = hl;
%     elseif (TimehrsUTC10(ifile) >= ss+1) && (TimehrsUTC10(ifile) >= sr+3)
    elseif (UTC10(ifile) >= ss+hours(1)) && (UTC10(ifile) <= srnext+hours(3)) || (UTC10(ifile) >= ssprev+hours(1)) && (UTC10(ifile) <= sr+hours(3))
        a = floor(amax/3);
        hl = floor(500/vres); %NSL 500m limit
        ma(ifile)=a; %to index transitionts
        detectlim(ifile) = hl;
    end
    
    %RL and Cloud settings
    rl = floor(mpbl(1,1)/vres); %residual layer height limit at zmax
    rli = floor(500/vres); %residual layer min limit above 500m  
    rla = amax*2; % rl dilation 

    %Haar pre-smoothing
    sgolaysmoothed=sgolayfilt(bsc10min(:,ifile),3,sgol); %additional smoothing
    clsmoothed=sgolayfilt(bsc10min(1:end,ifile),3,sgolclouds); %additional smoothing
    astep= 1; %step increment in dilations
    smoothed(:,ifile) = (sgolaysmoothed); %smoothed arrays to save
    smoothedcl(:,ifile) = (clsmoothed); %smoothed arrays to save
    
    %Apply Haar 
    c1 = cwtnew(sgolaysmoothed(1:b), astep:a, 'haar');
    crl = cwtnew(sgolaysmoothed(1:b), astep:rla, 'haar');
    ccl = cwtnew(clsmoothed(1:b), astep:cla, 'haar');
    dilations{ifile}(1:a,1:b) = c1; %save all dilations for standard dev calculations
    dilations_rl{ifile}(1:rla,1:b) = crl; %save all dilations for standard dev calculations
    
    % Calculate mean transforms and arrays
    c4mean = nanmean(c1);
    c4RLmean = nanmean(crl);
    c4CLmean = nanmean(ccl);
    c4t = transpose(c4mean); 
    c4(1:b,ifile) = (c4t);
    c4trl = transpose(c4RLmean); 
    c4rl(1:b,ifile) = (c4trl);
    c4tcl = transpose(c4CLmean); 
    c4cl(1:b,ifile) = (c4tcl);
        
    % Setup local min detection arrays and min settings 
    min_index=[];
    min_index_rl=[];
    i=ind;% start at index defined above
    max_cloud_index=[];
    min_cloud_index=[];
    c=3; %start of cloud retrievals 
    % Flag high cloud signals in profiles 
    xC4 = max(c4cl(min_cloud_threshold_height:end,:)); %min/max for cloud top/base
    mC4 = min(c4cl(min_cloud_threshold_height:end,:)); %min/max for cloud top/base
    if xC4(ifile) > cloud_threshold && mC4(ifile)< -(cloud_threshold)
          flag(ifile) = 1;
        else
          flag(ifile) = 0;
    end
    
    % Find cloud heights - up to three clouds layers 
    if flag(ifile) == 1 %&& rain(ifile) == 0; %no rain and yes clouds 
            while (c < b-1)
            %Cloud bases
            if ((c4CLmean(c) > c4CLmean(c-1))&&(c4CLmean(c) > c4CLmean(c-2)))    % Check to see if this pt is less than the prev 2 pts
            if ((c4CLmean(c) > c4CLmean(c+1))&&(c4CLmean(c) > c4CLmean(c+2)))   % Check to see if this pt is less than the next 2 pts
            if c4CLmean(c)>cloud_threshold
                max_cloud_index = [max_cloud_index c];   % Update the index array
                [max_cloud_index_val] = c4CLmean(max_cloud_index);
            end 
            end
            end
            %Cloud tops 
            if ((c4CLmean(c) < c4CLmean(c-1))&&(c4CLmean(c) < c4CLmean(c-2)))    % Check to see if this pt is less than the prev 2 pts
            if ((c4CLmean(c) < c4CLmean(c+1))&&(c4CLmean(c) < c4CLmean(c+2)))   % Check to see if this pt is less than the next 2 pts
            if c4CLmean(c)< -cloud_threshold
                min_cloud_index = [min_cloud_index c];   % Update the index array
                [min_cloud_index_val] = c4CLmean(min_cloud_index);
            end 
            end
            end
            c = c + 1;
            end
            
            % Find all Cloud base max
            [Ccmax Icmax] = max(c4CLmean(max_cloud_index(1:end)));
            [cbi,cbj] = sort(c4CLmean(max_cloud_index(1:end)),'descend');
            % Find all Cloud top mins
            [Ccmin Icmin] = min(c4CLmean(min_cloud_index(1:end)));
            [cti,ctj] = sort(c4CLmean(min_cloud_index(1:end)));
            
            % Identify cloud layer #1 
            if length(max_cloud_index)>=1 && length(min_cloud_index)>=1
                Cloud(1,ifile) = alt(max_cloud_index(Icmax)); %base
                Cloud(2,ifile) = alt(min_cloud_index(Icmin)); %top
                Cloud(3,ifile) = Cloud(2,ifile) - Cloud(1,ifile); %depth
                if Cloud(2,ifile) < Cloud(1,ifile) %remove if cloud top is cut off 
                    Cloud(2,ifile) = nan; 
                end 
                if Cloud(3,ifile) > 500 && Cloud(1,ifile) <= 2000 %remove if cloud depth>500 and below 2000m for precip
                    Cloud(1,ifile) = nan;
                    Cloud(2,ifile) = nan;
                    Cloud(3,ifile) = nan;
                end 
            else 
                Cloud(1,ifile) = nan;
                Cloud(2,ifile) = nan;
                Cloud(3,ifile) = nan;
            end 

            % Identify cloud layer #2 
            if length(max_cloud_index)>1 && length(min_cloud_index)>1 && max_cloud_index_val(cbj(1,2)) > ((max_cloud_index_val(cbj(1,1)))/cloud2fraction) && (max_cloud_index(cbj(1,2)))>ind && abs(max_cloud_index(cbj(1,2))-max_cloud_index(cbj(1,1)))>floor(500/vres) %second cloud layer is present 
                Cloud(4,ifile) = alt(max_cloud_index(cbj(1,2))); %Cloud2 base
                Cloud(5,ifile) = alt(min_cloud_index(ctj(1,2))); %Cloud2 top
                Cloud(6,ifile) = Cloud(5,ifile) - Cloud(4,ifile); %depth
            else 
                Cloud(4,ifile) = nan;
                Cloud(5,ifile) = nan;
                Cloud(6,ifile) = nan;
            end 
             if Cloud(6,ifile) > 500 && Cloud(4,ifile) <= 2000 %remove if cloud depth>500 and below 2000m for precip
                Cloud(4,ifile) = nan;
                Cloud(5,ifile) = nan;
                Cloud(6,ifile) = nan;
             end 
             
             % Identify cloud layer #3 
             if length(max_cloud_index)>2 && length(min_cloud_index)>2 && max_cloud_index_val(cbj(1,3)) >((max_cloud_index_val(cbj(1,1)))/cloud2fraction) && (max_cloud_index(cbj(1,3)))>ind && abs(max_cloud_index(cbj(1,3))-max_cloud_index(cbj(1,2)))>floor(1000/vres) 
                Cloud(7,ifile) = alt(max_cloud_index(cbj(1,3))); %Cloud3 base
                Cloud(8,ifile) = alt(min_cloud_index(ctj(1,3))); %Cloud3 top
                Cloud(9,ifile) = Cloud(8,ifile) - Cloud(7,ifile); %Cloud 3 depth 
            else 
                Cloud(7,ifile) = nan;
                Cloud(8,ifile) = nan;
                Cloud(9,ifile) = nan;
            end
            if Cloud(9,ifile) > 500 && Cloud(7,ifile) <= 2000 %remove if cloud depth>500 and below 2000m for precip
                Cloud(7,ifile) = nan;
                Cloud(8,ifile) = nan;
                Cloud(9,ifile) = nan;
            end 
            
    else %if no clouds found 
            Cloud(1,ifile) = nan; %Cloud1 base
            Cloud(2,ifile) = nan; %Cloud1 top
            Cloud(3,ifile) = nan; %Cloud1 depth
            Cloud(4,ifile) = nan; %Cloud2 base
            Cloud(5,ifile) = nan; %Cloud2 top
            Cloud(6,ifile) = nan; %Cloud2 depth
            Cloud(7,ifile) = nan; %Cloud3 base
            Cloud(8,ifile) = nan; %Cloud3 top
            Cloud(9,ifile) = nan; %Cloud3 depth
    end 
        
    %Identify Residual layer only local mins  
   
    for ii=rli:rl
        if ((c4RLmean(ii) < c4RLmean(ii-1))&&(c4RLmean(ii-1) < c4RLmean(ii-2))&&(c4RLmean(ii-2) < c4RLmean(ii-3)) && c4RLmean(ii)<0)    % Check to see if this pt is less than the prev 2 pts
        if ((c4RLmean(ii) < c4RLmean(ii+1))&&(c4RLmean(ii+1) < c4RLmean(ii+2))&&(c4RLmean(ii+2) < c4RLmean(ii+3)) && c4RLmean(ii)<0)  % Check to see if this pt is less than the next 2 pts
        min_index_rl = [min_index_rl ii];   % Update the index array
        [min_index_val_rl] = c4RLmean(min_index_rl);
        end 
        end
    end
    [Cminrl Iminrl] = min(c4RLmean(min_index_rl(1:end))); %all min rl index
    
    % Identify remaining PBL local mins  
    while (i < hl)
        if ((c4mean(i) < c4mean(i-1))&&(c4mean(i-1) < c4mean(i-2))&&(c4mean(i-2) < c4mean(i-3)) && c4mean(i)<0)    % Check to see if this pt is less than the prev 2 pts
        if ((c4mean(i) < c4mean(i+1))&&(c4mean(i+1) < c4mean(i+2))&&(c4mean(i+2) < c4mean(i+3)) && c4mean(i)<0)   % Check to see if this pt is less than the next 2 pts
        min_index = [min_index i];   % Update the index array
        [min_index_val] = c4mean(min_index);
        end 
        end
        i = i + 1;
    end

    [Cmin Imin] = min(c4mean(min_index(1:end)));
    [mi,mj] = sort(c4mean(min_index(1:end)));
    [mirl,mjrl] = sort(c4mean(min_index_rl(1:end)));
   
    clearvars ca

    %PBL attribution based on mean of last 20mins PBL heights with 3 local
    try 
        PBL(1,ifile) = alt(min_index(Imin)); %First try PBL height 
            if ifile > 2 && length(min_index)>=2 && ~all(isnan(PBL(1,ifile-3:ifile-1))) %test continuity 
                if (abs(PBL(1,ifile) - nanmean(PBL(1,ifile-3:ifile-1)))) > contlim 
                    Imin = mj(:,2);
                    PBL(1,ifile) = alt(min_index(Imin)); %new PBL if continuity failed 
                    if (abs(PBL(1,ifile) - nanmean(PBL(1,ifile-3:ifile-1)))) > contlim
                        Imin = mj(:,3);
                        PBL(1,ifile) = alt(min_index(Imin)); %new PBL if continuity failed 
                        if (abs(PBL(1,ifile) - nanmean(PBL(1,ifile-3:ifile-1)))) > contlim
                            PBL(1,ifile) = nan; %if all fail
                        end 
                    end 
                end 
%             elseif ifile > 2 && all(isnan(PBL(1,ifile-3:ifile-1)))%if no pbl was found in last 20 mins 
%                 PBL(1,ifile) = nan;
            elseif  ifile > 2 && length(min_index)==1%if only one local min found
                if (abs(PBL(1,ifile) - nanmean(PBL(1,ifile-3:ifile)))) > contlim %test continuation 
                    PBL(1,ifile) = nan; %if all fail
                end
            end 
    catch %continue if no local mins found
        warning('Improper assignment with rectangular empty matrix.'); 
        PBL(1,ifile) = nan;
    end     

     if isempty(min_index) || length(min_index)<5
        min_index(end+1:5)=1; 
     end 
    Grads(1:5,ifile) = transpose(alt(min_index(1,1:5)));
    
    %RL attribution based on mean of last 20mins PBL heights with 3 local
    try 
        if (UTC10(ifile) >= sr+hours(5)) && (UTC10(ifile) <= ss+hours(1)) || (UTC10(ifile) <= ss+hours(1)) && (UTC10(ifile) <= ssprev+hours(1)) %remove daytime RL finds 
            PBL(2,ifile) = nan;
        else 
            PBL(2,ifile) = alt(min_index_rl(Iminrl)); %First try RL height
            if ifile > 2 && length(min_index_rl)>=2 %test continuity 
                if (abs(PBL(2,ifile) - nanmean(PBL(2,ifile-3:ifile-1)))) > contlim %exclude current pbl from mean
                    Iminrl = mjrl(:,2);
                    PBL(2,ifile) = alt(min_index_rl(Iminrl)); %new PBL if continuity failed 
                    if (abs(PBL(2,ifile) - nanmean(PBL(2,ifile-3:ifile-1)))) > contlim
                        Iminrl = mjrl(:,3);
                        PBL(2,ifile) = alt(min_index_rl(Iminrl)); %new PBL if continuity failed 
                        if (abs(PBL(2,ifile) - nanmean(PBL(2,ifile-3:ifile-1)))) > contlim
                            PBL(2,ifile) = nan; %if all fail
                        end 
                    end 
                end
            end
            if ifile > 2 && length(min_index_rl)==1 %if only one local min found 
                if (abs(PBL(2,ifile) - nanmean(PBL(2,ifile-3:ifile)))) > contlim
                            PBL(2,ifile) = nan; %if all fail
                end 
            end 
        end
    catch %continue if no local mins found
        warning('Improper assignment with rectangular empty matrix.');
        PBL(2,ifile) = nan;
    end
    
    if isempty(min_index_rl) || length(min_index_rl)<5
        min_index_rl(end+1:5)=1; 
     end 
    Grads_rl(1:5,ifile) = transpose(alt(min_index_rl(1,1:5)));
    
end %finished PBL identification 

clearvars t2 ifile hl i Imin a astep b c1 c4mean c4t clclouds Cmin mi min_index min_index_val mj r2smoothed start stop A binaryVector c cbi Ccmax Ccmin Cminrl cti ctj Icmax Icmin ii Iminrl labeledVector ma mai max_cloud_index
clearvars max_cloud_index_val mirl mjrl numRegions p r threshold trans ts xC4 cbj mC4 measurements min_cloud_index min_cloud_index_val
%% Standard Deviation Estimates 
for ifile=1:n
    %Define time-dependant settings 
    if  (UTC10(ifile) >= sr+hours(5)) && (UTC10(ifile) <= ss+hours(1)) || (UTC10(ifile) <= ss+hours(1)) && (UTC10(ifile) <= ssprev+hours(1))
        a= amax;
        hl = floor(mpbl(1,1)/vres); 
        ma(ifile)=a; %to index transitionts 
        detectlim(ifile) = hl;
    elseif (UTC10(ifile) >= sr+hours(3))  && (UTC10(ifile) < sr+hours(5))
        a = floor(amax/2);
        hl = floor((mpbl(1,1)/1.5)/vres); 
        ma(ifile)=a; %to index transitionts
        detectlim(ifile) = hl;
    elseif (UTC10(ifile) >= ss+hours(1)) && (UTC10(ifile) <= srnext+hours(3)) || (UTC10(ifile) >= ssprev+hours(1)) && (UTC10(ifile) <= sr+hours(3))
        a = floor(amax/3);
        hl = floor(500/vres); %NSL 500m limit
        ma(ifile)=a; %to index transitionts
        detectlim(ifile) = hl;
    end
    
    %RL and Cloud settings
    rl = floor(mpbl(1,1)/vres); %residual layer height limit at zmax
    rli = floor(500/vres); %residual layer min limit above 500m  
    rla = amax*2; % rl dilation 
    cla = amax*2; % cloud dilation
    
    %create arrays for identification 
    min_index_d=[]; singleidx = {}; min_index_rl =[]; singleidx_rl = {};
    
  % Using calculated PBL transforms above, calculate each local min 
    for ii = 1:a 
        for i = ind:hl
                if ((dilations{ifile}(ii,i) < dilations{ifile}(ii,i-1)) && (dilations{ifile}(ii,i-1) < dilations{ifile}(ii,i-2)) && (dilations{ifile}(ii,i-2) < dilations{ifile}(ii,i-3)) && dilations{ifile}(ii,i)<0)    % Check to see if this pt is less than the prev 2 pts
                if ((dilations{ifile}(ii,i) < dilations{ifile}(ii,i+1)) && (dilations{ifile}(ii,i+1) < dilations{ifile}(ii,i+2)) && (dilations{ifile}(ii,i+2) < dilations{ifile}(ii,i+3)) && dilations{ifile}(ii,i)<0)   % Check to see if this pt is less than the next 2 pts
                min_index_d= [min_index_d i];   % Update the index array
                end 
                end 
        end
                singleidx{ii} = min_index_d;
                min_idx_all{ifile} = singleidx;
                min_index_d=[]; 
    end 
    clearvars ii i 
    
    %Calculate PBL from all transform local mins 
    for aa = 1:length((min_idx_all{1,ifile})) 
            [Cmin Imin] = min(dilations{1,ifile}(aa,(min_idx_all{1,ifile}{1,aa})));
            try
                single_PBLs{1,ifile}(1,aa) = alt(min_idx_all{1,ifile}{1,aa}(1,Imin));
            catch %if no local mins found 
                warning('Improper assignment with rectangular empty matrix.');
                single_PBLs{1,ifile}(1,aa) = nan;
            end
    end
    
    %calculate stats for all transforms 
    PBL_stats(ifile,1) = nanmean(single_PBLs{1,ifile}(:));
    PBL_stats(ifile,2) = nanstd(single_PBLs{1,ifile}(:));
    PBL_stats(ifile,3) = nanmax(single_PBLs{1,ifile}(:));
    PBL_stats(ifile,4) = nanmin(single_PBLs{1,ifile}(:));
    
    %Using calculated RL transforms above, calculate each local min 
    for rr = 1:rla 
        for r=rli:rl
                if ((dilations_rl{ifile}(rr,r) < dilations_rl{ifile}(rr,r-1)) && (dilations_rl{ifile}(rr,r-1) < dilations_rl{ifile}(rr,r-2)) && (dilations_rl{ifile}(rr,r-2) < dilations_rl{ifile}(rr,r-3)) && (dilations_rl{ifile}(rr,r)<0))    % Check to see if this pt is less than the prev 2 pts
                if ((dilations_rl{ifile}(rr,r) < dilations_rl{ifile}(rr,r+1)) && (dilations_rl{ifile}(rr,r+1) < dilations_rl{ifile}(rr,r+2)) && (dilations_rl{ifile}(rr,r+2) < dilations_rl{ifile}(rr,r+3)) && (dilations_rl{ifile}(rr,r)<0))   % Check to see if this pt is less than the next 2 pts
                min_index_rl= [min_index_rl r];   % Update the index array
                end 
                end
        end
        singleidx_rl{rr} = min_index_rl;
        min_idx_all_rl{ifile} = singleidx_rl;
        min_index_rl =[];
    end 
    clear r rr

    %Calculate RL heights from all transform local mins 
    for bb = 1:length((min_idx_all_rl{1,ifile})) 
            [Cmin_rl Imin_rl] = min(dilations_rl{1,ifile}(bb,(min_idx_all_rl{1,ifile}{1,bb})));
            try
                single_RLs{1,ifile}(1,bb) = alt(min_idx_all_rl{1,ifile}{1,bb}(1,Imin_rl));
            catch
                warning('Improper assignment with rectangular empty matrix.');
                single_RLs{1,ifile}(1,bb) = nan;
            end
    end

    %calculate stats for all RL transforms 
    RL_stats(ifile,1) = nanmean(single_RLs{1,ifile}(:));
    RL_stats(ifile,2) = nanstd(single_RLs{1,ifile}(:));
    RL_stats(ifile,3) = nanmax(single_RLs{1,ifile}(:));
    RL_stats(ifile,4) = nanmin(single_RLs{1,ifile}(:));
end
CODE=zeros(3,n);

%% Clean up based on time, std, clouds, and rain flags 
for p=1:n 
    %check again if continuation was not successful
    if p>=2 && abs(PBL(1,p) - PBL(1,p-1)) >=300
        PBL(1,p) = nan;
    end
    if p>=2 && abs(PBL(2,p) - PBL(2,p-1)) >=300
        PBL(2,p) = nan;
    end
    %Check for PBL not found 
    if isnan(PBL(1,p)) 
        CODE(1,p) = 5;
    end 
    if isnan(PBL(2,p)) 
        CODE(2,p) = 5;
    end 
    %Remove PBL/RL during precip (CODE 1)
    if precip(p) == 1
        PBL_stats(p,2) = nan; 
        PBL_stats(p,1) = nan;
        PBL_stats(p,3) = nan;
        PBL(1,p) = nan; 
        PBL(2,p) = nan;
        PBL(3,p) = nan;
        RL_stats(p,2) = nan; 
        RL_stats(p,1) = nan;
        CODE(1:3,p) = 1; 
    end
    
    %Remove PBL/RL if it is too close to first cloud (for att, cPBL, etc)
    %(CODE 2)
    
    %define sr in 10min bins 
    if (UTC10(p) >= sr+hours(4)) && (UTC10(p) <= ss+hours(1)) || (UTC10(p) <= ss+hours(1)) && (UTC10(p) <= ssprev+hours(1))
          %define cloud removals 
        clrmv=1500; %cloud removal distance 
        if Cloud(1,p) < 2500
            cPBL(1,p) = Cloud(1,p);
        elseif Cloud(4,p) < 2500
            cPBL(1,p) = Cloud(4,p);
        elseif Cloud(7,p) < 2500
            cPBL(1,p) = Cloud(7,p);
        else 
            cPBL(1,p) = nan;
        end 
    else
        clrmv=300; %cloud removal distance 
        cPBL(1,p) = nan;
    end 
        
    if (Cloud(1,p)-PBL(1,p)) < clrmv
        if (1<p) && (p<144) %except outter bounds 
            PBL_stats(p,2) = nan; PBL_stats(p+1,2) = nan; PBL_stats(p-1,2) = nan;
            PBL_stats(p,1) = nan; PBL_stats(p+1,1) = nan; PBL_stats(p-1,1) = nan;
            PBL(1,p) = nan; PBL(1,p+1) = nan; PBL(1,p-1) = nan;
            PBL(2,p) = nan; PBL(2,p+1) = nan; PBL(2,p-1) = nan;
            PBL(3,p) = nan; PBL(3,p+1) = nan; PBL(3,p-1) = nan;
            RL_stats(p,2) = nan; RL_stats(p+1,2) = nan; RL_stats(p-1,2) = nan;
            RL_stats(p,1) = nan; RL_stats(p+1,1) = nan; RL_stats(p-1,1) = nan;
            CODE(1:3,p) = 2;
        elseif p==1 %if first retrieval of run
            PBL_stats(p,2) = nan; PBL_stats(p+1,2) = nan; 
            PBL_stats(p,1) = nan; PBL_stats(p+1,1) = nan; 
            PBL(1,p) = nan; PBL(1,p+1) = nan; 
            PBL(2,p) = nan; PBL(2,p+1) = nan; 
            PBL(3,p) = nan; PBL(3,p+1) = nan; 
            RL_stats(p,2) = nan; RL_stats(p+1,2) = nan; 
            RL_stats(p,1) = nan; RL_stats(p+1,1) = nan; 
            CODE(1:3,p) = 2;
        elseif p==144 %if last retrieval of run
            PBL_stats(p,2) = nan; PBL_stats(p-1,2) = nan; 
            PBL_stats(p,1) = nan; PBL_stats(p-1,1) = nan; 
            PBL(1,p) = nan; PBL(1,p-1) = nan; 
            PBL(2,p) = nan; PBL(2,p-1) = nan; 
            PBL(3,p) = nan; PBL(3,p-1) = nan; 
            RL_stats(p,2) = nan; RL_stats(p-1,2) = nan; 
            RL_stats(p,1) = nan; RL_stats(p-1,1) = nan; 
            CODE(1:3,p) = 2;
        end 
    end 

    %Remove PBL/RL if second cloud layer is within 100 of time-dependent height limit (CODE 2) 
    if Cloud(4,p) <= alt(detectlim(p))+100 %(for att, cPBL, etc) 
        PBL_stats(p,2) = nan; 
        PBL_stats(p,1) = nan;
        PBL(1,p) = nan; 
        PBL(2,p) = nan;
        PBL(3,p) = nan;
        RL_stats(p,2) = nan; 
        RL_stats(p,1) = nan;
        CODE(1:3,p) = 2;
    end 
    
    %Remove PBL/RL if third cloud layer is within 100 of time-dependent height limit (CODE 2)
    if Cloud(7,p) <= alt(detectlim(p)) + 100 %for attenuation 
        PBL_stats(p,2) = nan; 
        PBL_stats(p,1) = nan;
        PBL(1,p) = nan; 
        PBL(2,p) = nan;
        PBL(3,p) = nan;
        RL_stats(p,2) = nan; 
        RL_stats(p,1) = nan;
        CODE(1:3,p) = 2;
    end
    
    %is this repetitive? 
    %Remove if cloud1 is within 200m of RL height limit (CODE 2)
    if Cloud(1,p) <= alt(rl) +200 %for residual layer 
        PBL(2,p) = nan;
        RL_stats(p,2) = nan; 
        RL_stats(p,1) = nan;
        CODE(2,p) = 2;
    end 
    %Remove if cloud2 is within 200m of RL height limit (CODE 2)
    if Cloud(4,p) <= alt(rl) +200 %for residual layer 
        PBL(2,p) = nan;
        RL_stats(p,2) = nan; 
        RL_stats(p,1) = nan;
        CODE(2,p) = 2;
    end
    
     %NSL based on time limits 
   if (UTC10(p) >= ss+hours(1)) && (UTC10(p) >= srnext+hours(3)) || (UTC10(p) >= ssprev+hours(1)) && (UTC10(p) <= sr+hours(3)) %remove nsl times 
        PBL(3,p) = nan;
        CODE(3,p) = 6; 
    else 
        PBL(3,p) = PBL(1,p); %make nsl = pbl
        PBL(3,p) = nan; %then remove nsl 
        CODE(3,p) = 6;
    end
    
    if UTC10(p)>= sr+hours(3) && abs(PBL(1,p)-PBL(2,p))<100 %transition times when RL and ML are too close
        PBL(2,p) = nan; %remove RL 
        CODE(2,p) = 6;
    end 
%     if TimehrsUTC10(p)<= sr %nighttime keep RL only 
%         PBL(3,p) = nan; %remove nsl
%         PBL(1,p) = nan; %remove ml 
%         CODE(1,p) = 6;
%         CODE(3,p) = 6;
%     end 
%     
    %Flag for stds larger than 200m 
    if PBL_stats(p,2)>stdvlim
        PBL_flag(p,1) = 1;
    else 
        PBL_flag(p,1) = 0;
    end 
    if RL_stats(p,2)>stdvlim
        RL_flag(p,1) = 1;
    else 
        RL_flag(p,1) = 0;
    end 
    
    %Remove stds if continuity did not find a PBL 
    if isnan(PBL(1,p)) && isnan(PBL(3,p))
        PBL_stats(p,:) = nan;
    end 
    if isnan(PBL(2,p))
        RL_stats(p,:) = nan;
    end 
    
    %Remove outliers (CODE 3) 
    if p>2 && p<141 && isnan(PBL(1,p+1)) && isnan(PBL(1,p-1)) && isnan(PBL(1,p+2)) && isnan(PBL(1,p-2))
        PBL(1,p)= nan;
        CODE(1,p) = 31;
    end
    if p>2 && p<141 && isnan(PBL(2,p+1)) && isnan(PBL(2,p-1)) && isnan(PBL(2,p+2)) && isnan(PBL(2,p-2))
        PBL(2,p)= nan;
        CODE(2,p) = 31;
    end
    if p>2 && p<141 && isnan(PBL(3,p+1)) && isnan(PBL(3,p-1)) && isnan(PBL(3,p+2)) && isnan(PBL(3,p-2))
        PBL(3,p)= nan;
        CODE(3,p) = 31;
    end 
end 
clear p 

%Remove outliers again 
for p=1:n 
    if p>2 && p<141 && isnan(PBL(1,p+1)) && isnan(PBL(1,p-1)) && isnan(PBL(1,p+2)) && isnan(PBL(1,p-2))
        PBL(1,p)= nan;
        CODE(1,p) = 32;
    end
    if p>2 && p<141 && isnan(PBL(2,p+1)) && isnan(PBL(2,p-1)) && isnan(PBL(2,p+2)) && isnan(PBL(2,p-2))
        PBL(2,p)= nan;
        CODE(2,p) = 32;
    end
    if p>2 && p<141 && isnan(PBL(3,p+1)) && isnan(PBL(3,p-1)) && isnan(PBL(3,p+2)) && isnan(PBL(3,p-2))
        PBL(3,p)= nan;
        CODE(3,p) = 32;
    end
end 
clear p 

%Remove outliers again 
for p=1:n 
    if p>2 && p<141 && isnan(PBL(1,p+1)) && isnan(PBL(1,p-1)) && isnan(PBL(1,p+2)) && isnan(PBL(1,p-2))
        PBL(1,p)= nan;
        CODE(1,p) = 33;
    end
    if p>2 && p<141 && isnan(PBL(2,p+1)) && isnan(PBL(2,p-1)) && isnan(PBL(2,p+2)) && isnan(PBL(2,p-2))
        PBL(2,p)= nan;
        CODE(2,p) = 33;
    end
    if p>2 && p<141 && isnan(PBL(3,p+1)) && isnan(PBL(3,p-1)) && isnan(PBL(3,p+2)) && isnan(PBL(3,p-2))
        PBL(3,p)= nan;
        CODE(3,p) = 33;
    end
end 


%duplicate for PBL removal with std flags 
PBL_rmvd = PBL;
PBL_stats_rmvd = PBL_stats;
RL_stats_rmvd =RL_stats;

%remove all retrievals with std above 200m (CODE 4) 
for s=1:length(PBL_stats)
    if PBL_stats(s,2)>stdvlim
       PBL_rmvd(1,s)=nan; 
       PBL_rmvd(3,s)=nan; 
       PBL_stats_rmvd(s,2)= nan;
       CODE(1,s) = 4;
       CODE(3,s) = 4;
    end 
    if RL_stats(s,2) > stdvlim 
        PBL_rmvd(2,s) = nan;
        RL_stats_rmvd(s,2)= nan;
        CODE(2,s) = 4;
    end 
end 

clear p 

%Perform outlier cleanup 
for p=1:n 
    if p>2 && p<141 && isnan(PBL_rmvd(1,p+1)) && isnan(PBL_rmvd(1,p-1)) && isnan(PBL_rmvd(1,p+2)) && isnan(PBL_rmvd(1,p-2))
        PBL_rmvd(1,p)= nan;
        CODE(1,p) = 34;
    end
    if p>2 && p<141 && isnan(PBL_rmvd(2,p+1)) && isnan(PBL_rmvd(2,p-1)) && isnan(PBL_rmvd(2,p+2)) && isnan(PBL_rmvd(2,p-2))
        PBL_rmvd(2,p)= nan;
        CODE(2,p) = 34;
    end
    if p>2 && p<141 && isnan(PBL_rmvd(3,p+1)) && isnan(PBL_rmvd(3,p-1)) && isnan(PBL_rmvd(3,p+2)) && isnan(PBL_rmvd(3,p-2))
        PBL_rmvd(3,p)= nan;
        CODE(3,p) = 34;
    end
end 

clear p 

%Perform outlier cleanup 
for p=1:n 
    if p>2 && p<141 && isnan(PBL_rmvd(1,p+1)) && isnan(PBL_rmvd(1,p-1)) && isnan(PBL_rmvd(1,p+2)) && isnan(PBL_rmvd(1,p-2))
        PBL_rmvd(1,p)= nan;
        CODE(1,p) = 35;
    end
    if p>2 && p<141 && isnan(PBL_rmvd(2,p+1)) && isnan(PBL_rmvd(2,p-1)) && isnan(PBL_rmvd(2,p+2)) && isnan(PBL_rmvd(2,p-2))
        PBL_rmvd(2,p)= nan;
        CODE(2,p) = 35;
    end
    if p>2 && p<141 && isnan(PBL_rmvd(3,p+1)) && isnan(PBL_rmvd(3,p-1)) && isnan(PBL_rmvd(3,p+2)) && isnan(PBL_rmvd(3,p-2))
        PBL_rmvd(3,p)= nan;
        CODE(3,p) = 35;
    end   
end 

%Perform last outlier cleanup that cleans up group of remaining outliers 
for p=1:n 
    %remove grouped outliers of two 
    if p>2 && p<141 && ~isnan(PBL_rmvd(1,p+1)) && isnan(PBL_rmvd(1,p-1)) && isnan(PBL_rmvd(1,p-2)) && isnan(PBL_rmvd(1,p+2)) && isnan(PBL_rmvd(1,p+3))
        PBL_rmvd(1,p)= nan;
        CODE(1,p) = 36;
    end
    if p>2 && p<141 && ~isnan(PBL_rmvd(2,p+1)) && isnan(PBL_rmvd(2,p-1)) && isnan(PBL_rmvd(2,p-2)) && isnan(PBL_rmvd(2,p+2)) && isnan(PBL_rmvd(2,p+3))
        PBL_rmvd(2,p)= nan;
        CODE(2,p) = 36;
    end
    if p>2 && p<141 && ~isnan(PBL_rmvd(3,p+1)) && isnan(PBL_rmvd(3,p-1)) && isnan(PBL_rmvd(3,p-2)) && isnan(PBL_rmvd(3,p+2)) && isnan(PBL_rmvd(3,p+3))
        PBL_rmvd(3,p)= nan;
        CODE(3,p) = 36;
    end
    %remove grouped outliers of three 
    if p>2 && p<141 && ~isnan(PBL_rmvd(1,p+1)) && ~isnan(PBL_rmvd(1,p+2)) && isnan(PBL_rmvd(1,p-1)) && isnan(PBL_rmvd(1,p-2)) && isnan(PBL_rmvd(1,p+3)) && isnan(PBL_rmvd(1,p+4))
        PBL_rmvd(1,p)= nan;
        CODE(1,p) = 37;
    end
    if p>2 && p<141 && ~isnan(PBL_rmvd(2,p+1)) && ~isnan(PBL_rmvd(2,p+2)) && isnan(PBL_rmvd(2,p-1)) && isnan(PBL_rmvd(2,p-2)) && isnan(PBL_rmvd(2,p+3)) && isnan(PBL_rmvd(2,p+4))
        PBL_rmvd(2,p)= nan;
        CODE(2,p) = 37;
    end
    if p>2 && p<141 && ~isnan(PBL_rmvd(3,p+1)) && ~isnan(PBL_rmvd(3,p+2)) && isnan(PBL_rmvd(3,p-1)) && isnan(PBL_rmvd(3,p-2)) && isnan(PBL_rmvd(3,p+3)) && isnan(PBL_rmvd(3,p+4))
        PBL_rmvd(3,p)= nan;
        CODE(3,p) = 37;
    end
    %now remove sigle outliers 
    if p>3 && p<141 && all(isnan(PBL_rmvd(1,p+1:p+3))) && all(isnan(PBL_rmvd(1,p-1:p-3)))
        PBL_rmvd(1,p)= nan;
        CODE(1,p) = 38;
    end
    if p>3 && p<141 && all(isnan(PBL_rmvd(2,p+1:p+3))) && all(isnan(PBL_rmvd(2,p-1:p-3)))
        PBL_rmvd(2,p)= nan;
        CODE(2,p) = 38;
    end
    if p>3 && p<141 && all(isnan(PBL_rmvd(3,p+1:p+3))) && all(isnan(PBL_rmvd(3,p-1:p-3)))
        PBL_rmvd(3,p)= nan;
        CODE(3,p) = 38;
    end
    %Remove stds if continuity did not find a PBL 
    if isnan(PBL_rmvd(1,p)) && isnan(PBL_rmvd(3,p))
        PBL_stats_rmvd(p,:) = nan;
    end 
    if isnan(PBL_rmvd(2,p))
        RL_stats_rmvd(p,:) = nan;
    end
end 

%Clear gradients that were filled in with 10 
Grads_rl(Grads_rl==10)=nan;
Grads(Grads==10)=nan;
 
% clearvars DayC New_Overlap Olxd_Overlap PBL_stats RL_stats APA TimehrsEST c4trl dstart newMap timeexEST TimehrsEST10 camp dstart10 timeexUTC Cloud TimehrsUTC caxis dstop rain timenumEST EST10 TimehrsUTC10 ccl dstop10 rl timevecEST EST10num UTC cei flag rla timevecEST10 InFile UTC10 cla h rli timevecUTC clsmoothed i s timevecUTC10 Lft_clouds crl l1 sd timewaveEST PBL bsc10min l2 sgolaysmoothed timewaveUTC Profile c4 detectlim l3 smoothed tt Profile_TEx c4CLmean di l4 startex Time c4RLmean di10 l5 stopex TimeExESTnum c4cl dilations min_index_rl te TimeExUTCnum c4rl dilations_rl min_index_val_rl texfile TimeExUTCvec c4tcl n tf 
% clearvars  a p n aa b bb camp caxis cei Cmin Cmin_rl detectlim di di10 dstart dstart10 dstop dstop10 figname h hl i ifile Imin Imin_rl ind l1 l2 l3 l4 l5 l6 lgnd ma mai n newMap p min_index_d min_index_val_rl min_index_rl rl rla rli s sd startex stopex sy te texfile tf tfile Tout trans ts tt dilations dilations_cl dilations_rl single_PBLs single_RLs singleidx singleidx_rl min_idx_all min_idx_all_rl
if (strcmp(instrument,'CL51'))
    PlotCL51
elseif (strcmp(instrument,'15k'))
    PlotLufft
end 
end 
end 