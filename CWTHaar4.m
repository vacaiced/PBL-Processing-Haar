function [exitcode] = CWTHaar4(instrument,location)
% exitcode = instrument;
% cloud1and3 ID in cl51 may need && max_cloud_index(cbj(1,1))<400 line?
% Instrument
% added gradient record 
%adding new ways to detect continuity 
% tStart = cputime;
rwpfiles = dir('*.bln');
files = dir('*.h5');
pfiles = dir('*v4.mat');
nof = length(files);

for K =32:nof
contlim = 200;    
stdvlim = 400; 
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
    bsc_variable_name = 'Profile_bsc';
    time_variable_name = 'UnixTime_UTC';
    precip_threshold1 = 2000;
    precip_threshold2 = 8000;
    cloud_threshold = 25;% 10;
    min_cloud_threshold_height = 5; %50meters
    cloud2fraction = 5; 
    cloud3fraction = 1.3;
    vres = 10;
    amax= 150; %max a dilations
    ind = 12; %min vertical height
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
    endname = [location 'rolling.nc']; %removed for historical 
%     endname = '20161204_Catonsville, MD_CHM160112_000.nc'; 
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
  Braw = double(h5read(InFile,'/Profile_bsc'));
  Time = double(h5read(InFile,'/UnixTime_UTC'));
  dateoffile = InFile(1,11:18);
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
    filedatetime = datetime(dateoffile,'InputFormat', 'yyyyMMdd');
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
    timewaveUTC = datetime(Time, 'ConvertFrom','posixtime');
    timevecUTC = datevec(timewaveUTC);
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
tfile = (floor(abs(seconds(t4-t3)/600)))*10;
% tfile = (floor(abs(seconds(timewaveUTC(1)-timewaveUTC(tf)))/600))*10;
UTC10 = transpose(t3 + minutes(0:10:tfile));
UTC10num = datenum(UTC10);
timevecUTC10 = datevec(UTC10);
 
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
        if all(isnan(bsc10min(:,ii)))
            emptyflag(ii) = 1; 
        else 
            emptyflag(ii) = 0;
        end 
end 

%%Load previous day
curr_day_vec = sprintf('%02.0f',timevecUTC(1,1:3)); %current day
prev_day_vec = datevec(timeexUTC(1,1)-days(1)); 
prev_date = sprintf('%02.0f',prev_day_vec(1,1:3)); %previous day
% prev_in = [instrument '_' location '_' prev_date '_v4']; 

for pd = 1:length(pfiles)
    if str2double((pfiles(pd).name(10:17))) == str2double(prev_date)
        previn = pfiles(pd).name; 
        %load needed variables
        UTC10_prev = load(previn,'UTC10');   
        UTC10_prev = UTC10_prev.('UTC10');
        PBL_rmvd_prev = load(previn,'PBL_rmvd');   
        PBL_rmvd_prev = PBL_rmvd_prev.('PBL_rmvd');
    end 
end


%%  Processing: Precipitation Flag
%Precip array 1=precip
n= length(UTC10num);

for r=1:n
    % Define sample data and a threshold value.
    A=transpose(bsc10min(:,r));

    % Find logical vector where bsc > threshold
    binaryVector = A > precip_threshold1;
    binaryVector2 = A > precip_threshold2;
    % Label each region with a label - an "ID" number.
    [labeledVector, numRegions] = bwlabel(binaryVector);
    [labeledVector2, numRegions2] = bwlabel(binaryVector2);
    % Measure lengths of each region and the indexes
    measurements = regionprops(labeledVector, A, 'Area', 'PixelValues','PixelIdxList');
    measurements2 = regionprops(labeledVector2, A, 'Area', 'PixelValues','PixelIdxList');
    % Find regions where the area (length) are 3 or greater
    [c,e]=max([measurements.Area]);[cc,ee]=max([measurements2.Area]);
    Area1min = floor(250/vres); %300meters
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
    if (UTC10(ifile) >= sr+hours(5)) && (UTC10(ifile) <= ss+hours(1)) || (UTC10(ifile) <= ss+hours(1)) && (UTC10(ifile) <= ssprev+hours(1))
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
        
    % Setup min settings 
    i=ind;% start at index defined above
    c=3; %start of cloud retrievals 
    
    % Flag high cloud signals in profiles 
    xC4 = max(c4cl(min_cloud_threshold_height:end-100,:)); %min/max for cloud top/base
    mC4 = min(c4cl(min_cloud_threshold_height:end-100,:)); %min/max for cloud top/base
    if xC4(ifile) > cloud_threshold && mC4(ifile)< -(cloud_threshold)
          flag(ifile) = 1;
        else
          flag(ifile) = 0;
    end
    
    if flag(ifile) == 1 
        %find local min and max in cloud transforms 
        CFmin = islocalmin(c4CLmean,'MinProminence',20);
        CFmax = islocalmax(c4CLmean,'MinProminence',20);

        if islogical(CFmin)
            Cclmin = c4CLmean(CFmin);
            Iclmin = alt(CFmin)/vres;
        end 

        if islogical(CFmax)
            Cclmax = c4CLmean(CFmax);
            Iclmax = alt(CFmax)/vres;
        end 

        %sort by strength
        [climin,cljmin] = sort(Cclmin); 
        [climax,cljmax] = sort(Cclmax,'descend'); 
    
    %find cloud heigts 
        if ~isempty(Iclmax) && ~isempty(cljmin)
            Cloud(1,ifile) = alt(Iclmax(cljmax(1))); %base
            Cloud(2,ifile) = alt(Iclmin(cljmin(1))); %top
            Cloud(3,ifile) = Cloud(2,ifile) - Cloud(1,ifile); %depth

            if Cloud(2,ifile) < Cloud(1,ifile) %remove if cloud top is cut off 
                    Cloud(2,ifile) = nan; 
            end 
            if Cloud(3,ifile) > 600 && Cloud(1,ifile) <= 2000 %remove if cloud depth>500 and below 2000m for precip
                Cloud(1,ifile) = nan;
                Cloud(2,ifile) = nan;
                Cloud(3,ifile) = nan;
            end

            if length(Iclmax)>1 && length(Iclmin)>1 && Cclmax(cljmax(1,2)) > ((Cclmax(cljmax(1,1)))/cloud2fraction) && (Cclmax(cljmax(1,2)))>ind && abs(Cclmax(cljmax(1,2))-Cclmax(cljmax(1,1)))>floor(500/vres) %second cloud layer is present
                Cloud(4,ifile) = alt(Iclmax(cljmax(2))); %base
                Cloud(5,ifile) = alt(Iclmin(cljmin(2))); %top
                Cloud(6,ifile) = Cloud(5,ifile) - Cloud(4,ifile); %depth
            else 
                Cloud(4,ifile) = nan;
                Cloud(5,ifile) = nan;
                Cloud(6,ifile) = nan;
            end 

            if Cloud(6,ifile) > 600 && Cloud(4,ifile) <= 2000 %remove if cloud depth>500 and below 2000m for precip
                Cloud(4,ifile) = nan;
                Cloud(5,ifile) = nan;
                Cloud(6,ifile) = nan;
            end

            %third cloud layer 
            if length(Iclmax)>2 && length(Iclmin)>2 && Cclmax(cljmax(1,3)) >((Cclmax(cljmax(1,1)))/cloud2fraction) && (Cclmax(cljmax(1,3)))>ind && abs(Cclmax(cljmax(1,3))-Cclmax(cljmax(1,2)))>floor(1000/vres) 
                Cloud(7,ifile) = alt(Iclmax(cljmax(3))); %base
                Cloud(8,ifile) = alt(Iclmin(cljmin(3))); %top
                Cloud(9,ifile) = Cloud(8,ifile) - Cloud(7,ifile); %depth
            else 
                Cloud(7,ifile) = nan;
                Cloud(8,ifile) = nan;
                Cloud(9,ifile) = nan;
            end

            if Cloud(9,ifile) > 600 && Cloud(7,ifile) <= 2000 %remove if cloud depth>500 and below 2000m for precip
                Cloud(7,ifile) = nan;
                Cloud(8,ifile) = nan;
                Cloud(9,ifile) = nan;
            end 
        end 
        
    elseif flag(ifile) == 0 %if no clouds found 
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
  
    clearvars CFmin CFmax Cclmin Iclmin Cclmax Iclmax climin cljmin climax cljmax
    
    %Find mins for PBL dectection
    TF_rl = islocalmin(c4RLmean(1:rl));
    if islogical(TF_rl)  
        Cminrl1 = c4RLmean(TF_rl);
        Iminrl1 = alt(TF_rl)/vres;
    end 
     
    TF = islocalmin(c4mean(1:hl));
    if islogical(TF) %&& alt(TF)>ind
        Cmin1 = c4mean(TF);
        Imin1 = alt(TF)/vres;
    end 
    
    Cmin = Cmin1(Imin1>ind); %limit bottom of RL retrieval to above NSL 500m limit 
    Imin = Imin1(Imin1>ind); %limit bottom of RL retrieval to above NSL 500m limit 
    [mi,mj] = sort(Cmin);
    
    Cminrl = Cminrl1(Iminrl1>rli); %limit bottom of RL retrieval to above NSL 500m limit 
    Iminrl = Iminrl1(Iminrl1>rli); %limit bottom of RL retrieval to above NSL 500m limit 
    [mirl,mjrl] = sort(Cminrl);
    clearvars ca
    
    %PBL attribution based on mean of last 20mins PBL heights with 3 local
    try 
        PBL(1,ifile) = alt(Imin(mj(1)));
            if ifile < 4 && exist('PBL_rmvd_prev')
                PBL_conc = horzcat(PBL_rmvd_prev(1,:), PBL(1,:)); 
                ifilep = length(PBL_conc);
                if (abs(PBL_conc(1,ifilep) - nanmean(PBL_conc(1,ifilep-3:ifilep-1)))) > contlim 
                    PBL(1,ifile) = alt(Imin(mj(2)));
                    PBL_conc(1,ifilep) = alt(Imin(mj(2)));
                    if (abs(PBL_conc(1,ifilep) - nanmean(PBL_conc(1,ifilep-3:ifilep-1)))) > contlim
                        PBL(1,ifile) = alt(Imin(mj(3)));
                        PBL_conc(1,ifilep) = alt(Imin(mj(3)));
                        if (abs(PBL_conc(1,ifilep) - nanmean(PBL_conc(1,ifilep-3:ifilep-1)))) > contlim
                            PBL(1,ifile) = nan; %if all fail
                            PBL_conc(1,ifilep) = nan;
                            fail(1,ifile) = 100; 
                        end 
                    end 
                end 
            elseif ifile >= 4 && length(Imin)>=3 && ~all(isnan(PBL(1,ifile-3:ifile-1))) %test continuity 
                if (abs(PBL(1,ifile) - nanmean(PBL(1,ifile-3:ifile-1)))) > contlim && Cmin(mj(2)) < (Cmin(mj(1))/4) 
                    PBL(1,ifile) = alt(Imin(mj(2)));
                    if (abs(PBL(1,ifile) - nanmean(PBL(1,ifile-3:ifile-1)))) > contlim
                        PBL(1,ifile) = alt(Imin(mj(3)));
                        if (abs(PBL(1,ifile) - nanmean(PBL(1,ifile-3:ifile-1)))) > contlim
                            PBL(1,ifile) = nan; %if all fail
                            fail(1,ifile) = 100; 
                        end 
                    end 
                end 
            elseif  ifile >= 4 && length(Imin)<3 %if only one local min found
                if (abs(PBL(1,ifile) - nanmean(PBL(1,ifile-3:ifile-1)))) > contlim %test continuation 
                    PBL(1,ifile) = nan; %if all fail
                    fail(1,ifile) = 100; 
                end
            end 
    catch %continue if no local mins found
        warning('Improper assignment with rectangular empty matrix.'); 
        PBL(1,ifile) = nan;
    end     

    try 
        if (UTC10(ifile) >= sr+hours(5)) && (UTC10(ifile) <= ss+hours(1)) || (UTC10(ifile) <= ss+hours(1)) && (UTC10(ifile) <= ssprev+hours(1)) %remove daytime RL finds 
            PBL(2,ifile) = nan;
        else 
            PBL(2,ifile) = alt(Iminrl(mjrl(1))); %First try RL height
            if ifile < 4 && exist('PBL_rmvd_prev')
                RL_conc = horzcat(PBL_rmvd_prev(2,:), PBL(2,:)); 
                ifilep = length(RL_conc);
                if (abs(RL_conc(1,ifilep) - nanmean(RL_conc(1,ifilep-3:ifilep-1)))) > contlim
                    PBL(2,ifile) = alt(Iminrl(mjrl(2)));
                    RL_conc(1,ifilep) = alt(Iminrl(mjrl(2)));
                    if (abs(RL_conc(1,ifilep) - nanmean(RL_conc(1,ifilep-3:ifilep-1)))) > contlim
                        PBL(2,ifile) = alt(Iminrl(mjrl(3)));
                        RL_conc(1,ifilep) = alt(Iminrl(mjrl(3)));
                        if (abs(RL_conc(1,ifilep) - nanmean(RL_conc(1,ifilep-3:ifilep-1)))) > contlim
                            PBL(2,ifile) = nan; %if all fail
                            RL_conc(1,ifilep) = nan;
                            fail(2,ifile) = 100; 
                        end 
                    end 
                end 
            
            elseif ifile >= 4 && length(Iminrl)>=3 %test continuity 
                if (abs(PBL(2,ifile) - nanmean(PBL(2,ifile-3:ifile-1)))) > contlim %exclude current pbl from mean
                    PBL(2,ifile) = alt(Iminrl(mjrl(2)));
                    if (abs(PBL(2,ifile) - nanmean(PBL(2,ifile-3:ifile-1)))) > contlim
                        PBL(2,ifile) = alt(Iminrl(mjrl(3)));
                        if (abs(PBL(2,ifile) - nanmean(PBL(2,ifile-3:ifile-1)))) > contlim
                            PBL(2,ifile) = nan; %if all fail
                            fail(2,ifile) = 100; 
                        end 
                    end 
                end
           
        
            elseif ifile >= 4 && length(Iminrl)<3 %if only one local min found
                if (abs(PBL(2,ifile) - nanmean(PBL(2,ifile-3:ifile-1)))) > contlim
                        PBL(2,ifile) = nan; %if all fail
                        fail(2,ifile) = 100; 
                end 
            end 
        end
    catch %continue if no local mins found
        warning('Improper assignment with rectangular empty matrix.');
        PBL(2,ifile) = nan;
    end
    
    %save five gradients
    apg1 = vertcat(alt(Imin(mj)), nan, nan, nan, nan, nan);
    apg2 = vertcat(alt(Iminrl(mjrl)), nan, nan, nan, nan, nan);
    
    Gradients_NSLML(:,ifile) = apg1(1:5); 
    Gradients_RL(:,ifile) = apg2(1:5);
    clearvars apg1 apg2 
end %finished PBL identification 

clearvars TF TF_rl Cminrl1 Iminrl1 t2 ifile hl i Imin a astep b c1 c4mean c4t clclouds Cmin mi min_index min_index_val mj r2smoothed start stop A binaryVector c cbi Ccmax Ccmin Cminrl cti ctj Icmax Icmin ii Iminrl labeledVector ma mai max_cloud_index
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

    for ii = 1:a
        TF_all{ifile}(ii,1:hl) = islocalmin(dilations{ifile}(ii,1:hl));
        if islogical(TF_all{ifile}(ii,1:hl))
           val_single{ii} = dilations{ifile}(ii,TF_all{ifile}(ii,:)); %actual values of mins
           alt_single{ii} = alt(TF_all{ifile}(ii,:));
           [ai,aj] = min(val_single{ii}); 
           PBL_single{ii} = alt_single{ii}(aj);
           PBL_std_all{ifile} = PBL_single; 
        end 
    end 
    clearvars ii i     
    
    %calculate stats for all transforms 
    if emptyflag(ifile) == 0
        PBL_stats(ifile,1) = nanmean(cell2mat(PBL_std_all{1,ifile}(:)));
        PBL_stats(ifile,2) = nanstd(cell2mat(PBL_std_all{1,ifile}(:)));
        PBL_stats(ifile,3) = nanmax(cell2mat(PBL_std_all{1,ifile}(:)));
        PBL_stats(ifile,4) = nanmin(cell2mat(PBL_std_all{1,ifile}(:)));
    else 
        PBL_stats(ifile,1) = nan; 
        PBL_stats(ifile,2) = nan;
        PBL_stats(ifile,3) = nan;
        PBL_stats(ifile,4) = nan;
    end 
    
    for rr = 1:rla
        TF_rl_all{ifile}(rr,1:rl) = islocalmin(dilations_rl{ifile}(rr,1:rl));
        if islogical(TF_rl_all{ifile}(rr,1:rl))
           val_rl_single{rr} = dilations_rl{ifile}(rr,TF_rl_all{ifile}(rr,:)); %actual values of mins
           alt_rl_single{rr} = alt(TF_rl_all{ifile}(rr,:));
           [ri,rj] = min(val_rl_single{rr}); 
           RL_single{rr} = alt_rl_single{rr}(rj);
           RL_std_all{ifile} = RL_single; 
        end 
    end 
    
    if emptyflag(ifile) == 0
        RL_stats(ifile,1) = nanmean(cell2mat(RL_std_all{1,ifile}(:)));
        RL_stats(ifile,2) = nanstd(cell2mat(RL_std_all{1,ifile}(:)));
        RL_stats(ifile,3) = nanmax(cell2mat(RL_std_all{1,ifile}(:)));
        RL_stats(ifile,4) = nanmin(cell2mat(RL_std_all{1,ifile}(:)));
    else 
        RL_stats(ifile,1) = nan;
        RL_stats(ifile,2) = nan;
        RL_stats(ifile,3) = nan;
        RL_stats(ifile,4) = nan; 
    end 
end
CODE=zeros(3,n);

%% Clean up based on time, std, clouds, and rain flags
PBL_stats_rmvd = PBL_stats;
PBL_rmvd = PBL;
RL_stats_rmvd = RL_stats;

for p=1:n 
    %Remove PBL/RL during precip (CODE 1)
    if precip(p) == 1
        PBL_stats_rmvd(p,:) = nan; 
        PBL_rmvd(:,p) = nan; 
        RL_stats_rmvd(p,:) = nan; 
        CODE(:,p) = 1; 
    end
    
    %define sr in 10min bins 
    if (UTC10(p) >= sr+hours(4)) && (UTC10(p) <= ss+hours(1)) || (UTC10(p) <= ss+hours(1)) && (UTC10(p) <= ssprev+hours(1))
          %define cloud removals 
        clrmv=1000; %cloud removal distance 
        if Cloud(1,p) < mpbl(1,1) && precip(p)==0
            cPBL(1,p) = Cloud(1,p);
        elseif Cloud(4,p) < mpbl(1,1) && precip(p)==0
            cPBL(1,p) = Cloud(4,p);
        elseif Cloud(7,p) < mpbl(1,1) && precip(p)==0
            cPBL(1,p) = Cloud(7,p);
        else 
            cPBL(1,p) = nan;
        end 
    else
        clrmv=300; %cloud removal distance 
        cPBL(1,p) = nan;
    end 
        
    clrmv = 1000;
    if (Cloud(1,p)-PBL(1,p)) < clrmv
        if (1<p) && (p<144) %except outter bounds 
            PBL_stats_rmvd(p,1) = nan; PBL_stats_rmvd(p+1,1) = nan; PBL_stats_rmvd(p-1,1) = nan;
            PBL_rmvd(1,p) = nan; %PBL_rmvd(1,p+1) = nan; PBL_rmvd(1,p-1) = nan;
            CODE(1:3,p) = 2;
        elseif p==1 %if first retrieval of run
            PBL_stats_rmvd(p,1) = nan; PBL_stats_rmvd(p+1,1) = nan; 
            PBL_rmvd(1,p) = nan; %PBL_rmvd(1,p+1) = nan; 
            CODE(1:3,p) = 2;
        elseif p==144 %if last retrieval of run
            PBL_stats_rmvd(p,1) = nan; PBL_stats_rmvd(p-1,1) = nan; 
            PBL_rmvd(1,p) = nan; %PBL_rmvd(1,p-1) = nan; 
            CODE(1:3,p) = 2;
        end 
    end 

    if (Cloud(1,p)-PBL(2,p)) < clrmv
        if (1<p) && (p<144) %except outter bounds 
            PBL_stats_rmvd(p,2) = nan; PBL_stats_rmvd(p+1,2) = nan; PBL_stats_rmvd(p-1,2) = nan;
            PBL_rmvd(2,p) = nan; PBL_rmvd(2,p+1) = nan; PBL_rmvd(2,p-1) = nan;
            RL_stats_rmvd(p,2) = nan; RL_stats_rmvd(p+1,2) = nan; RL_stats_rmvd(p-1,2) = nan;
            RL_stats_rmvd(p,1) = nan; RL_stats_rmvd(p+1,1) = nan; RL_stats_rmvd(p-1,1) = nan;
            CODE(1:3,p) = 2;
        elseif p==1 %if first retrieval of run
            PBL_stats_rmvd(p,2) = nan; PBL_stats_rmvd(p+1,2) = nan; 
            PBL_rmvd(2,p) = nan; PBL_rmvd(2,p+1) = nan; 
            RL_stats_rmvd(p,2) = nan; RL_stats_rmvd(p+1,2) = nan; 
            RL_stats_rmvd(p,1) = nan; RL_stats_rmvd(p+1,1) = nan; 
            CODE(1:3,p) = 2;
        elseif p==144 %if last retrieval of run
            PBL_stats_rmvd(p,2) = nan; PBL_stats_rmvd(p-1,2) = nan; 
            PBL_rmvd(2,p) = nan; PBL_rmvd(2,p-1) = nan; 
            RL_stats_rmvd(p,2) = nan; RL_stats_rmvd(p-1,2) = nan; 
            RL_stats_rmvd(p,1) = nan; RL_stats_rmvd(p-1,1) = nan; 
            CODE(1:3,p) = 2;
        end 
    end
    %Remove PBL/RL if second cloud layer is within 100 of time-dependent height limit (CODE 2) 
    if Cloud(4,p) <= alt(detectlim(p))+100 %(for att, cPBL, etc) 
        PBL_stats_rmvd(p,2) = nan; 
        PBL_stats_rmvd(p,1) = nan;
        PBL_rmvd(1,p) = nan; 
        PBL_rmvd(2,p) = nan;
        PBL_rmvd(3,p) = nan;
        RL_stats_rmvd(p,2) = nan; 
        RL_stats_rmvd(p,1) = nan;
        CODE(1:3,p) = 2;
    end 
    
    %Remove PBL/RL if third cloud layer is within 100 of time-dependent height limit (CODE 2)
    if Cloud(7,p) <= alt(detectlim(p)) + 100 %for attenuation 
        PBL_stats_rmvd(p,2) = nan; 
        PBL_stats_rmvd(p,1) = nan;
        PBL_rmvd(1,p) = nan; 
        PBL_rmvd(2,p) = nan;
        PBL_rmvd(3,p) = nan;
        RL_stats_rmvd(p,2) = nan; 
        RL_stats_rmvd(p,1) = nan;
        CODE(1:3,p) = 2;
    end
    
    %is this repetitive? 
    %Remove if cloud1 is within 200m of RL height limit (CODE 2)
    if Cloud(1,p) <= alt(rl) +200 %for residual layer 
        PBL_rmvd(2,p) = nan;
        RL_stats_rmvd(p,2) = nan; 
        RL_stats_rmvd(p,1) = nan;
        CODE(2,p) = 2;
    end 
    %Remove if cloud2 is within 200m of RL height limit (CODE 2)
    if Cloud(4,p) <= alt(rl) +200 %for residual layer 
        PBL_rmvd(2,p) = nan;
        RL_stats_rmvd(p,2) = nan; 
        RL_stats_rmvd(p,1) = nan;
        CODE(2,p) = 2;
    end
    
     %NSL based on time limits 
   if (UTC10(p) >= ss+hours(1)) && (UTC10(p) >= srnext+hours(3)) || (UTC10(p) >= ssprev+hours(1)) && (UTC10(p) <= sr+hours(3)) %remove nsl times 
        PBL_rmvd(3,p) = nan;
        CODE(3,p) = 6; 
    else 
        PBL_rmvd(3,p) = PBL_rmvd(1,p); %make nsl = pbl
        PBL_rmvd(3,p) = nan; %then remove nsl 
        CODE(3,p) = 6;
    end
    
    if UTC10(p)>= sr+hours(3) && abs(PBL_rmvd(1,p)-PBL_rmvd(2,p))<100 %transition times when RL and ML are too close
        PBL_rmvd(2,p) = nan; %remove RL 
        CODE(2,p) = 6;
    end 
    
    %invalid: create 30min windows to remove any windows with 2+ nans 
    missing1 = isnan(PBL_rmvd(1,:));
    M1 = reshape(missing1,3,48);
    Mflag1 = sum(M1);
    
    missing2 = isnan(PBL_rmvd(2,:));
    M2 = reshape(missing2,3,48);
    Mflag2 = sum(M2);
    
    m0 = [1:3:144]; 
    m1 = [3:3:144];
    for s=1:length(Mflag1)
        if Mflag1(s)>=2
            PBL_rmvd(1,m0(s):m1(s)) = nan; 
        end
    end 
       
    for s=1:length(Mflag2)
        if Mflag2(s)>=2
            PBL_rmvd(2,m0(s):m1(s)) = nan; 
        end
    end
    
    %Remove outliers (CODE 3) 
    if p>2 && p<141 && isnan(PBL_rmvd(1,p+1)) && isnan(PBL_rmvd(1,p-1)) && isnan(PBL_rmvd(1,p-2))
        PBL_rmvd(1,p)= nan;
        CODE(1,p) = 31;
    end
    if p>2 && p<141 && isnan(PBL_rmvd(2,p+1)) && isnan(PBL_rmvd(2,p-1)) && isnan(PBL_rmvd(2,p-2))
        PBL_rmvd(2,p)= nan;
        CODE(2,p) = 31;
    end

    missing3 = isnan(PBL_rmvd(1,:));
    M3 = reshape(missing3,6,24);
    Mflag3 = sum(M3);
    
    missing4 = isnan(PBL_rmvd(2,:));
    M4 = reshape(missing4,6,24);
    Mflag4 = sum(M4);
    
    m2 = [1:6:144]; 
    m3 = [6:6:144];
    for s=1:length(Mflag3)
        if Mflag3(s)>=4
            PBL_rmvd(1,m2(s):m3(s)) = nan; 
        end
    end 
       
    for s=1:length(Mflag4)
        if Mflag4(s)>=4
            PBL_rmvd(2,m2(s):m3(s)) = nan; 
        end
    end
  
% clearvars DayC New_Overlap Olxd_Overlap PBL_stats RL_stats APA TimehrsEST c4trl dstart newMap timeexEST TimehrsEST10 camp dstart10 timeexUTC Cloud TimehrsUTC caxis dstop rain timenumEST EST10 TimehrsUTC10 ccl dstop10 rl timevecEST EST10num UTC cei flag rla timevecEST10 InFile UTC10 cla h rli timevecUTC clsmoothed i s timevecUTC10 Lft_clouds crl l1 sd timewaveEST PBL bsc10min l2 sgolaysmoothed timewaveUTC Profile c4 detectlim l3 smoothed tt Profile_TEx c4CLmean di l4 startex Time c4RLmean di10 l5 stopex TimeExESTnum c4cl dilations min_index_rl te TimeExUTCnum c4rl dilations_rl min_index_val_rl texfile TimeExUTCvec c4tcl n tf 
% clearvars  a p n aa b bb camp caxis cei Cmin Cmin_rl detectlim di di10 dstart dstart10 dstop dstop10 figname h hl i ifile Imin Imin_rl ind l1 l2 l3 l4 l5 l6 lgnd ma mai n newMap p min_index_d min_index_val_rl min_index_rl rl rla rli s sd startex stopex sy te texfile tf tfile Tout trans ts tt dilations dilations_cl dilations_rl single_PBLs single_RLs singleidx singleidx_rl min_idx_all min_idx_all_rl
end

if (strcmp(instrument,'CL51') || (strcmp(instrument,'CL31')))
    ContourCL51UTC
elseif (strcmp(instrument,'15k'))
    PlotLufft
end 
end
% cputime - tStart
