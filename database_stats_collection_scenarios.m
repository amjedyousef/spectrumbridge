
tic;
clear all;
close all;
clc;

%Switch which database you want to query
spectrumbridge_test=1; %Query spectrumBridge database

%%
%Create legend for the figures
legend_string={'SpectrumBridge'};
legend_flag=[spectrumbridge_test];
legend_string(find(~legend_flag))=[];

%%
%Select which scenario to test
message_size_distribution=1;
response_error_calculation=0;
delay_distribution_per_location=0;
delay_distribution_area=0;

%%
%Plot parameters
ftsz=16;

%%
%Path to save files (select your own)
my_path='/home/amjed/Documents/Gproject/workspace/data/WSDB_DATA/';

%%
%General querying parameters


%Global SpectrumBridge parameters (refer to WSDB_TVBD_Interface_v1.0.pdf [provided by Peter Stanforth])
AntennaHeight='30'; %In meters; Ignored for personal/portable devices
DeviceType='3'; %Examples: 8-Fixed, 3-40 mW Mode II personal/portable; 4-100 mW Mode II personal/portable

if message_size_distribution==1
    
    %Location of start and finish query
    %Query start location
    WSDB_data{1}.name='LA'; %Los Aneles, CA, USA (Wilshire Blvd 1) [downtown]
    WSDB_data{1}.latitude='34.047955';
    WSDB_data{1}.longitude='-118.256013';

    
    %Query finish location
    WSDB_data{2}.name='CB'; %Carolina Beach, NC, USA [ocean coast]
    WSDB_data{2}.latitude='34.047955';
    WSDB_data{2}.longitude='-77.885639';

    
    longitude_start=str2num(WSDB_data{1}.longitude); %Start of the spectrum scanning trajectory
    longitude_end=str2num(WSDB_data{2}.longitude); %End of spectrum scanning trajectory
    
    longitude_interval=100;
    longitude_step=(longitude_end-longitude_start)/longitude_interval;
    
    delay_spectrumbridge=[];
    
    ggl_cnt=0;
    in = 0 ;
    for xx=longitude_start:longitude_step:longitude_end
        in=in+1;
        fprintf('Query no.: %d\n',in)
        
        %Fetch location data
        latitude=WSDB_data{1}.latitude;
        longitude=num2str(xx);
        
        instant_clock=clock; %Save clock for file name (if both WSDBs are queried)

        if spectrumbridge_test==1
            %Query SpectrumBridge
            instant_clock=clock; %Start clock again if scanning only one database
            cd([my_path,'/spectrumbridge']);
            delay_spectrumbridge_tmp_r=0;
            if DeviceType=='3'
                [msg_spectrumbridge,delay_spectrumbridge_tmp_r]=database_connect_spectrumbridge_register(...
                    AntennaHeight,DeviceType,latitude,longitude,[my_path,'/spectrumbridge']);
            end
            [msg_spectrumbridge,delay_spectrumbridge_tmp,error_spectrumbridge_tmp]=...
                database_connect_spectrumbridge(DeviceType,latitude,longitude);
            var_name=(['spectrumbridge_',num2str(longitude),'_',datestr(instant_clock, 'DD_mmm_YYYY_HH_MM_SS')]);
            fprintf('SpectrumBridge\n')
            if error_spectrumbridge_tmp==0
                dlmwrite([var_name,'.txt'],msg_spectrumbridge,'');
                delay_spectrumbridge_tmp=delay_spectrumbridge_tmp+delay_spectrumbridge_tmp_r;
                delay_spectrumbridge=[delay_spectrumbridge,delay_spectrumbridge_tmp];
            end
        end
    end

    if spectrumbridge_test==1
        %Clear old query results
        cd([my_path,'/spectrumbridge']);
        
        %Message size distribution (SpectrumBridge)
        list_dir=dir;
        [rowb,colb]=size({list_dir.bytes});
        spectrumbridge_resp_size=[];
        for x=4:colb
            spectrumbridge_resp_size=[spectrumbridge_resp_size,list_dir(x).bytes];
        end
        %system('rm *');
        
    end
    
    %%
    %Plot figure

    if spectrumbridge_test==1
        %figure('Position',[440 378 560 420/2]);
        [fs,xs]=ksdensity(spectrumbridge_resp_size,'support','positive');
        fs=fs./sum(fs);
        plot(xs,fs,'k-.');
        grid on;
        box on;
        hold on;
        set(gca,'FontSize',ftsz);
        xlabel('Message size (bytes)','FontSize',ftsz);
        ylabel('Probability','FontSize',ftsz);
    end
    %Add common legend
    legend(legend_string);
    
    %%
    %Calculate statistics of message sizes for each WSDB
    
    %Mean
    mean_spectrumbridge_resp_size=mean(spectrumbridge_resp_size)

    %Variance
    var_spectrumbridge_resp_size=var(spectrumbridge_resp_size)

end
if response_error_calculation==1
    
    %Location of start and finish query
    %Query start location
    WSDB_data{1}.name='LA'; %Los Aneles, CA, USA (Wilshire Blvd 1) [downtown]
    WSDB_data{1}.latitude='34.047955';
    WSDB_data{1}.longitude='-118.256013';

    
    %Query finish location
    WSDB_data{2}.name='CB'; %Carolina Beach, NC, USA [ocean coast]
    WSDB_data{2}.latitude='34.047955';
    WSDB_data{2}.longitude='-77.885639';

    
    number_queries=100;
    number_batches=20;
    
    %Initialize error counter vectors
    error_spectrumbridge_vec=[];
    
    %Initialize Google API request counter [important: it needs initliazed
    %manually every time as limit of 1e3 queries per API is enforced. Check
    %your Google API console to check how many queries are used already]
    ggl_cnt=0;
    
    for bb=1:number_batches
        %Initialize error counters

        error_spectrumbridge=0;
        %Initialize request number counter
        in=0;
        for xx=1:number_queries
            in=in+1;
            fprintf('[Batch no., Query no.]: %d, %d\n',bb,xx)
            
            %Fetch location data
            latitude=WSDB_data{1}.latitude;
            %Generate random longitude for one query
            a=str2num(WSDB_data{1}.longitude);
            b=str2num(WSDB_data{2}.longitude);
            longitude=num2str((b-a)*rand+a);
            
            instant_clock=clock; %Save clock for file name (if both WSDBs are queried)

            if spectrumbridge_test==1
                %Query SpectrumBridge
                instant_clock=clock; %Start clock again if scanning only one database
                cd([my_path,'/spectrumbridge']);
                delay_spectrumbridge_tmp_r=0;
                if DeviceType=='8'
                    [msg_spectrumbridge,delay_spectrumbridge_tmp_r]=database_connect_spectrumbridge_register(...
                        AntennaHeight,DeviceType,Latitude,Longitude,[my_path,'/spectrumbridge']);
                end
                delay_spectrumbridge_tmp=delay_spectrumbridge_tmp+delay_spectrumbridge_tmp_r;
                [msg_spectrumbridge,delay_spectrumbridge_tmp,error_spectrumbridge_tmp]=...
                    database_connect_spectrumbridge(DeviceType,latitude,longitude);
                if error_spectrumbridge_tmp==1
                    error_spectrumbridge=error_spectrumbridge+1;
                end
            end

        if spectrumbridge_test==1
            %Clear old query results
            cd([my_path,'/spectrumbridge']);
            error_spectrumbridge_vec=[error_spectrumbridge_vec,error_spectrumbridge/number_queries];
        end
    end

    if spectrumbridge_test==1
        er_spectrumbridge=mean(error_spectrumbridge_vec)*100
        var_spectrumbridge=var(error_spectrumbridge_vec)*100
    end
end
end 
if delay_distribution_per_location==1
    
    no_queries=50; %Select how many queries per location
    
    %Location data
    WSDB_data{1}.name='LA'; %Los Aneles, CA, USA (Wilshire Blvd 1) [downtown]
    WSDB_data{1}.latitude='34.047955';
    WSDB_data{1}.longitude='-118.256013';

    
    WSDB_data{2}.name='WV'; %West Village (Manhattan), NY, USA [urban canyon]
    WSDB_data{2}.latitude='40.729655';
    WSDB_data{2}.longitude='-74.002854';

    
    WSDB_data{3}.name='SC'; %Scipio, OH, USA [flatland]
    WSDB_data{3}.latitude='41.102884';
    WSDB_data{3}.longitude='-82.957361';

    
    WSDB_data{4}.name='LE'; %Cleveland (Lake Erie), USA [lake coast]
    WSDB_data{4}.latitude='41.575416';
    WSDB_data{4}.longitude='-81.585442';

    
    WSDB_data{5}.name='CB'; %Carolina Beach, NC, USA [ocean coast]
    WSDB_data{5}.latitude='34.047955';
    WSDB_data{5}.longitude='-77.885639';

    
    [wsbx,wsby]=size(WSDB_data); %Get location data size
  
    delay_spectrumbridge_vector=[];

    legend_label_spectrumbridge=[];
    
    %Initialize Google API request counter [important: it needs initliazed
    %manually every time as limit of 1e3 queries per API is enforced. Check
    %your Google API console to check how many queries are used already]
    ggl_cnt=250;
    
    for ln=1:wsby

        delay_spectrumbridge=[];
        for xx=1:no_queries
            fprintf('[Query no., Location no.]: %d, %d\n',xx,ln)
            
            %Fetch location data
            latitude=WSDB_data{ln}.latitude;
            longitude=WSDB_data{ln}.longitude;
            
            instant_clock=clock; %Save clock for file name (if both WSDBs are queried)

            if spectrumbridge_test==1
                %Query SpectrumBridge
                instant_clock=clock; %Start clock again if scanning only one database
                cd([my_path,'/spectrumbridge']);
                delay_spectrumbridge_tmp_r=0;
                if DeviceType=='8'
                    [msg_spectrumbridge,delay_spectrumbridge_tmp_r]=database_connect_spectrumbridge_register(...
                        AntennaHeight,DeviceType,Latitude,Longitude,[my_path,'/spectrumbridge']);
                end
                [msg_spectrumbridge,delay_spectrumbridge_tmp,error_spectrumbridge_tmp]=...
                    database_connect_spectrumbridge(DeviceType,latitude,longitude);
                var_name=(['spectrumbridge_',num2str(longitude),'_',datestr(instant_clock, 'DD_mmm_YYYY_HH_MM_SS')]);
                if error_spectrumbridge_tmp==0
                    dlmwrite([var_name,'.txt'],msg_spectrumbridge,'');
                    delay_spectrumbridge_tmp=delay_spectrumbridge_tmp+delay_spectrumbridge_tmp_r;
                    delay_spectrumbridge=[delay_spectrumbridge,delay_spectrumbridge_tmp];
                end
            end
        end

        if spectrumbridge_test==1
            %Clear old query results
            cd([my_path,'/spectrumbridge']);
            %system('rm *');
            
            %Save delay data per location
            WSDB_data{ln}.delay_spectrumbridge=delay_spectrumbridge;
            legend_label_spectrumbridge=[legend_label_spectrumbridge,...
                repmat(ln,1,length(delay_spectrumbridge))]; %Label items for boxplot
            delay_spectrumbridge_vector=[delay_spectrumbridge_vector,delay_spectrumbridge];
            labels_spectrumbridge(ln)={WSDB_data{ln}.name};
        end
    end
    
    %Query general web services for comparison
    delay_spectrumbridge_web=[];
    for xx=1:no_queries
        fprintf('Query no.: %d\n',xx)

        if spectrumbridge_test==1
            ds=connect_webserver(3);
            delay_spectrumbridge_web=[delay_spectrumbridge_web,ds];
        end
    end
    if spectrumbridge_test==1
        legend_label_spectrumbridge=[legend_label_spectrumbridge,...
            repmat(ln+1,1,length(delay_spectrumbridge_web))]; %Label items for boxplot
        delay_spectrumbridge_vector=[delay_spectrumbridge_vector,delay_spectrumbridge_web];
        labels_spectrumbridge(ln+1)={'[SB]'};
    end
    
    %%
    %Plot figure: Box plots for delay per location
    
    %Select maximum Y axis
    max_el=max(delay_spectrumbridge_vector(1:end));

    if spectrumbridge_test==1
        figure('Position',[440 378 560/2.5 420/2]);

        boxplot(delay_spectrumbridge_vector,legend_label_spectrumbridge,...
            'labels',labels_spectrumbridge,'symbol','k+','jitter',0,'notch','on',...
            'factorseparator',1);
        ylim([0 max_el]);
        set(gca,'FontSize',ftsz);
        ylabel('Response delay (sec)','FontSize',ftsz);
        set(findobj(gca,'Type','text'),'FontSize',ftsz); %Boxplot labels size
        %Move boxplot labels below to avoid overlap with x axis
        txt=findobj(gca,'Type','text');
        set(txt,'VerticalAlignment','Top');
    end
        
    %Plot figure: plot delay request PDF estimates per location
    Markers={'k-','r--','g.-','b-.','mx-','cv-'};
    
    %Reserve axex properties for all figures
    fm=[];
    xm=[];
    fs=[];
    xs=[];
    fg=[];
    xg=[];
    

    if spectrumbridge_test==1
        figure('Position',[440 378 560 420/3]);
        name_location_vector=[];
        for ln=1:wsby
            delay_spectrumbridge=WSDB_data{ln}.delay_spectrumbridge;
            
            %Outlier removal (SpectrumBridge delay)
            outliers_pos=abs(delay_spectrumbridge-median(delay_spectrumbridge))>3*std(delay_spectrumbridge);
            delay_spectrumbridge(outliers_pos)=[];
            
            [fs,xs]=ksdensity(delay_spectrumbridge,'support','positive');
            fs=fs./sum(fs);
            plot(xs,fs,Markers{ln});
            hold on;
            name_location=WSDB_data{ln}.name;
            name_location_vector=[name_location_vector,{name_location}];
        end
        %Add plot for general webservice
        
        %Outlier removal (SpectrumBridge delay)
        outliers_pos=abs(delay_spectrumbridge_web-median(delay_spectrumbridge_web))>3*std(delay_spectrumbridge_web);
        delay_spectrumbridge_web(outliers_pos)=[];
        
        name_location_vector=[name_location_vector,'[SB]'];
        
        [fs,xs]=ksdensity(delay_spectrumbridge_web,'support','positive');
        fs=fs./sum(fs);
        plot(xs,fs,Markers{wsby+1});
        
        box on;
        grid on;
        set(gca,'FontSize',ftsz);
        xlabel('Response delay (sec)','FontSize',ftsz);
        ylabel('Probability','FontSize',ftsz);
        legend(name_location_vector,'Location','Best');
    end
    
%Set y axis limit manually at the end of plot
ylim([0 max([fg,fm,fs])]);    
end

if delay_distribution_area==1
    
    %Location of start and finish query
    %Query start location
    WSDB_data{1}.name='LA'; %Los Aneles, CA, USA (Wilshire Blvd 1) [downtown]
    WSDB_data{1}.latitude='34.047955';
    WSDB_data{1}.longitude='-118.256013';

    %Query finish location
    WSDB_data{2}.name='CB'; %Carolina Beach, NC, USA [ocean coast]
    WSDB_data{2}.latitude='34.047955';
    WSDB_data{2}.longitude='-77.885639';

    
    longitude_start=str2num(WSDB_data{1}.longitude); %Start of the spectrum scanning trajectory
    longitude_end=str2num(WSDB_data{2}.longitude); %End of spectrum scanning trajectory
    
    longitude_interval=100;
    longitude_step=(longitude_end-longitude_start)/longitude_interval;
    no_queries=20; %Number of queries per individual location

    delay_spectrumbridge=[];
    
    inx=0; %Initialize position counter
    
    %Initialize Google API request counter [important: it needs initliazed
    %manually every time as limit of 1e3 queries per API is enforced. Check
    %your Google API console to check how many queries are used already]
    ggl_cnt=2029;
    
    for xx=longitude_start:longitude_step:longitude_end
        inx=inx+1;
        iny=0; %Initialize query counter
        for yy=1:no_queries
            iny=iny+1;
            fprintf('[Query no., Location no.]: %d, %d\n',iny,inx);
            
            %Fetch location data
            latitude=WSDB_data{1}.latitude;
            longitude=num2str(xx);
            
            instant_clock=clock; %Save clock for file name (if both WSDBs are queried)


            if spectrumbridge_test==1
                %Query SpectrumBridge
                fprintf('SpectrumBridge\n')
                instant_clock=clock; %Start clock again if scanning only one database
                cd([my_path,'/spectrumbridge']);
                delay_spectrumbridge_tmp_r=0;
                if DeviceType=='8'
                    [msg_spectrumbridge,delay_spectrumbridge_tmp_r]=database_connect_spectrumbridge_register(...
                        AntennaHeight,DeviceType,Latitude,Longitude,[my_path,'/spectrumbridge']);
                end
                [msg_spectrumbridge,delay_spectrumbridge_tmp,error_spectrumbridge_tmp]=...
                    database_connect_spectrumbridge(DeviceType,latitude,longitude);
                var_name=(['spectrumbridge_',num2str(longitude),'_',datestr(instant_clock, 'DD_mmm_YYYY_HH_MM_SS')]);
                if error_spectrumbridge_tmp==0
                    dlmwrite([var_name,'.txt'],msg_spectrumbridge,'');
                    delay_spectrumbridge_tmp=delay_spectrumbridge_tmp+delay_spectrumbridge_tmp_r;
                    delay_spectrumbridge=[delay_spectrumbridge,delay_spectrumbridge_tmp];
                end
            end
        end
        %%
        %Assign delay per location per WSDB to a new variable
        if spectrumbridge_test==1
            delay_spectrumbridge_loc{inx}=delay_spectrumbridge;
            delay_spectrumbridge=[];
        end
    end
    
    %%
    %Get elavation data
    Elev=[];
    for xx=longitude_start:longitude_step:longitude_end
        pause(0.5); %Google imposes cap on number of queries - delay query
        elevation=elevation_google_maps(str2num(latitude),xx);
        Elev=[Elev,elevation];
    end
    
    %%
    %Compute means of queries per location

    if spectrumbridge_test==1
        Vm_spectrumbridge=[];
        for xx=1:inx
            mtmp_spectrumbridge=delay_spectrumbridge_loc{xx};
            Vm_spectrumbridge=[Vm_spectrumbridge,mean(mtmp_spectrumbridge)];
        end
        %Clear old query results
        cd([my_path,'/spectrumbridge']);
        %system('rm *');
    end
    
    %%
    %Plot distribution curves
    Markers={'g-','b--','k-.'};
    %Plot figures
    if spectrumbridge_test==1
        %figure('Position',[440 378 560 420/3]);
        [fs,xs]=ksdensity(Vm_spectrumbridge,'support','positive');
        fs=fs./sum(fs);
        plot(xs,fs,Markers{3});
        hold on;
    end
    
    box on;
    grid on;
    set(gca,'FontSize',ftsz);
    xlabel('Response delay (sec)','FontSize',ftsz);
    ylabel('Probability','FontSize',ftsz);
    legend(legend_string,'Location','Best');
    
    %Plot delay per location curves
    if spectrumbridge_test==1
        %figure('Position',[440 378 560 420/2]);
        hold on;
        plot(1:longitude_interval+1,Vm_spectrumbridge./sum(Vm_spectrumbridge),Markers{3});
    end
    
    %Plot elevation
    %plot(Elev./sum(Elev),Markers{4});
    
    box on;
    grid on;
    set(gca,'FontSize',ftsz);
    xlim([0 longitude_interval+1]);
    xlabel('Location number','FontSize',ftsz);
    ylabel('Response delay (sec)','FontSize',ftsz);
    legend([legend_string],'Location','Best');
    %legend([legend_string,'Normalized elevation'],'Location','Best');
    
end

%%
['Elapsed time: ',num2str(toc/60),' min']