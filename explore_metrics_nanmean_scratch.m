%
% A script to explore Robin's data of graph metrics against AV binding
% In AD, PSP and Controls
% By Thomas Cope, July 2016
% Specify the group of interest with 'group' (controls, PSP, ADMCI or leave blank for all)


%function explore_metrics(plottype,group)

clear all

%% Read in data

% filteringtype = 'Filtered'; % Filtered or Unfiltered
% metrictype = 'Binary'; % Weighted or Binary
filteringtype = 'Relative'; % Filtered or Unfiltered
metrictype = 'Binary'; % Weighted or Binary
inputtype = 'txt'; % txt or csv (csv must be numeric only in an ROIxsubject table)

if ~exist('group','var') || strcmp(group,'all')==1
    groupstodraw = [1,2,3];
elseif strcmp(group,'controls')==1
    groupstodraw = 1;
elseif strcmp(group,'ADMCI')==1
    groupstodraw = 2;
elseif strcmp(group,'PSP')==1
    groupstodraw = 3;
else
    error('That group does not exist, try again!')
end

if strcmp(inputtype, 'csv')
    Metrics = ls([pwd '\' filteringtype '_' metrictype '_Metrics']); % List the contents of the data directory
    
    for filenumber = 3:size(Metrics,1)
        Thisfilepath = [pwd '\' filteringtype '_' metrictype '_Metrics\' deblank(Metrics(filenumber,:))]; % Get the filepath
        filename_deblanked = deblank(Metrics(filenumber,:)); %Get the filename without trailing blanks
        if strncmpi(Metrics(filenumber,:),'PSP',3) %Check if it's from one of these groups
            datatype = filename_deblanked(5:end-4); %Define what data is there
            eval(['PSPdata.' datatype '= xlsread(Thisfilepath);']) %Note that xlsread excludes the header row
        elseif strncmpi(Metrics(filenumber,:),'ADMCI',5)
            datatype = filename_deblanked(7:end-4);
            eval(['ADMCIdata.' datatype '= xlsread(Thisfilepath);'])
        elseif strncmpi(Metrics(filenumber,:),'Controls',8)
            datatype = filename_deblanked(10:end-4);
            eval(['controldata.' datatype '= xlsread(Thisfilepath);'])
        elseif strncmpi(Metrics(filenumber,:),'Atlas',5)
            datatype = 'Atlas';
            [~,atlasdata] = xlsread(Thisfilepath);
        end
    end
    
    
    %% Define what data we have
    
    controlfields = fieldnames(controldata);
    PSPfields = fieldnames(PSPdata);
    ADMCIfields = fieldnames(ADMCIdata);
    
    if isequal(controlfields,PSPfields,ADMCIfields)
        %     fprintf(['The fields are all the same, they are: \n\n'])
        %     controlfields
        %     fprintf(['\n\n'])
        %     input('Is this OK? Press any key to continue, or ctrl+c to terminate');
    else
        fprintf(['The fields are not all the same. For the controls: \n\n'])
        controlfields
        fprintf(['for the PSP:\n\n'])
        PSPfields
        fprintf(['for the ADMCI\n\n'])
        ADMCIfields
        fprintf(['\n\n'])
        input('Is this OK? Press any key to continue, or ctrl+c to terminate');
    end
    
    
    %% Now, re-extract tables from raw text files (test to avoid transcription errors)
    Metrics = ls([pwd '\results_' filteringtype '_' metrictype]);
    for filenumber = 1:size(Metrics,1)
        Thisfilepath = [pwd '\results_' filteringtype '_' metrictype '\' deblank(Metrics(filenumber,:))]; % Get the filepath
        filename_deblanked = deblank(Metrics(filenumber,:)); %Get the filename without trailing blanks
        
        try
            if strcmp(filename_deblanked(end-5:end-4),'05')
                network_threshold = 5;
            elseif strcmp(filename_deblanked(end-5:end-4),'.1')
                network_threshold = 10;
            elseif strcmp(filename_deblanked(end-5:end-4),'15')
                network_threshold = 15;
            elseif strcmp(filename_deblanked(end-5:end-4),'.2')
                network_threshold = 20;
            end
            if strcmp(filename_deblanked(end-3:end),'.txt')
                testdata = readtable(Thisfilepath,'FileType','text','Delimiter','tab','HeaderLines',0,'ReadRowNames',1,'ReadVariableNames',1);
                datatype = filename_deblanked(1:end-4);
                eval(['ADMCIdata_test.' datatype(1:5) '_' num2str(network_threshold) ' = testdata{1:size(ADMCIdata.AV_binding,2),:}'';']);
                eval(['controldata_test.' datatype(1:5) '_' num2str(network_threshold) ' = testdata{size(ADMCIdata.AV_binding,2)+1:size(ADMCIdata.AV_binding,2)+size(controldata.AV_binding,2),:}'';']);
                eval(['PSPdata_test.' datatype(1:5) '_' num2str(network_threshold) ' = testdata{size(ADMCIdata.AV_binding,2)+size(controldata.AV_binding,2)+1:end,:}'';']);
            end
        end
    end
    
    for density_threshold = 5:5:20; %density_threshold = 5:5:15;
        for y_variable = {'LocalEfficiency','BetCent','PartCoeff','Degree'};
            if strcmp(y_variable,'LocalEfficiency')
                if isequal(eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['ADMCIdata_test.effic_' num2str(density_threshold)]));
                    disp(['for ADMCI, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for ADMCI, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]), eval(['controldata_test.effic_' num2str(density_threshold)]));
                    disp(['for control, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for control, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['PSPdata_test.effic_' num2str(density_threshold)]));
                    disp(['for PSP, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for PSP, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
            elseif strcmp(y_variable,'BetCent')
                if isequal(eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['ADMCIdata_test.betwe_' num2str(density_threshold)]));
                    disp(['for ADMCI, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for ADMCI, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]), eval(['controldata_test.betwe_' num2str(density_threshold)]));
                    disp(['for control, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for control, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['PSPdata_test.betwe_' num2str(density_threshold)]));
                    disp(['for PSP, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for PSP, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
            elseif strcmp(y_variable,'PartCoeff')
                if isequal(eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['ADMCIdata_test.parti_' num2str(density_threshold)]));
                    disp(['for ADMCI, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for ADMCI, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]), eval(['controldata_test.parti_' num2str(density_threshold)]));
                    disp(['for control, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for control, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['PSPdata_test.parti_' num2str(density_threshold)]));
                    disp(['for PSP, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for PSP, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
            elseif strcmp(y_variable,'Degree')
                %             if strcmp(filteringtype, 'Filtered') && strcmp(metrictype, 'Weighted');  %GraphVar outputted an array of 1s for all others
                eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold) '=ADMCIdata_test.degre_' num2str(density_threshold) ';'])
                eval(['controldata.' y_variable{1} '_' num2str(density_threshold) '=controldata_test.degre_' num2str(density_threshold)])
                eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold) '=PSPdata_test.degre_' num2str(density_threshold)])
                %             end
            end
            
            
        end
    end
    input(['Do you want to continue (check whether data match) - ctrl+c for no, anything else for yes']);
    
else
    ADMCIdata.AV_binding = xlsread('ADMCI_AV_binding.csv');
    controldata.AV_binding = xlsread('Controls_AV_binding.csv');
    PSPdata.AV_binding = xlsread('PSP_AV_binding.csv');
    [~,atlasdata] = xlsread('Atlas_ROIs.csv');
    
    Metrics = ls([pwd '\relative_' metrictype '_thresholds']);
    
    for filenumber = 3:size(Metrics,1)
        Thisfilepath = [pwd '\relative_' metrictype '_thresholds\' deblank(Metrics(filenumber,:))]; % Get the filepath
        filename_deblanked = deblank(Metrics(filenumber,:)); %Get the filename without trailing blanks
        
        try
            if strcmp(filename_deblanked(end-5:end-4),'05')
                network_threshold = 5;
            elseif strcmp(filename_deblanked(end-5:end-4),'.1')
                network_threshold = 10;
            elseif strcmp(filename_deblanked(end-5:end-4),'15')
                network_threshold = 15;
            elseif strcmp(filename_deblanked(end-5:end-4),'.2')
                network_threshold = 20;
            end
            if strcmp(filename_deblanked(end-3:end),'.txt')
                testdata = readtable(Thisfilepath,'FileType','text','Delimiter','tab','HeaderLines',0,'ReadRowNames',1,'ReadVariableNames',1);
                datatype = filename_deblanked(1:end-4);
                eval(['ADMCIdata_test.' datatype(1:5) '_' num2str(network_threshold) ' = testdata{1:size(ADMCIdata.AV_binding,2),:}'';']);
                eval(['controldata_test.' datatype(1:5) '_' num2str(network_threshold) ' = testdata{size(ADMCIdata.AV_binding,2)+1:size(ADMCIdata.AV_binding,2)+size(controldata.AV_binding,2),:}'';']);
                eval(['PSPdata_test.' datatype(1:5) '_' num2str(network_threshold) ' = testdata{size(ADMCIdata.AV_binding,2)+size(controldata.AV_binding,2)+1:end,:}'';']);
            end
        end
    end
    
    for density_threshold = 5:5:20; %density_threshold = 5:5:15;
        try
        for y_variable = {'LocalEfficiency'}; %y_variable = {'LocalEfficiency','BetCent','PartCoeff','Degree'};
            if strcmp(y_variable,'LocalEfficiency')
                eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold) '=ADMCIdata_test.effic_' num2str(density_threshold) ';'])
                eval(['controldata.' y_variable{1} '_' num2str(density_threshold) '=controldata_test.effic_' num2str(density_threshold) ';'])
                eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold) '=PSPdata_test.effic_' num2str(density_threshold) ';'])
            elseif strcmp(y_variable,'BetCent')
                if isequal(eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['ADMCIdata_test.betwe_' num2str(density_threshold)]));
                    disp(['for ADMCI, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for ADMCI, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]), eval(['controldata_test.betwe_' num2str(density_threshold)]));
                    disp(['for control, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for control, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['PSPdata_test.betwe_' num2str(density_threshold)]));
                    disp(['for PSP, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for PSP, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
            elseif strcmp(y_variable,'PartCoeff')
                if isequal(eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['ADMCIdata_test.parti_' num2str(density_threshold)]));
                    disp(['for ADMCI, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for ADMCI, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]), eval(['controldata_test.parti_' num2str(density_threshold)]));
                    disp(['for control, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for control, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
                if isequal(eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]), eval(['PSPdata_test.parti_' num2str(density_threshold)]));
                    disp(['for PSP, data match for ' y_variable{1} '_' num2str(density_threshold)]);
                else
                    disp(['for PSP, DATA DO NOT MATCH FOR ' y_variable{1} '_' num2str(density_threshold)]);
                end
            elseif strcmp(y_variable,'Degree')
                %             if strcmp(filteringtype, 'Filtered') && strcmp(metrictype, 'Weighted');  %GraphVar outputted an array of 1s for all others
                eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold) '=ADMCIdata_test.degre_' num2str(density_threshold) ';'])
                eval(['controldata.' y_variable{1} '_' num2str(density_threshold) '=controldata_test.degre_' num2str(density_threshold) ';'])
                eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold) '=PSPdata_test.degre_' num2str(density_threshold) ';'])
                %             end
            end
            
            
        end
        end
    end
    
end

%% First draw a simple graph of each metric by group, ordered by region number.

if ~exist([pwd '\' filteringtype '_' metrictype '_figures\'] ,'dir')
    mkdir([pwd '\' filteringtype '_' metrictype '_figures\']);
end

for density_threshold = 5:5:20; %density_threshold = 5:5:15;
    try
    for y_variable = {'LocalEfficiency'}; %y_variable = {'LocalEfficiency','BetCent','PartCoeff','Degree'};
        
        h = figure;
        hold on
        
        for i = 1:length(groupstodraw)
            if groupstodraw(i) == 1
                y_data = eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]);
                h1 = plot(1:length(y_data),nanmean(y_data'),'k.');
            elseif groupstodraw(i) == 2
                y_data = eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]);
                h2 = plot(1:length(y_data),nanmean(y_data'),'r.');
            elseif groupstodraw(i) == 3
                y_data = eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]);
                h3 = plot(1:length(y_data),nanmean(y_data'),'g.');
            end
        end
        
        legend([h1, h2, h3], 'Controls', 'ADMCI', 'PSP');
        titletext = [y_variable{1} ' against region number at a density threshold of ' num2str(density_threshold)];
        
        xlabel('Region Number', 'Interpreter', 'none')
        ylabel(y_variable{1}, 'Interpreter', 'none')
        
        title(titletext, 'Interpreter', 'none')
        
        %now fill by brain region
        cols = hsv(numel(unique(atlasdata(:,2))));
        region_names = unique(atlasdata(:,2),'stable');
        
        region_order = atlasdata(:,2);
        for j = 1:length(region_names);
            region_order = strrep(region_order,region_names(j),num2str(j));
        end
        
        colourkey = cols(str2num(char(region_order))',:);
        
        currentlimits = ylim;
        ylim manual
        
        for i = 1:size(colourkey,1)
            thisfill = fill([i-0.5,i-0.5,i+0.5,i+0.5],[currentlimits,fliplr(currentlimits)],colourkey(i,:));
            set(thisfill,'facealpha',.2,'LineStyle','none')
        end
        
        xlim([0,length(atlasdata(:,2))+(length(atlasdata(:,2))/5)]);
        colourorder = unique(colourkey,'stable','rows');
        ydivided = range(ylim)/(length(unique(colourkey,'stable','rows'))+2);
        
        for i = 2:size(colourorder,1)+1
            thisfill = fill([length(atlasdata(:,2))+(length(atlasdata(:,2))/20), length(atlasdata(:,2))+(length(atlasdata(:,2))/20), length(atlasdata(:,2))+(length(atlasdata(:,2))/10), length(atlasdata(:,2))+(length(atlasdata(:,2))/10)],[currentlimits(2)-(i*ydivided), currentlimits(2)-((i+1)*ydivided), currentlimits(2)-((i+1)*ydivided), currentlimits(2)-(i*ydivided)], colourorder(i-1,:));
            set(thisfill,'facealpha',.2,'LineStyle','none')
            text(length(atlasdata(:,2))+(3*length(atlasdata(:,2))/20),currentlimits(2)-((i+0.5)*ydivided), region_names(i-1));
        end
        
        saveas(h,[filteringtype '_' metrictype '_figures\' titletext '.jpg']);
        close(h)
        
    end
    end
end

for y_variable = {'AV_binding'};
    
    h = figure;
    hold on
    
    for i = 1:length(groupstodraw)
        if groupstodraw(i) == 1
            y_data = eval(['controldata.' y_variable{1}]);
            h1 = plot(1:length(y_data),nanmean(y_data'),'k.');
        elseif groupstodraw(i) == 2
            y_data = eval(['ADMCIdata.' y_variable{1}]);
            h2 = plot(1:length(y_data),nanmean(y_data'),'r.');
        elseif groupstodraw(i) == 3
            y_data = eval(['PSPdata.' y_variable{1}]);
            h3 = plot(1:length(y_data),nanmean(y_data'),'g.');
        end
    end
    legend([h1, h2, h3], 'Controls', 'ADMCI', 'PSP');
    titletext = [y_variable{1} ' against region number'];
    
    xlabel('Region Number', 'Interpreter', 'none')
    ylabel(y_variable{1}, 'Interpreter', 'none')
    
    title(titletext, 'Interpreter', 'none')
    
    %now fill by brain region
    cols = hsv(numel(unique(atlasdata(:,2))));
    region_names = unique(atlasdata(:,2),'stable');
    
    region_order = atlasdata(:,2);
    for j = 1:length(region_names);
        region_order = strrep(region_order,region_names(j),num2str(j));
    end
    
    colourkey = cols(str2num(char(region_order))',:);
    currentlimits = ylim;
    ylim manual
    
    for i = 1:size(colourkey,1)
        thisfill = fill([i-0.5,i-0.5,i+0.5,i+0.5],[currentlimits,fliplr(currentlimits)],colourkey(i,:));
        set(thisfill,'facealpha',.2,'LineStyle','none')
    end
    
    xlim([0,length(atlasdata(:,2))+(length(atlasdata(:,2))/5)]);
    colourorder = unique(colourkey,'stable','rows');
    ydivided = range(ylim)/(length(unique(colourkey,'stable','rows'))+2);
    
    for i = 2:size(colourorder,1)+1
        thisfill = fill([length(atlasdata(:,2))+(length(atlasdata(:,2))/20), length(atlasdata(:,2))+(length(atlasdata(:,2))/20), length(atlasdata(:,2))+(length(atlasdata(:,2))/10), length(atlasdata(:,2))+(length(atlasdata(:,2))/10)],[currentlimits(2)-(i*ydivided), currentlimits(2)-((i+1)*ydivided), currentlimits(2)-((i+1)*ydivided), currentlimits(2)-(i*ydivided)], colourorder(i-1,:));
        set(thisfill,'facealpha',.2,'LineStyle','none')
        text(length(atlasdata(:,2))+(3*length(atlasdata(:,2))/20),currentlimits(2)-((i+0.5)*ydivided), region_names(i-1));
    end
    
    saveas(h,[filteringtype '_' metrictype '_figures\' titletext '.jpg']);
    close(h)

end


%% Next draw a simple graph of specified variable against AV for specified groups (controls = 1, ADMCI = 2, PSP = 3)

if ~exist([pwd '\' filteringtype '_' metrictype '_figures\'] ,'dir')
    mkdir([pwd '\' filteringtype '_' metrictype '_figures\']);
end

% Set what you want to plot. Assumes always AV on x axis, but y variable can be changed
for density_threshold = 5:5:20; %density_threshold = 5:5:15;
    try
    for y_variable = {'LocalEfficiency'}; %y_variable = {'LocalEfficiency','BetCent','PartCoeff','Degree'};
        x_variable = 'AV_binding';
        
        %First for all groups together
%         ncontrols = size(controldata.LocalEfficiency_10,2);
%         nADMCI = size(ADMCIdata.LocalEfficiency_10,2);
%         nPSP = size(PSPdata.LocalEfficiency_10,2);
        
%         x_data = eval(['[controldata.' x_variable ',ADMCIdata.' x_variable ',PSPdata.' x_variable ']']);
%         y_data = eval(['[controldata.' y_variable{1} '_' num2str(density_threshold) ',ADMCIdata.' y_variable{1} '_' num2str(density_threshold) ',PSPdata.' y_variable{1} '_' num2str(density_threshold) ']']);
%         groups = [ones(1,ncontrols), 2*ones(1,nADMCI), 3*ones(1,nPSP)];

        this_control_x_data = eval(['nanmean(controldata.' x_variable ',2)']);
        this_ADMCI_x_data = eval(['nanmean(ADMCIdata.' x_variable ',2)']);
        this_PSP_x_data = eval(['nanmean(PSPdata.' x_variable ',2)']);
        this_x_data = [this_control_x_data;this_ADMCI_x_data;this_PSP_x_data];
        this_control_y_data = eval(['nanmean(controldata.' y_variable{1} '_' num2str(density_threshold) ',2)']);
        this_ADMCI_y_data = eval(['nanmean(ADMCIdata.' y_variable{1} '_' num2str(density_threshold) ',2)']);
        this_PSP_y_data = eval(['nanmean(PSPdata.' y_variable{1} '_' num2str(density_threshold) ',2)']);
        this_y_data = [this_control_y_data;this_ADMCI_y_data;this_PSP_y_data];
        titletext = ['All subjects ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
 
        groups = [repmat({'Control'},1,length(this_control_x_data)), repmat({'ADMCI'},1,length(this_ADMCI_x_data)), repmat({'PSP'},1,length(this_PSP_x_data))]';

        h = figure;
        gscatter(this_x_data',this_y_data',groups,[0, 0, 0 ; 1, 0, 1 ; 0, 1, 0])
%         regres_coeffs = [this_x_data,ones(size(this_x_data))] \ this_y_data;
        regres_coeffs = [this_control_x_data,ones(size(this_control_x_data))] \ this_control_y_data;
        xs = min(this_x_data):(range(this_x_data)/100):max(this_x_data);
        hold on       
        plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)), 'k');
        regres_coeffs = [this_ADMCI_x_data,ones(size(this_ADMCI_x_data))] \ this_ADMCI_y_data;
        plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)), 'm');
        regres_coeffs = [this_PSP_x_data,ones(size(this_PSP_x_data))] \ this_PSP_y_data;
        plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)), 'g');
        
        xlabel(x_variable, 'Interpreter', 'none')
        ylabel(y_variable{1}, 'Interpreter', 'none')
        title(titletext, 'Interpreter', 'none')
        saveas(h,[filteringtype '_' metrictype '_figures\' titletext '.jpg']);
        close(h)
        
        
        for i = 1:length(groupstodraw)
            if groupstodraw(i) == 1
                x_data = eval(['controldata.' x_variable]);
                y_data = eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]);
                titletext = ['Controls ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
            elseif groupstodraw(i) == 2
                x_data = eval(['ADMCIdata.' x_variable]);
                y_data = eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]);
                titletext = ['ADMCI ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
            elseif groupstodraw(i) == 3
                x_data = eval(['PSPdata.' x_variable]);
                y_data = eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]);
                titletext = ['PSP ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
            end
            
            %Simple scatter plot
            h = figure;
            gscatter(nanmean(x_data'),nanmean(y_data'),atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))))
            regres_coeffs = [nanmean(x_data')',ones(size(nanmean(x_data')'))] \ nanmean(y_data')';
            hold on
            xs = min(nanmean(x_data')):(range(nanmean(x_data'))/100):max(nanmean(x_data'));
            plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)));
            xlabel(x_variable, 'Interpreter', 'none')
            ylabel(y_variable{1}, 'Interpreter', 'none')
            title(titletext, 'Interpreter', 'none')
            saveas(h,[filteringtype '_' metrictype '_figures\' titletext '.jpg']);
            close(h)
            
            %Multigradient plot
            h = figure;
            hold on
            all_gradients = [];
            cols = hsv(numel(unique(atlasdata(:,2))));
            regionorder = unique(atlasdata(:,2),'stable');
            all_gradients = zeros(1,size(y_data,1));
            R_gradients = zeros(1,size(y_data,1));
            P_gradients = zeros(1,size(y_data,1));
            for j = 1:size(y_data,1)
                regres_coeffs = [x_data(j,:)',ones(size(x_data(j,:)'))] \ y_data(j,:)';
                xs = min(x_data(j,:)):(range((x_data(j,:)))/100):max(x_data(j,:));
                plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)),'Color',cols(strcmp(atlasdata(j,2),regionorder),:));
                all_gradients(j) = regres_coeffs(1);
                [R,P] = corrcoef(x_data(j,:),y_data(j,:));
                R_gradients(j) = R(2,1);
                P_gradients(j) = P(2,1);
            end
            xlabel(x_variable, 'Interpreter', 'none')
            ylabel(y_variable{1}, 'Interpreter', 'none')
            title(['Multigradient plot for ' titletext], 'Interpreter', 'none')
            saveas(h,[filteringtype '_' metrictype '_figures\Multigradient_' titletext '.jpg']);
            close(h)
            
            %XXX TEST - plot only statistically significant gradients
            h = figure;
            hold on
            for j = 1:size(y_data,1)
                regres_coeffs = [x_data(j,:)',ones(size(x_data(j,:)'))] \ y_data(j,:)';
                xs = min(x_data(j,:)):(range((x_data(j,:)))/100):max(x_data(j,:));
                if P_gradients(j) < 0.05
                    plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)),'Color',cols(strcmp(atlasdata(j,2),regionorder),:));
                end
            end
            xlabel(x_variable, 'Interpreter', 'none')
            ylabel(y_variable{1}, 'Interpreter', 'none')
            title(['Significant multigradient plot for ' titletext], 'Interpreter', 'none')
            saveas(h,[filteringtype '_' metrictype '_figures\Significant_multigradient_' titletext '.jpg']);
            close(h)
            
            if groupstodraw(i) == 1
                controlgradients = all_gradients;
            elseif groupstodraw(i) == 2
                ADMCIgradients = all_gradients;
            elseif groupstodraw(i) == 3
                PSPgradients = all_gradients;
            end
            
            % Now plot a difference graph for each group
            if groupstodraw(i) == 2 || groupstodraw(i) == 3
                control_x_data = eval(['controldata.' x_variable]);
                control_y_data = eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]);
                if groupstodraw(i) == 2
                    x_data = eval(['ADMCIdata.' x_variable]);
                    y_data = eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]);
                    titletext = ['ADMCI normalised ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
                elseif groupstodraw(i) == 3
                    x_data = eval(['PSPdata.' x_variable]);
                    y_data = eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]);
                    titletext = ['PSP normalised ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
                end
                
                %Simple scatter plot
                h = figure;
                gscatter(nanmean(x_data')-nanmean(control_x_data'),nanmean(y_data')-nanmean(control_y_data'),atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))))
                regres_coeffs = [(nanmean(x_data')-nanmean(control_x_data'))',ones(size((nanmean(x_data')-nanmean(control_x_data'))'))] \ (nanmean(y_data')-nanmean(control_y_data'))';
                hold on
                xs = min(nanmean(x_data')-nanmean(control_x_data')):(range(nanmean(x_data')-nanmean(control_x_data'))/100):max(nanmean(x_data')-nanmean(control_x_data'));
                plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)));
                xlabel([x_variable ' difference'], 'Interpreter', 'none')
                ylabel([y_variable{1} ' difference'], 'Interpreter', 'none')
                title(titletext, 'Interpreter', 'none')
                saveas(h,[filteringtype '_' metrictype '_figures\' titletext '.jpg']);
                close(h)
                
            end
            
        end
        
        %Now plot the gradients and work out if significant
        if ~exist('group','var') || strcmp(group,'all')==1
            
            %First ordered by brain region
            [controlp,controlh] = signtest(controlgradients);
            [ADMCIp,ADMCIh] = signtest(ADMCIgradients);
            [PSPp,PSPh] = signtest(PSPgradients);
            
            xs = rand(length(controlgradients),1);
            xs = (xs-0.5)/16;
                     
            regionorder = unique(atlasdata(:,2),'stable');
            for j = 1:length(xs)
               xs(j) = xs(j)-0.33+((find(strcmp(atlasdata(j,2),regionorder))-1)/(1.5*(length(regionorder)-1)));                
            end

            h=figure;
            hold on
            gscatter(1+xs,controlgradients,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if controlh == 1
                %scatter(xs,controlgradients,'r')
                if nanmean(controlgradients)>0
                    plot(1,max([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                else
                    plot(1,min([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                end
            else
                %scatter(xs,controlgradients,'k')
            end
            gscatter(2+xs,ADMCIgradients,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if ADMCIh == 1
                %scatter(2*xs,ADMCIgradients,'r')
                if nanmean(ADMCIgradients)>0
                    plot(2,max([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                else
                    plot(2,min([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                end
            else
                %scatter(2*xs,ADMCIgradients,'k')
            end
            gscatter(3+xs,PSPgradients,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if PSPh == 1
                %scatter(3*xs,PSPgradients,'r')
                if nanmean(PSPgradients)>0
                    plot(3,max([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                else
                    plot(3,min([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                end
            else
                %scatter(3*xs,PSPgradients,'k')
            end
            
            plot([0,4],[0,0],'k--');
            xlabel(['Group'])
            ax = gca;
            ax.XTick = [1,2,3];
            ax.XLim = [0,4];
            ax.XTickLabel = {'Control';'ADMCI';'PSP'};
            ylabel(['Gradient'])
            title(['Allgradients plot for ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)], 'Interpreter', 'none')
            saveas(h,[filteringtype '_' metrictype '_figures\Allgradients_' y_variable{1} '_against_' x_variable ' at a density threshold of ' num2str(density_threshold) '.jpg']);
            close(h)
            
            
            %Same thing, ordered by AV p-value
            ADMCIabnormalregions = zeros(1,length(controldata.AV_binding)); %Flag significantly abnormal regions cf control
            PSPabnormalregions = zeros(1,length(controldata.AV_binding));
            ADMCIpvalues = zeros(2,length(controldata.AV_binding));
            PSPpvalues = zeros(2,length(controldata.AV_binding));
            ADMCIpvalues(1,:) = 1:length(controldata.AV_binding);
            PSPpvalues(1,:) = 1:length(controldata.AV_binding);
            
            for region = 1:length(controldata.AV_binding)
                [ADMCIabnormalregions(region),ADMCIpvalues(2,region)] = ttest2(controldata.AV_binding(region,:),ADMCIdata.AV_binding(region,:),'vartype','unequal');
                [PSPabnormalregions(region),PSPpvalues(2,region)] = ttest2(controldata.AV_binding(region,:),PSPdata.AV_binding(region,:),'vartype','unequal');
            end
            ADMCIabnormalregions = logical(ADMCIabnormalregions);
            PSPabnormalregions = logical(PSPabnormalregions);
            
            %Now plot the gradients ordered by AV t-score
            [controlp,controlh] = signtest(controlgradients);
            [ADMCIp,ADMCIh] = signtest(ADMCIgradients);
            [PSPp,PSPh] = signtest(PSPgradients);
            
            [~,ADsigorder]=sort(ADMCIpvalues(2,:));
            [~,PSPsigorder]=sort(PSPpvalues(2,:));
            
            ADxs = -0.33+(ADsigorder-1)/(1.5*(max(ADsigorder)-1));
            PSPxs = -0.33+(PSPsigorder-1)/(1.5*(max(PSPsigorder)-1));
            
            regionorder = unique(atlasdata(:,2),'stable');
            
            h=figure;
            hold on
            gscatter(1+ADxs,controlgradients,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if controlh == 1
                %scatter(xs,controlgradients,'r')
                if nanmean(controlgradients)>0
                    plot(1,max([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                else
                    plot(1,min([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                end
            else
                %scatter(xs,controlgradients,'k')
            end
            gscatter(2+ADxs,ADMCIgradients,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if ADMCIh == 1
                %scatter(2*xs,ADMCIgradients,'r')
                if nanmean(ADMCIgradients)>0
                    plot(2,max([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                else
                    plot(2,min([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                end
            else
                %scatter(2*xs,ADMCIgradients,'k')
            end
            gscatter(3+PSPxs,controlgradients,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if controlh == 1
                %scatter(xs,controlgradients,'r')
                if nanmean(controlgradients)>0
                    plot(1,max([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                else
                    plot(1,min([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                end
            else
                %scatter(xs,controlgradients,'k')
            end
            gscatter(4+PSPxs,PSPgradients,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if PSPh == 1
                %scatter(3*xs,PSPgradients,'r')
                if nanmean(PSPgradients)>0
                    plot(3,max([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                else
                    plot(3,min([controlgradients,ADMCIgradients,PSPgradients]),'*k')
                end
            else
                %scatter(3*xs,PSPgradients,'k')
            end
            
            plot([0,5],[0,0],'k--');
            xlabel(['Group'])
            ax = gca;
            ax.XTick = [1,2,3,4];
            ax.XLim = [0,5];
            ax.XTickLabel = {'Control';'ADMCI';'Control';'PSP'};
            ylabel(['Gradient'])
            title(['Allgradients ordered by pvalue for ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)], 'Interpreter', 'none')
            saveas(h,[filteringtype '_' metrictype '_figures\Allgradients_pordered_' y_variable{1} '_against_' x_variable ' at a density threshold of ' num2str(density_threshold) '.jpg']);
            close(h)
            
            
            %Now again ordered by brain region, normalised to AV data range
            %within each group
            controlgradientsnormed = eval(['controlgradients * nanmean(range(controldata.' x_variable '''))']);
            ADMCIgradientsnormed = eval(['ADMCIgradients * nanmean(range(ADMCIdata.' x_variable '''))']);
            PSPgradientsnormed = eval(['PSPgradients * nanmean(range(PSPdata.' x_variable '''))']);
            
            [controlp,controlh] = signtest(controlgradientsnormed);
            [ADMCIp,ADMCIh] = signtest(ADMCIgradientsnormed);
            [PSPp,PSPh] = signtest(PSPgradientsnormed);
            
            xs = rand(length(controlgradientsnormed),1);
            xs = (xs-0.5)/16;
                     
            regionorder = unique(atlasdata(:,2),'stable');
            for j = 1:length(xs)
               xs(j) = xs(j)-0.33+((find(strcmp(atlasdata(j,2),regionorder))-1)/(1.5*(length(regionorder)-1)));                
            end

            h=figure;
            hold on
            gscatter(1+xs,controlgradientsnormed,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if controlh == 1
                %scatter(xs,controlgradientsnormed,'r')
                if nanmean(controlgradientsnormed)>0
                    plot(1,max([controlgradientsnormed,ADMCIgradientsnormed,PSPgradientsnormed]),'*k')
                else
                    plot(1,min([controlgradientsnormed,ADMCIgradientsnormed,PSPgradientsnormed]),'*k')
                end
            else
                %scatter(xs,controlgradientsnormed,'k')
            end
            gscatter(2+xs,ADMCIgradientsnormed,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if ADMCIh == 1
                %scatter(2*xs,ADMCIgradientsnormed,'r')
                if nanmean(ADMCIgradientsnormed)>0
                    plot(2,max([controlgradientsnormed,ADMCIgradientsnormed,PSPgradientsnormed]),'*k')
                else
                    plot(2,min([controlgradientsnormed,ADMCIgradientsnormed,PSPgradientsnormed]),'*k')
                end
            else
                %scatter(2*xs,ADMCIgradientsnormed,'k')
            end
            gscatter(3+xs,PSPgradientsnormed,atlasdata(:,2),hsv(numel(unique(atlasdata(:,2)))));
            if PSPh == 1
                %scatter(3*xs,PSPgradientsnormed,'r')
                if nanmean(PSPgradientsnormed)>0
                    plot(3,max([controlgradientsnormed,ADMCIgradientsnormed,PSPgradientsnormed]),'*k')
                else
                    plot(3,min([controlgradientsnormed,ADMCIgradientsnormed,PSPgradientsnormed]),'*k')
                end
            else
                %scatter(3*xs,PSPgradientsnormed,'k')
            end
            
            plot([0,4],[0,0],'k--');
            xlabel(['Group'])
            ax = gca;
            ax.XTick = [1,2,3];
            ax.XLim = [0,4];
            ax.XTickLabel = {'Control';'ADMCI';'PSP'};
            ylabel(['Gradient'])
            title(['Allgradients plot for ' y_variable{1} ' normalised against ' x_variable ' at a density threshold of ' num2str(density_threshold)], 'Interpreter', 'none')
            saveas(h,[filteringtype '_' metrictype '_figures\Allgradients_normalised_' y_variable{1} '_against_' x_variable ' at a density threshold of ' num2str(density_threshold) '.jpg']);
            close(h)
            
        end
    end
    end
end

%% Next step - analyse only abnormal areas.
for density_threshold = 5:5:20; %density_threshold = 5:5:15;
    try
    for y_variable = {'LocalEfficiency'}; %y_variable = {'LocalEfficiency','BetCent','PartCoeff','Degree'};
        % First detect abnormal areas by unpaired ttest vs controls with unequal variances
        ADMCIabnormalregions = zeros(1,length(controldata.AV_binding)); %Flag significantly abnormal regions cf control
        PSPabnormalregions = zeros(1,length(controldata.AV_binding));
        ADMCIpvalues = zeros(2,length(controldata.AV_binding));
        PSPpvalues = zeros(2,length(controldata.AV_binding));
        ADMCIpvalues(1,:) = 1:length(controldata.AV_binding);
        PSPpvalues(1,:) = 1:length(controldata.AV_binding);
        
        for region = 1:length(controldata.AV_binding)
            [ADMCIabnormalregions(region),ADMCIpvalues(2,region)] = ttest2(controldata.AV_binding(region,:),ADMCIdata.AV_binding(region,:),'vartype','unequal');
            [PSPabnormalregions(region),PSPpvalues(2,region)] = ttest2(controldata.AV_binding(region,:),PSPdata.AV_binding(region,:),'vartype','unequal');
        end
        ADMCIabnormalregions = logical(ADMCIabnormalregions);
        PSPabnormalregions = logical(PSPabnormalregions);
        
        % Then replot graphs with only unequal values
        % XXX
        
        %
        for i = 1:length(groupstodraw)
            if groupstodraw(i) == 1
                x_data = eval(['controldata.' x_variable]);
                y_data = eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]);
                titletext = ['Controls ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
            elseif groupstodraw(i) == 2
                x_data = eval(['ADMCIdata.' x_variable]);
                y_data = eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]);
                titletext = ['Coloured ADMCI ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
            elseif groupstodraw(i) == 3
                x_data = eval(['PSPdata.' x_variable]);
                y_data = eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]);
                titletext = ['Coloured PSP ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
            end
            
            bothabnormal = and(ADMCIabnormalregions,PSPabnormalregions);
            neitherabnormal = not(ADMCIabnormalregions) & not(PSPabnormalregions);
            onlyADabnormal = ADMCIabnormalregions & not(PSPabnormalregions);
            onlyPSPabnormal = not(ADMCIabnormalregions) & PSPabnormalregions;
            
            allxnanmeans = nanmean(x_data');
            allynanmeans = nanmean(y_data');
            
            %Simple scatter plot
            h = figure;
            hold on
            h_black = scatter(allxnanmeans(neitherabnormal),allynanmeans(neitherabnormal),'k');
            h_red = scatter(allxnanmeans(bothabnormal),allynanmeans(bothabnormal),'r');
            h_green = scatter(allxnanmeans(onlyADabnormal),allynanmeans(onlyADabnormal),'g');
            h_cyan = scatter(allxnanmeans(onlyPSPabnormal),allynanmeans(onlyPSPabnormal),'c');
            legend([h_black, h_red, h_green, h_cyan],{'Neither black','Both red','OnlyAD green','OnlyPSP cyan'})
            
            regres_coeffs = [allxnanmeans(neitherabnormal)',ones(size(allxnanmeans(neitherabnormal)'))] \ allynanmeans(neitherabnormal)';
            xs = min(allxnanmeans(neitherabnormal)):(range(allxnanmeans(neitherabnormal))/100):max(allxnanmeans(neitherabnormal));
            plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)), 'k');
            regres_coeffs = [allxnanmeans(bothabnormal)',ones(size(allxnanmeans(bothabnormal)'))] \ allynanmeans(bothabnormal)';
            xs = min(allxnanmeans(bothabnormal)):(range(allxnanmeans(bothabnormal))/100):max(allxnanmeans(bothabnormal));
            plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)), 'r');
            regres_coeffs = [allxnanmeans(onlyADabnormal)',ones(size(allxnanmeans(onlyADabnormal)'))] \ allynanmeans(onlyADabnormal)';
            xs = min(allxnanmeans(onlyADabnormal)):(range(allxnanmeans(onlyADabnormal))/100):max(allxnanmeans(onlyADabnormal));
            plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)), 'g');
            regres_coeffs = [allxnanmeans(onlyPSPabnormal)',ones(size(allxnanmeans(onlyPSPabnormal)'))] \ allynanmeans(onlyPSPabnormal)';
            xs = min(allxnanmeans(onlyPSPabnormal)):(range(allxnanmeans(onlyPSPabnormal))/100):max(allxnanmeans(onlyPSPabnormal));
            plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)), 'c');
            
            xlabel(x_variable, 'Interpreter', 'none')
            ylabel(y_variable{1}, 'Interpreter', 'none')
            title(titletext, 'Interpreter', 'none')
            saveas(h,[filteringtype '_' metrictype '_figures\Coloured_' titletext '.jpg']);
            close(h)
            
%                     %Multigradient plot
%                     h = figure
%                     hold on
%                     all_sig_gradients = [];
%                     all_nonsig_gradients = [];
%                     for j = 1:size(y_data,1)
%                         regres_coeffs = [x_data(j,:)',ones(size(x_data(j,:)'))] \ y_data(j,:)';
%                         xs = min(x_data(j,:)):(range((x_data(j,:)))/100):max(x_data(j,:));
%                         plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)));
%                         all_gradients(j) = regres_coeffs(1);
%             
%                     end
%                     xlabel(x_variable, 'Interpreter', 'none')
%                     ylabel(y_variable{1}, 'Interpreter', 'none')
%                     title(['Multigradient plot for' titletext], 'Interpreter', 'none')
%                     saveas(h,[filteringtype '_' metrictype '_figures\Multigradient_' titletext '.jpg']);
%                     close(h)
%             
%                     if groupstodraw(i) == 1
%                         controlgradients = all_gradients;
%                     elseif groupstodraw(i) == 2
%                         ADMCIgradients = all_gradients;
%                     elseif groupstodraw(i) == 3
%                         PSPgradients = all_gradients;
%                     end
            
            % Now plot a difference graph for each group
            if groupstodraw(i) == 2 || groupstodraw(i) == 3
                control_x_data = eval(['controldata.' x_variable]);
                control_y_data = eval(['controldata.' y_variable{1} '_' num2str(density_threshold)]);
                if groupstodraw(i) == 2
                    x_data = eval(['ADMCIdata.' x_variable]);
                    y_data = eval(['ADMCIdata.' y_variable{1} '_' num2str(density_threshold)]);
                    titletext = ['Coloured ADMCI normalised ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
                elseif groupstodraw(i) == 3
                    x_data = eval(['PSPdata.' x_variable]);
                    y_data = eval(['PSPdata.' y_variable{1} '_' num2str(density_threshold)]);
                    titletext = ['Coloured PSP normalised ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
                end
                
                %Simple scatter plot
                h = figure;
                if groupstodraw(i)==2
                    gscatter(nanmean(x_data')-nanmean(control_x_data'),nanmean(y_data')-nanmean(control_y_data'),ADMCIabnormalregions)
                else
                    gscatter(nanmean(x_data')-nanmean(control_x_data'),nanmean(y_data')-nanmean(control_y_data'),PSPabnormalregions)
                end
                regres_coeffs = [(nanmean(x_data')-nanmean(control_x_data'))',ones(size((nanmean(x_data')-nanmean(control_x_data'))'))] \ (nanmean(y_data')-nanmean(control_y_data'))';
                hold on
                xs = min(nanmean(x_data')-nanmean(control_x_data')):(range(nanmean(x_data')-nanmean(control_x_data'))/100):max(nanmean(x_data')-nanmean(control_x_data'));
                plot(xs, (regres_coeffs(1)*xs + regres_coeffs(2)));
                xlabel([x_variable ' difference'], 'Interpreter', 'none')
                ylabel([y_variable{1} ' difference'], 'Interpreter', 'none')
                title(titletext, 'Interpreter', 'none')
                saveas(h,[filteringtype '_' metrictype '_figures\Coloured_' titletext '.jpg']);
                close(h)
            end
        end
        
    end
        %         end
    end
end

%% Next step - perform representational similarity analyses
Correl_method = 'Pearson'; %Specify Spearman (non-parametric) or Pearson (parametric)
linkage_method = 'average'; %Specify linkage method - average, complete, nearest
mds_method = 'metricstress'; %Specify multidimensional scaling method

ncontrols = size(controldata.AV_binding,2);
nADMCI = size(ADMCIdata.AV_binding,2);
nPSP = size(PSPdata.AV_binding,2);

patientgroups = [repmat({'HC'},1,ncontrols),repmat({'AD'},1,nADMCI),repmat({'PSP'},1,nPSP)];
patientlabels = [1:ncontrols,1:nADMCI,1:nPSP];

all_variables = {'AV_binding'};
for density_threshold = 5:5:20; %density_threshold = 5:5:15;
    for y_variable = {'LocalEfficiency'}; %y_variable = {'LocalEfficiency','BetCent','PartCoeff','Degree'};
        all_variables{end+1} = [y_variable{1} '_' num2str(density_threshold)];
    end
end

if ~exist([pwd '\' filteringtype '_' metrictype '_figures\RSA\'] ,'dir')
    mkdir([pwd '\' filteringtype '_' metrictype '_figures\RSA\']);
end

for this_variable = all_variables
    eval(['full_array = [controldata.' this_variable{1} ', ADMCIdata.' this_variable{1} ', PSPdata.' this_variable{1} '];']);
    
    if sum(sum(isnan(full_array))) == 1 %ZZZ TEMPORARY FIX FOR THE CONTROL WITH A SINGLE MISSING VALUE DUE TO TOP OF AV PET BEING CUT OFF
        [row_nan,col_nan] = find(isnan(full_array));
        full_array(isnan(full_array)) = full_array(row_nan-1,col_nan);
    elseif sum(sum(isnan(full_array))) > 1
        error('More than one missing value in the array - something has gone wrong')
    end
    corr_full_array = corr(full_array,'type',Correl_method);
    
    h = figure;
    imagesc(corr_full_array);
    colormap('jet');colorbar;
    title(['Correlation matrix of ' this_variable{1} ' for all subjects - 1-' num2str(ncontrols) ' = Controls, ' num2str(ncontrols+1) '-' num2str(ncontrols+nADMCI) ' = ADMCI, ' num2str(ncontrols+nADMCI+1) '-' num2str(ncontrols+nADMCI+nPSP) ' = PSP', ], 'Interpreter', 'none')
    saveas(h,[filteringtype '_' metrictype '_figures\RSA\Correlation_self_' this_variable{1} '.jpg']);
    close(h)
    
    tri_corr = tril(corr_full_array,-1);
    dissimilarity = 1 - tri_corr(find(tri_corr))';
    cutoff = 0.7;
    Z = linkage(dissimilarity,linkage_method);
    groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
    c = cophenet(Z,dissimilarity);
    h = figure;
    dendrogram(Z,0,'colorthreshold',cutoff)
    title(['Dendrogram for all subjects based on ' this_variable{1}], 'Interpreter', 'none')
    
    group_subj = cell(size(patientgroups));
    for i = 1:length(patientgroups)
        group_subj(i) = {[patientgroups{i} '\newline' num2str(patientlabels(i))]};
    end
    
    perm=str2num(get(gca,'XtickLabel'));
    set(gca,'XTickLabel',group_subj(perm),'FontSize',8)
    saveas(h,[filteringtype '_' metrictype '_figures\RSA\Dendrogram_' this_variable{1} '.jpg']);
    close(h)
    
    %Now plot the mds graph
    points = mdscale([1-corr_full_array],2,'criterion',mds_method);
    ADcolor = [1 0 1];
    %MCIcolor = [0 1 1];
    PSPcolor = [0 1 0];
    Controlcolor = [0 0 0];
    allcolors = [repmat(Controlcolor,ncontrols,1);repmat(ADcolor,nADMCI,1);repmat(PSPcolor,nPSP,1)];
    h = figure;
    scatter(points(:,1),points(:,2),48,allcolors,'filled')
    title(['Multi dimensional scaling of ' this_variable{1} ' distribution dissimilarities (Control black, AD pink, PSP green)'], 'Interpreter', 'none')
    saveas(h,[filteringtype '_' metrictype '_figures\RSA\Multi_Dimensional_Scaling_' this_variable{1} '.jpg']);
    close(h)
    eval(['full_' this_variable{1} '_array = full_array;']);
    eval(['corr_full_' this_variable{1} '_array = corr_full_array;']);
    eval(['mds_' this_variable{1} '_points = points;']);
    
end

for first_variable = all_variables
    for second_variable = all_variables
        
        if strcmp(first_variable{1},second_variable{1}) == 0
            eval(['corr_full_array = corr(full_' first_variable{1} '_array,full_' second_variable{1} '_array,''type'',Correl_method);'])
            h = figure;
            imagesc(corr_full_array);
            colormap('jet');colorbar;
            ylabel(first_variable{1}, 'Interpreter', 'none')
            xlabel(second_variable{1}, 'Interpreter', 'none')
            title(['Correlation matrix of ' first_variable{1} ' against ' second_variable{1} ' for all subjects - 1-' num2str(ncontrols) ' = Controls, ' num2str(ncontrols+1) '-' num2str(ncontrols+nADMCI) ' = ADMCI, ' num2str(ncontrols+nADMCI+1) '-' num2str(ncontrols+nADMCI+nPSP) ' = PSP', ], 'Interpreter', 'none')
            saveas(h,[filteringtype '_' metrictype '_figures\RSA\Correlation_matrix_' first_variable{1} '_against_' second_variable{1} '.jpg']);
            close(h)
            
        end
        
    end
end




