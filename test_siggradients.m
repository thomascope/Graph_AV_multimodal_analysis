for groupstodraw = 1:3
    
    for density_threshold = 5:5:15;
        for y_variable = {'LocalEfficiency','BetCent','PartCoeff'};
            x_variable = 'AV_binding';
            
            this_control_x_data = eval(['mean(controldata.' x_variable ',2)']);
            this_ADMCI_x_data = eval(['mean(ADMCIdata.' x_variable ',2)']);
            this_PSP_x_data = eval(['mean(PSPdata.' x_variable ',2)']);
            this_x_data = [this_control_x_data;this_ADMCI_x_data;this_PSP_x_data];
            this_control_y_data = eval(['mean(controldata.' y_variable{1} '_' num2str(density_threshold) ',2)']);
            this_ADMCI_y_data = eval(['mean(ADMCIdata.' y_variable{1} '_' num2str(density_threshold) ',2)']);
            this_PSP_y_data = eval(['mean(PSPdata.' y_variable{1} '_' num2str(density_threshold) ',2)']);
            this_y_data = [this_control_y_data;this_ADMCI_y_data;this_PSP_y_data];
            titletext = ['All subjects ' y_variable{1} ' against ' x_variable ' at a density threshold of ' num2str(density_threshold)];
            
            groups = [repmat({'Control'},1,length(this_control_x_data)), repmat({'ADMCI'},1,length(this_ADMCI_x_data)), repmat({'PSP'},1,length(this_PSP_x_data))]';
            
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
                
                %Multigradient plot
                all_gradients = [];
                cols = hsv(numel(unique(atlasdata(:,2))));
                regionorder = unique(atlasdata(:,2),'stable');
                all_gradients = zeros(1,size(y_data,1));
                R_gradients = zeros(1,size(y_data,1));
                P_gradients = zeros(1,size(y_data,1));
                for j = 1:size(y_data,1)
                    regres_coeffs = [x_data(j,:)',ones(size(x_data(j,:)'))] \ y_data(j,:)';
                    xs = min(x_data(j,:)):(range((x_data(j,:)))/100):max(x_data(j,:));
                    all_gradients(j) = regres_coeffs(1);
                    [R,P] = corrcoef(x_data(j,:),y_data(j,:));
                    R_gradients(j) = R(2,1);
                    P_gradients(j) = P(2,1);
                end
                
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
                
            end
        end
    end
    
end