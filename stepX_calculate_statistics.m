%% -------------------- Analyze all screens -----------------------------------------
close all, clear all

%log_out = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/data.out';
root_path = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/';

%logfile_ID = fopen(log_out,'a');
%fprintf(logfile_ID,'%s\n', '------');
%fclose(logfile_ID);

%logfile_ID = fopen(log_out,'a');


folders = dir(root_path);

discard = {'.', '..', 'AAA_TEMPLATE', 'ZZZfolder_of_shame_aka_missing_data', 'ZZZnon_standard_folding_screens', ...
    'ZZZ_output'};

i_discard = [];
for i=1:length(folders)
    if any(strcmp(discard, folders(i).name))
        disp(['Discarding folder ' folders(i).name])
        i_discard = [i_discard, i];
    end
    if ~folders(i).isdir
        i_discard = [i_discard, i];
    end
end
folders(i_discard) = [];

metrics = cell(length(folders),1);
for i=1:length(folders)
    if folders(i).isdir
        disp('---')
        %fprintf(logfile_ID,'%s\n', '------');
        txt_files = dir([folders(i).folder filesep folders(i).name filesep '*.txt']);
        tif_files = dir([folders(i).folder filesep folders(i).name filesep '*.tif']);
        mat_files = dir([folders(i).folder filesep folders(i).name filesep '*.mat']);
        
        if length(txt_files)==1
           % [tmp, warnings] = parse_gel_info([txt_files(1).folder filesep txt_files(1).name], log_out);
           %     check_parsed_gel_info(tmp);
        else
            disp(['More or less than one txt file found in ' folders(i).name ' (' num2str(length(txt_files)) ' found).' ])
            %fprintf(logfile_ID,'%s\n', ['More or less than one txt file found in ' folders(i).name ' (' num2str(length(txt_files)) ' found).' ]);
        end
        
        if length(tif_files)==1
            %[tmp, warnings] = parse_gel_info_simple([txt_files(1).folder filesep txt_files(1).name]);
            %    check_parsed_gel_info(tmp);
        else
            disp(['More or less than one tif file found in ' folders(i).name ' (' num2str(length(tif_files)) ' found).' ])
            %fprintf(logfile_ID,'%s\n', ['More or less than one tif file found in ' folders(i).name ' (' num2str(length(tif_files)) ' found).' ]);
        end
        
        if length(mat_files)==1
            %[tmp, warnings] = parse_gel_info_simple([txt_files(1).folder filesep txt_files(1).name]);
            %    check_parsed_gel_info(tmp);
            cur_data = load([mat_files(1).folder filesep mat_files(1).name]);
            if isfield(cur_data.profileData, 'monomerFits')
                metrics{i} = get_best_folding(cur_data.profileData, cur_data.gelInfo, cur_data.gelData, false);
            end
            
        else
            disp(['More or less than one .mat file found in ' folders(i).name ' (' num2str(length(mat_files)) ' found).' ])
            %fprintf(logfile_ID,'%s\n', ['More or less than one tif file found in ' folders(i).name ' (' num2str(length(tif_files)) ' found).' ]);
        end
       
       
       
    end
end
%fclose(logfile_ID);

%%
close all
path_out = '/Users/jonasfunke/Dropbox (DIETZ LAB)/FOLDINGSCREENS/ZZZ_output/';

cur_fig = figure(1);
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 20 15], 'PaperSize', [ 20 15]);
d = [];
for i=1:length(metrics)
    if ~isempty(metrics{i})
        d = [d; metrics{i}.fractionMonomer(metrics{i}.bestFoldingIndex)];
        if metrics{i}.fractionMonomer(metrics{i}.bestFoldingIndex) > 1 
            disp(folders(i).name)
        end
    end
end
histogram(d, 20)

xlabel('Monomer fraction of best folding condition')
ylabel('Frequency')

print(cur_fig, '-dpdf', [path_out 'monomer_fraction_of_best_folding.pdf']); %save figure

%%
cur_fig = figure(2);
set(gcf,'Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters', ...
    'PaperPosition', [0 0 30 15], 'PaperSize', [ 30 15]);
best_folding_condition = {};
best_T = [];
best_Mg = [];

for i=1:length(metrics)
    if ~isempty(metrics{i})
        best_folding_condition = [best_folding_condition; metrics{i}.bestFolding];
        if strcmpi(metrics{i}.bestTscrn, 'T1')
            cur_T = 48;
        end
        if strcmpi(metrics{i}.bestTscrn, 'T2')
            cur_T = 50;
        end
        if strcmpi(metrics{i}.bestTscrn, 'T3')
            cur_T = 52;
        end
        if strcmpi(metrics{i}.bestTscrn, 'T4')
            cur_T = 54;
        end
        if strcmpi(metrics{i}.bestTscrn, 'T5')
            cur_T = 56;
        end
        if strcmpi(metrics{i}.bestTscrn, 'T6')
            cur_T = 58;
        end
        if strcmpi(metrics{i}.bestTscrn, 'T7')
            cur_T = 60;
        end
        if strcmpi(metrics{i}.bestTscrn, 'T8')
            cur_T = 62;
        end
        best_T = [best_T; cur_T];
        
        
        
        
        
        best_Mg = [best_Mg; str2double(metrics{i}.bestMgscrn(2:end))];
       
    end
end
subplot(1, 3, 1)
tmp = categorical(best_folding_condition);
histogram(tmp)
ylabel('Frequency')
title('Best folding condition')

subplot(1, 3, 2)
histogram(best_T, [46.5:1:63])
ylabel('Frequency')
xlabel('Folding temperature')
title('Best folding temperature')
set(gca, 'XTick', [48:2:62])

subplot(1, 3, 3)
histogram(best_Mg, [2.5:5:32.5])
xlabel('Mg concentration')
ylabel('Frequency')
title('Best Mg concentration')

print(cur_fig, '-dpdf', [path_out 'best_folding_T_Mg.pdf']); %save figure

%%
figure(3)
plot(best_T, d, '.')
%%

