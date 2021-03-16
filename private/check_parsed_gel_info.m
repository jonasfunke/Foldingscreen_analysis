function [parsed] = check_parsed_gel_info(parsed_data)
    % this script checks if the parsed data from the gel-info.txt is ok
    %% check if all important objects are there
    parsed = false;
    fields_required = {'user', 'project', 'design_name', 'date', 'scaffold_type', ...
        'scaffold_concentration', 'staple_concentration', 'comment', ...
        'lanes', 'gelsize', 'agarose_concentration', 'staining', 'mg_concentration', ...
        'voltage', 'running_time', 'cooling', 'tem_verified', 'comment'};

    missing_fields = {};

    fileID = fopen(parsed_data.log_file,'a');
    
    for i=1:length(fields_required)
        if ~isfield(parsed_data, fields_required{i})
            disp(['Warning: ' fields_required{i} ' not found. Check the gel_info.txt file.'])
            fprintf(fileID,'%s\n', ['Warning: ' fields_required{i} ' not found. Check the gel_info.txt file.']);
            missing_fields = [missing_fields fields_required{i}];

        end
    end

    if isempty(missing_fields)
        disp(['Data OK: ' parsed_data.filepath])
        fprintf(fileID,'%s\n', ['Data OK: ' parsed_data.filepath]);
        parsed = true;
    else
        tmp = join(missing_fields(:), ', ');
        disp(['Warning: Parsed data has missing fields: ' tmp{1}] )
        fprintf(fileID,'%s\n', ['Warning: Parsed data has missing fields: ' tmp{1}] );
    end

    if isfield(parsed_data, 'lanes')
        tmp = join(parsed_data.lanes(:), ', ');
        disp([num2str(length(parsed_data.lanes)) ' lanes detected: ' tmp{1}])
        fprintf(fileID,'%s\n', [num2str(length(parsed_data.lanes)) ' lanes detected: ' tmp{1}] );
    else
        disp(['WARNING: NO lanes detected. Fix gel_info_file'])
        fprintf(fileID,'%s\n', ['WARNING: NO lanes detected. Fix gel_info_file']);
        parsed = false;
    end
    
    fclose(fileID);

end

