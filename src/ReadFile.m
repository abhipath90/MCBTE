function file_data = ReadFile(file_name)
    
% return filedata if file exists and is not empty
% if file is empty it retruns empty variable
% if file does not exists it return file name
    try
        file_data = load(file_name);
    catch exception
        file_data = file_name;
    end
    

    %{
    if isfile(file_name) % if file exists
        fid = fopen(file_name);
        if fseek(fid,1,'bof') == -1
            % empty file
            file_data = [];
        else 
            file_data = load(file_name);
        end % if fseek(fid,
        
    else
        % file does not exists
        file_data = file_name;
    end % if isfile(file_name) % if file exists
    %}
    
end

