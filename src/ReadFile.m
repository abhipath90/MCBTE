function file_data = ReadFile(file_name)
    
% return filedata if file exists and is not empty
% if file is empty it retruns empty variable
% if file does not exists it return file name
    try
        file_data = load(file_name);
    catch exception
        file_data = file_name;
    end
        
end

