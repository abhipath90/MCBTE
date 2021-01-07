function f=Distributed_computing()
    %% Submission script for using the MDCS at CCR, API-2
    startupSlurm %need this to load API-2 settings
    Dir = 'fullpath_to_current_directory';
    src = 'fullpath_to_source_code_directory';
    StorageLocation=Dir; %<== enter YOUR dir
    ppn=32;
    time='1:00:00';
    %%<== enter YOUR e-mail followed by any extra SBATCH directives
    email='abc@xyz.mno --mail-type=FAIL --partition=partition_name --qos=qos --nodes=1 --output=matlab.out --error=matlab.err'; 

    %% pass to CommunicatingSubmitFcn using env. vars 
    setenv('MBATCH_PPN',int2str(ppn))
    setenv('MBATCH_TIME',time)
    setenv('MBATCH_EMAIL',email)

    %%Don't modify this line'%% changing storage location to the current folder
    set(u2,'JobStorageLocation',StorageLocation);

    %% for spmd code
    pjob = createCommunicatingJob(u2,'Type', 'SPMD');
    pjob.NumWorkersRange = [1 32]; %%Set Maximum Workers (cores);
    
    %% auto attach all the file in the directory assigned in the path
    pjob.AutoAttachFiles = true;

    %%*****************

    %% Example for submitting a function (My_MDCS_Script.m) w/ 1 output and 1 input:   
    set(pjob, 'AdditionalPaths', {Dir,src}); %
    createTask(pjob,@BTE_solution_3D,1,{300}); % The code to run, 300 is a dummy number
    
    
    submit(pjob);
    f=pjob;

    % save job information to file, to use in conjunction with 
    % mbatchGetOutput
    diary(['JobInfo' num2str(f.ID) '.txt'])
    disp(f);
    disp(u2);
    diary off;

end


