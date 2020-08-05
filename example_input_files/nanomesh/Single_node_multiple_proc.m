c = parcluster(); % creating parallel cluster
job = createCommunicatingJob(c,'Type','SPMD'); % creating communicating job
job.AutoAttachFiles = true;
Dir =  'Full_path_to_current_folder'; % location of input files and where the output will be written
src = 'full_path_to_src_folder'; % Location of the source files

set(job,'AdditionalPaths',{Dir,src});

job.NumWorkersRange = 2; % Number of processors to be used.

task = createTask(job,@BTE_solution_3D,1,{300}); % creating taks. 300 is a dummy number.

submit(job);