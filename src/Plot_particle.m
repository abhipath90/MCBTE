% This one plots the trajectories between two scattering events with
% straight lines and events as scatter plot

% internal boundaries
In_bnd = load('In_bnd.txt');
In_bnd = In_bnd(:,1:4);
[len, ~] = size(In_bnd);

movieLength = 5; % duration of movie in seconds
N =10; % number of trajectories
for ii=1:N
    filename = [num2str(ii) '.txt'];
    Traj = load(filename);
    
    figure();
    hold on
    % Plotting internal boundaries
    for jj=1:len
        plot(In_bnd(jj,1:2:3),In_bnd(jj,2:2:4));
    end
    %removing zero rows
    %Traj(~any(Traj,2),:) = [];
    [row, col] = size(Traj);
    for i=1:row
        if(i~=1 && Traj(i,8)~=4)
            plot(Traj((i-1):i,2),Traj((i-1):i,3));
        end
        plot(Traj(i,2),Traj(i,3), 'ro');
        xlim([0 34e-9]);
        ylim([0 34e-9])
        pause(movieLength/length(Traj));
        title(['Particle ID ' num2str(ii)])
    end
    disp(ii);
    %pause;
    close all
end
