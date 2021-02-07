
%grad_file = load('Thermal_gradient.txt');
%grad = norm(grad_file(1,3:5));
detect = load('detector_location.txt');

%detect = load(['det' num2str(len(ii)) '.txt']);
% Temp = sum(load(['T' num2str(len(ii)) '.txt']),2);
%Qx = sum(load(['Qx' num2str(len(ii)) '.txt']),2);
%Qy = sum(load(['Qy' num2str(len(ii)) '.txt']),2);
% Qz = sum(load(['Qz' num2str(len(ii)) '.txt']),2);

Temp_modal = load('T300.txt');
[row,col] = size(Temp_modal);
Temp = sum(Temp_modal,2);

Qz_modal = load('Qz300.txt');
Qz = sum(Qz_modal,2);

%Kappa = mean(Qz)/grad;

centroid = [(detect(:,1)+detect(:,2))/2 (detect(:,3)+detect(:,4))/2 (detect(:,5)+detect(:,6))/2];
zpts = unique(centroid(:,3));

index_z = cell(length(zpts),1);

for ii=1:length(index_z)
    index_z{ii} = find(centroid(:,3)==zpts(ii));
end

Tpts = zeros(length(zpts),1);
Tpts_modal = zeros(length(zpts),col);

% Qxpts = zeros(length(zpts),1);
% Qypts = zeros(length(zpts),1);
Qzpts = zeros(length(zpts),1);
for ii=1:length(zpts)
    Tpts(ii) = mean(Temp(index_z{ii}));
%     Qxpts(ii) = mean(Qx(index_z{ii}));
%     Qypts(ii) = mean(Qy(index_z{ii}));
    Qzpts(ii) = mean(Qz(index_z{ii}));
    
    for jj=1:col
        Tpts_modal(ii,jj) = mean(Temp_modal(index_z{ii},jj));
    end
end

font_size = 16;
figure(1);
plot(zpts*1e9,Tpts+300,'*');
set(gca, 'FontSize', font_size);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
xlabel('Z-coordinate (nm)');
ylabel('Temperature (K)');
saveas(gcf,'temp.fig');

%figure(2);
%plot(zpts*1e9,Tpts_modal(:,1)+300);
%hold on;
%for jj=2:col
%    plot(zpts*1e9,Tpts_modal(:,jj)+300);
%end
%set(gca, 'FontSize', font_size);
%set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
%xlabel('Z-coordinate (nm)');
%ylabel('Temperature (K)');
%saveas(gcf,'temp_modal.fig');


% 
% figure(2);
% plot(zpts*1e9,Qxpts,'*');
% set(gca, 'FontSize', font_size);
% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
% xlabel('Z-coordinate (nm)');
% ylabel('q_x (W/m^2)');
% saveas(gcf, 'Qx.fig');
% 
% figure(3);
% plot(zpts*1e9,Qypts,'*');
% set(gca, 'FontSize', font_size);
% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
% xlabel('Z-coordinate (nm)');
% ylabel('q_y (W/m^2)');
% saveas(gcf,'Qy.fig');
% 
% figure(4);
% plot(zpts*1e9,Qzpts*6,'*');
% set(gca, 'FontSize', font_size);
% set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
% xlabel('Z-coordinate (nm)');
% ylabel('q_z (W/m^2)');
% saveas(gcf,'Qz.fig');
