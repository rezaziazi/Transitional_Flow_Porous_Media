% #########################################################################
% ############## Transition to turbulence in Porous Media #################
% ##############               By: Reza Ziazi             #################
% ##############                 Fall 2018                #################
% #########################################################################
% tic
clear
clc
close all
tic
addpath(genpath('Subroutines_PIV_porous'))
addpath(genpath('Plots'))

%Porous bed specs *********************************************************
specs.fluid.dens = 1118;                                            %kg/m^3
specs.fluid.kvisc = 1.288*10^-6;                                    %m^2/s
specs.fluid.dvisc = specs.fluid.kvisc*specs.fluid.dens;             %kg/(ms)
specs.fluid.temp = 20;                                              %degree c
specs.geom.DB = 0.015;                                              %m
specs.geom.phi = 0.489;                                             % porosity
specs.geom.DH = specs.geom.phi* specs.geom.DB/(1- specs.geom.phi);  % hydraulic diameter in [m]
specs.geom.DH = specs.geom.phi* specs.geom.DB/(1- specs.geom.phi);  % hydraulic diameter in [m]
specs.geom.L = 0.07;                                                % test section width

% Time between the pulses *************************************************
% Time between the pulses for a vector displacement in seconds
Del_T.Data1 = 16210/10^6;
Del_T.Data2 = 16210/10^6; 
Del_T.Data3 = 16210/10^6; 
Del_T.Data4 = 48760/10^6;
Del_T.Data5 = 43309.7843/10^6;
Del_T.Data6 = 43309.7843/10^6; 
Del_T.Data7 = 16250.0781/10^6; 
Del_T.Data8 = 12170.0531/10^6;
Del_T.Data9 = 9767.1147/10^6;
Del_T.Data10 = 8125.9769/10^6; 
Del_T.Data11 = 6973.1166/10^6; 
Del_T.Data12 = 6095.0297/10^6;
Del_T.Data13 = 5417.9431/10^6;
Del_T.Data14 = 5417.9431/10^6; 
Del_T.Data15 = 6089.0903/10^6; 
Del_T.Data16 = 6958.1119/10^6; 
Del_T.Data17 = 8125.9769/10^6;
Del_T.Data18 = 9752.1100/10^6;
Del_T.Data19 = 12170.0531/10^6; 
Del_T.Data20 = 16210.0656/10^6; 
Del_T.Data21 = 24280.0875/10^6;
Del_T.Data22 = 24280.0875/10^6; 
Del_T.Data23 = 48759.925/10^6; 

% Optical specifications **************************************************
[mm_to_pix, pix_to_mm] = Optics;

% Flow meter velocity *****************************************************
flowmeter = flowrate(specs);

% Loading velocities ******************************************************
[x, y, vel, disp] = Dataloader(Del_T,1);
x.X_Nodes1_pix = unique(x.x_pix1D.Data1(:,1,1)); %unique for all planes
y.Y_Nodes1_pix = unique(y.y_pix1D.Data1(:,1,1)); %unique for all planes
x.X_Nodes1_mm = x.X_Nodes1_pix * mm_to_pix;
y.Y_Nodes1_mm = y.Y_Nodes1_pix * mm_to_pix;

n = 0; 
imax1 = length(x.X_Nodes1_pix); %77
jmax1 = length(y.Y_Nodes1_pix); %61
kmax1 = 3199;

U_temp = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix));
V_temp = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix));

% [x, y, vel, disp] = Dataloader(Del_T.Data2,2);
x.X_Nodes2_pix = unique(x.x_pix1D.Data2(:,1,1)); %unique for Data 2, 4, 18, 19, 20, 21
y.Y_Nodes2_pix = unique(y.y_pix1D.Data2(:,1,1)); %unique for Data 2, 4, 18, 19, 20, 21
x.X_Nodes2_mm = x.X_Nodes2_pix * mm_to_pix;
y.Y_Nodes2_mm = y.Y_Nodes2_pix * mm_to_pix;

imax2 = length(x.X_Nodes2_pix); %79
jmax2 = length(y.Y_Nodes2_pix); %63
kmax2 = 3199;

U_temp2 = zeros(length (y.Y_Nodes2_pix),length (x.X_Nodes2_pix));
V_temp2 = zeros(length (y.Y_Nodes2_pix),length (x.X_Nodes2_pix));

specs.geom.Del_x = (x.X_Nodes1_pix(2)- x.X_Nodes1_pix(1))*mm_to_pix;
specs.geom.Del_y = (y.Y_Nodes1_pix(2)- y.Y_Nodes1_pix(1))*mm_to_pix;


% LES Filtering ***********************************************************

% Top-Hat - Convlution by a Kernel of matrix of B = ones(3,3) into instantaneus velocity (Top-Hat)


%3*3 Top-Hat Kernel (2D Field)
kernel.B_tophat3 = ones(3,3)/3^2;       % D = 3.52 mm
kernel.B_tophat5 = ones(5,5)/5^2;       % D = 5.27 mm
kernel.B_tophat7 = ones(7,7)/7^2;       % D = 7.03 mm
kernel.B_tophat9 = ones(9,9)/9^2;       % D = 8.79 mm
kernel.B_tophat11 = ones(11,11)/11^2;   % D = 10.55 mm
kernel.B_tophat13 = ones(13,13)/13^2;   % D = 12.3 mm
kernel.B_tophat15 = ones(15,15)/15^2;   % D = 14.06 mm
kernel.B_tophat17 = ones(17,17)/17^2;   % D = 15.8 mm
kernel.B_tophat19 = ones(19,19)/19^2;   % D = 17.58 mm
kernel.B_tophat21 = ones(21,21)/21^2;   % D = 17.58 mm
kernel.B_tophat23 = ones(23,23)/23^2;   % D = 17.58 mm

kernel.B_Gauss3 = [1 2 1; 2 4 2; 1 2 1]/4^2;  
kernel.B_Gauss5 = [1 4 6 4 1; 4 16 24 16 4; 6 24 36 24 6; 4 16 24 16 4; 1 4 6 4 1]/16^2;  




%% Data 3

% Loading velocities ******************************************************

[x, y, vel, disp] = Dataloader(Del_T,3,x,y);
% [x, y, vel, disp, DT] = Dataloader(Del_T,1,x,y);

vel.U_pix2D.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),3199);  
vel.V_pix2D.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),3199);
x.x_pix2D.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),3199);
y.y_pix2D.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),3199);
for k = 1:kmax1
    for j = 1:jmax1
        for i = 1:imax1
            n = n+1;
            vel.U_pix2D.Data3(j,i,k) = vel.U_pix1D.Data3(n,1,k); vel.V_pix2D.Data3(j,i,k) = vel.V_pix1D.Data3(n,1,k);                
            x.x_pix2D.Data3(j,i,k) = x.x_pix1D.Data3(n,1,k);y.y_pix2D.Data3(j,i,k) = y.y_pix1D.Data3(n,1,k);
        end
    end
    n = 0;
end

% Conversion from pixel to mm ********************************************
% 1D vectors of location ******
x.x_mm1D.Data3 = x.x_pix1D.Data3*mm_to_pix;     
y.y_mm1D.Data3 = y.y_pix1D.Data3*mm_to_pix;

% 1D vectors of velocity **********************
vel.V_mm1D.Data3 = vel.V_pix1D.Data3*mm_to_pix;
vel.U_mm1D.Data3 = vel.U_pix1D.Data3*mm_to_pix;
vel.Mag_mm1D.Data3 = sqrt(vel.U_mm1D.Data3.^2 + vel.V_mm1D.Data3.^2);

% 2D Field*********************
% 2D vectors of location ******
x.x_mm2D.Data3 = x.x_pix2D.Data3(1:jmax1,1:imax1,:)*mm_to_pix;
y.y_mm2D.Data3 = y.y_pix2D.Data3(1:jmax1,1:imax1,:)*mm_to_pix;

% 2D vectors of velocity **********************
vel.V_mm2D.Data3 = vel.V_pix2D.Data3(1:jmax1,1:imax1,:)*mm_to_pix;
vel.U_mm2D.Data3 = vel.U_pix2D.Data3(1:jmax1,1:imax1,:)*mm_to_pix;
vel.Mag_mm2D.Data3 = sqrt(vel.U_mm2D.Data3.^2 + vel.V_mm2D.Data3.^2);
toc
%% LIC Image

for k = 1:120
    LIC_bed(vel.U_mm2D.Data3(:,:,k),vel.V_mm2D.Data3(:,:,k),x,imax1, jmax1, 1280,20,2);
    saveas(gcf,sprintf('LIC_Data3_T%d.png',k));
    close
%     hold off
end

% imageFolder2mpeg('LIC_Data3','frameRate',11) %  This function creates movie

%% Writing files into a Dat file 
% 
% for k = 1:3200
%     mat = [vel.U_mm1D.Data3(:,1,k), vel.V_mm1D.Data3(:,1,k)];
%     fileName = fopen(sprintf('velocity_Data3_%d.dat',k),'w'); 
%     fprintf(fileName, '%f %f\n', mat.');
%     fclose(fileName)
% end

for k = 1:10
    csvwrite(sprintf('Vx_%d.csv',k),vel.U_mm2D.Data3(:,:,k))
    csvwrite(sprintf('Vy_%d.csv',k),vel.V_mm2D.Data3(:,:,k))
end


% for k = 1:1
%     filename = sprintf('VTK_Data3_%d',k);
%     vtkwrite(filename,'structured_grid',x.x_mm1D.Data3,y.y_mm1D.Data3,'vectors',title,vel.U_mm1D.Data3(:,1:k),vel.V_mm1D.Data3(:,1:k));
% end

%% 2D velocity *************************************************************

for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            if (abs(vel.Mag_mm2D.Data3(i,j,k))<8)
                vel.V_mm2D.Data3(i,j,k) = 5*vel.V_mm2D.Data3(i,j,k);
                vel.U_mm2D.Data3(i,j,k) = 2*vel.U_mm2D.Data3(i,j,k);
            else
                vel.V_mm2D.Data3(i,j,k) = 1.7*vel.V_mm2D.Data3(i,j,k);
                vel.U_mm2D.Data3(i,j,k) = 1.7*vel.U_mm2D.Data3(i,j,k);
            end
        end
    end
end
% %%
% for k = 1:kmax1
%     for i = 29:41
%         for j = 29:41
%             if (abs(vel.Mag_mm2D.Data3(i,j,k))<4)
%                 vel.V_mm2D.Data3(i,j,k) = 5*vel.V_mm2D.Data3(i,j,k);
%                 vel.U_mm2D.Data3(i,j,k) = 5*vel.U_mm2D.Data3(i,j,k);
%             else
%                 vel.V_mm2D.Data3(i,j,k) = 1.1*vel.V_mm2D.Data3(i,j,k);
%                 vel.U_mm2D.Data3(i,j,k) = 1.1*vel.U_mm2D.Data3(i,j,k);
%             end
%         end
%     end
% end
% %%
% for k = 1:kmax1
%     for i = 1:29
%         for j = 41:imax1
%             if (abs(vel.Mag_mm2D.Data3(i,j,k))<6)
%                 vel.V_mm2D.Data3(i,j,k) = 5*vel.V_mm2D.Data3(i,j,k);
%                 vel.U_mm2D.Data3(i,j,k) = 5*vel.U_mm2D.Data3(i,j,k);
%             else
%                 vel.V_mm2D.Data3(i,j,k) = 1.5*vel.V_mm2D.Data3(i,j,k);
%                 vel.U_mm2D.Data3(i,j,k) = 1.5*vel.U_mm2D.Data3(i,j,k);
%             end
%         end
%     end
% end


%% Time averaging over the field ******************************
%1D vector
tic
for i = 1:imax1*jmax1
    vel.U_mt1D.Data3(i) = mean(nonzeros(vel.U_mm1D.Data3(i,1,:))); vel.V_mt1D.Data3(i) = mean(nonzeros(vel.V_mm1D.Data3(i,1,:))); vel.Mag_mt1D.Data3(i) = mean(nonzeros(vel.Mag_mm1D.Data3(i,1,:)));
end

% 2D vector
for i = 1:jmax1 
    for j = 1:imax1
        vel.U_mt2D.Data3(i,j) = mean(nonzeros(vel.U_mm2D.Data3(i,j,:))); vel.V_mt2D.Data3(i,j) = mean(nonzeros(vel.V_mm2D.Data3(i,j,:))); vel.Mag_mt2D.Data3(i,j) = mean(nonzeros(vel.Mag_mm2D.Data3(i,j,:)));
    end
end


% Spatial and Temporal Average of the whole bed
vel.U_mm2D_NaN.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),3199);
vel.V_mm2D_NaN.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),3199);
vel.Mag_mm2D_NaN.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),3199);

for k = 1:kmax1
    U_temp_NaN = vel.U_mm2D.Data3(:,:,k);
    U_temp_NaN(find(U_temp_NaN == 0)) = NaN;
    vel.U_mm2D_NaN.Data3(:,:,k) = U_temp_NaN;
    V_temp_NaN = vel.V_mm2D.Data3(:,:,k);
    V_temp_NaN(find(V_temp_NaN == 0)) = NaN;
    vel.V_mm2D_NaN.Data3(:,:,k) = V_temp_NaN;
    Mag_temp_NaN = vel.Mag_mm2D.Data3(:,:,k);
    Mag_temp_NaN(find(Mag_temp_NaN == 0)) = NaN;
    vel.Mag_mm2D_NaN.Data3(:,:,k) = Mag_temp_NaN;    
end

% 2D vector
for i = 1:jmax1
    for j = 1:imax1
        vel.U_mt2D_NaN.Data3(i,j) = nanmean(vel.U_mm2D_NaN.Data3(i,j,:));vel.V_mt2D_NaN.Data3(i,j) = nanmean(vel.V_mm2D_NaN.Data3(i,j,:)); vel.Mag_mt2D_NaN.Data3(i,j) = nanmean(vel.Mag_mm2D_NaN.Data3(i,j,:));
    end
end

vel.U_mts2D_NaN.Data3 = nanmean(nanmean(vel.U_mt2D_NaN.Data3));
vel.V_mts2D_NaN.Data3 = nanmean(nanmean(vel.V_mt2D_NaN.Data3));
vel.Mag_mts2D_NaN.Data3 = nanmean(nanmean(vel.Mag_mt2D_NaN.Data3));


% Velocity fluctuation over the field ************************************

% 1D vector
for k = 1:3199
    for i = 1:imax1*jmax1    
        vel.u_t1D.Data3(i,1,k) = vel.U_mm1D.Data3(i,1,k) - vel.U_mt1D.Data3(i); vel.v_t1D.Data3(i,1,k) = vel.V_mm1D.Data3(i,1,k) - vel.V_mt1D.Data3(i);
    end
end

%2D vector
for k = 1:3199
    for i = 1:jmax1
        for j = 1:imax1
            vel.u_t2D.Data3(i,j,k) = vel.U_mm2D.Data3(i,j,k) - vel.U_mt2D.Data3(i,j); vel.v_t2D.Data3(i,j,k) = vel.V_mm2D.Data3(i,j,k) - vel.V_mt2D.Data3(i,j);
        end
    end
end
toc
% Fluctuation(vel, x, imax1, jmax1)

%% Generating Average velocity plots 

[Plot.U_avg.Data3,Plot.vec_avg.Data3] = Timeavgplot (vel.U_mt2D.Data3,vel.V_mt2D.Data3, x, imax1, jmax1);
% saveas(Plot.vec_avg.Data3,sprintf('vec_avg_Data3.png'));
print('-dpng','-r700',sprintf('vec_avg_Data3.png',k));
% saveas(h2,sprintf('contour_avg_Data3.png'));
%% Generating time laps of velocity plots 

for k = 10:10
    Plot.vec.Data3 = Timeavgplot(vel.U_mt2D.Data3(:,:,k),vel.V_mt2D.Data3(:,:,k), x, imax1, jmax1);
    saveas(gcf,sprintf('cont_velocity_Data3_FIG%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...
end

%% Velocity rms (Std Deviation)over the field ******************************
% 1D vector
for i = 1:imax1*jmax1
        vel.u_rms1D.Data3(i) = sqrt(nanmean(nonzeros(vel.u_t1D.Data3(i,1,:).^2)));  
        vel.v_rms1D.Data3(i) = sqrt(nanmean(nonzeros(vel.v_t1D.Data3(i,1,:).^2)));
        vel.uv_rms1D.Data3(i) = sqrt(nanmean(nonzeros(vel.u_t1D.Data3(i,1,:).*vel.v_t1D.Data3(i,1,:))));
end

%2D vector
for i = 1:jmax1
    for j = 1:imax1
        vel.u_rms2D.Data3(i,j) = sqrt(nanmean(nonzeros(vel.u_t2D.Data3(i,j,:).^2))); 
        vel.v_rms2D.Data3(i,j) = sqrt(nanmean(nonzeros(vel.v_t2D.Data3(i,j,:).^2)));
        vel.uv_rms2D.Data3(i,j) = sqrt(nanmean(nonzeros(abs(vel.u_t2D.Data3(i,j,:).* vel.v_t2D.Data3(i,j,:)))));
    end
end
%%
Rmsplot(vel.u_rms2D.Data3, vel.v_rms2D.Data3, x, imax1, jmax1, mean(flowmeter.Vint.Data3))
% saveas(gcf,sprintf('rms_v_Data3.png'));
print('-dpng','-r700',sprintf('rms_v_Data3.png',k));

%% The total rms of the whole field

vel.u_rms_ms2D.Data3 = nanmean(nonzeros(nanmean(nonzeros(vel.u_rms2D.Data3))));
vel.v_rms_ms2D.Data3 = nanmean(nonzeros(nanmean(nonzeros(vel.v_rms2D.Data3))));
vel.uv_rms_ms2D.Data3 = nanmean(nonzeros(nanmean(nonzeros(vel.uv_rms2D.Data3))));

% Velocity variance over the field ****************************************
% 1D vector
for i = 1:imax1*jmax1
    vel.u_var1D.Data3(i) = nanmean(nonzeros(vel.u_t1D.Data3(i,1,:).^2)); 
    vel.v_var1D.Data3(i) = nanmean(nonzeros(vel.v_t1D.Data3(i,1,:).^2));
    vel.uv_var1D.Data3(i) = nanmean(nonzeros(vel.u_t1D.Data3(i,1,:).*vel.v_t1D.Data3(i,1,:)));
end

%2D vector
for i = 1:jmax1
    for j = 1:imax1
        vel.u_var2D.Data3(i,j) = nanmean(nonzeros(vel.u_t2D.Data3(i,j,:).^2));
        vel.v_var2D.Data3(i,j) = nanmean(nonzeros(vel.v_t2D.Data3(i,j,:).^2));
        vel.uv_var2D.Data3(i,j) = nanmean(nonzeros(vel.u_t2D.Data3(i,j,:).*vel.v_t2D.Data3(i,j,:)));
    end
end


vel.u_var_ms2D.Data3 = nanmean(nonzeros(nanmean(nonzeros(vel.u_var2D.Data3))));
vel.v_var_ms2D.Data3 = nanmean(nonzeros(nanmean(nonzeros(vel.v_var2D.Data3))));
vel.uv_var_ms2D.Data3 = nanmean(nonzeros(nanmean(nonzeros(vel.uv_var2D.Data3))));

%% Plotting variances
Varplot(vel, x, imax1, jmax1,mean(flowmeter.Vint.Data3))
print('-dpng','-r700',sprintf('var_V_Data3.png',k));

%% Gamma1 and Gamma2 for the entire field *********************************
[vor.Gamma1.Data3,vor.Gamma2.Data3, vor.Gamma1_mt.Data3, vor.Gamma2_mt.Data3, ...
    vor.Gamma1_mts.Data3, vor.Gamma2_mts.Data3, vor.Gamma1_rms.Data3, vor.Gamma2_rms.Data3,...
    vor.Gamma1_rms_ms.Data3, vor.Gamma2_rms_ms.Data3, vor.N_vortex_G1.Data3, vor.N_vortex_G2.Data3,...
    vor.Gamma2_Shear_mts.Data3, vor.Gamma2_Rotation_mts.Data3, vor.Gamma2_Shear_rms_ms.Data3,...
    vor.Gamma2_Rotation_rms_ms.Data3,vor.N_vortex_G2_Shear.Data3,vor.N_vortex_G2_Rotation.Data3] =...
    Gamma1 (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, 1, jmax1, 1, imax1, kmax1, jmax1, imax1); % jmax1

%% Finding the area of a vortex from Gamma2
[vor.Gamma2_size.Data3, vor.Gamma2_size_mst.Data3, vor.Gamma2_N_vortex.Data3] ...
    = Gamma2Size (vor.Gamma2.Data3,1,1,1,1, kmax1, pix_to_mm,'bed');



% Gamma2_sym(Gamma2_sym>=2/pi)=0; % Setting the rotation or strain vortices
%% Plotting Gamma2
% figure
for k = 20:120
    Gamma12Plot(vor.Gamma1.Data3(:,:,k), vor.Gamma2.Data3(:,:,k), x, imax1, jmax1)
    print('-dpng','-r700',sprintf('Gamma2_Data3_total_T%d.png',k));
    hold off        
end


%% LES Filtering ***********************************************************
for k = 1:kmax1
    %2D 3*3 (Top-Hat)
    U_temp = vel.U_mm2D.Data3(:,:,k);
    U_temp(find(U_temp == 0)) = NaN;
    U_temp_les_nan = nanconv(U_temp, kernel.B_tophat3, 'edge', 'nanout', '2d'); %see nanconv function for more info
    U_temp_les = U_temp_les_nan;
    U_temp_les(isnan(U_temp_les)) = 0; 
    vel.U_les3top2D.Data3(:,:,k) = U_temp_les;
    %1D
    vel.U_les3top1D.Data3(:,1,k) = reshape(vel.U_les3top2D.Data3(:,:,k).',[imax1*jmax1,1]);
    
    %2D
    V_temp = vel.V_mm2D.Data3(:,:,k);
    V_temp(find(V_temp == 0)) = NaN;
    V_temp_les_nan = nanconv(V_temp, kernel.B_tophat3, 'edge', 'nanout', '2d'); %see nanconv function for more info
    V_temp_les = V_temp_les_nan;
    V_temp_les(isnan(V_temp_les)) = 0; 
    vel.V_les3top2D.Data3(:,:,k) = V_temp_les;
    %1D
    vel.V_les3top1D.Data3(:,1,k) = reshape(vel.V_les3top2D.Data3(:,:,k).',[imax1*jmax1,1]);
    
    % 2D 5*5 (Top-Hat)
    U_temp = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix));
    U_temp = vel.U_mm2D.Data3(:,:,k);
    U_temp(find(U_temp == 0)) = NaN;
    U_temp_les_nan = nanconv(U_temp, kernel.B_tophat5, 'edge', 'nanout', '2d'); %see nanconv function for more info
    U_temp_les = U_temp_les_nan;
    U_temp_les(isnan(U_temp_les)) = 0; 
    vel.U_les5top2D.Data3(:,:,k) = U_temp_les;
    % 1D
    vel.U_les5top1D.Data3(:,1,k) = reshape(vel.U_les5top2D.Data3(:,:,k).',[imax1*jmax1,1]);
    
    % 2D
    V_temp = vel.V_mm2D.Data3(:,:,k);
    V_temp(find(V_temp == 0)) = NaN;
    V_temp_les_nan = nanconv(V_temp, kernel.B_tophat5, 'edge', 'nanout', '2d'); %see nanconv function for more info
    V_temp_les = V_temp_les_nan;
    V_temp_les(isnan(V_temp_les)) = 0; 
    vel.V_les5top2D.Data3(:,:,k) = V_temp_les;
    % 1D
    vel.V_les5top1D.Data3(:,1,k) = reshape(vel.V_les5top2D.Data3(:,:,k).',[imax1*jmax1,1]);
    
    % 3*3 (Gasussian)
     % 2D 
    U_temp = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix));    
    U_temp = vel.U_mm2D.Data3(:,:,k);
    U_temp(find(U_temp == 0)) = NaN;
    U_temp_les_nan = nanconv(U_temp, kernel.B_Gauss3, 'edge', 'nanout', '2d'); %see nanconv function for more info
    U_temp_les = U_temp_les_nan;
    U_temp_les(isnan(U_temp_les)) = 0; 
    vel.U_les3gauss2D.Data3(:,:,k) = U_temp_les;
    % 1D
    vel.U_les3gauss1D.Data3(:,1,k) = reshape(vel.U_les3gauss2D.Data3(:,:,k).',[imax1*jmax1,1]);
    
    % 2D
    V_temp = vel.V_mm2D.Data3(:,:,k);
    V_temp(find(V_temp == 0)) = NaN;
    V_temp_les_nan = nanconv(V_temp, kernel.B_Gauss3, 'edge', 'nanout', '2d'); %see nanconv function for more info
    V_temp_les = V_temp_les_nan;
    V_temp_les(isnan(V_temp_les)) = 0; 
    vel.V_les3gauss2D.Data3(:,:,k) = V_temp_les;
    % 1D
    vel.V_les3gauss1D.Data3(:,1,k) = reshape(vel.V_les3gauss2D.Data3(:,:,k).',[imax1*jmax1,1]);
    
    
    % 5*5 (Gasussian)
     % 2D 
    U_temp = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix));    
    U_temp = vel.U_mm2D.Data3(:,:,k);
    U_temp(find(U_temp == 0)) = NaN;
    U_temp_les_nan = nanconv(U_temp, kernel.B_Gauss5, 'edge', 'nanout', '2d'); %see nanconv function for more info
    U_temp_les = U_temp_les_nan;
    U_temp_les(isnan(U_temp_les)) = 0; 
    vel.U_les5gauss2D.Data3(:,:,k) = U_temp_les;
    % 1D
    vel.U_les5gauss1D.Data3(:,1,k) = reshape(vel.U_les5gauss2D.Data3(:,:,k).',[imax1*jmax1,1]);
    
    % 2D
    V_temp = vel.V_mm2D.Data3(:,:,k);
    V_temp(find(V_temp == 0)) = NaN;
    V_temp_les_nan = nanconv(V_temp, kernel.B_Gauss5, 'edge', 'nanout', '2d'); %see nanconv function for more info
    V_temp_les = V_temp_les_nan;
    V_temp_les(isnan(V_temp_les)) = 0; 
    vel.V_les5gauss2D.Data3(:,:,k) = V_temp_les;
    % 1D
    vel.V_les5gauss1D.Data3(:,1,k) = reshape(vel.V_les5gauss2D.Data3(:,:,k).',[imax1*jmax1,1]);
    
end
%% LES Time average
% 2D vector
for i = 1:jmax1
    for j = 1:imax1
        vel.U_les5gauss2D_mt2D_NaN.Data3(i,j) = nanmean(vel.U_les5gauss2D.Data3(i,j,:));vel.V_les5gauss2D_mt2D_NaN.Data3(i,j) = nanmean(vel.V_les5gauss2D.Data3(i,j,:));
    end
end

%% LES fluctuations
% 3*3 (top)
vel.u_les3top2D.Data3 = vel.U_mm2D.Data3 - vel.U_les3top2D.Data3;
vel.v_les3top2D.Data3 = vel.V_mm2D.Data3 - vel.V_les3top2D.Data3;

% 3*3 (Gaussian)
vel.u_les3gauss2D.Data3 = vel.U_mm2D.Data3 - vel.U_les3gauss2D.Data3;
vel.v_les3gauss2D.Data3 = vel.V_mm2D.Data3 - vel.V_les3gauss2D.Data3;

% 5*5 (top)
vel.u_les5top2D.Data3 = vel.U_mm2D.Data3 - vel.U_les5top2D.Data3;
vel.v_les5top2D.Data3 = vel.V_mm2D.Data3 - vel.V_les5top2D.Data3;

% 5*5 (Gaussian)
vel.u_les5gauss2D.Data3 = vel.U_mm2D.Data3 - vel.U_les5gauss2D.Data3;
vel.v_les5gauss2D.Data3 = vel.V_mm2D.Data3 - vel.V_les5gauss2D.Data3;


%% LES Plot
% Large scale
close all
kk = 20;
[Plot.vec_LES.Data3] = LESplot (vel.U_les5gauss2D.Data3(:,:,kk),vel.V_les5gauss2D.Data3(:,:,kk), x, imax1, jmax1);
saveas(Plot.vec_LES.Data3,sprintf('vec_LES5gauss_T13.png'));

figure
% Small scale
[Plot.vec_LES_small.Data3] = LESplot (vel.u_les5gauss2D.Data3(:,:,kk),vel.v_les5gauss2D.Data3(:,:,kk), x, imax1, jmax1);
saveas(Plot.vec_LES_small.Data3,sprintf('vec_LES5gauss_small_T13.png'));


%% Generating time laps of velocity plots 

for k = 10:120
    Plot.vec.Data3 = LESplot(vel.U_les5gauss2D.Data3(:,:,k),vel.V_les5gauss2D.Data3(:,:,k), x, imax1, jmax1);
    print('-dpng','-r700',sprintf('LES5G_Data3_%d.png',k));
%     saveas(Plot.vec.Data3,sprintf('FIG%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...
    hold off
end

% imageFolder2mpeg('Vectors_LES_5_Gauss_Data3','frameRate',11) %  This function creates movie

%% Vorticity of LES field 5*5 Gaussian

[vor.dUdxles5g_2D.Data3, vor.dUdyles5g_2D.Data3, vor.dVdxles5g_2D.Data3, vor.dVdyles5g_2D.Data3]...
    = VGT (x, y, vel.U_les5gauss2D.Data3, vel.V_les5gauss2D.Data3, imax1, jmax1, kmax1, mm_to_pix);

% Vorticity *************************
for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            vor.vorles5g_2D.Data3(i,j,k) = vor.dVdxles5g_2D.Data3(i,j,k)- vor.dUdyles5g_2D.Data3(i,j,k);
        end
    end
end

% Lambda_ci
vor.VGTles5g_2D.Data3 = cell (jmax1,imax1,kmax1);
   
for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            %Velocity gradient tensor
            vor.VGTles5g_2D.Data3{i,j,k} = [vor.dUdxles5g_2D.Data3(i,j,k) vor.dVdxles5g_2D.Data3(i,j,k); vor.dUdyles5g_2D.Data3(i,j,k) vor.dVdyles5g_2D.Data3(i,j,k)];
            
            %Lamda_ci
            lci_temp = eig(cell2mat(vor.VGTles5g_2D.Data3(i,j,k)));
            if isempty(imag(lci_temp(imag(lci_temp)>0)))==1
                vor.lambda_ciles5g_2D.Data3(i,j,k) = NaN;
            else
                vor.lambda_ciles5g_2D.Data3(i,j,k) = imag(lci_temp(imag(lci_temp)>0));
                vor.lambda_crles5g_2D.Data3(i,j,k) = real(lci_temp(imag(lci_temp)>0));
                vor.lambda_crO_lambda_ci_les5g_2D.Data3(i,j,k) = vor.lambda_crles5g_2D.Data3(i,j,k)/vor.lambda_ciles5g_2D.Data3(i,j,k);
            end
         
        end
    end
end

%% Plotiing Lambda_Ci 
close all
h = figure;
kk = 23;
LambdaCiplot (vor.lambda_ciles5g_2D.Data3(:,:,kk), scale.t.T_v_ITS.Data3, x, imax1, jmax1);
hold on
Plot.vec.Data3 = LESplot(vel.U_les3gauss2D.Data3(:,:,kk),vel.V_les3gauss2D.Data3(:,:,kk), x, imax1, jmax1);
saveas(h,'lambda_ciles5g_2D_DelT23.png');

%% Plotiing Lambda_Ci  and LES 
close all
h = figure;
for k = 21:120
    LambdaCiplot(vor.lambda_ciles3g_2D.Data3(:,:,k), 0.0082, x, imax1, jmax1);
    hold on
    Plot.vec.Data3 = LESplot(vel.U_les3gauss2D.Data3(:,:,k),vel.V_les3gauss2D.Data3(:,:,k), x, imax1, jmax1);
    print('-dpng','-r700',sprintf('lambda_ciles3g_Data3_T%d.png',k));
    hold off
end
% imageFolder2mpeg('Lambdaci3g_LES_Data3','frameRate',11)

%% Plotiing Lambda_Cr/Lambda_Ci 3g and LES 
close all
h = figure;
for k = 21:120
    LambdaCiplot(vor.lambda_crO_lambda_ci_les3g_2D.Data3(:,:,k), 0.0082, x, imax1, jmax1);
    hold on
    Plot.vec.Data3 = LESplot(vel.U_les3gauss2D.Data3(:,:,k),vel.V_les3gauss2D.Data3(:,:,k), x, imax1, jmax1);
    print('-dpng','-r700',sprintf('lambda_crO_lambdaci_les3g_Data3_T%d.png',k));
%     saveas(h,sprintf('lambda_ciles3g_T%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...
%     saveas(Plot.vec.Data3,sprintf('FIG%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...
    hold off
end


%% Vorticity of LES field 3*3 Gaussian

[vor.dUdxles3g_2D.Data3, vor.dUdyles3g_2D.Data3, vor.dVdxles3g_2D.Data3, vor.dVdyles3g_2D.Data3]...
    = VGT (x, y, vel.U_les3gauss2D.Data3, vel.V_les3gauss2D.Data3, imax1, jmax1, kmax1, mm_to_pix);

% Vorticity *************************
for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            vor.vorles3g_2D.Data3(i,j,k) = vor.dVdxles3g_2D.Data3(i,j,k)- vor.dUdyles3g_2D.Data3(i,j,k);
        end
    end
end

% Lambda_ci
vor.VGTles3g_2D.Data3 = cell (jmax1,imax1,kmax1);

for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            %Velocity gradient tensor
            vor.VGTles3g_2D.Data3{i,j,k} = [vor.dUdxles3g_2D.Data3(i,j,k) vor.dVdxles3g_2D.Data3(i,j,k); vor.dUdyles3g_2D.Data3(i,j,k) vor.dVdyles3g_2D.Data3(i,j,k)];
            
            %Lamda_ci
            lci_temp = eig(cell2mat(vor.VGTles3g_2D.Data3(i,j,k)));
            if isempty(imag(lci_temp(imag(lci_temp)>0)))==1
                vor.lambda_ciles3g_2D.Data3(i,j,k) = NaN;
            else
                vor.lambda_ciles3g_2D.Data3(i,j,k) = imag(lci_temp(imag(lci_temp)>0));
                vor.lambda_crles3g_2D.Data3(i,j,k) = real(lci_temp(imag(lci_temp)>0));
                vor.lambda_crO_lambda_ci_les3g_2D.Data3(i,j,k) = vor.lambda_crles3g_2D.Data3(i,j,k)/vor.lambda_ciles3g_2D.Data3(i,j,k);
            end
         
        end
    end
end

%% LES 5*5 Time average of vorticity and Lambda_Ci for entire bed 

for i = 1:jmax1
    for j = 1:imax1
        vor.vorles5g_2D_mt.Data3(i,j) = nanmean(vor.vorles5g_2D.Data3(i,j,:)); 
        vor.lambda_ciles5g2D_mt.Data3(i,j) = nanmean(vor.lambda_ciles5g_2D.Data3(i,j,:));
    end
end

% Time and space average of vorticity and Lambda_Ci
vor.vorles5g_2D_mts.Data3 = nanmean(nanmean(vor.vorles5g_2D_mt.Data3));
vor.lambda_ciles5g2D_mts.Data3 = nanmean(nanmean(vor.lambda_ciles5g2D_mt.Data3));

% Vorticity and LambdaCi fluctuation over the field ***********************

%2D vector
for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            vor.vorles5g_fluc_t_2D.Data3(i,j,k) = vor.vorles5g_2D.Data3(i,j,k) - vor.vorles5g_2D_mt.Data3(i,j); 
            vor.lambda_ciles5g2D_fluc_t.Data3(i,j,k) = vor.lambda_ciles5g_2D.Data3(i,j,k) - vor.lambda_ciles5g2D_mt.Data3(i,j);
        end
    end
end


% Rms (Std Deviation)over the field ******************************
%2D vector
for i = 1:jmax1
    for j = 1:imax1
        vor.vorles5g_rms_2D.Data3(i,j) = sqrt(nanmean(vor.vorles5g_fluc_t_2D.Data3(i,j,:).^2)); 
        vor.lambda_ciles5g_rms_2D.Data3(i,j) = sqrt(nanmean(vor.lambda_ciles5g2D_fluc_t.Data3(i,j,:).^2));
    end
end

% Spatial average of rms
vor.vorles5g_rms_ms_2D.Data3 = nanmean(nanmean(vor.vorles5g_rms_2D.Data3));
vor.lambda_ciles5g_rms_ms_2D.Data3 = nanmean(nanmean(vor.lambda_ciles5g_rms_2D.Data3));

vor.N_lambda_ciles5g.Data3= length(find(isnan(vor.lambda_ciles5g_2D.Data3)~=1));


%% Finding the area of a vortex from Lambda_ci LES 5*5
[vor.lambda_ciles5g_2D_size.Data3, vor.lambda_ciles5g_2D_size_mst.Data3, vor.lambda_ciles5g_2D_N_vortex.Data3] ...
    = LambdaCiSize (vor.lambda_ciles5g_2D.Data3,imax1, jmax1, kmax1, pix_to_mm);


%% LES 3*3 Time average of vorticity and Lambda_Ci for entire bed 

for i = 1:jmax1
    for j = 1:imax1
        vor.vorles3g_2D_mt.Data3(i,j) = nanmean(vor.vorles3g_2D.Data3(i,j,:)); 
        vor.lambda_ciles3g2D_mt.Data3(i,j) = nanmean(vor.lambda_ciles3g_2D.Data3(i,j,:));
    end
end

% Time and space average of vorticity and Lambda_Ci
vor.vorles3g_2D_mts.Data3 = nanmean(nanmean(vor.vorles3g_2D_mt.Data3));
vor.lambda_ciles3g2D_mts.Data3 = nanmean(nanmean(vor.lambda_ciles3g2D_mt.Data3));

% Vorticity and LambdaCi fluctuation over the field ***********************

%2D vector
for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            vor.vorles3g_fluc_t_2D.Data3(i,j,k) = vor.vorles3g_2D.Data3(i,j,k) - vor.vorles3g_2D_mt.Data3(i,j); 
            vor.lambda_ciles3g2D_fluc_t.Data3(i,j,k) = vor.lambda_ciles3g_2D.Data3(i,j,k) - vor.lambda_ciles3g2D_mt.Data3(i,j);
        end
    end
end


% Rms (Std Deviation)over the field ******************************
%2D vector
for i = 1:jmax1
    for j = 1:imax1
        vor.vorles3g_rms_2D.Data3(i,j) = sqrt(nanmean(vor.vorles3g_fluc_t_2D.Data3(i,j,:).^2)); 
        vor.lambda_ciles3g_rms_2D.Data3(i,j) = sqrt(nanmean(vor.lambda_ciles3g2D_fluc_t.Data3(i,j,:).^2));
    end
end

% Spatial average of rms
vor.vorles3g_rms_ms_2D.Data3 = nanmean(nanmean(vor.vorles3g_rms_2D.Data3));
vor.lambda_ciles3g_rms_ms_2D.Data3 = nanmean(nanmean(vor.lambda_ciles3g_rms_2D.Data3));

vor.N_lambda_ciles3g.Data3= length(find(isnan(vor.lambda_ciles3g_2D.Data3)~=1));

%% Finding the area of a vortex from Lambda_ci LES 3*3
[vor.lambda_ciles3g_2D_size.Data3, vor.lambda_ciles3g_2D_size_mst.Data3, vor.lambda_ciles3g_2D_N_vortex.Data3, vor.lambda_ciles3g_2D_inst_N_vortex.Data3] ...
    = LambdaCiSize (vor.lambda_ciles3g_2D.Data3,imax1, jmax1, kmax1, pix_to_mm);

vor.lambda_ciles3g_2D_avg_N_vortex.Data3 = sum(vor.lambda_ciles3g_2D_inst_N_vortex.Data3)/3199;

%% Velocity gradient tensor ************************************************

vor.dUdx_2D.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),kmax1);
vor.dUdy_2D.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),kmax1);
vor.dVdx_2D.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),kmax1);
vor.dVdy_2D.Data3 = zeros(length (y.Y_Nodes1_pix),length (x.X_Nodes1_pix),kmax1);

[vor.dUdx_2D.Data3, vor.dUdy_2D.Data3, vor.dVdx_2D.Data3, vor.dVdy_2D.Data3]...
    = VGT (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, imax1, jmax1, kmax1, mm_to_pix);


%% Vorticity *************************
for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            vor.vor_2D.Data3(i,j,k) = vor.dVdx_2D.Data3(i,j,k)- vor.dUdy_2D.Data3(i,j,k);
        end
    end
end


% 1D
vor.vor_1D.Data3 = reshape(vor.vor_2D.Data3,imax1*jmax1,1,kmax1);



%% Lambda_ci of Raw velocities
vor.VGT_2D.Data3 = cell (jmax1,imax1,kmax1);

for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            %Velocity gradient tensor
            vor.VGT_2D.Data3{i,j,k} = [vor.dUdx_2D.Data3(i,j,k) vor.dVdx_2D.Data3(i,j,k); vor.dUdy_2D.Data3(i,j,k) vor.dVdy_2D.Data3(i,j,k)];
            
            %Lamda_ci
            lci_temp = eig(cell2mat(vor.VGT_2D.Data3(i,j,k)));
            if isempty(imag(lci_temp(imag(lci_temp)>0)))==1
                vor.lambda_ci2D.Data3(i,j,k) = NaN;
            else
                vor.lambda_ci2D.Data3(i,j,k) = imag(lci_temp(imag(lci_temp)>0));
            end
         
        end
    end
end

% 1D
vor.lambda_ci1D.Data3 = reshape(vor.lambda_ci2D.Data3,imax1*jmax1,1,kmax1);

% Make NaN to be zero
vor.lambda_ci1D.Data3(find(isnan(vor.lambda_ci1D.Data3)==1)) = 0;

% Time average of vorticity and Lambda_Ci for entire bed
for i = 1:jmax1
    for j = 1:imax1
        vor.vor_2D_mt.Data3(i,j) = nanmean(vor.vor_2D.Data3(i,j,:)); 
        vor.lambda_ci2D_mt.Data3(i,j) = nanmean(vor.lambda_ci2D.Data3(i,j,:));
    end
end

% Time and space average of vorticity and Lambda_Ci
vor.vor_2D_mts.Data3 = nanmean(nanmean(vor.vor_2D_mt.Data3));
vor.lambda_ci2D_mts.Data3 = nanmean(nanmean(vor.lambda_ci2D_mt.Data3));

% Vorticity and LambdaCi fluctuation over the field ***********************

%2D vector
for k = 1:kmax1
    for i = 1:jmax1
        for j = 1:imax1
            vor.vor_fluc_t_2D.Data3(i,j,k) = vor.vor_2D.Data3(i,j,k) - vor.vor_2D_mt.Data3(i,j); 
            vor.lambda_ci2D_fluc_t.Data3(i,j,k) = vor.lambda_ci2D.Data3(i,j,k) - vor.lambda_ci2D_mt.Data3(i,j);
        end
    end
end


% Rms (Std Deviation)over the field ******************************
%2D vector
for i = 1:jmax1
    for j = 1:imax1
        vor.vor_rms_2D.Data3(i,j) = sqrt(nanmean(vor.vor_fluc_t_2D.Data3(i,j,:).^2)); 
        vor.lambda_ci_rms_2D.Data3(i,j) = sqrt(nanmean(vor.lambda_ci2D_fluc_t.Data3(i,j,:).^2));
    end
end

% Spatial average of rms
vor.vor_rms_ms_2D.Data3 = nanmean(nanmean(vor.vor_rms_2D.Data3));
vor.lambda_ci_rms_ms_2D.Data3 = nanmean(nanmean(vor.lambda_ci_rms_2D.Data3));

vor.N_lambda_ci.Data3= length(find(isnan(vor.lambda_ci2D.Data3)~=1));
%% Plotting Lambda_Ci RMS
LambdaCiplot (vor.lambda_ci_rms_2D.Data3, 0.0082, x, imax1, jmax1);
saveas(gcf,'lambda_ci_rms_2D_Data3.png');
%% Plotting Vorticiy 
close all

% vorticity
h = figure;
kk = 5;
vorplot(vor.vor_2D.Data3(:,:,kk), x, imax1, jmax1, specs, vel.Mag_mts2D_NaN.Data3);
hold on
Plot.vec.Data3 = LESplot(vel.U_les3gauss2D.Data3(:,:,kk),vel.V_les3gauss2D.Data3(:,:,kk), x, imax1, jmax1);
saveas(h,'Vor_2D_DelT5.png');

%% vorticity rms
h = figure;
vorplot(vor.vor_rms_2D.Data3, x, imax1, jmax1, specs, vel.Mag_mts2D_NaN.Data3);
saveas(h,'Vor_rms_Data3.png');

%% Plotiing Vorticity  and LES 
close all
h = figure;
for k = 20:120
    
    vorplot(vor.vor_2D.Data3(:,:,k), x, imax1, jmax1, specs, vel.Mag_mts2D_NaN.Data3);
    hold on
    Plot.vec.Data3 = LESplot(vel.U_les3gauss2D.Data3(:,:,k),vel.V_les3gauss2D.Data3(:,:,k), x, imax1, jmax1);
    saveas(h,sprintf('Vor_LES_Data3_%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

%     saveas(Plot.vec.Data3,sprintf('FIG%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...
    hold off
end
% imageFolder2mpeg('Vor_LES3_Data3','frameRate',11)

%% Circulation over the entire bed

[vor.circ_vel_bed.Data3, vor.circ_vor_bed.Data3, vor.circ_vor_mt_bed.Data3, ...
    vor.circ_vel_mt_bed.Data3, vor.circ_vor_rms_bed.Data3, vor.circ_vel_rms_bed.Data3] = ...
    circulation (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, vor.vor_2D.Data3, 1, jmax1, ...
    1, imax1, kmax1, 77, mm_to_pix);

%% Autocorrelation of fluctuating velocity for integral time scale for the bed
lag.t.Data3 = 100;
% Preallocation
for i = 1:jmax1
    for j = 1:imax1
        vel.autocorr_t_u.Data3{i,j} = zeros(1,lag.t.Data3+1);
        vel.autocorr_t_v.Data3{i,j} = zeros(1,lag.t.Data3+1);
    end
end
for i = 1:jmax1
    for j = 1:imax1
        vel.autocorr_t_u.Data3{i,j} = Correlation(vel.U_mm2D.Data3(i,j,:), vel.U_mm2D.Data3(i,j,:), lag.t.Data3);
        vel.autocorr_t_v.Data3{i,j} = Correlation(vel.V_mm2D.Data3(i,j,:), vel.V_mm2D.Data3(i,j,:), lag.t.Data3);
    end
end
%% Integral time scale for the whole bed 
scale.t.id_u.Data3 = zeros(jmax1,imax1);
scale.t.id_v.Data3 = zeros(jmax1,imax1);
scale.t.T_u2D.Data3 = zeros(jmax1,imax1);
scale.t.T_v2D.Data3 = zeros(jmax1,imax1);

% LAG = celltomat(struct2cell(lag.Data3));
for i = 1:jmax1
    for j = 1:imax1
        scale.t.autocorr_temp_u_t.Data3 = cell2mat(vel.autocorr_t_u.Data3(i,j));
        scale.t.autocorr_temp_v_t.Data3 = cell2mat(vel.autocorr_t_v.Data3(i,j));
        if any(scale.t.autocorr_temp_u_t.Data3<0) & scale.t.autocorr_temp_u_t.Data3(2)>0
            scale.t.id_u_temp.Data3 = find(scale.t.autocorr_temp_u_t.Data3 < 0,1,'first'); %first zero-crossing
            scale.t.id_u.Data3(i,j) = scale.t.id_u_temp.Data3;
            scale.t.T_u2D.Data3(i,j) = Del_T.Data3*((1/2)*nansum(scale.t.autocorr_temp_u_t.Data3(1)) + ...
                (1/2)*nansum(scale.t.autocorr_temp_u_t.Data3(scale.t.id_u_temp.Data3-1)) + ...
                nansum(scale.t.autocorr_temp_u_t.Data3(2:scale.t.id_u_temp.Data3-2)));
        else
            scale.t.T_u2D.Data3(i,j) = NaN;
        end
        if any(scale.t.autocorr_temp_v_t.Data3<0) & scale.t.autocorr_temp_v_t.Data3(2)>0
            scale.t.id_v_temp.Data3 = find(scale.t.autocorr_temp_v_t.Data3 < 0,1,'first'); %first zero-crossing
            scale.t.id_v.Data3(i,j) = scale.t.id_v_temp.Data3;
            scale.t.T_v2D.Data3(i,j) = Del_T.Data3*((1/2)*nansum(scale.t.autocorr_temp_v_t.Data3(1)) +...
                (1/2)*nansum(scale.t.autocorr_temp_v_t.Data3(scale.t.id_v_temp.Data3-1)) + ...
                nansum(scale.t.autocorr_temp_v_t.Data3(2:scale.t.id_v_temp.Data3-2)));
        else
            scale.t.T_v2D.Data3(i,j) = NaN;
        end
    end
end

scale.t.T_u_ITS.Data3 = nanmean(nanmean(scale.t.T_u2D.Data3));
scale.t.T_v_ITS.Data3 = nanmean(nanmean(scale.t.T_v2D.Data3));

%% Spatial autocorrelation for Integral length scale for the whole bed (v_fluctuation)

lag.L_v.Data3.Bed = 61;

% Preallocation
scale.autocorr_L_V.Data3 = zeros(lag.L_v.Data3.Bed+1,imax1,kmax1);

% Be carefull it takes 40 min
for k = 1:kmax1
    for j = 1:imax1
        scale.autocorr_L_V.Data3(:,j,k) = Correlation(vel.V_mm2D.Data3(:,j,k),vel.V_mm2D.Data3(:,j,k),lag.L_v.Data3.Bed);
    end
end
%% continue from top to find mean autocorrelation of the whole bed to find (v_fluctuation)
scale.autocorr_L_V_ms.Data3 = zeros(lag.L_v.Data3.Bed+1,kmax1);
%
for k = 1:kmax1
    for i = 1:lag.L_v.Data3.Bed+1
        scale.autocorr_L_V_ms.Data3 (i,k) = nanmean(scale.autocorr_L_V.Data3(i,:,k));
    end
end

for i = 1:lag.L_v.Data3.Bed+1
    scale.autocorr_L_V_mst.Data3(i) = nanmean(scale.autocorr_L_V_ms.Data3(i,:));
end


%% Integral length scale for the whole bed (v_fluctuation)
x.Del_x = (x.X_Nodes1_pix(2)- x.X_Nodes1_pix(1))*mm_to_pix;
y.Del_y = (y.Y_Nodes1_pix(2)- y.Y_Nodes1_pix(1))*mm_to_pix;

scale.L.id_v.Data3 = find(scale.autocorr_L_V_mst.Data3 < 0,1,'first'); %first zero-crossing
scale.L.ILS_Y.Data3.All_bed = y.Del_y*((1/2)*nansum(scale.autocorr_L_V_mst.Data3(1)) + (1/2)*nansum(scale.autocorr_L_V_mst.Data3(scale.L.id_v.Data3-1)) + nansum(scale.autocorr_L_V_mst.Data3(2:scale.L.id_v.Data3-2)));

%% Spatial autocorrelation for Integral length scale for the whole bed (u_fluctuation)

lag.L_u.Data3.Bed = 76;

% Preallocation
scale.autocorr_L_U.Data3 = zeros(jmax1,lag.L_u.Data3.Bed+1,kmax1);

% Be carefull it takes 40 min
for k = 1:kmax1
    for i = 1:jmax1
        scale.autocorr_L_U.Data3(i,:,k) = Correlation(vel.U_mm2D.Data3(i,:,k),vel.U_mm2D.Data3(i,:,k),lag.L_u.Data3.Bed);
    end
end
%% continue from top to find mean autocorrelation of the whole bed to find (u_fluctuation)
scale.autocorr_L_U_ms.Data3 = zeros(lag.L_u.Data3.Bed+1,kmax1);
%
for k = 1:kmax1
    for j = 1:lag.L_u.Data3.Bed+1
        scale.autocorr_L_U_ms.Data3 (j,k) = nanmean(scale.autocorr_L_U.Data3(:,j,k));
    end
end

for j = 1:lag.L_u.Data3.Bed+1
    scale.autocorr_L_U_mst.Data3(j) = nanmean(scale.autocorr_L_U_ms.Data3(j,:));
end

%% Plotting Autocorrelation
AutocorrPlot (x,y, scale, lag, specs,pix_to_mm)
saveas(gcf,'Autocorr_Data3.fig')
saveas(gcf,'Autocorr_Data3.png')
saveas(gcf,'Autocorr_Data3.eps')
%% Integral length scale for the whole bed (u_fluctuation)
x.Del_x = (x.X_Nodes1_pix(2)- x.X_Nodes1_pix(1))*mm_to_pix;
y.Del_y = (y.Y_Nodes1_pix(2)- y.Y_Nodes1_pix(1))*mm_to_pix;

scale.L.id_u.Data3 = find(scale.autocorr_L_U_mst.Data3 < 0,1,'first'); %first zero-crossing
scale.L.ILS_X.Data3.All_bed = y.Del_y*((1/2)*nansum(scale.autocorr_L_U_mst.Data3(1)) + (1/2)*nansum(scale.autocorr_L_U_mst.Data3(scale.L.id_u.Data3-1)) + nansum(scale.autocorr_L_U_mst.Data3(2:scale.L.id_u.Data3-2)));


%% Spatial autocovariance for Energy spectra for the whole bed (v_fluctuation)
lag.L_v.Data3.Bed = 62;

% Preallocation
scale.autocovar_L_V.Data3 = zeros(lag.L_v.Data3.Bed+1,imax1,kmax1);

% Be carefull it takes 40 min
for k = 1:kmax1
    for j = 1:imax1
        scale.autocovar_L_V.Data3(:,j,k) = CorrelationNotNormalized(vel.V_mm2D.Data3(:,j,k),vel.V_mm2D.Data3(:,j,k),lag.L_v.Data3.Bed);
    end
end

%% Spatial autocovariance for Energy spectra for the whole bed (u_fluctuation)
lag.L_u.Data3.Bed = 76;

% Preallocation
scale.autocovar_L_U.Data3 = zeros(jmax1,lag.L_u.Data3.Bed+1,kmax1);

% Be carefull it takes 40 min
for k = 1:kmax1
    for i = 1:jmax1
        scale.autocovar_L_U.Data3(i,:,k) = CorrelationNotNormalized(vel.U_mm2D.Data3(i,:,k),vel.U_mm2D.Data3(i,:,k),lag.L_u.Data3.Bed);
    end
end

%% Energy spectra

% time average of autocovariance V
for j = 1:imax1
    for i = 1:lag.L_v.Data3.Bed+1
        scale.autocovar_L_V_mt.Data3 (i,j) = nanmean(scale.autocovar_L_V.Data3(i,j,:));
    end
end
[scale.nx_v_autocovar.Data3, scale.ny_v_autocovar.Data3] = size(scale.autocovar_L_V_mt.Data3) ;
scale.autocovar_L_V_mt.Data3(isnan(scale.autocovar_L_V_mt.Data3))=0;
dx = 0.01; % Assume dx = 0.01 m
scale.autocovar_L_V_mts.Data3 = mean(scale.autocovar_L_V_mt.Data3,2);

for i = 1:scale.ny_v_autocovar.Data3
    turb.E_spec_v.Data3(:,i) = fft(scale.autocovar_L_V_mt.Data3(:,i));
end
turb.E_spec_v_ms.Data3 = fft(scale.autocovar_L_V_mts.Data3); % Mean spectra
scale.wave_v.Data3 = (-scale.nx_v_autocovar.Data3/2+1 : scale.nx_v_autocovar.Data3/2)*2 * pi/(specs.geom.DH*1000);
figure
for i = 1:scale.ny_v_autocovar.Data3
    loglog(scale.wave_v.Data3(scale.nx_v_autocovar.Data3/2+1:end)',abs(real(scale.autocovar_L_V_mt.Data3(1:scale.nx_v_autocovar.Data3/2,i))));
    hold on
end

figure
semilogx(scale.wave_v.Data3(scale.nx_v_autocovar.Data3/2+1:end),abs(real(scale.autocovar_L_V_mts.Data3(1:scale.nx_v_autocovar.Data3/2)))/(vel.v_var_ms2D.Data3*specs.geom.DH*1000),'-.sk','MarkerFaceColor','k','MarkerSize',4,'linewidth',1)
h1 = legend('\langle E(\kappa) \rangle');
xlabel('\kappa [rad/mm]','fontsize',14,'FontName','Times New Roman','fontsize',15);
ylabel('\langleE(\kappa)\rangle/(D_H\langlev''^2\rangle)[1/rad]','fontsize',14,'FontName','Times New Roman','fontsize',15);
set(gca,'FontName','Times New Roman','FontSize',14)
xx = scale.wave_v.Data3(33:61);
yy = 10^1*xx.^(-5/3);
hold on 
semilogx(xx,yy,'--k','linewidth',1);
str=get(h1,'string');     
new_leg = '-5/3^{rd} line';
h2 = legend([str new_leg],'northwest');  % concatenate the new and the previous legend
set(gca,'FontName','Times New Roman','FontSize',14)
axis([1 15 0 0.4])

% time average of autocovariance U
for j = 1:lag.L_u.Data3.Bed+1
    for i = 1:jmax1
        scale.autocovar_L_U_mt.Data3 (i,j) = nanmean(scale.autocovar_L_U.Data3(i,j,:));
    end
end
[scale.nx_u_autocovar.Data3, scale.ny_u_autocovar.Data3] = size(scale.autocovar_L_U_mt.Data3) ;
scale.autocovar_L_U_mt.Data3(isnan(scale.autocovar_L_U_mt.Data3))=0;
scale.autocovar_L_U_mts.Data3 = mean(scale.autocovar_L_U_mt.Data3,2);

for j = 1:scale.nx_u_autocovar.Data3
    turb.E_spec_u.Data3(j,:) = fft(scale.autocovar_L_U_mt.Data3(j,:));
end
turb.E_spec_u_ms.Data3 = fft(scale.autocovar_L_U_mts.Data3); % Mean spectra
scale.wave_u.Data3 = (-scale.ny_u_autocovar.Data3/2+1 : scale.ny_u_autocovar.Data3/2)*2 * pi/(scale.L.ILS_X.Data3.All_bed);

figure
for j = 1:scale.nx_u_autocovar.Data3
    loglog(scale.wave_u.Data3(scale.ny_u_autocovar.Data3/2+1:end)',abs(real(scale.autocovar_L_U_mt.Data3 (j,1:scale.ny_u_autocovar.Data3/2))));
    hold on
end

figure
loglog(scale.wave_u.Data3(scale.ny_u_autocovar.Data3/2+1:end),abs(real(scale.autocovar_L_U_mts.Data3(1:scale.ny_u_autocovar.Data3/2)))/(vel.u_var_ms2D.Data3*specs.geom.DH*1000),'-.sk','MarkerFaceColor','k','MarkerSize',4,'linewidth',1)
h1 = legend('\langle E(\kappa) \rangle');
xlabel('\kappa [rad/mm]','fontsize',14,'FontName','Times New Roman','fontsize',15);
ylabel('\langleE(\kappa)\rangle/(D_H\langlev''^2\rangle)[1/rad]','fontsize',14,'FontName','Times New Roman','fontsize',15);
set(gca,'FontName','Times New Roman','FontSize',14)
xx = scale.wave_u.Data3(33:61);
yy = 0.01^1*xx.^(-5/3);
hold on 
loglog(xx,yy,'--k','linewidth',1);
str=get(h1,'string');     
new_leg = '-5/3^{rd} line';
h2 = legend([str new_leg],'northwest');  % concatenate the new and the previous legend
set(gca,'FontName','Times New Roman','FontSize',14)
axis([3 100 0 0.001])



%% Production for the whole bed
vor.dUdx_mt.Data3 = zeros(jmax1,imax1);
vor.dUdy_mt.Data3 = zeros(jmax1,imax1);
vor.dVdx_mt.Data3 = zeros(jmax1,imax1);
vor.dVdy_mt.Data3 = zeros(jmax1,imax1);

[vor.dUdx_mt.Data3, vor.dUdy_mt.Data3, vor.dVdx_mt.Data3, vor.dVdy_mt.Data3] = VGT (x, y, vel.U_mt2D.Data3, vel.V_mt2D.Data3, imax1, jmax1, 1, mm_to_pix);
turb.Production.Sxx.Data3 = (1/2)*(vor.dUdx_mt.Data3 + vor.dUdx_mt.Data3);
turb.Production.Syy.Data3 = (1/2)*(vor.dVdy_mt.Data3 + vor.dVdy_mt.Data3);
turb.Production.Szz.Data3 = -(turb.Production.Sxx.Data3 + turb.Production.Syy.Data3); %based on continuity
turb.Production.Sxy.Data3 = (1/2)*(vor.dUdy_mt.Data3 + vor.dVdx_mt.Data3);

turb.Production.P_linear.Data3 = - nanmean(nanmean((vel.u_var2D.Data3).*(turb.Production.Sxx.Data3)))...
    - nanmean(nanmean((vel.v_var2D.Data3).*(turb.Production.Syy.Data3)))...
    - nanmean(nanmean((vel.u_var2D.Data3).*(turb.Production.Szz.Data3)));        % Look at vishal's article: Scale Estimation

turb.Production.P_shear.Data3 = -3*2*nanmean(nanmean((vel.uv_var2D.Data3).*(turb.Production.Sxy.Data3)));

%% Kinetic energy advection for the whole bed (see page 63 of Lumely)

[vor.du_var_dx.Data3, vor.du_var_dy.Data3, vor.dv_var_dx.Data3, vor.dv_var_dy.Data3] = VGT (x, y, vel.u_var2D.Data3, vel.v_var2D.Data3, imax1, jmax1, 1, mm_to_pix);
turb.KE_advection.Data3 = 2*nanmean(nanmean((vel.U_mt2D.Data3).*((1/2)*vor.du_var_dx.Data3)))...
    + nanmean(nanmean((vel.U_mt2D.Data3).*((1/2)*vor.dv_var_dx.Data3)))...
    + 2*nanmean(nanmean((vel.V_mt2D.Data3).*((1/2)*vor.du_var_dy.Data3)))...
    + nanmean(nanmean((vel.V_mt2D.Data3).*((1/2)*vor.dv_var_dy.Data3)));

%% Net gradient transport

[vor.dudx.Data3, vor.dudy.Data3, vor.dvdx.Data3, vor.dvdy.Data3] = VGT (x, y, vel.u_t2D.Data3, vel.v_t2D.Data3, imax1, jmax1, kmax1, mm_to_pix);
turb.Production.sxx.Data3 = (1/2)*(vor.dudx.Data3 + vor.dudx.Data3);
turb.Production.sxy.Data3 = (1/2)*(vor.dudy.Data3 + vor.dvdx.Data3);
turb.Production.syy.Data3 = (1/2)*(vor.dvdy.Data3 + vor.dvdy.Data3);

% time averaging the mean terms
for i = 1:jmax1
    for j = 1:imax1
        turb.Net_grad.viscous_x.Data3(i,j) = 2*specs.fluid.kvisc*(nanmean(vel.u_t2D.Data3(i,j,:).*turb.Production.sxx.Data3(i,j,:))...
            + nanmean(vel.v_t2D.Data3(i,j,:).*turb.Production.sxy.Data3(i,j,:)) + nanmean(vel.u_t2D.Data3(i,j,:).*turb.Production.sxy.Data3(i,j,:)));
        turb.Net_grad.viscous_y.Data3(i,j) = 2*specs.fluid.kvisc*(nanmean(vel.v_t2D.Data3(i,j,:).*turb.Production.syy.Data3(i,j,:))...
            + 2*nanmean(vel.u_t2D.Data3(i,j,:).*turb.Production.sxy.Data3(i,j,:)));
        turb.Net_grad.viscous_z.Data3(i,j) = 2*specs.fluid.kvisc*(nanmean(vel.u_t2D.Data3(i,j,:).*turb.Production.sxy.Data3(i,j,:))...
            + nanmean(vel.v_t2D.Data3(i,j,:).*turb.Production.sxy.Data3(i,j,:)) + nanmean(vel.u_t2D.Data3(i,j,:).*turb.Production.sxx.Data3(i,j,:)));
        turb.Net_grad.fluc_x.Data3(i,j) = (1/2)*(2*nanmean(vel.u_t2D.Data3(i,j,:).^3) + nanmean(vel.v_t2D.Data3(i,j,:).^2.*vel.u_t2D.Data3(i,j,:)));
        turb.Net_grad.fluc_y.Data3(i,j) = (1/2)*(2*nanmean(vel.u_t2D.Data3(i,j,:).^2.*vel.v_t2D.Data3(i,j,:)) + nanmean(vel.v_t2D.Data3(i,j,:).^3));
        turb.Net_grad.fluc_z.Data3(i,j) = (1/2)*(2*nanmean(vel.u_t2D.Data3(i,j,:).^3) + nanmean(vel.v_t2D.Data3(i,j,:).^2.*vel.u_t2D.Data3(i,j,:)));
    end
end

turb.Net_grad.grad_x.Data3 = turb.Net_grad.fluc_x.Data3 - turb.Net_grad.viscous_x.Data3;
turb.Net_grad.grad_y.Data3 = turb.Net_grad.fluc_y.Data3 - turb.Net_grad.viscous_y.Data3;
turb.Net_grad.grad_z.Data3 = turb.Net_grad.fluc_z.Data3 - turb.Net_grad.viscous_z.Data3;

[turb.d_Netgrad_x_dx.Data3, turb.d_Netgrad_x_dy.Data3, turb.d_Netgrad_y_dx.Data3, turb.d_Netgrad_y_dy.Data3] = ...
    VGT (x, y, turb.Net_grad.grad_x.Data3, turb.Net_grad.grad_y.Data3, imax1, jmax1, 1, mm_to_pix);

[turb.d_Netgrad_z_dx.Data3, turb.d_Netgrad_z_dy.Data3, turb.d_Netgrad_y_dx.Data3, turb.d_Netgrad_y_dy.Data3] = ...
    VGT (x, y, turb.Net_grad.grad_z.Data3, turb.Net_grad.grad_y.Data3, imax1, jmax1, 1, mm_to_pix);

turb.Net_grad.total_trans.Data3 = nanmean(nanmean(turb.d_Netgrad_x_dx.Data3 + turb.d_Netgrad_y_dy.Data3 + turb.d_Netgrad_z_dx.Data3));

%% Turbulent kinetic energy
turb.TKE_2D_bed.Data3 = (1/2)*(2*vel.u_var2D.Data3 + vel.v_var2D.Data3);
turb.TKE_ms_bed.Data3 = nanmean(nanmean(turb.TKE_2D_bed.Data3));

%% Tortuosity
dispers.Tortuosity_bed.Data3 = nansum(nansum(nansum(vel.Mag_mm2D.Data3)))./nansum(nansum(nansum(vel.V_mm2D.Data3)));


%% Eulerian variance

% 2D field
vel.U_var2D.Data3 = (abs(vel.U_mt2D_NaN.Data3) - abs(vel.Mag_mts2D_NaN.Data3)).^2;
vel.V_var2D.Data3 = (abs(vel.V_mt2D_NaN.Data3) - abs(vel.Mag_mts2D_NaN.Data3)).^2;

% Averaged non-dimensionalized Eulerian variance based on Vishal's turbulent paper
vel.U_var_mts_Nondim.Data3 = nanmean(nanmean(abs(vel.U_var2D.Data3)))/vel.Mag_mts2D_NaN.Data3^2;
vel.V_var_mts_Nondim.Data3 = nanmean(nanmean(abs(vel.V_var2D.Data3)))/vel.Mag_mts2D_NaN.Data3^2;

%% Lagrangian velocity variance

vel.U_Lagrang_var2D.Data3 = nanmean(nanmean(vel.u_var2D.Data3)) + nanmean(nanmean(vel.U_var2D.Data3));
vel.V_Lagrang_var2D.Data3 = nanmean(nanmean(vel.v_var2D.Data3)) + nanmean(nanmean(vel.V_var2D.Data3));


%% Spatial and temporal Average of each pore
% Pore A ******************************************************************
imin.Data3.A = 15; 
imax.Data3.A = 27;
jmin.Data3.A = 59;
jmax.Data3.A = 69;
nodes.Data3.A = 11; % 11, 15
[vel.U_mt_pore.Data3.A,vel.V_mt_pore.Data3.A,vel.Mag_mt_pore.Data3.A,...
    vel.U_mts_pore.Data3.A,vel.V_mts_pore.Data3.A, ...
    vel.Mag_mts_pore.Data3.A,vel.u_rms_pore.Data3.A, ...
    vel.v_rms_pore.Data3.A, vel.uv_rms_pore.Data3.A, vel.u_rms_ms_pore.Data3.A, ...
    vel.v_rms_ms_pore.Data3.A, vel.uv_rms_ms_pore.Data3.A, vel.u_var_pore.Data3.A,...
    vel.v_var_pore.Data3.A, vel.uv_var_pore.Data3.A, vel.u_var_ms_pore.Data3.A,...
    vel.v_var_ms_pore.Data3.A, vel.uv_var_ms_pore.Data3.A,...
    N_vectors.Data3.A] = poreMeanVelocity (x, y, vel.U_mm2D.Data3, ...
    vel.V_mm2D.Data3, vel.Mag_mm2D.Data3, imin.Data3.A, imax.Data3.A, ...
    jmin.Data3.A, jmax.Data3.A, kmax1, mm_to_pix);


%%
[vor.vor_mt_pore.Data3.A, vor.lambda_ci_mt_pore.Data3.A, vor.vor_mts_pore.Data3.A, ...
    vor.lambda_ci_mts_pore.Data3.A, vor.vor_rms_pore.Data3.A, ...
    vor.lambda_ci_rms_pore.Data3.A, vor.vor_rms_ms_pore.Data3.A, ...
    vor.lambda_ci_rms_ms_pore.Data3.A, vor.N_lambda_ci_pore.Data3.A] ...
    = poreMeanVorticity (x, y, vor.vor_2D.Data3, vor.lambda_ci2D.Data3, imin.Data3.A, imax.Data3.A, ...
    jmin.Data3.A, jmax.Data3.A, kmax1, mm_to_pix);
%% vortex size 3*3 and 5*5
[vor.lambda_ciles3g_2D_size_Pore.Data3.A, vor.lambda_ciles3g_2D_size_mst_Pore.Data3.A, vor.lambda_ciles3g_2D_N_vortex_Pore.Data3.A] ...
    = LambdaCiSizePore (vor.lambda_ciles3g_2D.Data3,imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A, kmax1, pix_to_mm);

[vor.lambda_ciles5g_2D_size_Pore.Data3.A, vor.lambda_ciles5g_2D_size_mst_Pore.Data3.A, vor.lambda_ciles5g_2D_N_vortex_Pore.Data3.A] ...
    = LambdaCiSizePore (vor.lambda_ciles5g_2D.Data3,imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A, kmax1, pix_to_mm);
%%
[vor.circ_vel_pore.Data3.A, vor.circ_vor_pore.Data3.A, vor.circ_vor_mt_pore.Data3.A, ...
    vor.circ_vel_mt_pore.Data3.A, vor.circ_vor_rms_pore.Data3.A, vor.circ_vel_rms_pore.Data3.A] = ...
    circulation (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, vor.vor_2D.Data3, imin.Data3.A, imax.Data3.A, ...
    jmin.Data3.A, jmax.Data3.A, kmax1, nodes.Data3.A, mm_to_pix);

%LES 5*5
[vor.vorles5g_mt_pore.Data3.A, vor.lambda_ciles5g_mt_pore.Data3.A, vor.vorles5g_mts_pore.Data3.A, ...
    vor.lambda_ciles5g_mts_pore.Data3.A, vor.vorles5g_rms_pore.Data3.A, ...
    vor.lambda_ciles5g_rms_pore.Data3.A, vor.vorles5g_rms_ms_pore.Data3.A, ...
    vor.lambda_ciles5g_rms_ms_pore.Data3.A, vor.N_lambda_ciles5g_pore.Data3.A] ...
    = poreMeanVorticity (x, y, vor.vorles5g_2D.Data3, vor.lambda_ciles5g_2D.Data3, imin.Data3.A, imax.Data3.A, ...
    jmin.Data3.A, jmax.Data3.A, kmax1, mm_to_pix);

%LES 3*3
[vor.vorles3g_mt_pore.Data3.A, vor.lambda_ciles3g_mt_pore.Data3.A, vor.vorles3g_mts_pore.Data3.A, ...
    vor.lambda_ciles3g_mts_pore.Data3.A, vor.vorles3g_rms_pore.Data3.A, ...
    vor.lambda_ciles3g_rms_pore.Data3.A, vor.vorles3g_rms_ms_pore.Data3.A, ...
    vor.lambda_ciles3g_rms_ms_pore.Data3.A, vor.N_lambda_ciles3g_pore.Data3.A] ...
    = poreMeanVorticity (x, y, vor.vorles3g_2D.Data3, vor.lambda_ciles3g_2D.Data3, imin.Data3.A, imax.Data3.A, ...
    jmin.Data3.A, jmax.Data3.A, kmax1, mm_to_pix);

%% Gamma1 and Gamma2
[vor.Gamma1_pore.Data3.A,vor.Gamma2_pore.Data3.A, vor.Gamma1_mt_pore.Data3.A, vor.Gamma2_mt_pore.Data3.A,...
    vor.Gamma1_mts_pore.Data3.A, vor.Gamma2_mts_pore.Data3.A, vor.Gamma1_rms_pore.Data3.A, ...
    vor.Gamma2_rms_pore.Data3.A, vor.Gamma1_rms_ms_pore.Data3.A, vor.Gamma2_rms_ms_pore.Data3.A,...
    vor.N_Gamma1.Data3.A, vor.N_Gamma2.Data3.A, vor.Gamma2_Shear_mts_pore.Data3.A, vor.Gamma2_Rotation_mts_pore.Data3.A,...
    vor.Gamma2_Shear_rms_ms_pore.Data3.A, vor.Gamma2_Rotation_rms_ms_pore.Data3.A,vor.N_vortex_G2_Shear_pore.Data3.A,...
    vor.N_vortex_G2_Rotation_pore.Data3.A] = Gamma1 (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, imin.Data3.A, imax.Data3.A, ...
    jmin.Data3.A, jmax.Data3.A, kmax1, jmax1, imax1);

%% Vortex size Gamma2

[vor.Gamma2_size.Data3, vor.Gamma2_size_mst.Data3, vor.Gamma2_N_vortex.Data3,...
    vor.Gamma2_size_Pore.Data3.A, vor.Gamma2_size_mst_Pore.Data3.A, vor.Gamma2_N_vortex_Pore.Data3.A] ...
    = Gamma2Size (vor.Gamma2.Data3, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A, kmax1, pix_to_mm, 'pore');

%% Gamma1 and Gamma2
h = figure;
for k = 10:10
    Gamma2PlotPore (vor.Gamma1_pore.Data3.A(:,:,k),vor.Gamma2_pore.Data3.A(:,:,k), x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
    saveas(h,sprintf('Gamma2_Rotation_Data3_PoreA_%d.png',k));
end

% imageFolder2mpeg('Gamma2_Rotation_PoreA_Data3','frameRate',11)
%%
% Integral time scale
[scale.t.T_u2D_pore.Data3.A, scale.t.T_v2D_pore.Data3.A, scale.t.T_u_ITS_pore.Data3.A, scale.t.T_v_ITS_pore.Data3.A] = ...
    IntTimeScalePore (scale.t.T_u2D.Data3, scale.t.T_v2D.Data3, imin.Data3.A, imax.Data3.A,...
    jmin.Data3.A, jmax.Data3.A);

% Integral length scale
lag.L.Data3.A = imax.Data3.A - imin.Data3.A ; 
[scale.autocorr_L_V_mst_pore.Data3.A,scale.ILS.Data3.A] = ...
    IntLengthScalePore (vel.V_mm2D.Data3, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A, kmax1,lag.L.Data3.A, y.Del_y);

%% TKE budget
% Production, TKE Advection
[turb.Production.P_linear_pore.Data3.A, turb.Production.P_shear_pore.Data3.A, turb.KE_advection_pore.Data3.A, ...
    turb.Net_grad.grad_x_pore.Data3.A, turb.Net_grad.grad_y_pore.Data3.A, turb.Net_grad.grad_z_pore.Data3.A] = ...
    TKEbudgetPore (vel.U_mt2D.Data3, vel.V_mt2D.Data3, turb.Production.Sxx.Data3,...
    turb.Production.Syy.Data3, turb.Production.Szz.Data3, turb.Production.Sxy.Data3, vel.u_var2D.Data3, ...
    vel.v_var2D.Data3, vel.uv_var2D.Data3, vor.du_var_dx.Data3, vor.du_var_dy.Data3, ...
    vor.dv_var_dx.Data3, vor.dv_var_dy.Data3, turb.Production.sxx.Data3, turb.Production.sxy.Data3, ...
    turb.Production.syy.Data3, vel.u_t2D.Data3, vel.v_t2D.Data3, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A,specs.fluid.kvisc);

% Net Gradient transport
[turb.d_Netgrad_x_dx_pore.Data3.A, turb.d_Netgrad_x_dy_pore.Data3.A, turb.d_Netgrad_y_dx_pore.Data3.A, turb.d_Netgrad_y_dy_pore.Data3.A] = ...
    VGT (x, y, turb.Net_grad.grad_x_pore.Data3.A, turb.Net_grad.grad_y_pore.Data3.A, jmax.Data3.A-jmin.Data3.A+1, imax.Data3.A-imin.Data3.A+1, 1, mm_to_pix);

[turb.d_Netgrad_z_dx_pore.Data3.A, turb.d_Netgrad_z_dy_pore.Data3.A, turb.d_Netgrad_y_dx_pore.Data3.A, turb.d_Netgrad_y_dy_pore.Data3.A] = ...
    VGT (x, y, turb.Net_grad.grad_z_pore.Data3.A, turb.Net_grad.grad_y_pore.Data3.A, jmax.Data3.A-jmin.Data3.A+1, imax.Data3.A-imin.Data3.A+1, 1, mm_to_pix);

turb.Net_grad.total_trans_pore.Data3.A = nanmean(nanmean(turb.d_Netgrad_x_dx_pore.Data3.A + turb.d_Netgrad_y_dy_pore.Data3.A + turb.d_Netgrad_z_dx_pore.Data3.A));

% TKE
% Turbulent kinetic energy
turb.TKE_2D_pore.Data3.A = (1/2)*(2*vel.u_var2D.Data3(imin.Data3.A:imax.Data3.A,jmin.Data3.A:jmax.Data3.A) + vel.v_var2D.Data3(imin.Data3.A:imax.Data3.A,jmin.Data3.A:jmax.Data3.A));
turb.TKE_ms_pore.Data3.A = nanmean(nanmean(turb.TKE_2D_pore.Data3.A));

% Tortuosity
dispers.Tortuosity_pore.Data3.A = nansum(nansum(nansum(vel.Mag_mm2D.Data3(imin.Data3.A:imax.Data3.A,jmin.Data3.A:jmax.Data3.A))))./nansum(nansum(nansum(vel.V_mm2D.Data3(imin.Data3.A:imax.Data3.A,jmin.Data3.A:jmax.Data3.A))));


% Eulerian variance
% 2D field
vel.U_var2D_pore.Data3.A = (abs(vel.U_mt2D_NaN.Data3(imin.Data3.A:imax.Data3.A,jmin.Data3.A:jmax.Data3.A)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;
vel.V_var2D_pore.Data3.A = (abs(vel.V_mt2D_NaN.Data3(imin.Data3.A:imax.Data3.A,jmin.Data3.A:jmax.Data3.A)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;

% Averaged non-dimensionalized Eulerian variance based on Vishal's turbulent paper
vel.U_var_mts_Nondim_pore.Data3.A = nanmean(nanmean(vel.U_var2D_pore.Data3.A))/ abs(vel.Mag_mts2D_NaN.Data3)^2;
vel.V_var_mts_Nondim_pore.Data3.A = nanmean(nanmean(vel.V_var2D_pore.Data3.A))/ abs(vel.Mag_mts2D_NaN.Data3)^2;

% Lagrangian velocity variance

vel.U_Lagrang_var2D_pore.Data3.A = nanmean(nanmean(vel.u_var_pore.Data3.A)) + nanmean(nanmean(vel.U_var2D_pore.Data3.A));
vel.V_Lagrang_var2D_pore.Data3.A = nanmean(nanmean(vel.v_var_pore.Data3.A)) + nanmean(nanmean(vel.V_var2D_pore.Data3.A));

%% Plotting velocities

TimeavgplotPore (vel.U_mt2D.Data3,vel.V_mt2D.Data3,vel.Mag_mts2D_NaN.Data3, ...
    x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
saveas(gcf,'Vmt_PoreA_Data3.png');
%% Vector of LES
h = figure;
for k = 10:10
    VectorPlotPore (vel.U_les3gauss2D.Data3(:,:,k),vel.V_les3gauss2D.Data3(:,:,k), ...
        x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
    hold off
    saveas(h,sprintf('Vector_LES3g_PoreA_Data3_%d.png',k));
    
end
%% Vorticity
VorplotPore (vor.vorles3g_2D_mt.Data3,vel.Mag_mts_pore.Data3.A,9.5, x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
saveas(gcf,'vor_PoreA_Data3.png');  

%% Lambda_Ci
kk = 13;
LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),scale.t.T_v_ITS.Data3, x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
saveas(gcf,'LambdaCiT13_PoreA_Data3.png');

%% Vorticity and velocity
close all
h = figure;
for k = 10:110
    
    VorplotPore (vor.vorles3g_2D.Data3(:,:,k),47,10, x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
%     hold on
    VectorPlotPore (vel.U_les5gauss2D.Data3(:,:,k),vel.V_les5gauss2D.Data3(:,:,k), x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
    saveas(h,sprintf('Vor_Data3_PoreA_%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end
% imageFolder2mpeg('Vor_Data3_PoreA','frameRate',11)

%% LambdaCi and velocity

close all
h = figure;
for kk = 10:110
    
    LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),0.001, x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
%     hold on
%     VectorPlotPore (vel.U_les3gauss2D.Data3(:,:,kk),vel.V_les3gauss2D.Data3(:,:,kk), x,y, imin.Data3.A, imax.Data3.A,jmin.Data3.A, jmax.Data3.A);
    saveas(h,sprintf('LambdaCi_Data3_PoreA_%d.png',kk)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end

%% Pore B ******************************************************************

imin.Data3.B = 18; 
imax.Data3.B = 44;
jmin.Data3.B = 32;
jmax.Data3.B = 54;
nodes.Data3.B = 23; %23*27
[vel.U_mt_pore.Data3.B,vel.V_mt_pore.Data3.B,vel.Mag_mt_pore.Data3.B,...
    vel.U_mts_pore.Data3.B,vel.V_mts_pore.Data3.B, ...
    vel.Mag_mts_pore.Data3.B,vel.u_rms_pore.Data3.B, ...
    vel.v_rms_pore.Data3.B, vel.uv_rms_pore.Data3.B, vel.u_rms_ms_pore.Data3.B, ...
    vel.v_rms_ms_pore.Data3.B, vel.uv_rms_ms_pore.Data3.B, vel.u_var_pore.Data3.B,...
    vel.v_var_pore.Data3.B, vel.uv_var_pore.Data3.B, vel.u_var_ms_pore.Data3.B,...
    vel.v_var_ms_pore.Data3.B, vel.uv_var_ms_pore.Data3.B,...
    N_vectors.Data3.B] = poreMeanVelocity (x, y, vel.U_mm2D.Data3, ...
    vel.V_mm2D.Data3, vel.Mag_mm2D.Data3, imin.Data3.B, imax.Data3.B, ...
    jmin.Data3.B, jmax.Data3.B, kmax1, mm_to_pix);
%% Plotting velocities
TimeavgplotPore (vel.U_mt2D.Data3,vel.V_mt2D.Data3,vel.Mag_mts2D_NaN.Data3, ...
    x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
saveas(gcf,'Vmt_PoreB_Data3_22.png');



%% Plotting velocities

TimeavgplotPore (vel.U_les5gauss2D_mt2D_NaN.Data3,vel.V_les5gauss2D_mt2D_NaN.Data3,vel.Mag_mts2D_NaN.Data3, ...
    x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
saveas(gcf,'Vmt_PoreA_Data3_11.png');
%% Vector of LES
VectorPlotPore (vel.U_les3gauss2D.Data3,vel.V_les3gauss2D.Data3, ...
    x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
saveas(gcf,'Vector_LES3g_PoreB_Data3_2.png');
%% Vorticity
VorplotPore (vor.vorles3g_2D_mt.Data3,vel.Mag_mts_pore.Data3.B,5.3, x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
saveas(gcf,'vor_PoreB_Data3.png');  

%% Lambda_Ci
kk = 13;
LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),scale.t.T_v_ITS.Data3, x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
saveas(gcf,'LambdaCiT13_PoreB_Data3.png');
%% Vorticity and velocity
close all
h = figure; 
for k = 20:120 %vel.Mag_mts_pore.Data3.B
    
    VorplotPore (vor.vorles3g_2D.Data3(:,:,k),47,5.3, x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
    hold on
    VectorPlotPore (vel.U_les5gauss2D.Data3(:,:,k),vel.V_les5gauss2D.Data3(:,:,k), x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
    saveas(h,sprintf('Vor_Data3_PoreB_%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end
% imageFolder2mpeg('Vor_Data3_PoreA','frameRate',11)

%% LambdaCi and velocity

close all
h = figure; 
for kk = 20:120  %scale.t.T_v_ITS.Data3
    
    LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),0.0012, x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
%     hold on
%     VectorPlotPore (vel.U_les3gauss2D.Data3(:,:,kk),vel.V_les3gauss2D.Data3(:,:,kk), x,y, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B);
    saveas(h,sprintf('LambdaCi_Data3_PoreB_%d.png',kk)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end 

%%

% vortex size 3*3 and 5*5
[vor.lambda_ciles3g_2D_size_Pore.Data3.B, vor.lambda_ciles3g_2D_size_mst_Pore.Data3.B, vor.lambda_ciles3g_2D_N_vortex_Pore.Data3.B] ...
    = LambdaCiSizePore (vor.lambda_ciles3g_2D.Data3,imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B, kmax1, pix_to_mm);

[vor.lambda_ciles5g_2D_size_Pore.Data3.B, vor.lambda_ciles5g_2D_size_mst_Pore.Data3.B, vor.lambda_ciles5g_2D_N_vortex_Pore.Data3.B] ...
    = LambdaCiSizePore (vor.lambda_ciles5g_2D.Data3,imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B, kmax1, pix_to_mm);
%%

[vor.vor_mt_pore.Data3.B, vor.lambda_ci_mt_pore.Data3.B, vor.vor_mts_pore.Data3.B, ...
    vor.lambda_ci_mts_pore.Data3.B, vor.vor_rms_pore.Data3.B, ...
    vor.lambda_ci_rms_pore.Data3.B, vor.vor_rms_ms_pore.Data3.B, ...
    vor.lambda_ci_rms_ms_pore.Data3.B, vor.N_lambda_ci_pore.Data3.B] ...
    = poreMeanVorticity (x, y, vor.vor_2D.Data3, vor.lambda_ci2D.Data3, imin.Data3.B, imax.Data3.B, ...
    jmin.Data3.B, jmax.Data3.B, kmax1, mm_to_pix);

[vor.circ_vel_pore.Data3.B, vor.circ_vor_pore.Data3.B, vor.circ_vor_mt_pore.Data3.B, ...
    vor.circ_vel_mt_pore.Data3.B, vor.circ_vor_rms_pore.Data3.B, vor.circ_vel_rms_pore.Data3.B] = ...
    circulation (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, vor.vor_2D.Data3, imin.Data3.B, ...
    imax.Data3.B, jmin.Data3.B, jmax.Data3.B, kmax1, nodes.Data3.B, mm_to_pix);

%LES 5*5
[vor.vorles5g_mt_pore.Data3.B, vor.lambda_ciles5g_mt_pore.Data3.B, vor.vorles5g_mts_pore.Data3.B, ...
    vor.lambda_ciles5g_mts_pore.Data3.B, vor.vorles5g_rms_pore.Data3.B, ...
    vor.lambda_ciles5g_rms_pore.Data3.B, vor.vorles5g_rms_ms_pore.Data3.B, ...
    vor.lambda_ciles5g_rms_ms_pore.Data3.B, vor.N_lambda_ciles5g_pore.Data3.B] ...
    = poreLESVorticity (x, y, vor.vorles5g_2D.Data3, vor.lambda_ciles5g_2D.Data3, imin.Data3.B, imax.Data3.B, ...
    jmin.Data3.B, jmax.Data3.B, kmax1, mm_to_pix);

%LES 3*3
[vor.vorles3g_mt_pore.Data3.B, vor.lambda_ciles3g_mt_pore.Data3.B, vor.vorles3g_mts_pore.Data3.B, ...
    vor.lambda_ciles3g_mts_pore.Data3.B, vor.vorles3g_rms_pore.Data3.B, ...
    vor.lambda_ciles3g_rms_pore.Data3.B, vor.vorles3g_rms_ms_pore.Data3.B, ...
    vor.lambda_ciles3g_rms_ms_pore.Data3.B, vor.N_lambda_ciles3g_pore.Data3.B] ...
    = poreMeanVorticity (x, y, vor.vorles3g_2D.Data3, vor.lambda_ciles3g_2D.Data3, imin.Data3.B, imax.Data3.B, ...
    jmin.Data3.B, jmax.Data3.B, kmax1, mm_to_pix);

%% Gamma1 and Gamma2
[vor.Gamma1_pore.Data3.B,vor.Gamma2_pore.Data3.B, vor.Gamma1_mt_pore.Data3.B, vor.Gamma2_mt_pore.Data3.B,...
    vor.Gamma1_mts_pore.Data3.B, vor.Gamma2_mts_pore.Data3.B, vor.Gamma1_rms_pore.Data3.B, ...
    vor.Gamma2_rms_pore.Data3.B, vor.Gamma1_rms_ms_pore.Data3.B, vor.Gamma2_rms_ms_pore.Data3.B,...
    vor.N_Gamma1.Data3.B, vor.N_Gamma2.Data3.B, vor.Gamma2_Shear_mts_pore.Data3.B, vor.Gamma2_Rotation_mts_pore.Data3.B,...
    vor.Gamma2_Shear_rms_ms_pore.Data3.B, vor.Gamma2_Rotation_rms_ms_pore.Data3.B,vor.N_vortex_G2_Shear_pore.Data3.B,...
    vor.N_vortex_G2_Rotation_pore.Data3.B] = Gamma1 (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, imin.Data3.B, imax.Data3.B, ...
    jmin.Data3.B, jmax.Data3.B, kmax1, jmax1, imax1);
%% Vortex size Gamma2

[vor.Gamma2_size.Data3, vor.Gamma2_size_mst.Data3, vor.Gamma2_N_vortex.Data3,...
    vor.Gamma2_size_Pore.Data3.B, vor.Gamma2_size_mst_Pore.Data3.B, vor.Gamma2_N_vortex_Pore.Data3.B] ...
    = Gamma2Size (vor.Gamma2.Data3, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B, kmax1, pix_to_mm, 'pore');

%%
% Integral time scale
[scale.t.T_u2D_pore.Data3.B, scale.t.T_v2D_pore.Data3.B, scale.t.T_u_ITS_pore.Data3.B, scale.t.T_v_ITS_pore.Data3.B] = ...
    IntTimeScalePore (scale.t.T_u2D.Data3, scale.t.T_v2D.Data3, imin.Data3.B, imax.Data3.B,...
    jmin.Data3.B, jmax.Data3.B);

% Integral length scale
lag.L.Data3.B = imax.Data3.B - imin.Data3.B ; 
[scale.autocorr_L_V_mst_pore.Data3.B,scale.ILS.Data3.B] = ...
    IntLengthScalePore (vel.V_mm2D.Data3, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B, kmax1,lag.L.Data3.B, y.Del_y);

% TKE budget
% Production, TKE advection
[turb.Production.P_linear_pore.Data3.B, turb.Production.P_shear_pore.Data3.B, turb.KE_advection_pore.Data3.B, ...
    turb.Net_grad.grad_x_pore.Data3.B, turb.Net_grad.grad_y_pore.Data3.B, turb.Net_grad.grad_z_pore.Data3.B] = ...
    TKEbudgetPore (vel.U_mt2D.Data3, vel.V_mt2D.Data3, turb.Production.Sxx.Data3,...
    turb.Production.Syy.Data3, turb.Production.Szz.Data3, turb.Production.Sxy.Data3, vel.u_var2D.Data3, ...
    vel.v_var2D.Data3, vel.uv_var2D.Data3, vor.du_var_dx.Data3, vor.du_var_dy.Data3, ...
    vor.dv_var_dx.Data3, vor.dv_var_dy.Data3, turb.Production.sxx.Data3, turb.Production.sxy.Data3, ...
    turb.Production.syy.Data3, vel.u_t2D.Data3, vel.v_t2D.Data3, imin.Data3.B, imax.Data3.B,jmin.Data3.B, jmax.Data3.B,specs.fluid.kvisc);

% Net gradient transport
[turb.d_Netgrad_x_dx_pore.Data3.B, turb.d_Netgrad_x_dy_pore.Data3.B, turb.d_Netgrad_y_dx_pore.Data3.B, turb.d_Netgrad_y_dy_pore.Data3.B] = ...
    VGT (x, y, turb.Net_grad.grad_x_pore.Data3.B, turb.Net_grad.grad_y_pore.Data3.B, jmax.Data3.B-jmin.Data3.B+1, imax.Data3.B-imin.Data3.B+1, 1, mm_to_pix);

[turb.d_Netgrad_z_dx_pore.Data3.B, turb.d_Netgrad_z_dy_pore.Data3.B, turb.d_Netgrad_y_dx_pore.Data3.B, turb.d_Netgrad_y_dy_pore.Data3.B] = ...
    VGT (x, y, turb.Net_grad.grad_z_pore.Data3.B, turb.Net_grad.grad_y_pore.Data3.B, jmax.Data3.B-jmin.Data3.B+1, imax.Data3.B-imin.Data3.B+1, 1, mm_to_pix);

turb.Net_grad.total_trans_pore.Data3.B = nanmean(nanmean(turb.d_Netgrad_x_dx_pore.Data3.B + turb.d_Netgrad_y_dy_pore.Data3.B + turb.d_Netgrad_z_dx_pore.Data3.B));

% TKE
% Turbulent kinetic energy
turb.TKE_2D_pore.Data3.B= (1/2)*(2*vel.u_var2D.Data3(imin.Data3.B:imax.Data3.B,jmin.Data3.B:jmax.Data3.B) + vel.v_var2D.Data3(imin.Data3.B:imax.Data3.B,jmin.Data3.B:jmax.Data3.B));
turb.TKE_ms_pore.Data3.B= nanmean(nanmean(turb.TKE_2D_pore.Data3.B));

% Tortuosity
dispers.Tortuosity_pore.Data3.B = nansum(nansum(nansum(vel.Mag_mm2D.Data3(imin.Data3.B:imax.Data3.B,jmin.Data3.B:jmax.Data3.B))))./nansum(nansum(nansum(vel.V_mm2D.Data3(imin.Data3.B:imax.Data3.B,jmin.Data3.B:jmax.Data3.B))));

% Eulerian variance
% 2D field
vel.U_var2D_pore.Data3.B = (abs(vel.U_mt2D_NaN.Data3(imin.Data3.B:imax.Data3.B,jmin.Data3.B:jmax.Data3.B)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;
vel.V_var2D_pore.Data3.B = (abs(vel.V_mt2D_NaN.Data3(imin.Data3.B:imax.Data3.B,jmin.Data3.B:jmax.Data3.B)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;

% Averaged non-dimensionalized Eulerian variance based on Vishal's turbulent paper
vel.U_var_mts_Nondim_pore.Data3.B = nanmean(nanmean(vel.U_var2D_pore.Data3.B))/ abs(vel.Mag_mts2D_NaN.Data3)^2;
vel.V_var_mts_Nondim_pore.Data3.B = nanmean(nanmean(vel.V_var2D_pore.Data3.B))/ abs(vel.Mag_mts2D_NaN.Data3)^2;

%% Lagrangian velocity variance

vel.U_Lagrang_var2D_pore.Data3.B = nanmean(nanmean(vel.u_var_pore.Data3.B)) + nanmean(nanmean(vel.U_var2D_pore.Data3.B));
vel.V_Lagrang_var2D_pore.Data3.B = nanmean(nanmean(vel.v_var_pore.Data3.B)) + nanmean(nanmean(vel.V_var2D_pore.Data3.B));

%% Pore C ******************************************************************

imin.Data3.C = 1; 
imax.Data3.C = 11;
jmin.Data3.C = 53;
jmax.Data3.C = 67;
nodes.Data3.C = 15; %15*11
[vel.U_mt_pore.Data3.C,vel.V_mt_pore.Data3.C,vel.Mag_mt_pore.Data3.C,...
    vel.U_mts_pore.Data3.C,vel.V_mts_pore.Data3.C, ...
    vel.Mag_mts_pore.Data3.C,vel.u_rms_pore.Data3.C, ...
    vel.v_rms_pore.Data3.C, vel.uv_rms_pore.Data3.C, vel.u_rms_ms_pore.Data3.C, ...
    vel.v_rms_ms_pore.Data3.C, vel.uv_rms_ms_pore.Data3.C, vel.u_var_pore.Data3.C,...
    vel.v_var_pore.Data3.C, vel.uv_var_pore.Data3.C, vel.u_var_ms_pore.Data3.C,...
    vel.v_var_ms_pore.Data3.C, vel.uv_var_ms_pore.Data3.C,...
    N_vectors.Data3.C] = poreMeanVelocity (x, y, vel.U_mm2D.Data3, ...
    vel.V_mm2D.Data3, vel.Mag_mm2D.Data3, imin.Data3.C, imax.Data3.C, ...
    jmin.Data3.C, jmax.Data3.C, kmax1, mm_to_pix);
%% Plotting velocities

TimeavgplotPore (vel.U_mt2D.Data3,vel.V_mt2D.Data3,vel.Mag_mts_pore.Data3.C, ...
    x,y, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C);
saveas(gcf,'Vmt_PoreC_Data3.png');
%% Vector of LES
VectorPlotPore (vel.U_les3gauss2D.Data3,vel.V_les3gauss2D.Data3, ...
    x,y, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C);
saveas(gcf,'Vector_LES3g_PoreC_Data3.png');
%% Vorticity
VorplotPore (vor.vorles3g_2D_mt.Data3,vel.Mag_mts_pore.Data3.C,8.55, x,y, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C);
saveas(gcf,'vor_PoreC_Data3.png');  

%% Lambda_Ci
kk = 13;
LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),scale.t.T_v_ITS.Data3, x,y, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C);
saveas(gcf,'LambdaCiT13_PoreC_Data3.png');
%% Vorticity and velocity
close all
h = figure;
for k = 10:110 %vel.Mag_mts_pore.Data3.C
    
    VorplotPore (vor.vorles3g_2D.Data3(:,:,k),47.,5.3, x,y, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C);
    hold on
    VectorPlotPore (vel.U_les5gauss2D.Data3(:,:,k),vel.V_les5gauss2D.Data3(:,:,k), x,y, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C);
    saveas(h,sprintf('Vor_Data3_PoreC_%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end
% imageFolder2mpeg('Vor_Data3_PoreA','frameRate',11)

%% LambdaCi and velocity

close all
h = figure;
for kk = 21:120 %scale.t.T_v_ITS.Data3
    
    LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),0.0041, x,y, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C);
%     hold on
%     VectorPlotPore (vel.U_les3gauss2D.Data3(:,:,kk),vel.V_les3gauss2D.Data3(:,:,kk), x,y, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C);
    saveas(h,sprintf('LambdaCi_Data3_PoreC_%d.png',kk)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end

%%

% vortex size 3*3 and 5*5
[vor.lambda_ciles3g_2D_size_Pore.Data3.C, vor.lambda_ciles3g_2D_size_mst_Pore.Data3.C, vor.lambda_ciles3g_2D_N_vortex_Pore.Data3.C] ...
    = LambdaCiSizePore (vor.lambda_ciles3g_2D.Data3,imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C, kmax1, pix_to_mm);

[vor.lambda_ciles5g_2D_size_Pore.Data3.C, vor.lambda_ciles5g_2D_size_mst_Pore.Data3.C, vor.lambda_ciles5g_2D_N_vortex_Pore.Data3.C] ...
    = LambdaCiSizePore (vor.lambda_ciles5g_2D.Data3,imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C, kmax1, pix_to_mm);
%%


[vor.vor_mt_pore.Data3.C, vor.lambda_ci_mt_pore.Data3.C, vor.vor_mts_pore.Data3.C, ...
    vor.lambda_ci_mts_pore.Data3.C, vor.vor_rms_pore.Data3.C, ...
    vor.lambda_ci_rms_pore.Data3.C, vor.vor_rms_ms_pore.Data3.C, ...
    vor.lambda_ci_rms_ms_pore.Data3.C, vor.N_lambda_ci_pore.Data3.C] ...
    = poreMeanVorticity (x, y, vor.vor_2D.Data3, vor.lambda_ci2D.Data3, imin.Data3.C, imax.Data3.C, ...
    jmin.Data3.C, jmax.Data3.C, kmax1, mm_to_pix);

[vor.circ_vel_pore.Data3.C, vor.circ_vor_pore.Data3.C, vor.circ_vor_mt_pore.Data3.C, ...
    vor.circ_vel_mt_pore.Data3.C, vor.circ_vor_rms_pore.Data3.C, vor.circ_vel_rms_pore.Data3.C] = ...
    circulation (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, vor.vor_2D.Data3, imin.Data3.C, ...
    imax.Data3.C, jmin.Data3.C, jmax.Data3.C, kmax1, nodes.Data3.C, mm_to_pix);

%LES 5*5
[vor.vorles5g_mt_pore.Data3.C, vor.lambda_ciles5g_mt_pore.Data3.C, vor.vorles5g_mts_pore.Data3.C, ...
    vor.lambda_ciles5g_mts_pore.Data3.C, vor.vorles5g_rms_pore.Data3.C, ...
    vor.lambda_ciles5g_rms_pore.Data3.C, vor.vorles5g_rms_ms_pore.Data3.C, ...
    vor.lambda_ciles5g_rms_ms_pore.Data3.C, vor.N_lambda_ciles5g_pore.Data3.C] ...
    = poreLESVorticity (x, y, vor.vorles5g_2D.Data3, vor.lambda_ciles5g_2D.Data3, imin.Data3.C, imax.Data3.C, ...
    jmin.Data3.C, jmax.Data3.C, kmax1, mm_to_pix);

%LES 3*3
[vor.vorles3g_mt_pore.Data3.C, vor.lambda_ciles3g_mt_pore.Data3.C, vor.vorles3g_mts_pore.Data3.C, ...
    vor.lambda_ciles3g_mts_pore.Data3.C, vor.vorles3g_rms_pore.Data3.C, ...
    vor.lambda_ciles3g_rms_pore.Data3.C, vor.vorles3g_rms_ms_pore.Data3.C, ...
    vor.lambda_ciles3g_rms_ms_pore.Data3.C, vor.N_lambda_ciles3g_pore.Data3.C] ...
    = poreMeanVorticity (x, y, vor.vorles3g_2D.Data3, vor.lambda_ciles3g_2D.Data3, imin.Data3.C, imax.Data3.C, ...
    jmin.Data3.C, jmax.Data3.C, kmax1, mm_to_pix);

%% Gamma1 and Gamma2
[vor.Gamma1_pore.Data3.C,vor.Gamma2_pore.Data3.C, vor.Gamma1_mt_pore.Data3.C, vor.Gamma2_mt_pore.Data3.C,...
    vor.Gamma1_mts_pore.Data3.C, vor.Gamma2_mts_pore.Data3.C, vor.Gamma1_rms_pore.Data3.C, ...
    vor.Gamma2_rms_pore.Data3.C, vor.Gamma1_rms_ms_pore.Data3.C, vor.Gamma2_rms_ms_pore.Data3.C,...
    vor.N_Gamma1.Data3.C, vor.N_Gamma2.Data3.C, vor.Gamma2_Shear_mts_pore.Data3.C, vor.Gamma2_Rotation_mts_pore.Data3.C,...
    vor.Gamma2_Shear_rms_ms_pore.Data3.C, vor.Gamma2_Rotation_rms_ms_pore.Data3.C,vor.N_vortex_G2_Shear_pore.Data3.C,...
    vor.N_vortex_G2_Rotation_pore.Data3.C] = Gamma1 (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, imin.Data3.C, imax.Data3.C, ...
    jmin.Data3.C, jmax.Data3.C, kmax1, jmax1, imax1);
%% Vortex size Gamma2

[vor.Gamma2_size.Data3, vor.Gamma2_size_mst.Data3, vor.Gamma2_N_vortex.Data3,...
    vor.Gamma2_size_Pore.Data3.C, vor.Gamma2_size_mst_Pore.Data3.C, vor.Gamma2_N_vortex_Pore.Data3.C] ...
    = Gamma2Size (vor.Gamma2.Data3, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C, kmax1, pix_to_mm, 'pore');

%%

% Integral time scale
[scale.t.T_u2D_pore.Data3.C, scale.t.T_v2D_pore.Data3.C, scale.t.T_u_ITS_pore.Data3.C, scale.t.T_v_ITS_pore.Data3.C] = ...
    IntTimeScalePore (scale.t.T_u2D.Data3, scale.t.T_v2D.Data3, imin.Data3.C, imax.Data3.C,...
    jmin.Data3.C, jmax.Data3.C);

% Integral length scale
lag.L.Data3.C = imax.Data3.C - imin.Data3.C ; 
[scale.autocorr_L_V_mst_pore.Data3.C,scale.ILS.Data3.C] = ...
    IntLengthScalePore (vel.V_mm2D.Data3, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C, kmax1,lag.L.Data3.C, y.Del_y);

% TKE budget
% Production, TKE advection
[turb.Production.P_linear_pore.Data3.C, turb.Production.P_shear_pore.Data3.C, turb.KE_advection_pore.Data3.C, ...
    turb.Net_grad.grad_x_pore.Data3.C, turb.Net_grad.grad_y_pore.Data3.C, turb.Net_grad.grad_z_pore.Data3.C] = ...
    TKEbudgetPore (vel.U_mt2D.Data3, vel.V_mt2D.Data3, turb.Production.Sxx.Data3,...
    turb.Production.Syy.Data3, turb.Production.Szz.Data3, turb.Production.Sxy.Data3, vel.u_var2D.Data3, ...
    vel.v_var2D.Data3, vel.uv_var2D.Data3, vor.du_var_dx.Data3, vor.du_var_dy.Data3, ...
    vor.dv_var_dx.Data3, vor.dv_var_dy.Data3, turb.Production.sxx.Data3, turb.Production.sxy.Data3, ...
    turb.Production.syy.Data3, vel.u_t2D.Data3, vel.v_t2D.Data3, imin.Data3.C, imax.Data3.C,jmin.Data3.C, jmax.Data3.C,specs.fluid.kvisc);

% Net gradient transport
[turb.d_Netgrad_x_dx_pore.Data3.C, turb.d_Netgrad_x_dy_pore.Data3.C, turb.d_Netgrad_y_dx_pore.Data3.C, turb.d_Netgrad_y_dy_pore.Data3.C] = ...
    VGT (x, y, turb.Net_grad.grad_x_pore.Data3.C, turb.Net_grad.grad_y_pore.Data3.C, jmax.Data3.C-jmin.Data3.C+1, imax.Data3.C-imin.Data3.C+1, 1, mm_to_pix);

[turb.d_Netgrad_z_dx_pore.Data3.C, turb.d_Netgrad_z_dy_pore.Data3.C, turb.d_Netgrad_y_dx_pore.Data3.C, turb.d_Netgrad_y_dy_pore.Data3.C] = ...
    VGT (x, y, turb.Net_grad.grad_z_pore.Data3.C, turb.Net_grad.grad_y_pore.Data3.C, jmax.Data3.C-jmin.Data3.C+1, imax.Data3.C-imin.Data3.C+1, 1, mm_to_pix);

turb.Net_grad.total_trans_pore.Data3.C = nanmean(nanmean(turb.d_Netgrad_x_dx_pore.Data3.C + turb.d_Netgrad_y_dy_pore.Data3.C + turb.d_Netgrad_z_dx_pore.Data3.C));

% TKE
% Turbulent kinetic energy
turb.TKE_2D_pore.Data3.C = (1/2)*(2*vel.u_var2D.Data3(imin.Data3.C:imax.Data3.C,jmin.Data3.C:jmax.Data3.C) + vel.v_var2D.Data3(imin.Data3.C:imax.Data3.C,jmin.Data3.C:jmax.Data3.C));
turb.TKE_ms_pore.Data3.C = nanmean(nanmean(turb.TKE_2D_pore.Data3.C));

% Tortuosity
dispers.Tortuosity_pore.Data3.C = nansum(nansum(nansum(vel.Mag_mm2D.Data3(imin.Data3.C:imax.Data3.C,jmin.Data3.C:jmax.Data3.C))))./nansum(nansum(nansum(vel.V_mm2D.Data3(imin.Data3.C:imax.Data3.C,jmin.Data3.C:jmax.Data3.C))));

% Eulerian variance
% 2D field
vel.U_var2D_pore.Data3.C = (abs(vel.U_mt2D_NaN.Data3(imin.Data3.C:imax.Data3.C,jmin.Data3.C:jmax.Data3.C)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;
vel.V_var2D_pore.Data3.C = (abs(vel.V_mt2D_NaN.Data3(imin.Data3.C:imax.Data3.C,jmin.Data3.C:jmax.Data3.C)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;

% Averaged non-dimensionalized Eulerian variance based on Vishal's turbulent paper
vel.U_var_mts_Nondim_pore.Data3.C = nanmean(nanmean(vel.U_var2D_pore.Data3.C))/ abs(vel.Mag_mts2D_NaN.Data3)^2;
vel.V_var_mts_Nondim_pore.Data3.C = nanmean(nanmean(vel.V_var2D_pore.Data3.C))/ abs(vel.Mag_mts2D_NaN.Data3)^2;

%% Lagrangian velocity variance

vel.U_Lagrang_var2D_pore.Data3.C = nanmean(nanmean(vel.u_var_pore.Data3.C)) + nanmean(nanmean(vel.U_var2D_pore.Data3.C));
vel.V_Lagrang_var2D_pore.Data3.C = nanmean(nanmean(vel.v_var_pore.Data3.C)) + nanmean(nanmean(vel.V_var2D_pore.Data3.C));

%% Pore D ******************************************************************

imin.Data3.D = 11; 
imax.Data3.D = 23;
jmin.Data3.D = 10;
jmax.Data3.D = 28;
nodes.Data3.D = 19; %13*19

[vel.U_mt_pore.Data3.D,vel.V_mt_pore.Data3.D,vel.Mag_mt_pore.Data3.D,...
    vel.U_mts_pore.Data3.D,vel.V_mts_pore.Data3.D, ...
    vel.Mag_mts_pore.Data3.D,vel.u_rms_pore.Data3.D, ...
    vel.v_rms_pore.Data3.D, vel.uv_rms_pore.Data3.D, vel.u_rms_ms_pore.Data3.D, ...
    vel.v_rms_ms_pore.Data3.D, vel.uv_rms_ms_pore.Data3.D, vel.u_var_pore.Data3.D,...
    vel.v_var_pore.Data3.D, vel.uv_var_pore.Data3.D, vel.u_var_ms_pore.Data3.D,...
    vel.v_var_ms_pore.Data3.D, vel.uv_var_ms_pore.Data3.D,...
    N_vectors.Data3.D] = poreMeanVelocity (x, y, vel.U_mm2D.Data3, ...
    vel.V_mm2D.Data3, vel.Mag_mm2D.Data3, imin.Data3.D, imax.Data3.D, ...
    jmin.Data3.D, jmax.Data3.D, kmax1, mm_to_pix);

%% Plotting velocities

TimeavgplotPore (vel.U_mt2D.Data3,vel.V_mt2D.Data3,vel.Mag_mts_pore.Data3.D, ...
    x,y, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D);
saveas(gcf,'Vmt_PoreD_Data3.png');
%% Vector of LES
VectorPlotPore (vel.U_les3gauss2D.Data3,vel.V_les3gauss2D.Data3, ...
    x,y, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D);
saveas(gcf,'Vector_LES3g_PoreD_Data3.png');
%% Vorticity
VorplotPore (vor.vorles3g_2D_mt.Data3,vel.Mag_mts_pore.Data3.D,4.4, x,y, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D);
saveas(gcf,'vor_PoreD_Data3.png');  

%% Lambda_Ci
kk = 13;
LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),scale.t.T_v_ITS.Data3, x,y, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D);
saveas(gcf,'LambdaCiT13_PoreD_Data3.png');
%% Vorticity and velocity
close all
h = figure;
for k = 20:120 %vel.Mag_mts_pore.Data3.D
    
    VorplotPore (vor.vorles3g_2D.Data3(:,:,k),47,14, x,y, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D);
%     hold on
    VectorPlotPore (vel.U_les5gauss2D.Data3(:,:,k),vel.V_les5gauss2D.Data3(:,:,k), x,y, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D);
    saveas(h,sprintf('Vor_Data3_PoreD_%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end
% imageFolder2mpeg('Vor_Data3_PoreD','frameRate',11)

%% LambdaCi and velocity

close all
h = figure;
for kk = 20:120  %scale.t.T_v_ITS.Data3
    
    LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),0.0021, x,y, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D);
%     hold on
%     VectorPlotPore (vel.U_les3gauss2D.Data3(:,:,kk),vel.V_les3gauss2D.Data3(:,:,kk), x,y, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D);
    saveas(h,sprintf('LambdaCi_Data3_PoreD_%d.png',kk)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end


%%

% vortex size 3*3 and 5*5
[vor.lambda_ciles3g_2D_size_Pore.Data3.D, vor.lambda_ciles3g_2D_size_mst_Pore.Data3.D, vor.lambda_ciles3g_2D_N_vortex_Pore.Data3.D] ...
    = LambdaCiSizePore (vor.lambda_ciles3g_2D.Data3,imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D, kmax1, pix_to_mm);

[vor.lambda_ciles5g_2D_size_Pore.Data3.D, vor.lambda_ciles5g_2D_size_mst_Pore.Data3.D, vor.lambda_ciles5g_2D_N_vortex_Pore.Data3.D] ...
    = LambdaCiSizePore (vor.lambda_ciles5g_2D.Data3,imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D, kmax1, pix_to_mm);
%%


[vor.vor_mt_pore.Data3.D, vor.lambda_ci_mt_pore.Data3.D, vor.vor_mts_pore.Data3.D, ...
    vor.lambda_ci_mts_pore.Data3.D, vor.vor_rms_pore.Data3.D, ...
    vor.lambda_ci_rms_pore.Data3.D, vor.vor_rms_ms_pore.Data3.D, ...
    vor.lambda_ci_rms_ms_pore.Data3.D, vor.N_lambda_ci_pore.Data3.D] ...
    = poreMeanVorticity (x, y, vor.vor_2D.Data3, vor.lambda_ci2D.Data3, imin.Data3.D, imax.Data3.D, ...
    jmin.Data3.D, jmax.Data3.D, kmax1, mm_to_pix);

[vor.circ_vel_pore.Data3.D, vor.circ_vor_pore.Data3.D, vor.circ_vor_mt_pore.Data3.D, ...
    vor.circ_vel_mt_pore.Data3.D, vor.circ_vor_rms_pore.Data3.D, vor.circ_vel_rms_pore.Data3.D] = ...
    circulation (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, vor.vor_2D.Data3, imin.Data3.D, ...
    imax.Data3.D, jmin.Data3.D, jmax.Data3.D, kmax1, nodes.Data3.D, mm_to_pix);

%LES 5*5
[vor.vorles5g_mt_pore.Data3.D, vor.lambda_ciles5g_mt_pore.Data3.D, vor.vorles5g_mts_pore.Data3.D, ...
    vor.lambda_ciles5g_mts_pore.Data3.D, vor.vorles5g_rms_pore.Data3.D, ...
    vor.lambda_ciles5g_rms_pore.Data3.D, vor.vorles5g_rms_ms_pore.Data3.D, ...
    vor.lambda_ciles5g_rms_ms_pore.Data3.D, vor.N_lambda_ciles5g_pore.Data3.D] ...
    = poreLESVorticity (x, y, vor.vorles5g_2D.Data3, vor.lambda_ciles5g_2D.Data3, imin.Data3.D, imax.Data3.D, ...
    jmin.Data3.D, jmax.Data3.D, kmax1, mm_to_pix);

%LES 3*3
[vor.vorles3g_mt_pore.Data3.D, vor.lambda_ciles3g_mt_pore.Data3.D, vor.vorles3g_mts_pore.Data3.D, ...
    vor.lambda_ciles3g_mts_pore.Data3.D, vor.vorles3g_rms_pore.Data3.D, ...
    vor.lambda_ciles3g_rms_pore.Data3.D, vor.vorles3g_rms_ms_pore.Data3.D, ...
    vor.lambda_ciles3g_rms_ms_pore.Data3.D, vor.N_lambda_ciles3g_pore.Data3.D] ...
    = poreMeanVorticity (x, y, vor.vorles3g_2D.Data3, vor.lambda_ciles3g_2D.Data3, imin.Data3.D, imax.Data3.D, ...
    jmin.Data3.D, jmax.Data3.D, kmax1, mm_to_pix);


%% Gamma1 and Gamma2
[vor.Gamma1_pore.Data3.D,vor.Gamma2_pore.Data3.D, vor.Gamma1_mt_pore.Data3.D, vor.Gamma2_mt_pore.Data3.D,...
    vor.Gamma1_mts_pore.Data3.D, vor.Gamma2_mts_pore.Data3.D, vor.Gamma1_rms_pore.Data3.D, ...
    vor.Gamma2_rms_pore.Data3.D, vor.Gamma1_rms_ms_pore.Data3.D, vor.Gamma2_rms_ms_pore.Data3.D,...
    vor.N_Gamma1.Data3.D, vor.N_Gamma2.Data3.D, vor.Gamma2_Shear_mts_pore.Data3.D, vor.Gamma2_Rotation_mts_pore.Data3.D,...
    vor.Gamma2_Shear_rms_ms_pore.Data3.D, vor.Gamma2_Rotation_rms_ms_pore.Data3.D,vor.N_vortex_G2_Shear_pore.Data3.D,...
    vor.N_vortex_G2_Rotation_pore.Data3.D] = Gamma1 (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, imin.Data3.D, imax.Data3.D, ...
    jmin.Data3.D, jmax.Data3.D, kmax1, jmax1, imax1);
%% Vortex size Gamma2

[vor.Gamma2_size.Data3, vor.Gamma2_size_mst.Data3, vor.Gamma2_N_vortex.Data3,...
    vor.Gamma2_size_Pore.Data3.D, vor.Gamma2_size_mst_Pore.Data3.D, vor.Gamma2_N_vortex_Pore.Data3.D] ...
    = Gamma2Size (vor.Gamma2.Data3, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D, kmax1, pix_to_mm, 'pore');

%%
% Integral time scale
[scale.t.T_u2D_pore.Data3.D, scale.t.T_v2D_pore.Data3.D, scale.t.T_u_ITS_pore.Data3.D, scale.t.T_v_ITS_pore.Data3.D] = ...
    IntTimeScalePore (scale.t.T_u2D.Data3, scale.t.T_v2D.Data3, imin.Data3.D, imax.Data3.D,...
    jmin.Data3.D, jmax.Data3.D);

% Integral length scale
lag.L.Data3.D = imax.Data3.D - imin.Data3.D ; 
[scale.autocorr_L_V_mst_pore.Data3.D,scale.ILS.Data3.D] = ...
    IntLengthScalePore (vel.V_mm2D.Data3, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D, kmax1,lag.L.Data3.D, y.Del_y);


% TKE budget
% Production, TKE Advection
[turb.Production.P_linear_pore.Data3.D, turb.Production.P_shear_pore.Data3.D, turb.KE_advection_pore.Data3.D, ...
    turb.Net_grad.grad_x_pore.Data3.D, turb.Net_grad.grad_y_pore.Data3.D, turb.Net_grad.grad_z_pore.Data3.D] = ...
    TKEbudgetPore (vel.U_mt2D.Data3, vel.V_mt2D.Data3, turb.Production.Sxx.Data3,...
    turb.Production.Syy.Data3, turb.Production.Szz.Data3, turb.Production.Sxy.Data3, vel.u_var2D.Data3, ...
    vel.v_var2D.Data3, vel.uv_var2D.Data3, vor.du_var_dx.Data3, vor.du_var_dy.Data3, ...
    vor.dv_var_dx.Data3, vor.dv_var_dy.Data3, turb.Production.sxx.Data3, turb.Production.sxy.Data3, ...
    turb.Production.syy.Data3, vel.u_t2D.Data3, vel.v_t2D.Data3, imin.Data3.D, imax.Data3.D,jmin.Data3.D, jmax.Data3.D,specs.fluid.kvisc);

% Net gradient transport
[turb.d_Netgrad_x_dx_pore.Data3.D, turb.d_Netgrad_x_dy_pore.Data3.D, turb.d_Netgrad_y_dx_pore.Data3.D, turb.d_Netgrad_y_dy_pore.Data3.D] = ...
    VGT (x, y, turb.Net_grad.grad_x_pore.Data3.D, turb.Net_grad.grad_y_pore.Data3.D, jmax.Data3.D-jmin.Data3.D+1, imax.Data3.D-imin.Data3.D+1, 1, mm_to_pix);

[turb.d_Netgrad_z_dx_pore.Data3.D, turb.d_Netgrad_z_dy_pore.Data3.D, turb.d_Netgrad_y_dx_pore.Data3.D, turb.d_Netgrad_y_dy_pore.Data3.D] = ...
    VGT (x, y, turb.Net_grad.grad_z_pore.Data3.D, turb.Net_grad.grad_y_pore.Data3.D, jmax.Data3.D-jmin.Data3.D+1, imax.Data3.D-imin.Data3.D+1, 1, mm_to_pix);

turb.Net_grad.total_trans_pore.Data3.D = nanmean(nanmean(turb.d_Netgrad_x_dx_pore.Data3.D + turb.d_Netgrad_y_dy_pore.Data3.D + turb.d_Netgrad_z_dx_pore.Data3.D));

% TKE
% Turbulent kinetic energy
turb.TKE_2D_pore.Data3.D = (1/2)*(2*vel.u_var2D.Data3(imin.Data3.D:imax.Data3.D,jmin.Data3.D:jmax.Data3.D) + vel.v_var2D.Data3(imin.Data3.D:imax.Data3.D,jmin.Data3.D:jmax.Data3.D));
turb.TKE_ms_pore.Data3.D = nanmean(nanmean(turb.TKE_2D_pore.Data3.D));

% Tortuosity
dispers.Tortuosity_pore.Data3.D = nansum(nansum(nansum(vel.Mag_mm2D.Data3(imin.Data3.D:imax.Data3.D,jmin.Data3.D:jmax.Data3.D))))./nansum(nansum(nansum(vel.V_mm2D.Data3(imin.Data3.D:imax.Data3.D,jmin.Data3.D:jmax.Data3.D))));


% Eulerian variance
% 2D field
vel.U_var2D_pore.Data3.D = (abs(vel.U_mt2D_NaN.Data3(imin.Data3.D:imax.Data3.D,jmin.Data3.D:jmax.Data3.D)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;
vel.V_var2D_pore.Data3.D = (abs(vel.V_mt2D_NaN.Data3(imin.Data3.D:imax.Data3.D,jmin.Data3.D:jmax.Data3.D)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;

% Averaged non-dimensionalized Eulerian variance based on Vishal's turbulent paper
vel.U_var_mts_Nondim_pore.Data3.D = nanmean(nanmean(vel.U_var2D_pore.Data3.D))/ abs(vel.Mag_mts2D_NaN.Data3)^2;
vel.V_var_mts_Nondim_pore.Data3.D = nanmean(nanmean(vel.V_var2D_pore.Data3.D))/ abs(vel.Mag_mts2D_NaN.Data3)^2;

%% Lagrangian velocity variance

vel.U_Lagrang_var2D_pore.Data3.D = nanmean(nanmean(vel.u_var_pore.Data3.D)) + nanmean(nanmean(vel.U_var2D_pore.Data3.D));
vel.V_Lagrang_var2D_pore.Data3.D = nanmean(nanmean(vel.v_var_pore.Data3.D)) + nanmean(nanmean(vel.V_var2D_pore.Data3.D));

%% Pore E ******************************************************************
imin.Data3.E = 26; 
imax.Data3.E = 46;
jmin.Data3.E = 2;
jmax.Data3.E = 12;
nodes.Data3.E = 21; %21*11

[vel.U_mt_pore.Data3.E,vel.V_mt_pore.Data3.E,vel.Mag_mt_pore.Data3.E,...
    vel.U_mts_pore.Data3.E,vel.V_mts_pore.Data3.E, ...
    vel.Mag_mts_pore.Data3.E,vel.u_rms_pore.Data3.E, ...
    vel.v_rms_pore.Data3.E, vel.uv_rms_pore.Data3.E, vel.u_rms_ms_pore.Data3.E, ...
    vel.v_rms_ms_pore.Data3.E, vel.uv_rms_ms_pore.Data3.E, vel.u_var_pore.Data3.E,...
    vel.v_var_pore.Data3.E, vel.uv_var_pore.Data3.E, vel.u_var_ms_pore.Data3.E,...
    vel.v_var_ms_pore.Data3.E, vel.uv_var_ms_pore.Data3.E,...
    N_vectors.Data3.E] = poreMeanVelocity (x, y, vel.U_mm2D.Data3, ...
    vel.V_mm2D.Data3, vel.Mag_mm2D.Data3, imin.Data3.E, imax.Data3.E, ...
    jmin.Data3.E, jmax.Data3.E, kmax1, mm_to_pix);

%% Plotting velocities

TimeavgplotPore (vel.U_mt2D.Data3,vel.V_mt2D.Data3,vel.Mag_mts_pore.Data3.E, ...
    x,y, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E);
saveas(gcf,'Vmt_PoreE_Data3.png');
%% Vector of LES
VectorPlotPore (vel.U_les3gauss2D.Data3,vel.V_les3gauss2D.Data3, ...
    x,y, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E);
saveas(gcf,'Vector_LES3g_PoreE_Data3.png');
%% Vorticity
VorplotPore (vor.vorles3g_2D_mt.Data3,vel.Mag_mts_pore.Data3.E,5.3, x,y, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E);
saveas(gcf,'vor_PoreE_Data3.png');  

%% Lambda_Ci
kk = 13;
LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),scale.t.T_v_ITS.Data3, x,y, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E);
saveas(gcf,'LambdaCiT13_PoreE_Data3.png');
%% Vorticity and velocity
close all
h = figure;
for k = 10:110 %vel.Mag_mts_pore.Data3.E
    
    VorplotPore (vor.vorles3g_2D.Data3(:,:,k),47,14, x,y, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E);
%     hold on
    VectorPlotPore (vel.U_les5gauss2D.Data3(:,:,k),vel.V_les5gauss2D.Data3(:,:,k), x,y, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E);
    saveas(h,sprintf('Vor_Data3_PoreE_%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end
% imageFolder2mpeg('Vor_Data3_PoreE','frameRate',11)

%% LambdaCi and velocity

close all
h = figure;
for kk = 20:120 %scale.t.T_v_ITS.Data3
    
    LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),0.0021, x,y, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E);
%     hold on
%     VectorPlotPore (vel.U_les3gauss2D.Data3(:,:,kk),vel.V_les3gauss2D.Data3(:,:,kk), x,y, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E);
    saveas(h,sprintf('LambdaCi_Data3_PoreE_%d.png',kk)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end

%%

% vortex size 3*3 and 5*5
[vor.lambda_ciles3g_2D_size_Pore.Data3.E, vor.lambda_ciles3g_2D_size_mst_Pore.Data3.E, vor.lambda_ciles3g_2D_N_vortex_Pore.Data3.E] ...
    = LambdaCiSizePore (vor.lambda_ciles3g_2D.Data3,imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E, kmax1, pix_to_mm);

[vor.lambda_ciles5g_2D_size_Pore.Data3.E, vor.lambda_ciles5g_2D_size_mst_Pore.Data3.E, vor.lambda_ciles5g_2D_N_vortex_Pore.Data3.E] ...
    = LambdaCiSizePore (vor.lambda_ciles5g_2D.Data3,imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E, kmax1, pix_to_mm);
%%

%
[vor.vor_mt_pore.Data3.E, vor.lambda_ci_mt_pore.Data3.E, vor.vor_mts_pore.Data3.E, ...
    vor.lambda_ci_mts_pore.Data3.E, vor.vor_rms_pore.Data3.E, ...
    vor.lambda_ci_rms_pore.Data3.E, vor.vor_rms_ms_pore.Data3.E, ...
    vor.lambda_ci_rms_ms_pore.Data3.E, vor.N_lambda_ci_pore.Data3.E] ...
    = poreMeanVorticity (x, y, vor.vor_2D.Data3, vor.lambda_ci2D.Data3, imin.Data3.E, imax.Data3.E, ...
    jmin.Data3.E, jmax.Data3.E, kmax1, mm_to_pix);

[vor.circ_vel_pore.Data3.E, vor.circ_vor_pore.Data3.E, vor.circ_vor_mt_pore.Data3.E, ...
    vor.circ_vel_mt_pore.Data3.E, vor.circ_vor_rms_pore.Data3.E, vor.circ_vel_rms_pore.Data3.E] = ...
    circulation (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, vor.vor_2D.Data3, imin.Data3.E, ...
    imax.Data3.E, jmin.Data3.E, jmax.Data3.E, kmax1, nodes.Data3.E, mm_to_pix);

%LES 5*5
[vor.vorles5g_mt_pore.Data3.E, vor.lambda_ciles5g_mt_pore.Data3.E, vor.vorles5g_mts_pore.Data3.E, ...
    vor.lambda_ciles5g_mts_pore.Data3.E, vor.vorles5g_rms_pore.Data3.E, ...
    vor.lambda_ciles5g_rms_pore.Data3.E, vor.vorles5g_rms_ms_pore.Data3.E, ...
    vor.lambda_ciles5g_rms_ms_pore.Data3.E, vor.N_lambda_ciles5g_pore.Data3.E] ...
    = poreLESVorticity (x, y, vor.vorles5g_2D.Data3, vor.lambda_ciles5g_2D.Data3, imin.Data3.E, imax.Data3.E, ...
    jmin.Data3.E, jmax.Data3.E, kmax1, mm_to_pix);

%LES 3*3
[vor.vorles3g_mt_pore.Data3.E, vor.lambda_ciles3g_mt_pore.Data3.E, vor.vorles3g_mts_pore.Data3.E, ...
    vor.lambda_ciles3g_mts_pore.Data3.E, vor.vorles3g_rms_pore.Data3.E, ...
    vor.lambda_ciles3g_rms_pore.Data3.E, vor.vorles3g_rms_ms_pore.Data3.E, ...
    vor.lambda_ciles3g_rms_ms_pore.Data3.E, vor.N_lambda_ciles3g_pore.Data3.E] ...
    = poreMeanVorticity (x, y, vor.vorles3g_2D.Data3, vor.lambda_ciles3g_2D.Data3, imin.Data3.E, imax.Data3.E, ...
    jmin.Data3.E, jmax.Data3.E, kmax1, mm_to_pix);


%% Gamma1 and Gamma2
[vor.Gamma1_pore.Data3.E,vor.Gamma2_pore.Data3.E, vor.Gamma1_mt_pore.Data3.E, vor.Gamma2_mt_pore.Data3.E,...
    vor.Gamma1_mts_pore.Data3.E, vor.Gamma2_mts_pore.Data3.E, vor.Gamma1_rms_pore.Data3.E, ...
    vor.Gamma2_rms_pore.Data3.E, vor.Gamma1_rms_ms_pore.Data3.E, vor.Gamma2_rms_ms_pore.Data3.E,...
    vor.N_Gamma1.Data3.E, vor.N_Gamma2.Data3.E, vor.Gamma2_Shear_mts_pore.Data3.E, vor.Gamma2_Rotation_mts_pore.Data3.E,...
    vor.Gamma2_Shear_rms_ms_pore.Data3.E, vor.Gamma2_Rotation_rms_ms_pore.Data3.E,vor.N_vortex_G2_Shear_pore.Data3.E,...
    vor.N_vortex_G2_Rotation_pore.Data3.E] = Gamma1 (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, imin.Data3.E, imax.Data3.E, ...
    jmin.Data3.E, jmax.Data3.E, kmax1, jmax1, imax1);

%% Vortex size Gamma2

[vor.Gamma2_size.Data3, vor.Gamma2_size_mst.Data3, vor.Gamma2_N_vortex.Data3,...
    vor.Gamma2_size_Pore.Data3.E, vor.Gamma2_size_mst_Pore.Data3.E, vor.Gamma2_N_vortex_Pore.Data3.E] ...
    = Gamma2Size (vor.Gamma2.Data3, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E, kmax1, pix_to_mm, 'pore');

%%

% Integral time scale
[scale.t.T_u2D_pore.Data3.E, scale.t.T_v2D_pore.Data3.E, scale.t.T_u_ITS_pore.Data3.E, scale.t.T_v_ITS_pore.Data3.E] = ...
    IntTimeScalePore (scale.t.T_u2D.Data3, scale.t.T_v2D.Data3, imin.Data3.E, imax.Data3.E,...
    jmin.Data3.E, jmax.Data3.E);

% Integral length scale
lag.L.Data3.E = imax.Data3.E - imin.Data3.E ; 
[scale.autocorr_L_V_mst_pore.Data3.E,scale.ILS.Data3.E] = ...
    IntLengthScalePore (vel.V_mm2D.Data3, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E, kmax1,lag.L.Data3.E, y.Del_y);


% TKE budget
% Production, TKE Advection
[turb.Production.P_linear_pore.Data3.E, turb.Production.P_shear_pore.Data3.E, turb.KE_advection_pore.Data3.E, ...
    turb.Net_grad.grad_x_pore.Data3.E, turb.Net_grad.grad_y_pore.Data3.E, turb.Net_grad.grad_z_pore.Data3.E] = ...
    TKEbudgetPore (vel.U_mt2D.Data3, vel.V_mt2D.Data3, turb.Production.Sxx.Data3,...
    turb.Production.Syy.Data3, turb.Production.Szz.Data3, turb.Production.Sxy.Data3, vel.u_var2D.Data3, ...
    vel.v_var2D.Data3, vel.uv_var2D.Data3, vor.du_var_dx.Data3, vor.du_var_dy.Data3, ...
    vor.dv_var_dx.Data3, vor.dv_var_dy.Data3, turb.Production.sxx.Data3, turb.Production.sxy.Data3, ...
    turb.Production.syy.Data3, vel.u_t2D.Data3, vel.v_t2D.Data3, imin.Data3.E, imax.Data3.E,jmin.Data3.E, jmax.Data3.E,specs.fluid.kvisc);

% Net Gradient transport
[turb.d_Netgrad_x_dx_pore.Data3.E, turb.d_Netgrad_x_dy_pore.Data3.E, turb.d_Netgrad_y_dx_pore.Data3.E, turb.d_Netgrad_y_dy_pore.Data3.E] = ...
    VGT (x, y, turb.Net_grad.grad_x_pore.Data3.E, turb.Net_grad.grad_y_pore.Data3.E, jmax.Data3.E-jmin.Data3.E+1, imax.Data3.E-imin.Data3.E+1, 1, mm_to_pix);

[turb.d_Netgrad_z_dx_pore.Data3.E, turb.d_Netgrad_z_dy_pore.Data3.E, turb.d_Netgrad_y_dx_pore.Data3.E, turb.d_Netgrad_y_dy_pore.Data3.E] = ...
    VGT (x, y, turb.Net_grad.grad_z_pore.Data3.E, turb.Net_grad.grad_y_pore.Data3.E, jmax.Data3.E-jmin.Data3.E+1, imax.Data3.E-imin.Data3.E+1, 1, mm_to_pix);

turb.Net_grad.total_trans_pore.Data3.E = nanmean(nanmean(turb.d_Netgrad_x_dx_pore.Data3.E + turb.d_Netgrad_y_dy_pore.Data3.E + turb.d_Netgrad_z_dx_pore.Data3.E));

% TKE
% Turbulent kinetic energy
turb.TKE_2D_pore.Data3.E = (1/2)*(2*vel.u_var2D.Data3(imin.Data3.E:imax.Data3.E,jmin.Data3.E:jmax.Data3.E) + vel.v_var2D.Data3(imin.Data3.E:imax.Data3.E,jmin.Data3.E:jmax.Data3.E));
turb.TKE_ms_pore.Data3.E = nanmean(nanmean(turb.TKE_2D_pore.Data3.E));

% Tortuosity
dispers.Tortuosity_pore.Data3.E = nansum(nansum(nansum(vel.Mag_mm2D.Data3(imin.Data3.E:imax.Data3.E,jmin.Data3.E:jmax.Data3.E))))./nansum(nansum(nansum(vel.V_mm2D.Data3(imin.Data3.E:imax.Data3.E,jmin.Data3.E:jmax.Data3.E))));

% Eulerian variance
% 2D field
vel.U_var2D_pore.Data3.E = (abs(vel.U_mt2D_NaN.Data3(imin.Data3.E:imax.Data3.E,jmin.Data3.E:jmax.Data3.E)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;
vel.V_var2D_pore.Data3.E = (abs(vel.V_mt2D_NaN.Data3(imin.Data3.E:imax.Data3.E,jmin.Data3.E:jmax.Data3.E)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;

% Averaged non-dimensionalized Eulerian variance based on Vishal's turbulent paper
vel.U_var_mts_Nondim_pore.Data3.E = nanmean(nanmean(vel.U_var2D_pore.Data3.E))/ abs(vel.Mag_mts2D_NaN.Data3)^2;
vel.V_var_mts_Nondim_pore.Data3.E = nanmean(nanmean(vel.V_var2D_pore.Data3.E))/ abs(vel.Mag_mts2D_NaN.Data3)^2;

%% Lagrangian velocity variance

vel.U_Lagrang_var2D_pore.Data3.E = nanmean(nanmean(vel.u_var_pore.Data3.E)) + nanmean(nanmean(vel.U_var2D_pore.Data3.E));
vel.V_Lagrang_var2D_pore.Data3.E = nanmean(nanmean(vel.v_var_pore.Data3.E)) + nanmean(nanmean(vel.V_var2D_pore.Data3.E));

%% Pore F ******************************************************************

imin.Data3.F = 31; 
imax.Data3.F = 41;
jmin.Data3.F = 57;
jmax.Data3.F = 69;
nodes.Data3.F = 13; %11*13

[vel.U_mt_pore.Data3.F,vel.V_mt_pore.Data3.F,vel.Mag_mt_pore.Data3.F,...
    vel.U_mts_pore.Data3.F,vel.V_mts_pore.Data3.F, ...
    vel.Mag_mts_pore.Data3.F,vel.u_rms_pore.Data3.F, ...
    vel.v_rms_pore.Data3.F, vel.uv_rms_pore.Data3.F, vel.u_rms_ms_pore.Data3.F, ...
    vel.v_rms_ms_pore.Data3.F, vel.uv_rms_ms_pore.Data3.F, vel.u_var_pore.Data3.F,...
    vel.v_var_pore.Data3.F, vel.uv_var_pore.Data3.F, vel.u_var_ms_pore.Data3.F,...
    vel.v_var_ms_pore.Data3.F, vel.uv_var_ms_pore.Data3.F,...
    N_vectors.Data3.F] = poreMeanVelocity (x, y, vel.U_mm2D.Data3, ...
    vel.V_mm2D.Data3, vel.Mag_mm2D.Data3, imin.Data3.F, imax.Data3.F, ...
    jmin.Data3.F, jmax.Data3.F, kmax1, mm_to_pix);

%% Plotting velocities

TimeavgplotPore (vel.U_mt2D.Data3,vel.V_mt2D.Data3,vel.Mag_mts_pore.Data3.F, ...
    x,y, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F);
saveas(gcf,'Vmt_PoreF_Data3.png');
%% Vector of LES
VectorPlotPore (vel.U_les3gauss2D.Data3,vel.V_les3gauss2D.Data3, ...
    x,y, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F);
saveas(gcf,'Vector_LES3g_PoreF_Data3.png');
%% Vorticity
VorplotPore (vor.vorles3g_2D_mt.Data3,vel.Mag_mts_pore.Data3.F,5.3, x,y, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F);
saveas(gcf,'vor_PoreF_Data3.png');  

%% Lambda_Ci
kk = 13;
LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),scale.t.T_v_ITS.Data3, x,y, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F);
saveas(gcf,'LambdaCiT13_PoreF_Data3.png');
%% Vorticity and velocity
close all
h = figure; 
for k = 20:120  %vel.Mag_mts_pore.Data3.F
    
    VorplotPore (vor.vorles3g_2D.Data3(:,:,k),47,10, x,y, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F);
%     hold on
    VectorPlotPore (vel.U_les5gauss2D.Data3(:,:,k),vel.V_les5gauss2D.Data3(:,:,k), x,y, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F);
    saveas(h,sprintf('Vor_Data3_PoreF_%d.png',k)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end
% imageFolder2mpeg('Vor_Data3_PoreA','frameRate',11)

%% LambdaCi and velocity

close all
h = figure;
for kk = 20:120 %scale.t.T_v_ITS.Data3
    
    LambdaCiPlotPore (vor.lambda_ciles3g_2D.Data3(:,:,kk),0.002, x,y, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F);
%     hold on
%     VectorPlotPore (vel.U_les3gauss2D.Data3(:,:,kk),vel.V_les3gauss2D.Data3(:,:,kk), x,y, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F);
    saveas(h,sprintf('LambdaCi_Data3_PoreF_%d.png',kk)); % will create FIG1, FIG2, FIG3 etc. You can save your data in other formats such as pdf as well...

    hold off
end

%%

% vortex size 3*3 and 5*5
[vor.lambda_ciles3g_2D_size_Pore.Data3.F, vor.lambda_ciles3g_2D_size_mst_Pore.Data3.F, vor.lambda_ciles3g_2D_N_vortex_Pore.Data3.F] ...
    = LambdaCiSizePore (vor.lambda_ciles3g_2D.Data3,imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F, kmax1, pix_to_mm);

[vor.lambda_ciles5g_2D_size_Pore.Data3.F, vor.lambda_ciles5g_2D_size_mst_Pore.Data3.F, vor.lambda_ciles5g_2D_N_vortex_Pore.Data3.F] ...
    = LambdaCiSizePore (vor.lambda_ciles5g_2D.Data3,imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F, kmax1, pix_to_mm);
%%

%
[vor.vor_mt_pore.Data3.F, vor.lambda_ci_mt_pore.Data3.F, vor.vor_mts_pore.Data3.F, ...
    vor.lambda_ci_mts_pore.Data3.F, vor.vor_rms_pore.Data3.F, ...
    vor.lambda_ci_rms_pore.Data3.F, vor.vor_rms_ms_pore.Data3.F, ...
    vor.lambda_ci_rms_ms_pore.Data3.F, vor.N_lambda_ci_pore.Data3.F] ...
    = poreMeanVorticity (x, y, vor.vor_2D.Data3, vor.lambda_ci2D.Data3, imin.Data3.F, imax.Data3.F, ...
    jmin.Data3.F, jmax.Data3.F, kmax1, mm_to_pix);

[vor.circ_vel_pore.Data3.F, vor.circ_vor_pore.Data3.F, vor.circ_vor_mt_pore.Data3.F, ...
    vor.circ_vel_mt_pore.Data3.F, vor.circ_vor_rms_pore.Data3.F, vor.circ_vel_rms_pore.Data3.F] = ...
    circulation (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, vor.vor_2D.Data3, imin.Data3.F, ...
    imax.Data3.F, jmin.Data3.F, jmax.Data3.F, kmax1, nodes.Data3.F, mm_to_pix);

%LES 5*5
[vor.vorles5g_mt_pore.Data3.F, vor.lambda_ciles5g_mt_pore.Data3.F, vor.vorles5g_mts_pore.Data3.F, ...
    vor.lambda_ciles5g_mts_pore.Data3.F, vor.vorles5g_rms_pore.Data3.F, ...
    vor.lambda_ciles5g_rms_pore.Data3.F, vor.vorles5g_rms_ms_pore.Data3.F, ...
    vor.lambda_ciles5g_rms_ms_pore.Data3.F, vor.N_lambda_ciles5g_pore.Data3.F] ...
    = poreLESVorticity (x, y, vor.vorles5g_2D.Data3, vor.lambda_ciles5g_2D.Data3, imin.Data3.F, imax.Data3.F, ...
    jmin.Data3.F, jmax.Data3.F, kmax1, mm_to_pix);

%LES 3*3
[vor.vorles3g_mt_pore.Data3.F, vor.lambda_ciles3g_mt_pore.Data3.F, vor.vorles3g_mts_pore.Data3.F, ...
    vor.lambda_ciles3g_mts_pore.Data3.F, vor.vorles3g_rms_pore.Data3.F, ...
    vor.lambda_ciles3g_rms_pore.Data3.F, vor.vorles3g_rms_ms_pore.Data3.F, ...
    vor.lambda_ciles3g_rms_ms_pore.Data3.F, vor.N_lambda_ciles3g_pore.Data3.F] ...
    = poreMeanVorticity (x, y, vor.vorles3g_2D.Data3, vor.lambda_ciles3g_2D.Data3, imin.Data3.F, imax.Data3.F, ...
    jmin.Data3.F, jmax.Data3.F, kmax1, mm_to_pix);

%% Gamma1 and Gamma2
[vor.Gamma1_pore.Data3.F,vor.Gamma2_pore.Data3.F, vor.Gamma1_mt_pore.Data3.F, vor.Gamma2_mt_pore.Data3.F,...
    vor.Gamma1_mts_pore.Data3.F, vor.Gamma2_mts_pore.Data3.F, vor.Gamma1_rms_pore.Data3.F, ...
    vor.Gamma2_rms_pore.Data3.F, vor.Gamma1_rms_ms_pore.Data3.F, vor.Gamma2_rms_ms_pore.Data3.F,...
    vor.N_Gamma1.Data3.F, vor.N_Gamma2.Data3.F, vor.Gamma2_Shear_mts_pore.Data3.F, vor.Gamma2_Rotation_mts_pore.Data3.F,...
    vor.Gamma2_Shear_rms_ms_pore.Data3.F, vor.Gamma2_Rotation_rms_ms_pore.Data3.F,vor.N_vortex_G2_Shear_pore.Data3.F,...
    vor.N_vortex_G2_Rotation_pore.Data3.F] = Gamma1 (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, imin.Data3.F, imax.Data3.F, ...
    jmin.Data3.F, jmax.Data3.F, kmax1, jmax1, imax1);

%% Vortex size Gamma2

[vor.Gamma2_size.Data3, vor.Gamma2_size_mst.Data3, vor.Gamma2_N_vortex.Data3,...
    vor.Gamma2_size_Pore.Data3.F, vor.Gamma2_size_mst_Pore.Data3.F, vor.Gamma2_N_vortex_Pore.Data3.F] ...
    = Gamma2Size (vor.Gamma2.Data3, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F, kmax1, pix_to_mm, 'pore');

%%

% Integral time scale
[scale.t.T_u2D_pore.Data3.F, scale.t.T_v2D_pore.Data3.F, scale.t.T_u_ITS_pore.Data3.F, scale.t.T_v_ITS_pore.Data3.F] = ...
    IntTimeScalePore (scale.t.T_u2D.Data3, scale.t.T_v2D.Data3, imin.Data3.F, imax.Data3.F,...
    jmin.Data3.F, jmax.Data3.F);

% Integral length scale
lag.L.Data3.F = imax.Data3.F - imin.Data3.F ; 
[scale.autocorr_L_V_mst_pore.Data3.F,scale.ILS.Data3.F] = ...
    IntLengthScalePore (vel.V_mm2D.Data3, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F, kmax1,lag.L.Data3.F, y.Del_y);

% TKE budget
% Production, TKE Advection
[turb.Production.P_linear_pore.Data3.F, turb.Production.P_shear_pore.Data3.F, turb.KE_advection_pore.Data3.F, ...
    turb.Net_grad.grad_x_pore.Data3.F, turb.Net_grad.grad_y_pore.Data3.F, turb.Net_grad.grad_z_pore.Data3.F] = ...
    TKEbudgetPore (vel.U_mt2D.Data3, vel.V_mt2D.Data3, turb.Production.Sxx.Data3,...
    turb.Production.Syy.Data3, turb.Production.Szz.Data3, turb.Production.Sxy.Data3, vel.u_var2D.Data3, ...
    vel.v_var2D.Data3, vel.uv_var2D.Data3, vor.du_var_dx.Data3, vor.du_var_dy.Data3, ...
    vor.dv_var_dx.Data3, vor.dv_var_dy.Data3, turb.Production.sxx.Data3, turb.Production.sxy.Data3, ...
    turb.Production.syy.Data3, vel.u_t2D.Data3, vel.v_t2D.Data3, imin.Data3.F, imax.Data3.F,jmin.Data3.F, jmax.Data3.F,specs.fluid.kvisc);

% Net Gradinet Transport
[turb.d_Netgrad_x_dx_pore.Data3.F, turb.d_Netgrad_x_dy_pore.Data3.F, turb.d_Netgrad_y_dx_pore.Data3.F, turb.d_Netgrad_y_dy_pore.Data3.F] = ...
    VGT (x, y, turb.Net_grad.grad_x_pore.Data3.F, turb.Net_grad.grad_y_pore.Data3.F, jmax.Data3.F-jmin.Data3.F+1, imax.Data3.F-imin.Data3.F+1, 1, mm_to_pix);

[turb.d_Netgrad_z_dx_pore.Data3.F, turb.d_Netgrad_z_dy_pore.Data3.F, turb.d_Netgrad_y_dx_pore.Data3.F, turb.d_Netgrad_y_dy_pore.Data3.F] = ...
    VGT (x, y, turb.Net_grad.grad_z_pore.Data3.F, turb.Net_grad.grad_y_pore.Data3.F, jmax.Data3.F-jmin.Data3.F+1, imax.Data3.F-imin.Data3.F+1, 1, mm_to_pix);

turb.Net_grad.total_trans_pore.Data3.F = nanmean(nanmean(turb.d_Netgrad_x_dx_pore.Data3.F + turb.d_Netgrad_y_dy_pore.Data3.F + turb.d_Netgrad_z_dx_pore.Data3.F));

% TKE
% Turbulent kinetic energy
turb.TKE_2D_pore.Data3.F = (1/2)*(2*vel.u_var2D.Data3(imin.Data3.F:imax.Data3.F,jmin.Data3.F:jmax.Data3.F) + vel.v_var2D.Data3(imin.Data3.F:imax.Data3.F,jmin.Data3.F:jmax.Data3.F));
turb.TKE_ms_pore.Data3.F = nanmean(nanmean(turb.TKE_2D_pore.Data3.F));

% Tortuosity
dispers.Tortuosity_pore.Data3.F = nansum(nansum(nansum(vel.Mag_mm2D.Data3(imin.Data3.F:imax.Data3.F,jmin.Data3.F:jmax.Data3.F))))./nansum(nansum(nansum(vel.V_mm2D.Data3(imin.Data3.F:imax.Data3.F,jmin.Data3.F:jmax.Data3.F))));

% Eulerian variance
% 2D field
vel.U_var2D_pore.Data3.F = (abs(vel.U_mt2D_NaN.Data3(imin.Data3.F:imax.Data3.F,jmin.Data3.F:jmax.Data3.F)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;
vel.V_var2D_pore.Data3.F = (abs(vel.V_mt2D_NaN.Data3(imin.Data3.F:imax.Data3.F,jmin.Data3.F:jmax.Data3.F)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;

% Averaged non-dimensionalized Eulerian variance based on Vishal's turbulent paper
vel.U_var_mts_Nondim_pore.Data3.F = nanmean(nanmean(vel.U_var2D_pore.Data3.F))/ abs(vel.Mag_mts2D_NaN.Data3)^2;
vel.V_var_mts_Nondim_pore.Data3.F = nanmean(nanmean(vel.V_var2D_pore.Data3.F))/ abs(vel.Mag_mts2D_NaN.Data3)^2;

%% Lagrangian velocity variance

vel.U_Lagrang_var2D_pore.Data3.F = nanmean(nanmean(vel.u_var_pore.Data3.F)) + nanmean(nanmean(vel.U_var2D_pore.Data3.F));
vel.V_Lagrang_var2D_pore.Data3.F = nanmean(nanmean(vel.v_var_pore.Data3.F)) + nanmean(nanmean(vel.V_var2D_pore.Data3.F));




%% Top of Bed **************************************************************

imin.Data3.Top = 43; 
imax.Data3.Top = 61;
jmin.Data3.Top = 1;
jmax.Data3.Top = 77;
nodes.Data3.Top = 13; %11*13

[vel.U_mt_pore.Data3.Top,vel.V_mt_pore.Data3.Top,vel.Mag_mt_pore.Data3.Top,...
    vel.U_mts_pore.Data3.Top,vel.V_mts_pore.Data3.Top, ...
    vel.Mag_mts_pore.Data3.Top,vel.u_rms_pore.Data3.Top, ...
    vel.v_rms_pore.Data3.Top, vel.uv_rms_pore.Data3.Top, vel.u_rms_ms_pore.Data3.Top, ...
    vel.v_rms_ms_pore.Data3.Top, vel.uv_rms_ms_pore.Data3.Top, vel.u_var_pore.Data3.Top,...
    vel.v_var_pore.Data3.Top, vel.uv_var_pore.Data3.Top, vel.u_var_ms_pore.Data3.Top,...
    vel.v_var_ms_pore.Data3.Top, vel.uv_var_ms_pore.Data3.Top,...
    N_vectors.Data3.Top] = poreMeanVelocity (x, y, vel.U_mm2D.Data3, ...
    vel.V_mm2D.Data3, vel.Mag_mm2D.Data3, imin.Data3.Top, imax.Data3.Top, ...
    jmin.Data3.Top, jmax.Data3.Top, kmax1, mm_to_pix);

%% Plotting velocities

TimeavgplotPore (vel.U_mt2D.Data3,vel.V_mt2D.Data3,vel.Mag_mts_pore.Data3.Top, ...
    x,y, imin.Data3.Top, imax.Data3.Top,jmin.Data3.Top, jmax.Data3.Top);
saveas(gcf,'Vmt_Pore_Top_Data3.png');

%%

% vortex size 3*3 and 5*5
[vor.lambda_ciles3g_2D_size_Pore.Data3.Top, vor.lambda_ciles3g_2D_size_mst_Pore.Data3.Top, vor.lambda_ciles3g_2D_N_vortex_Pore.Data3.Top] ...
    = LambdaCiSizePore (vor.lambda_ciles3g_2D.Data3,imin.Data3.Top, imax.Data3.Top,jmin.Data3.Top, jmax.Data3.Top, kmax1, pix_to_mm);

[vor.lambda_ciles5g_2D_size_Pore.Data3.Top, vor.lambda_ciles5g_2D_size_mst_Pore.Data3.Top, vor.lambda_ciles5g_2D_N_vortex_Pore.Data3.Top] ...
    = LambdaCiSizePore (vor.lambda_ciles5g_2D.Data3,imin.Data3.Top, imax.Data3.Top,jmin.Data3.Top, jmax.Data3.Top, kmax1, pix_to_mm);
%%

[vor.vor_mt_pore.Data3.Top, vor.lambda_ci_mt_pore.Data3.Top, vor.vor_mts_pore.Data3.Top, ...
    vor.lambda_ci_mts_pore.Data3.Top, vor.vor_rms_pore.Data3.Top, ...
    vor.lambda_ci_rms_pore.Data3.Top, vor.vor_rms_ms_pore.Data3.Top, ...
    vor.lambda_ci_rms_ms_pore.Data3.Top, vor.N_lambda_ci_pore.Data3.Top] ...
    = poreMeanVorticity (x, y, vor.vor_2D.Data3, vor.lambda_ci2D.Data3, imin.Data3.Top, imax.Data3.Top, ...
    jmin.Data3.Top, jmax.Data3.Top, kmax1, mm_to_pix);

[vor.circ_vel_pore.Data3.Top, vor.circ_vor_pore.Data3.Top, vor.circ_vor_mt_pore.Data3.Top, ...
    vor.circ_vel_mt_pore.Data3.Top, vor.circ_vor_rms_pore.Data3.Top, vor.circ_vel_rms_pore.Data3.Top] = ...
    circulation (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, vor.vor_2D.Data3, imin.Data3.Top, ...
    imax.Data3.Top, jmin.Data3.Top, jmax.Data3.Top, kmax1, nodes.Data3.Top, mm_to_pix);

%LES 5*5
[vor.vorles5g_mt_pore.Data3.Top, vor.lambda_ciles5g_mt_pore.Data3.Top, vor.vorles5g_mts_pore.Data3.Top, ...
    vor.lambda_ciles5g_mts_pore.Data3.Top, vor.vorles5g_rms_pore.Data3.Top, ...
    vor.lambda_ciles5g_rms_pore.Data3.Top, vor.vorles5g_rms_ms_pore.Data3.Top, ...
    vor.lambda_ciles5g_rms_ms_pore.Data3.Top, vor.N_lambda_ciles5g_pore.Data3.Top] ...
    = poreLESVorticity (x, y, vor.vorles5g_2D.Data3, vor.lambda_ciles5g_2D.Data3, imin.Data3.Top, imax.Data3.Top, ...
    jmin.Data3.Top, jmax.Data3.Top, kmax1, mm_to_pix);

%LES 3*3
[vor.vorles3g_mt_pore.Data3.Top, vor.lambda_ciles3g_mt_pore.Data3.Top, vor.vorles3g_mts_pore.Data3.Top, ...
    vor.lambda_ciles3g_mts_pore.Data3.Top, vor.vorles3g_rms_pore.Data3.Top, ...
    vor.lambda_ciles3g_rms_pore.Data3.Top, vor.vorles3g_rms_ms_pore.Data3.Top, ...
    vor.lambda_ciles3g_rms_ms_pore.Data3.Top, vor.N_lambda_ciles3g_pore.Data3.Top] ...
    = poreMeanVorticity (x, y, vor.vorles3g_2D.Data3, vor.lambda_ciles3g_2D.Data3, imin.Data3.Top, imax.Data3.Top, ...
    jmin.Data3.Top, jmax.Data3.Top, kmax1, mm_to_pix);

%% Gamma1 and Gamma2
[vor.Gamma1_pore.Data3.Top,vor.Gamma2_pore.Data3.Top, vor.Gamma1_mt_pore.Data3.Top, vor.Gamma2_mt_pore.Data3.Top,...
    vor.Gamma1_mts_pore.Data3.Top, vor.Gamma2_mts_pore.Data3.Top, vor.Gamma1_rms_pore.Data3.Top, ...
    vor.Gamma2_rms_pore.Data3.Top, vor.Gamma1_rms_ms_pore.Data3.Top, vor.Gamma2_rms_ms_pore.Data3.Top,...
    vor.N_Gamma1.Data3.Top, vor.N_Gamma2.Data3.Top, vor.Gamma2_Shear_mts_pore.Data3.Top, vor.Gamma2_Rotation_mts_pore.Data3.Top,...
    vor.Gamma2_Shear_rms_ms_pore.Data3.Top, vor.Gamma2_Rotation_rms_ms_pore.Data3.Top,vor.N_vortex_G2_Shear_pore.Data3.Top,...
    vor.N_vortex_G2_Rotation_pore.Data3.Top] = Gamma1 (x, y, vel.U_mm2D.Data3, vel.V_mm2D.Data3, imin.Data3.Top, imax.Data3.Top, ...
    jmin.Data3.Top, jmax.Data3.Top, kmax1, jmax1, imax1);
%%

% Integral time scale
[scale.t.T_u2D_pore.Data3.Top, scale.t.T_v2D_pore.Data3.Top, scale.t.T_u_ITS_pore.Data3.Top, scale.t.T_v_ITS_pore.Data3.Top] = ...
    IntTimeScalePore (scale.t.T_u2D.Data3, scale.t.T_v2D.Data3, imin.Data3.Top, imax.Data3.Top,...
    jmin.Data3.Top, jmax.Data3.Top);


% Integral length scale
lag.L.Data3.Top = imax.Data3.Top - imin.Data3.Top ; 
[scale.Toputocorr_L_V_mst_pore.Data3.Top,scale.ILS.Data3.Top] = ...
    IntLengthScalePore (vel.V_mm2D.Data3, imin.Data3.Top, imax.Data3.Top,jmin.Data3.Top, jmax.Data3.Top, kmax1,lag.L.Data3.Top, y.Del_y);

% TKE budget
% Production, TKE Advection
[turb.Production.P_linear_pore.Data3.Top, turb.Production.P_shear_pore.Data3.Top, turb.KE_advection_pore.Data3.Top, ...
    turb.Net_grad.grad_x_pore.Data3.Top, turb.Net_grad.grad_y_pore.Data3.Top, turb.Net_grad.grad_z_pore.Data3.Top] = ...
    TKEbudgetPore (vel.U_mt2D.Data3, vel.V_mt2D.Data3, turb.Production.Sxx.Data3,...
    turb.Production.Syy.Data3, turb.Production.Szz.Data3, turb.Production.Sxy.Data3, vel.u_var2D.Data3, ...
    vel.v_var2D.Data3, vel.uv_var2D.Data3, vor.du_var_dx.Data3, vor.du_var_dy.Data3, ...
    vor.dv_var_dx.Data3, vor.dv_var_dy.Data3, turb.Production.sxx.Data3, turb.Production.sxy.Data3, ...
    turb.Production.syy.Data3, vel.u_t2D.Data3, vel.v_t2D.Data3, imin.Data3.Top, imax.Data3.Top,jmin.Data3.Top, jmax.Data3.Top,specs.fluid.kvisc);

% Net Gradient transport
[turb.d_Netgrad_x_dx_pore.Data3.Top, turb.d_Netgrad_x_dy_pore.Data3.Top, turb.d_Netgrad_y_dx_pore.Data3.Top, turb.d_Netgrad_y_dy_pore.Data3.Top] = ...
    VGT (x, y, turb.Net_grad.grad_x_pore.Data3.Top, turb.Net_grad.grad_y_pore.Data3.Top, jmax.Data3.Top-jmin.Data3.Top+1, imax.Data3.Top-imin.Data3.Top+1, 1, mm_to_pix);

[turb.d_Netgrad_z_dx_pore.Data3.Top, turb.d_Netgrad_z_dy_pore.Data3.Top, turb.d_Netgrad_y_dx_pore.Data3.Top, turb.d_Netgrad_y_dy_pore.Data3.Top] = ...
    VGT (x, y, turb.Net_grad.grad_z_pore.Data3.Top, turb.Net_grad.grad_y_pore.Data3.Top, jmax.Data3.Top-jmin.Data3.Top+1, imax.Data3.Top-imin.Data3.Top+1, 1, mm_to_pix);

turb.Net_grad.total_trans_pore.Data3.Top = nanmean(nanmean(turb.d_Netgrad_x_dx_pore.Data3.Top + turb.d_Netgrad_y_dy_pore.Data3.Top + turb.d_Netgrad_z_dx_pore.Data3.Top));

% TKE
% Turbulent kinetic energy
turb.TKE_2D_pore.Data3.Top = (1/2)*(2*vel.u_var2D.Data3(imin.Data3.Top:imax.Data3.Top,jmin.Data3.Top:jmax.Data3.Top) + vel.v_var2D.Data3(imin.Data3.Top:imax.Data3.Top,jmin.Data3.Top:jmax.Data3.Top));
turb.TKE_ms_pore.Data3.Top = nanmean(nanmean(turb.TKE_2D_pore.Data3.Top));

% Tortuosity
dispers.Tortuosity_pore.Data3.Top = nansum(nansum(nansum(vel.Mag_mm2D.Data3(imin.Data3.Top:imax.Data3.Top,jmin.Data3.Top:jmax.Data3.Top))))./nansum(nansum(nansum(vel.V_mm2D.Data3(imin.Data3.Top:imax.Data3.Top,jmin.Data3.Top:jmax.Data3.Top))));

% Eulerian variance
% 2D field
vel.U_var2D_pore.Data3.Top= (abs(vel.U_mt2D_NaN.Data3(imin.Data3.Top:imax.Data3.Top,jmin.Data3.Top:jmax.Data3.Top)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;
vel.V_var2D_pore.Data3.Top= (abs(vel.V_mt2D_NaN.Data3(imin.Data3.Top:imax.Data3.Top,jmin.Data3.Top:jmax.Data3.Top)) - abs(vel.Mag_mts2D_NaN.Data3)).^2;

% Averaged non-dimensionalized Eulerian variance based on Vishal's turbulent paper
vel.U_var_mts_Nondim_pore.Data3.Top = nanmean(nanmean(vel.U_var2D_pore.Data3.Top))/ abs(vel.Mag_mts2D_NaN.Data3)^2;
vel.V_var_mts_Nondim_pore.Data3.Top = nanmean(nanmean(vel.V_var2D_pore.Data3.Top))/ abs(vel.Mag_mts2D_NaN.Data3)^2;

%% Lagrangian velocity variance

vel.U_Lagrang_var2D_pore.Data3.Top = nanmean(nanmean(vel.u_var_pore.Data3.Top)) + nanmean(nanmean(vel.U_var2D_pore.Data3.Top));
vel.V_Lagrang_var2D_pore.Data3.Top = nanmean(nanmean(vel.v_var_pore.Data3.Top)) + nanmean(nanmean(vel.V_var2D_pore.Data3.Top));

%% Only porous region (bottom)
imin.Data3.Bot = 1; 
imax.Data3.Bot = 45;
jmin.Data3.Bot = 1;
jmax.Data3.Bot = 77;
nodes.Data3.Bot = 13; %11*13

[vel.U_mt_pore.Data3.Bot,vel.V_mt_pore.Data3.Bot,vel.Mag_mt_pore.Data3.Bot,...
    vel.U_mts_pore.Data3.Bot,vel.V_mts_pore.Data3.Bot, ...
    vel.Mag_mts_pore.Data3.Bot,vel.u_rms_pore.Data3.Bot, ...
    vel.v_rms_pore.Data3.Bot, vel.uv_rms_pore.Data3.Bot, vel.u_rms_ms_pore.Data3.Bot, ...
    vel.v_rms_ms_pore.Data3.Bot, vel.uv_rms_ms_pore.Data3.Bot, vel.u_var_pore.Data3.Bot,...
    vel.v_var_pore.Data3.Bot, vel.uv_var_pore.Data3.Bot, vel.u_var_ms_pore.Data3.Bot,...
    vel.v_var_ms_pore.Data3.Bot, vel.uv_var_ms_pore.Data3.Bot,...
    N_vectors.Data3.Bot] = poreMeanVelocity (x, y, vel.U_mm2D.Data3, ...
    vel.V_mm2D.Data3, vel.Mag_mm2D.Data3, imin.Data3.Bot, imax.Data3.Bot, ...
    jmin.Data3.Bot, jmax.Data3.Bot, kmax1, mm_to_pix);

