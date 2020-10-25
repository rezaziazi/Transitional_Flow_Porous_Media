% #########################################################################
% ##############        Circulation in porous media       #################
% ##############               By: Reza Ziazi             #################
% ##############                 Fall 2018                #################
% #########################################################################


% This algorithm uses the inputs of a 2D velocity field found from CFD or PIV
% to find the circulation of the flow over each pore. The circulation is
% computed by two methods of direct vorticity summation and by line integral
% over a closed path for each pore. 
% 
% INPUTS:
% 
% imin, imax, jmin, jmax:
% These re the coordinates of the rectangular area of the pore, and kmax is
% the total number of time steps. 
% 
% nodes:
% The total number of nodes for each rectangular pore
% 
% mm_to_pix:
% The conversion from milimiters to pixel specifically for PIV data

function [circ_vel, circ_vor, circ_vor_mt, circ_vel_mt, circ_vor_rms, circ_vel_rms] = ...
    circulation (x, y, U, V, vor, imin, imax, jmin, jmax, kmax, nodes, mm_to_pix)

% n is the number of nodes and (n-1) is the number of cells
Del_x = (x.X_Nodes2_pix(2)- x.X_Nodes2_pix(1))*mm_to_pix;
Del_y = (y.Y_Nodes2_pix(2)- y.Y_Nodes2_pix(1))*mm_to_pix;



% U = flip(U); V = flip(V); % flipping the U and V upside down to perform circulation 

% Circulation from velocity

if mod (nodes,2) ~= 0 % Simpson rule
    for k = 1:kmax
        % Bottom and top of the field
        if (isnan(U(imin,jmin,k)))~= 1 && isnan(U(imin,jmax,k))~=1 && isnan(U(imax,jmax,k))~=1 && isnan(U(imax,jmin,k))~=1 
            circ_bot(k) = (1/3)*Del_x*(nansum(U(imin,jmin,k)) + nansum(U(imin,jmax,k)) + 4*nansum(U(imin,jmin+1:2:jmax-1,k)) + 2*nansum(U(imin,jmin+2:2:jmax-2,k)));
            circ_top(k) = (1/3)*Del_x*(nansum(U(imax,jmin,k)) + nansum(U(imax,jmax,k)) + 4*nansum(U(imax,jmin+1:2:jmax-1,k)) + 2*nansum(U(imax,jmin+2:2:jmax-2,k)));
            
            % Left and right of the field
            circ_rig(k) = (1/3)*Del_y*(nansum(V(imin,jmax,k)) + nansum(V(imax,jmax,k)) + 4*nansum(V(imin+1:2:imax-1,jmax,k)) + 2*nansum(V(imin+2:2:imax-2,jmax,k)));
            circ_lef(k) = (1/3)*Del_y*(nansum(V(imin,jmin,k)) + nansum(V(imax,jmin,k)) + 4*nansum(V(imin+1:2:imax-1,jmin,k)) + 2*nansum(V(imin+2:2:imax-2,jmin,k)));
            
            circ_vel(k) = circ_bot(k) + circ_rig(k) - (circ_top(k) + circ_lef(k));
        else
            disp('There is NaN in one of the cells for circulation')
        end
    end
else            % Trapezoid rule
    for k = 1:kmax 
        if (isnan(U(imin,jmin,k)))~= 1 && isnan(U(imin,jmax,k))~=1 && isnan(U(imax,jmax,k))~=1 && isnan(U(imax,jmin,k))~=1
        % Bottom and top of the field
        circ_bot(k) = (1/2)*Del_x*(U(imin,jmin,k) + U(imin,jmax,k) + 2*sum(U(imin,jmin+1:jmax-1,k)));
        circ_top(k) = (1/2)*Del_x*(U(imax,jmin,k) + U(imax,jmax,k) + 2*sum(U(imax,jmin+1:jmax-1,k)));
        
        % Left and right of the field
        circ_rig(k) = (1/2)*Del_y*(V(imin,jmax,k) + V(imax,jmax,k) + 2*sum(V(imin+1:imax-1,jmax,k)));
        circ_lef(k) = (1/2)*Del_y*(V(imin,jmin,k) + V(imax,jmin,k) + 2*sum(V(imin+1:imax-1,jmin,k)));
        
        circ_vel(k) = circ_bot(k) + circ_rig(k) - (circ_top(k) + circ_lef(k));
        else
            disp('There is NaN in one of the cells for circulation')
        end
    end
end


vor(find(vor == 0)) = NaN;


% Circulation from vorticity summation
for k = 1:kmax
    circ_vor(k) = nansum(nansum(vor(imin:imax,jmin:jmax,k)));
end

% Time average of total circulation
circ_vor_mt = nanmean(circ_vor);
circ_vel_mt = nanmean(circ_vel);

% circulation fluctuation
circ_vor_fluc = circ_vor - circ_vor_mt;
circ_vel_fluc = circ_vel - circ_vel_mt;

% circulation rms
circ_vor_rms = sqrt(nanmean(circ_vor_fluc(:).^2));
circ_vel_rms = sqrt(nanmean(circ_vel_fluc(:).^2));

end