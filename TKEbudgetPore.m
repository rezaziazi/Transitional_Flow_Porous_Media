
function [P_linear, P_shear, TKE_advection, Net_grad_grad_x, Net_grad_grad_y, Net_grad_grad_z] = TKEbudgetPore (U_mt2D, V_mt2D, Sxx, Syy, Szz, Sxy, u_var2D, ...
    v_var2D, uv_var2D, du_var_dx, du_var_dy, dv_var_dx, dv_var_dy, sxx, sxy, syy, u_t2D,v_t2D, imin, imax, jmin, jmax,kvisc)

U_mt2D = U_mt2D(imin:imax,jmin:jmax);
V_mt2D = V_mt2D(imin:imax,jmin:jmax);
Sxx = Sxx(imin:imax,jmin:jmax);
Syy = Syy(imin:imax,jmin:jmax);
Szz = Szz(imin:imax,jmin:jmax);
Sxy = Sxy(imin:imax,jmin:jmax);
u_var2D = u_var2D(imin:imax,jmin:jmax);
v_var2D = v_var2D(imin:imax,jmin:jmax);
uv_var2D = uv_var2D(imin:imax,jmin:jmax);
du_var_dx = du_var_dx(imin:imax,jmin:jmax);
du_var_dy = du_var_dy(imin:imax,jmin:jmax);
dv_var_dx = dv_var_dx(imin:imax,jmin:jmax);
dv_var_dy = dv_var_dy(imin:imax,jmin:jmax);


% Production***************************************************************
P_linear = - nanmean(nanmean(u_var2D.*Sxx))...
    - nanmean(nanmean(v_var2D.*Syy))...
    - nanmean(nanmean(u_var2D.*Szz));        % Look at vishal's article: Scale Estimation

P_shear = -3*2*nanmean(nanmean(uv_var2D.*Sxy));



% Kinetic energy advection ************************************************
TKE_advection = 2*nanmean(nanmean(U_mt2D.*((1/2)*du_var_dx)))...
    + nanmean(nanmean(U_mt2D.*((1/2)*dv_var_dx)))...
    + 2*nanmean(nanmean(V_mt2D.*((1/2)*du_var_dy)))...
    + nanmean(nanmean(V_mt2D.*((1/2)*dv_var_dy)));

% Gradient transport ******************************************************

% time averaging the mean terms
for i = imin:imax
    for j = jmin:jmax
        Net_grad_viscous_x(i,j) = 2*kvisc*(nanmean(u_t2D(i,j,:).*sxx(i,j,:))...
            + nanmean(v_t2D(i,j,:).*sxy(i,j,:)) + nanmean(u_t2D(i,j,:).*sxy(i,j,:)));
        Net_grad_viscous_y(i,j) = 2*kvisc*(nanmean(v_t2D(i,j,:).*syy(i,j,:))...
            + 2*nanmean(u_t2D(i,j,:).*sxy(i,j,:)));
        Net_grad_viscous_z(i,j) = 2*kvisc*(nanmean(u_t2D(i,j,:).*sxy(i,j,:))...
            + nanmean(v_t2D(i,j,:).*sxy(i,j,:)) + nanmean(u_t2D(i,j,:).*sxx(i,j,:)));
        Net_grad_fluc_x(i,j) = (1/2)*(2*nanmean(u_t2D(i,j,:).^3) + nanmean(v_t2D(i,j,:).^2.*u_t2D(i,j,:)));
        Net_grad_fluc_y(i,j) = (1/2)*(2*nanmean(u_t2D(i,j,:).^2.*v_t2D(i,j,:)) + nanmean(v_t2D(i,j,:).^3));
        Net_grad_fluc_z(i,j) = (1/2)*(2*nanmean(u_t2D(i,j,:).^3) + nanmean(v_t2D(i,j,:).^2.*u_t2D(i,j,:)));
    end
end

Net_grad_grad_x = Net_grad_fluc_x - Net_grad_viscous_x;
Net_grad_grad_y = Net_grad_fluc_y - Net_grad_viscous_y;
Net_grad_grad_z = Net_grad_fluc_z - Net_grad_viscous_z;

Net_grad_grad_x = Net_grad_grad_x(imin:imax,jmin:jmax);
Net_grad_grad_y = Net_grad_grad_y(imin:imax,jmin:jmax);
Net_grad_grad_z = Net_grad_grad_z(imin:imax,jmin:jmax);


end