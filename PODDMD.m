function [q_U_Data1] = PODDMD (U, V, x, imax, jmax, mode) % jmax1

Del_x = (x.X_Nodes2_mm(2)- x.X_Nodes2_mm(1));
U(isnan(U))=0;
V(isnan(V))=0;

A = zeros((imax-jmax),imax);


r = size(U,1); 
c = size(U,2);
% Transform the data into column vectors.
data=[reshape(U,r*c,size(U,3));reshape(V,r*c,size(V,3))];


%%%%%%%%%%%%%%%%%%% POD %%%%%%%%%%%%%%%%%%%%%%
% Perform the POD - this is mean subtracted for POD not for DMD
[Phi ,~, C]=svd(data-repmat(mean(data,2),[1 size(data,2)]),'econ');
PHI_U_contour = reshape(Phi(1:r*c,mode),r,c);
PHI_V_contour = reshape(Phi(r*c+1:end,mode),r,c);


Phi_U_sym = [PHI_U_contour;A];
Phi_V_sym = [PHI_V_contour;A];
[x_mesh,y_mesh] = meshgrid(x.X_Nodes1_mm,x.X_Nodes1_mm);

q_U_Data1 = contourf(x_mesh,y_mesh,Phi_U_sym,50,'edgecolor','none');
axis([x_mesh(1,1) x_mesh(imax,imax) y_mesh(1,1) y_mesh(jmax,jmax)])
xlabel('X [mm]','FontSize',18,'FontName','Times New Roman');ylabel('Y [mm]','FontSize',18,'FontName','Times New Roman')
colorbar; 
shading interp
h = colorbar; 
colormap((jet))

figure
q_V_Data1 = contourf(x_mesh,y_mesh,Phi_V_sym,50,'edgecolor','none');
axis([x_mesh(1,1) x_mesh(imax,imax) y_mesh(1,1) y_mesh(jmax,jmax)])
xlabel('X [mm]','FontSize',18,'FontName','Times New Roman');ylabel('Y [mm]','FontSize',18,'FontName','Times New Roman')
colorbar; 
shading interp
h = colorbar; 
colormap((jet))


% % Plane3 ******************************************************************
% % Circle
% hold on
% filledCircle([8,-3],8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([8.2,11.4],6.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR) 
% hold on
% filledCircle([7.7,23.3],5.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([20.6,5.5],7.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([22,22],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([30.0,35],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([15.3,35],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([33.3,12.4],6.4,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([41,1],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([49,13.5],7.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([51,29.0],6.7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62,12.7],5.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62,-1],8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,24],6.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,37.5],7.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% 
% xlim([Del_x*1 Del_x*(imax-2)])
% ylim([Del_x*2 Del_x*(jmax-4)])
% daspect([1 1 1])
% set(gca,'FontSize',14,'linewidth',1.2)

% Data 9 ******************************************************************
% Circle
hold on
filledCircle([8.2,39],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([8.2,24],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR) 
hold on
filledCircle([22.5,30.5],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([22,14],6.6,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([15,2],6.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([31,2],6.6,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([52,-2],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([42,10.5],6.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([34,22],5.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([48,24],6.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([41,36],7.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([56.5,36.5],8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([62,24],3.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
hold on
filledCircle([60,12],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)

xlim([Del_x*3 Del_x*(imax)])
ylim([Del_x*2 Del_x*(jmax-4)])
daspect([1 1 1])
set(gca,'FontSize',14,'linewidth',1.2)

% % Plane6 ******************************************************************
% % Circle
% hold on
% filledCircle([8.2,39],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([8.2,24],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR) 
% hold on
% filledCircle([22.5,30.5],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([22,14],6.6,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([15,2],6.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([31,2],6.6,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([52,-2],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([42,10.5],6.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([34,22],5.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([48,24],6.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([41,36],6.9,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([56.5,36.5],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62,24],3.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([60,12],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% 
% xlim([Del_x*3 Del_x*(imax-8)])
% ylim([Del_x*2 Del_x*(jmax-4)])
% daspect([1 1 1])
% set(gca,'FontSize',14,'linewidth',1.2)

% % % Plane 5,6 ******************************************************************
% % % Circle
% hold on
% filledCircle([7,23],9,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([14,3],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR) 
% hold on
% filledCircle([25,13],7.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([24.5,28.5],7.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([39,7],7.9,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([39,22],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([42,35],4.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([52,-1],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([53,25],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([59,12],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([59,38],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% 
% 
% xlim([Del_x*2.9 Del_x*(imax-4)])
% ylim([Del_x*2 Del_x*(jmax-4)])
% daspect([1 1 1])
% set(gca,'FontSize',14,'linewidth',1.2)


% % Data 4 ******************************************************************
% % Circle
% hold on
% filledCircle([8.2,34],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([29.5,35.1],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR) 
% hold on
% filledCircle([7.4,22.5],4.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([22,22],6.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([42.2,39.2],4.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([7,11.5],6.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([19.1,4.3],7.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([5,-5],9.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([32.7,11],7.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([42,1],6.4,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,2],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([49,14.1],7.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([50,29],6.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,38.1],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62.5,14.7],5.3,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62.5,25.9],4.7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([17.5,40],3.7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% 
% 
% xlim([Del_x*2 Del_x*(imax-3)])
% ylim([Del_x*3 Del_x*(jmax-4)])
% daspect([1 1 1])
% set(gca,'FontSize',14,'linewidth',1.2)

% % Data1 *******************************************************************
% % Circle
% hold on
% filledCircle([14.75,36],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([29,35.1],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR) 
% hold on
% filledCircle([6,22],4,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([22,22],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([41.7,39.8],1.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([7,11.5],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([19,4.3],7.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([5,-5],9,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([32.7,11],6.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([42,1],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,2],7.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([49,14.1],6,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([50,29],6.3,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,38.1],6.7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62.5,14.7],5.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62.5,25.9],5.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% xlim([Del_x*2 Del_x*(imax-3)])
% ylim([Del_x*3 Del_x*(jmax-4)])
% daspect([1 1 1])
% set(gca,'FontSize',14,'linewidth',1.2)


% % Data 2 ******************************************************************
% % Circle
% hold on
% filledCircle([14.75,36],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([29.5,35.1],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR) 
% hold on
% filledCircle([6.4,22.5],4,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([22,22],6.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([42.2,36.2],3.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([7,11.5],6.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([19.1,4.3],7.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([5,-5],9.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([32.7,11],6.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([42,1],6.4,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,2],7.0,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([49,14.1],5.9,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([50,29],6.2,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,38.1],6.50,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62.5,14.7],5.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62.5,25.9],5.7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([34.1,23],2.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% 
% xlim([Del_x*2.0 Del_x*(imax-3)])
% ylim([Del_x Del_x*(jmax-6)])
% daspect([1 1 1])
% set(gca,'FontSize',14,'linewidth',1.2)


% % Circle
% hold on
% filledCircle([14.75,36],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([29,35.1],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR) 
% hold on
% filledCircle([6,22],4,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([22,22],7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([41.7,39.8],1.8,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([7,11.5],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([19,4.3],7.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([5,-5],9,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([32.7,11],6.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([42,1],6.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,2],7.1,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([49,14.1],6,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([50,29],6.3,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([63,38.1],6.7,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62.5,14.7],5.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% hold on
% filledCircle([62.5,25.9],5.5,1000,'w'); % filledCircle(CENTER,R,N,COLOR)
% 
% % Axis border lines 
% set(gca,'FontSize',14,'linewidth',1.2)
% xL = xlim;
% yL = ylim;
% line([0 0], yL,'Color','k','linewidth',1.2);  %x-axis
% line(xL, [0 0],'Color','k','linewidth',1.2);  %y-axis
% axis([x_mesh(1,1) x_mesh(imax,imax) y_mesh(1,1) y_mesh(jmax-5,jmax-5)],'equal')
% xlim([x_mesh(1,1) x_mesh(imax,imax)])
% ylim([y_mesh(1,1) y_mesh(jmax-5,jmax-5)])

% Generating the 
function h = filledCircle(center,r,N,color)
%---------------------------------------------------------------------------------------------
% FILLEDCIRCLE Filled circle drawing
% 
% filledCircle(CENTER,R,N,COLOR) draws a circle filled with COLOR that 
% has CENTER as its center and R as its radius, by using N points on the 
% periphery.
%
% Usage Examples,
%
% filledCircle([1,3],3,1000,'b'); 
% filledCircle([2,4],2,1000,'r');
%
% Sadik Hava <sadik.hava@gmail.com>
% May, 2010
%
% Inspired by: circle.m [Author: Zhenhai Wang]
%---------------------------------------------------------------------------------------------
THETA=linspace(0,2*pi,N);
RHO=ones(1,N)*r;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
h=fill(X,Y,color);
set(h,'EdgeColor','none');
set(h,'FaceColor',[0.45 0.45 0.45]);
axis square;
axx = gca;
axx.XColor = 'k'; % Red
axx.YColor = 'k'; % Blue
box on
% set(h,'facealpha',.5)
end

end

