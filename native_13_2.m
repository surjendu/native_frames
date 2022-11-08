
%% newton_diagram with respect to 2nd hit without the boundary 
% clear some mass etc are reading from previous run
clear 
close all 
clc

hatoev=  27.211396; % taken from ksu webpage/ twice of ionization potential of H
har_to_ev=  27.211396;
mpvsme=1822.888486424682;

frag={ {'Br81', 1, 1} {'Br81', 1, 1} {'C',1,'H',1,'Br81', 1, 1}     };
[frag_m,frag_m_z,frag_m_z_str]=fragments(frag);
m(1)=frag_m(1);
m(2)=frag_m(2);
m(3)=frag_m(3);
% bin_nw =40/5000*5; %0.008
% bin_nv_ke = .1*2;
% bin_nv_an = 0.5*2;

bin_nw =40/5000; %0.008
bin_nv_ke = .1;
bin_nv_an = 0.5;

 %% Newton plot wrt 2   1 up 3 down
 

clearvars  px1 py1 pz1 px2 py2 pz2 px3 py3 pz3;
clc
close all
% px1=px1(j_KE_all); py1=py1(j_KE_all);pz1=pz1(j_KE_all);
% px2=px2(j_KE_all); py2=py2(j_KE_all);pz2=pz2(j_KE_all);
% px3=px3(j_KE_all); py3=py3(j_KE_all);pz3=pz3(j_KE_all);
p_gated=dlmread('p_gated.csv');

px1=p_gated(:,1); py1=p_gated(:,2);pz1=p_gated(:,3);
px2=p_gated(:,4); py2=p_gated(:,5);pz2=p_gated(:,6);
px3=p_gated(:,7); py3=p_gated(:,8);pz3=p_gated(:,9);


%% lab frame to center of mass transformation
px1_cm= px1 - (m(1)/(m(1)+m(2)+m(3)))*(px1+px2+px3);
px2_cm= px2 - (m(2)/(m(1)+m(2)+m(3)))*(px1+px2+px3);
px3_cm= px3 - (m(3)/(m(1)+m(2)+m(3)))*(px1+px2+px3);

py1_cm= py1 - (m(1)/(m(1)+m(2)+m(3)))*(py1+py2+py3);
py2_cm= py2 - (m(2)/(m(1)+m(2)+m(3)))*(py1+py2+py3);
py3_cm= py3 - (m(3)/(m(1)+m(2)+m(3)))*(py1+py2+py3);

pz1_cm= pz1- (m(1)/(m(1)+m(2)+m(3)))*(pz1+pz2+pz3);
pz2_cm= pz2- (m(2)/(m(1)+m(2)+m(3)))*(pz1+pz2+pz3);
pz3_cm= pz3- (m(3)/(m(1)+m(2)+m(3)))*(pz1+pz2+pz3);
%% renonate 
px1 = px1_cm; px2 = px2_cm; px3 = px3_cm;
py1 = py1_cm; py2 = py2_cm; py3 = py3_cm;
pz1 = pz1_cm; pz2 = pz2_cm; pz3 = pz3_cm;

%%
p1_m=sqrt(px1.*px1 + py1.*py1 + pz1.*pz1);
p2_m=sqrt(px2.*px2 + py2.*py2 + pz2.*pz2);
p3_m=sqrt(px3.*px3 + py3.*py3 + pz3.*pz3);


p1_dot_p2 = (px1.*px2 + py1.*py2 + pz1.*pz2);
p1_dot_p3 = (px1.*px3 + py1.*py3 + pz1.*pz3); 
p2_dot_p3 = (px2.*px3 + py2.*py3 + pz2.*pz3); 

cos_theta_12=p1_dot_p2./(p1_m.*p2_m);
cos_theta_13=p1_dot_p3./(p1_m.*p3_m);
cos_theta_23=p2_dot_p3./(p2_m.*p3_m);

close all;
acute_theta_12=(pi - acos(cos_theta_12))*180/pi;
acute_theta_13=(pi - acos(cos_theta_13))*180/pi;
acute_theta_23=(pi - acos(cos_theta_23))*180/pi;

p3_m_np2=p3_m./p2_m;
p1_m_np2=p1_m./p2_m;

p3_m_np2_x = -(p3_m_np2).*cos(acute_theta_23.*pi/180); %momentum coordinate x
p3_m_np2_y = - (p3_m_np2).*sin(acute_theta_23.*pi/180); %momentum coordinate y

p1_m_np2_x = -(p1_m_np2).*cos(acute_theta_12.*pi/180);
p1_m_np2_y = (p1_m_np2).*sin(acute_theta_12.*pi/180);


p31_x=[p3_m_np2_x;p1_m_np2_x];
p31_y=[p3_m_np2_y;p1_m_np2_y];


%making 2d histogram
% Xedges = min(p31_x(:,1)):(max(p31_x(:,1)) - min(p31_x(:,1)) )/5000:max(p31_x(:,1));
% Yedges = min(p31_y(:,1)):(max(p31_y(:,1)) - min(p31_y(:,1)) )/5000:max(p31_y(:,1));

Xedges = -20: bin_nw: 20;
Yedges = -20: bin_nw: 20;


count_allevt_nt = histcounts2(p31_x(:,1),p31_y(:,1),Xedges,Yedges);

count_mod = max(count_allevt_nt,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([-2 2 -2 2])
% xticks([-2:1:2])
% yticks([-2:1:2])
hold on;
quiver(0,0,1.,0,'-r','LineWidth',2,'MaxHeadSize',3);
set(gca,'FontSize',50)
axis equal;
pbaspect([1 1 1]);

xlabel('rel. P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('rel. P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',50)
set(gca,'colorscale','log');

caxis([ 1 120]);


% annotation('textbox',[0.322 0.24 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.51 0.61 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.84 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');
% 
%% native frame intermediate 3rd and 1st
% close all;
figure

m_123=(m(1)+m(2)+m(3));
m_31=(m(3)+m(1));
mu_31=1/(1/m(3)+1/m(1));

p31_x=mu_31*(px1./m(1)-px3./m(3));
p31_y=mu_31*(py1./m(1)-py3./m(3));
p31_z=mu_31*(pz1./m(1)-pz3./m(3));

p31_m=sqrt(p31_x.^2 + p31_y.^2 + p31_z.^2);

p31_2_x=(m_31/m_123)*px2 - (m(2)/m_123)*(px3+px1);
p31_2_y=(m_31/m_123)*py2 - (m(2)/m_123)*(py3+py1);
p31_2_z=(m_31/m_123)*pz2 - (m(2)/m_123)*(pz3+pz1);
p31_2_m=sqrt(p31_2_x.^2 + p31_2_y.^2 + p31_2_z.^2);

p31_2_dot_p31=(p31_2_x.*p31_x + p31_2_y.*p31_y +p31_2_z.*p31_z );

ke_31 = p31_m.^2/(2*mu_31)*hatoev;
theta_31_2=acos(p31_2_dot_p31./(p31_2_m.*p31_m))*180/pi;

Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;

count_allevt_nv = histcounts2(ke_31(:,1),theta_31_2(:,1),Xedges,Yedges);

count_mod = max(count_allevt_nv,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);
hold on;


set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);


xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');
caxis([ 1 200]);

hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');


% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[20 45 45 20 20], 'k', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[135 160 160 135 135], 'b', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[92 114 114 92 92], '--k', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[66 88 88 66 66], '--b', 'Linewidth', 2)


%%

theta_min=20 ; theta_max=45;  
j_nat_ke31 = ke_31 >= 2.6 & ke_31 <= 4.1 & theta_31_2 >= theta_min & theta_31_2 <= theta_max ; % new value to make compare with dalitz plot



%%
j_native_gate= j_nat_ke31 ;
j_native_gate=j_native_gate;

ke_31_1 = ke_31(j_native_gate );
theta_31_2_1 = theta_31_2(j_native_gate);

px1_1 = px1(j_native_gate);   py1_1 = py1(j_native_gate);   pz1_1 = pz1(j_native_gate);
px2_1 = px2(j_native_gate);   py2_1 = py2(j_native_gate);   pz2_1 = pz2(j_native_gate);
px3_1 = px3(j_native_gate);   py3_1 = py3(j_native_gate);   pz3_1 = pz3(j_native_gate);

p31_x_1 = p31_x(j_native_gate); p31_y_1 = p31_y(j_native_gate); p31_z_1 = p31_z(j_native_gate);

p31_2_x_1 = p31_2_x(j_native_gate); p31_2_y_1 = p31_2_y(j_native_gate); p31_2_z_1 = p31_2_z(j_native_gate);



%%
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count = histcounts2(ke_31_1(:,1),theta_31_2_1(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);

% figure
% % dXedges_mod=(max(ke_23(:,1)) - min(ke_23(:,1)) )/200;
% % Xedges_mod = min(ke_23(:,1))+dXedges_mod/2:dXedges_mod:max(ke_23(:,1))-dXedges_mod/2;
% 
% dXedges_mod=0.1;
% Xedges_mod = 0+dXedges_mod/2:dXedges_mod:15-dXedges_mod/2;
% 
% plot(Xedges_mod,sum(count,2))
% hold on
plot([2.6 2.6 4.1 4.1 2.6],[20 45 45 20 20], 'k', 'Linewidth', 2)





%% reconstructed reflection
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count = histcounts2(ke_31_1(:,1),180 - theta_31_2_1(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);

% figure
% % dXedges_mod=(max(ke_23(:,1)) - min(ke_23(:,1)) )/200;
% % Xedges_mod = min(ke_23(:,1))+dXedges_mod/2:dXedges_mod:max(ke_23(:,1))-dXedges_mod/2;
% 
% dXedges_mod=0.1;
% Xedges_mod = 0+dXedges_mod/2:dXedges_mod:15-dXedges_mod/2;
% 
% plot(Xedges_mod,sum(count,2))
% hold on
plot([2.6 2.6 4.1 4.1 2.6],[135 160 160 135 135], 'b', 'Linewidth', 2)


%%  another cut

theta_min_1=92 ; theta_max_1=114;  
j_nat_ke31 = ke_31 >= 2.6 & ke_31 <= 4.1 & theta_31_2 > theta_min_1 & theta_31_2 < theta_max_1 ; % new value to make compare with dalitz plot



%%
j_native_gate= j_nat_ke31 ;
j_native_gate=j_native_gate;


ke_31_2 = ke_31(j_native_gate );
theta_31_2_2 = theta_31_2(j_native_gate);

px1_2 = px1(j_native_gate);   py1_2 = py1(j_native_gate);   pz1_2 = pz1(j_native_gate);
px2_2 = px2(j_native_gate);   py2_2 = py2(j_native_gate);   pz2_2 = pz2(j_native_gate);
px3_2 = px3(j_native_gate);   py3_2 = py3(j_native_gate);   pz3_2 = pz3(j_native_gate);

p31_x_2 = p31_x(j_native_gate); p31_y_2 = p31_y(j_native_gate); p31_z_2 = p31_z(j_native_gate);

p31_2_x_2 = p31_2_x(j_native_gate); p31_2_y_2 = p31_2_y(j_native_gate); p31_2_z_2 = p31_2_z(j_native_gate);



%%
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count = histcounts2(ke_31_2(:,1),theta_31_2_2(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);

% figure
% % dXedges_mod=(max(ke_23(:,1)) - min(ke_23(:,1)) )/200;
% % Xedges_mod = min(ke_23(:,1))+dXedges_mod/2:dXedges_mod:max(ke_23(:,1))-dXedges_mod/2;
% 
% dXedges_mod=0.1;
% Xedges_mod = 0+dXedges_mod/2:dXedges_mod:15-dXedges_mod/2;
% 
% plot(Xedges_mod,sum(count,2))
% hold on
plot([2.6 2.6 4.1 4.1 2.6],[92 114 114 92 92], '--k', 'Linewidth', 2)



%% reconstructed reflection
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count = histcounts2(ke_31_2(:,1),180 - theta_31_2_2(:,1),Xedges,Yedges);

count_mod = max(count,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);

% figure
% % dXedges_mod=(max(ke_23(:,1)) - min(ke_23(:,1)) )/200;
% % Xedges_mod = min(ke_23(:,1))+dXedges_mod/2:dXedges_mod:max(ke_23(:,1))-dXedges_mod/2;
% 
% dXedges_mod=0.1;
% Xedges_mod = 0+dXedges_mod/2:dXedges_mod:15-dXedges_mod/2;
% 
% plot(Xedges_mod,sum(count,2))
% hold on
plot([2.6 2.6 4.1 4.1 2.6],[66 88 88 66 66], '--b', 'Linewidth', 2)





%%  actual data

j_nat_ke31 = ke_31 >= 2.6 & ke_31 <= 4.1 & theta_31_2 >= 0  & theta_31_2 <= 66 | ...
             ke_31 >= 2.6 & ke_31 <= 4.1 & theta_31_2 >= 88 & theta_31_2 <= 135 | ...
             ke_31 >= 2.6 & ke_31 <= 4.1 & theta_31_2 >= 160 & theta_31_2 <= 180; % new value to make compare with dalitz plot




         %%
         clc
j_native_gate= j_nat_ke31 ;
j_native_gate=j_native_gate;


ke_31_3 = ke_31(j_native_gate );
theta_31_2_3 = theta_31_2(j_native_gate);

px1_3 = px1(j_native_gate);   py1_3 = py1(j_native_gate);   pz1_3 = pz1(j_native_gate);
px2_3 = px2(j_native_gate);   py2_3 = py2(j_native_gate);   pz2_3 = pz2(j_native_gate);
px3_3 = px3(j_native_gate);   py3_3 = py3(j_native_gate);   pz3_3 = pz3(j_native_gate);

p31_x_3 = p31_x(j_native_gate); p31_y_3 = p31_y(j_native_gate); p31_z_3 = p31_z(j_native_gate);

p31_2_x_3 = p31_2_x(j_native_gate); p31_2_y_3 = p31_2_y(j_native_gate); p31_2_z_3 = p31_2_z(j_native_gate);

%% one extra plot
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count_rec_f = histcounts2([ke_31_3(:,1)],[theta_31_2_3(:,1)],Xedges,Yedges);

count_mod = max(count_rec_f,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);

% figure
% % dXedges_mod=(max(ke_23(:,1)) - min(ke_23(:,1)) )/200;
% % Xedges_mod = min(ke_23(:,1))+dXedges_mod/2:dXedges_mod:max(ke_23(:,1))-dXedges_mod/2;
% 
% dXedges_mod=0.1;
% Xedges_mod = 0+dXedges_mod/2:dXedges_mod:15-dXedges_mod/2;
% 
% plot(Xedges_mod,sum(count,2))
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[92 114 114 92 92], 'k', 'Linewidth', 2)


hold on
plot([2.6 2.6 4.1 4.1 2.6],[20 45 45 20 20], 'k', 'Linewidth', 2)
hold on
plot([2.6 2.6 4.1 4.1 2.6],[135 160 160 135 135], 'b', 'Linewidth', 2)
hold on
plot([2.6 2.6 4.1 4.1 2.6],[92 114 114 92 92], '--k', 'Linewidth', 2)
hold on
plot([2.6 2.6 4.1 4.1 2.6],[66 88 88 66 66], '--b', 'Linewidth', 2)


%%
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count_rec_f = histcounts2([ke_31_1(:,1);ke_31_2(:,1);ke_31_3(:,1)],[(180-theta_31_2_1(:,1));(180-theta_31_2_2(:,1));theta_31_2_3(:,1)],Xedges,Yedges);

count_mod = max(count_rec_f,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([0 15 0 180])
xticks([0:5:15]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);

% figure
% % dXedges_mod=(max(ke_23(:,1)) - min(ke_23(:,1)) )/200;
% % Xedges_mod = min(ke_23(:,1))+dXedges_mod/2:dXedges_mod:max(ke_23(:,1))-dXedges_mod/2;
% 
% dXedges_mod=0.1;
% Xedges_mod = 0+dXedges_mod/2:dXedges_mod:15-dXedges_mod/2;
% 
% plot(Xedges_mod,sum(count,2))


% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[20 45 45 20 20], 'k', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[135 160 160 135 135], 'b', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[92 114 114 92 92], '--k', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[66 88 88 66 66], '--b', 'Linewidth', 2)
% 
% 
%%

close all
count_rec_sub_nv= count_allevt_nv - count_rec_f;

count_rec_f_mod = max(count_rec_sub_nv,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_rec_f_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
% axis([0 15 0 180])
% xticks([0:5:15]);
% yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');


caxis([ 1 200]);





%%  New momentum formation

px1_rec_f = [px1_1; px1_2; px1_3]; py1_rec_f = [py1_1; py1_2; py1_3]; pz1_rec_f = [pz1_1; pz1_2; pz1_3];
px2_rec_f = [px2_1; px2_2; px2_3]; py2_rec_f = [py2_1; py2_2; py2_3]; pz2_rec_f = [pz2_1; pz2_2; pz2_3];
px3_rec_f = [px3_1; px3_2; px3_3]; py3_rec_f = [py3_1; py3_2; py3_3]; pz3_rec_f = [pz3_1; pz3_2; pz3_3];


theta_rec_n = [(180-theta_31_2_1(:,1));(180-theta_31_2_2(:,1));theta_31_2_3(:,1)];
p31_2_x_n = [p31_2_x_1; p31_2_x_2; p31_2_x_3];
p31_2_y_n = [p31_2_y_1; p31_2_y_2; p31_2_y_3];
p31_2_z_n = [p31_2_z_1; p31_2_z_2; p31_2_z_3];

p31_x_n = [p31_x_1; p31_x_2; p31_x_3];
p31_y_n = [p31_y_1; p31_y_2; p31_y_3];
p31_z_n = [p31_z_1; p31_z_2; p31_z_3];


%%
z_fp = cross([p31_2_x_n  p31_2_y_n  p31_2_z_n], [p31_x_n p31_y_n p31_z_n]);
z_fp_m = sqrt(z_fp(:,1).^2 + z_fp(:,2).^2 + z_fp(:,3).^2);
z_fp_u = z_fp./z_fp_m;

y_fp =  [p31_2_x_n  p31_2_y_n p31_2_z_n];
y_fp_m = sqrt(y_fp(:,1).^2 + y_fp(:,2).^2 + y_fp(:,3).^2);
y_fp_u = y_fp./y_fp_m;


x_fp = cross(y_fp_u, z_fp_u);
x_fp_m = sqrt(x_fp(:,1).^2 + x_fp(:,2).^2 + x_fp(:,3).^2);
x_fp_u = x_fp./x_fp_m;

p31_x_fp = (p31_x_n.*x_fp_u(:,1) + p31_y_n.*x_fp_u(:,2) + p31_z_n.*x_fp_u(:,3));
p31_y_fp = (p31_x_n.*y_fp_u(:,1) + p31_y_n.*y_fp_u(:,2) + p31_z_n.*y_fp_u(:,3));
p31_z_fp = (p31_x_n.*z_fp_u(:,1) + p31_y_n.*z_fp_u(:,2) + p31_z_n.*z_fp_u(:,3));

%%
p31_m=sqrt(p31_x_fp.^2 + p31_y_fp.^2 + p31_z_fp.^2);
p31_x_fp_n = (p31_m).*sin(theta_rec_n.*pi/180);
p31_y_fp_n = -(p31_m).*cos(theta_rec_n.*pi/180);
p31_z_fp_n= p31_z_fp;

p31_lab = repmat(p31_x_fp_n,1,3).*x_fp_u + repmat(p31_y_fp_n,1,3).*y_fp_u + repmat(p31_z_fp_n,1,3).*z_fp_u;

p31_x_f = p31_lab(:,1);
p31_y_f = p31_lab(:,2);
p31_z_f = p31_lab(:,3);


%% total reconstructed sequential

pxs = px1_rec_f+px2_rec_f +px3_rec_f;
pys = py1_rec_f+py2_rec_f +py3_rec_f;
pzs = pz1_rec_f+pz2_rec_f +pz3_rec_f;
ke_31_f  = (p31_x_f.^2 + p31_y_f.^2 + p31_z_f.^2)./(2*mu_31)*hatoev;


%%
close all
count_rec_f = histcounts2(ke_31_f,theta_rec_n,Xedges,Yedges);  % this hitogram is exactly same as the histogram of previous  count_rec_f

count_rec_f_mod = max(count_rec_f,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_rec_f_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
% axis([0 15 0 180])
% xticks([0:5:15]);
% yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(1)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(1), Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);

%%
   
  

    PAx=  - (m(3)/(m(1)+m(3)))*(p31_2_x_n) +  (p31_x_f);
    PAy=  - (m(3)/(m(1)+m(3)))*(p31_2_y_n) + (p31_y_f);
    PAz=  - (m(3)/(m(1)+m(3)))*(p31_2_z_n) + (p31_z_f);


    PBx= -(m(1)/(m(1)+m(3)))*(p31_2_x_n) - (p31_x_f);
    PBy= -(m(1)/(m(1)+m(3)))*(p31_2_y_n) - (p31_y_f);
    PBz= -(m(1)/(m(1)+m(3)))*(p31_2_z_n) - (p31_z_f);
     

    PCx=(p31_2_x_n) ;
    PCy=(p31_2_y_n)  ;
    PCz=(p31_2_z_n)  ;
    


    %%  Newton diagram
   px3=PAx;py3=PAy;pz3=PAz; % because AB is the intermeduate i.e. 2nd and 3rd hit
   px1=PBx;py1=PBy;pz1=PBz;
   px2=PCx;py2=PCy;pz2=PCz; 

   % for combining all sequential events
   
   px3_xx=PAx;py3_xx=PAy;pz3_xx=PAz; % because AB is the intermeduate i.e. 2nd and 3rd hit
   px1_xx=PBx;py1_xx=PBy;pz1_xx=PBz;
   px2_xx=PCx;py2_xx=PCy;pz2_xx=PCz; 

%% newton_diagram with respect to 2 hit with reconstructed eevnts 1 up 2 down


p1_m=sqrt(px1.*px1 + py1.*py1 + pz1.*pz1);
p2_m=sqrt(px2.*px2 + py2.*py2 + pz2.*pz2);
p3_m=sqrt(px3.*px3 + py3.*py3 + pz3.*pz3);


p1_dot_p2 = (px1.*px2 + py1.*py2 + pz1.*pz2);
p1_dot_p3 = (px1.*px3 + py1.*py3 + pz1.*pz3); 
p2_dot_p3 = (px2.*px3 + py2.*py3 + pz2.*pz3); 

cos_theta_12=p1_dot_p2./(p1_m.*p2_m);
cos_theta_13=p1_dot_p3./(p1_m.*p3_m);
cos_theta_23=p2_dot_p3./(p2_m.*p3_m);

close all;
acute_theta_12=(pi - acos(cos_theta_12))*180/pi;
acute_theta_13=(pi - acos(cos_theta_13))*180/pi;
acute_theta_23=(pi - acos(cos_theta_23))*180/pi;


p3_m_np2=p3_m./p2_m;
p1_m_np2=p1_m./p2_m;

p3_m_np2_x = -(p3_m_np2).*cos(acute_theta_23.*pi/180); %momentum coordinate x
p3_m_np2_y = - (p3_m_np2).*sin(acute_theta_23.*pi/180); %momentum coordinate y

p1_m_np2_x = -(p1_m_np2).*cos(acute_theta_12.*pi/180);
p1_m_np2_y = (p1_m_np2).*sin(acute_theta_12.*pi/180);


p31_x=[p3_m_np2_x;p1_m_np2_x];
p31_y=[p3_m_np2_y;p1_m_np2_y];


%making 2d histogram
% Xedges = min(p31_x(:,1)):(max(p31_x(:,1)) - min(p31_x(:,1)) )/5000:max(p31_x(:,1));
% Yedges = min(p31_y(:,1)):(max(p31_y(:,1)) - min(p31_y(:,1)) )/5000:max(p31_y(:,1));

Xedges = -20: bin_nw: 20;
Yedges = -20: bin_nw: 20;

count_rec_nt = histcounts2(p31_x(:,1),p31_y(:,1),Xedges,Yedges);

count_mod = max(count_rec_nt,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;

imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([-2 2 -2 2])
% xticks([-2:1:2])
% yticks([-2:1:2])
axis([-2 2 -2 2])
% xticks([-2:1:2])
% yticks([-2:1:2])
hold on;
quiver(0,0,1.,0,'-r','LineWidth',2,'MaxHeadSize',3);
set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('rel. P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('rel. P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',50)
set(gca,'colorscale','log');
caxis([ 1 120]);


% annotation('textbox',[0.322 0.24 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.51 0.61 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.84 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');

    
%% subtratced newton's diagram
close all
count_sub_nt= count_allevt_nt - count_rec_nt ;
count_mod = max(count_sub_nt ,1); 
imagesc( Xedges,Yedges, ((count_mod')));
  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([-2 2 -2 2])
% xticks([-2:1:2])
% yticks([-2:1:2])
hold on;
quiver(0,0,1.,0,'-r','LineWidth',2,'MaxHeadSize',3);
set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('rel. P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('rel. P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',50)
set(gca,'colorscale','log');
caxis([ 1 120]);





































