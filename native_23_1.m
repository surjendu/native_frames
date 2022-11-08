close all 
clc
% clear

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

%% newton_diagram with respect to 1st hit without the boundary 
% clear some mass etc are reading from previous run
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

p2_m_np1=p2_m./p1_m;
p3_m_np1=p3_m./p1_m;

p2_m_np1_x = -(p2_m_np1).*cos(acute_theta_12.*pi/180); %momentum coordinate x
p2_m_np1_y =  (p2_m_np1).*sin(acute_theta_12.*pi/180); %momentum coordinate y

p3_m_np1_x = -(p3_m_np1).*cos(acute_theta_13.*pi/180);
p3_m_np1_y = -(p3_m_np1).*sin(acute_theta_13.*pi/180);

p23_x=[p2_m_np1_x;p3_m_np1_x];
p23_y=[p2_m_np1_y;p3_m_np1_y];

%making 2d histogram
% Xedges = min(p23_x(:,1)):(max(p23_x(:,1)) - min(p23_x(:,1)) )/5000:max(p23_x(:,1));
% Yedges = min(p23_y(:,1)):(max(p23_y(:,1)) - min(p23_y(:,1)) )/5000:max(p23_y(:,1));

Xedges = -20: bin_nw : 20;
Yedges = -20: bin_nw : 20;

count_allevt_nt = histcounts2(p23_x(:,1),p23_y(:,1),Xedges,Yedges);
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
set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('rel. P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('rel. P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',50)
set(gca,'colorscale','log');
% 
% annotation('textbox',[0.51 0.63 0 0.04],'EdgeColor','w','String', strcat(frag_m_z_str(1),'(',num2str(1),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.84 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(2),'(',num2str(2),')'),'FontSize',20,'FontWeight','normal','Color','k');
% annotation('textbox',[0.322 0.27 0 0.04],'EdgeColor','w','String',strcat(frag_m_z_str(3),'(',num2str(3),')'),'FontSize',20,'FontWeight','normal','Color','k');
caxis([ 1 120]);

scatter([-0.489675374 -0.510307908 ],[0.871919914 -0.871939817],300,'x','g','LineWidth',3) %concerted
p_wrt1_seq=dlmread('p_wrt1_seq.csv');
scatter(p_wrt1_seq(:,1),p_wrt1_seq(:,2),'.','b','LineWidth',3)

%%  native frame intermediate 2nd and 3rd
close all
clc


m_123=(m(1)+m(2)+m(3));
m_23=(m(2)+m(3));
mu_23=1/(1/m(2)+1/m(3));

p23_x=-mu_23*(px3./m(3)-px2./m(2));
p23_y=-mu_23*(py3./m(3)-py2./m(2));
p23_z=-mu_23*(pz3./m(3)-pz2./m(2));
p23_m=sqrt(p23_x.^2 + p23_y.^2 + p23_z.^2);

p23_1_x=(m_23/m_123)*px1 - (m(1)/m_123)*(px2+px3);
p23_1_y=(m_23/m_123)*py1 - (m(1)/m_123)*(py2+py3);
p23_1_z=(m_23/m_123)*pz1 - (m(1)/m_123)*(pz2+pz3);
p23_1_m=sqrt(p23_1_x.^2 + p23_1_y.^2 + p23_1_z.^2);

p23_1_dot_p23=(p23_1_x.*p23_x + p23_1_y.*p23_y +p23_1_z.*p23_z );

ke_23 = p23_m.^2/(2*mu_23)*hatoev;
theta_23_1=acos(p23_1_dot_p23./(p23_1_m.*p23_m))*180/pi;

% native frames plot without gating intermediate 2nd and 3rdS
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count_allevt_nv = histcounts2(ke_23(:,1),theta_23_1(:,1),Xedges,Yedges);

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

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)/^{\circ}}', 'FontWeight', 'normal','FontName', 'Arial');


hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');



% figure
% % dXedges_mod=(max(ke_23(:,1)) - min(ke_23(:,1)) )/200;
% % Xedges_mod = min(ke_23(:,1))+dXedges_mod/2:dXedges_mod:max(ke_23(:,1))-dXedges_mod/2;
% 
% dXedges_mod=0.1;
% Xedges_mod = 0+dXedges_mod/2:dXedges_mod:15-dXedges_mod/2;
% 
% plot(Xedges_mod,sum(count,2))
% j_nat_ke23 = ke_23 < 4.15; % old value

caxis([ 1 200]);

% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[20 45 45 20 20], 'k', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[135 160 160 135 135], 'b', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[92 114 114 92 92], '--k', 'Linewidth', 2)
% hold on
% plot([2.6 2.6 4.1 4.1 2.6],[66 88 88 66 66], '--b', 'Linewidth', 2)


hold on
plot([2.6 2.6 4.1 4.1 2.6],[20 45 45 20 20], 'b', 'Linewidth', 4)
hold on
plot([2.6 2.6 4.1 4.1 2.6],[135 160 160 135 135], 'r', 'Linewidth', 4)
hold on
plot([2.6 2.6 4.1 4.1 2.6],[92 114 114 92 92], 'b', 'Linewidth', 4)
hold on
plot([2.6 2.6 4.1 4.1 2.6],[66 88 88 66 66], 'r', 'Linewidth', 4)
hold on
plot([4.5 4.5 8 8 4.5],[80 100 100 80 80], 'm', 'Linewidth',4)


%% 1d projection
figure
% dXedges_mod=(max(ke_23(:,1)) - min(ke_23(:,1)) )/200;
% Xedges_mod = min(ke_23(:,1))+dXedges_mod/2:dXedges_mod:max(ke_23(:,1))-dXedges_mod/2;

dXedges_mod=bin_nv_ke;
Xedges_mod = 0+dXedges_mod/2:dXedges_mod:15-dXedges_mod/2;
axis xy;
plot(Xedges_mod,sum(count_allevt_nv,2), 'k', 'LineWidth', 2)
xticks([0:5:15]);
set(gca,'FontSize',30)
xlabel('KER_{CHBr-Br(2)}/eV','FontWeight', 'normal','FontName', 'Arial');
pbaspect([1 1 1]);


%%



theta_min_1=20 ; theta_max_1=45;  
j_nat_ke23 = ke_23 >= 2.6 & ke_23 <= 4.1 & theta_23_1 > theta_min_1 & theta_23_1 < theta_max_1 ; % new value to make compare with dalitz plot


%%
j_native_gate= j_nat_ke23 ;
j_native_gate=j_native_gate;


ke_23_1 = ke_23(j_native_gate );
theta_23_1_1 = theta_23_1(j_native_gate);

px1_1 = px1(j_native_gate);   py1_1 = py1(j_native_gate);   pz1_1 = pz1(j_native_gate);
px2_1 = px2(j_native_gate);   py2_1 = py2(j_native_gate);   pz2_1 = pz2(j_native_gate);
px3_1 = px3(j_native_gate);   py3_1 = py3(j_native_gate);   pz3_1 = pz3(j_native_gate);

p23_x_1 = p23_x(j_native_gate); p23_y_1 = p23_y(j_native_gate); p23_z_1 = p23_z(j_native_gate);
p23_1_x_1 = p23_1_x(j_native_gate); p23_1_y_1 = p23_1_y(j_native_gate); p23_1_z_1 = p23_1_z(j_native_gate);



%%
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count = histcounts2(ke_23_1(:,1),theta_23_1_1(:,1),Xedges,Yedges);

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

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)}', 'FontWeight', 'normal','FontName', 'Arial');


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
count = histcounts2(ke_23_1(:,1),180 - theta_23_1_1(:,1),Xedges,Yedges);

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

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)}', 'FontWeight', 'normal','FontName', 'Arial');


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
j_nat_ke23 = ke_23 >= 2.6 & ke_23 <= 4.1 & theta_23_1 > theta_min_1 & theta_23_1 < theta_max_1 ; % new value to make compare with dalitz plot

%%
j_native_gate= j_nat_ke23 ;
j_native_gate=j_native_gate;


ke_23_2 = ke_23(j_native_gate );
theta_23_1_2 = theta_23_1(j_native_gate);

px1_2 = px1(j_native_gate);   py1_2 = py1(j_native_gate);   pz1_2 = pz1(j_native_gate);
px2_2 = px2(j_native_gate);   py2_2 = py2(j_native_gate);   pz2_2 = pz2(j_native_gate);
px3_2 = px3(j_native_gate);   py3_2 = py3(j_native_gate);   pz3_2 = pz3(j_native_gate);

p23_x_2 = p23_x(j_native_gate); p23_y_2 = p23_y(j_native_gate); p23_z_2 = p23_z(j_native_gate);

p23_1_x_2 = p23_1_x(j_native_gate); p23_1_y_2 = p23_1_y(j_native_gate); p23_1_z_2 = p23_1_z(j_native_gate);



%%
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count = histcounts2(ke_23_2(:,1),theta_23_1_2(:,1),Xedges,Yedges);

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

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)}', 'FontWeight', 'normal','FontName', 'Arial');


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
count = histcounts2(ke_23_2(:,1),180 - theta_23_1_2(:,1),Xedges,Yedges);

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

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)}', 'FontWeight', 'normal','FontName', 'Arial');


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

j_nat_ke23 = ke_23 >= 2.6 & ke_23 <= 4.1 & theta_23_1 >= 0  & theta_23_1 <= 66 | ...
             ke_23 >= 2.6 & ke_23 <= 4.1 & theta_23_1 >= 88 & theta_23_1 <= 135 | ...
             ke_23 >= 2.6 & ke_23 <= 4.1 & theta_23_1 >= 160 & theta_23_1 <= 180; % new value to make compare with dalitz plot

%%
j_native_gate= j_nat_ke23 ;
% j_native_gate=j_native_gate;


ke_23_3 = ke_23(j_native_gate );
theta_23_1_3 = theta_23_1(j_native_gate);

px1_3 = px1(j_native_gate);   py1_3 = py1(j_native_gate);   pz1_3 = pz1(j_native_gate);
px2_3 = px2(j_native_gate);   py2_3 = py2(j_native_gate);   pz2_3 = pz2(j_native_gate);
px3_3 = px3(j_native_gate);   py3_3 = py3(j_native_gate);   pz3_3 = pz3(j_native_gate);

p23_x_3 = p23_x(j_native_gate); p23_y_3 = p23_y(j_native_gate); p23_z_3 = p23_z(j_native_gate);

p23_1_x_3 = p23_1_x(j_native_gate); p23_1_y_3 = p23_1_y(j_native_gate); p23_1_z_3 = p23_1_z(j_native_gate);

%%  one extra plot
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count_rec_f = histcounts2([ke_23_3(:,1)],[theta_23_1_3(:,1)],Xedges,Yedges);

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

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)}', 'FontWeight', 'normal','FontName', 'Arial');


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

% 

%%
close all
% Xedges = min(ke_23(:,1)):(max(ke_23(:,1)) - min(ke_23(:,1)) )/200:max(ke_23(:,1));
% Yedges = min(theta_23_1(:,1)):(max(theta_23_1(:,1)) - min(theta_23_1(:,1)) )/250:max(theta_23_1(:,1));
Xedges = 0:bin_nv_ke:15;
Yedges = 0:bin_nv_an:180;
count_rec_f = histcounts2([ke_23_1(:,1);ke_23_2(:,1);ke_23_3(:,1)],[(180-theta_23_1_1(:,1));(180-theta_23_1_2(:,1));theta_23_1_3(:,1)],Xedges,Yedges);

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

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)/^{\circ}}', 'FontWeight', 'normal','FontName', 'Arial');

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

%save('fck', 'count_rec_f', '-v7.3');


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

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)}', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');


caxis([ 1 200]);


%%  New momentum formation
px1_rec_f = [px1_1; px1_2; px1_3]; py1_rec_f = [py1_1; py1_2; py1_3]; pz1_rec_f = [pz1_1; pz1_2; pz1_3];
px2_rec_f = [px2_1; px2_2; px2_3]; py2_rec_f = [py2_1; py2_2; py2_3]; pz2_rec_f = [pz2_1; pz2_2; pz2_3];
px3_rec_f = [px3_1; px3_2; px3_3]; py3_rec_f = [py3_1; py3_2; py3_3]; pz3_rec_f = [pz3_1; pz3_2; pz3_3];


theta_rec_n = [(180-theta_23_1_1(:,1));(180-theta_23_1_2(:,1));theta_23_1_3(:,1)];
p23_1_x_n = [p23_1_x_1; p23_1_x_2; p23_1_x_3];
p23_1_y_n = [p23_1_y_1; p23_1_y_2; p23_1_y_3];
p23_1_z_n = [p23_1_z_1; p23_1_z_2; p23_1_z_3];

p23_x_n = [p23_x_1; p23_x_2; p23_x_3];
p23_y_n = [p23_y_1; p23_y_2; p23_y_3];
p23_z_n = [p23_z_1; p23_z_2; p23_z_3];

%%
z_fp = cross([p23_1_x_n  p23_1_y_n  p23_1_z_n], [p23_x_n p23_y_n p23_z_n]);
z_fp_m = sqrt(z_fp(:,1).^2 + z_fp(:,2).^2 + z_fp(:,3).^2);
z_fp_u = z_fp./z_fp_m;

y_fp =  [p23_1_x_n  p23_1_y_n p23_1_z_n];
y_fp_m = sqrt(y_fp(:,1).^2 + y_fp(:,2).^2 + y_fp(:,3).^2);
y_fp_u = y_fp./y_fp_m;


x_fp = cross(y_fp_u, z_fp_u);
x_fp_m = sqrt(x_fp(:,1).^2 + x_fp(:,2).^2 + x_fp(:,3).^2);
x_fp_u = x_fp./x_fp_m;

p23_x_fp = (p23_x_n.*x_fp_u(:,1) + p23_y_n.*x_fp_u(:,2) + p23_z_n.*x_fp_u(:,3));
p23_y_fp = (p23_x_n.*y_fp_u(:,1) + p23_y_n.*y_fp_u(:,2) + p23_z_n.*y_fp_u(:,3));
p23_z_fp = (p23_x_n.*z_fp_u(:,1) + p23_y_n.*z_fp_u(:,2) + p23_z_n.*z_fp_u(:,3));
%% line of nodes
x_lab = repmat([1 0 0],length(z_fp_u(:,1)),1);
y_lab = repmat([0 1 0],length(z_fp_u(:,1)),1);
z_lab = repmat([0 0 1],length(z_fp_u(:,1)),1);

n_lon = cross( [z_lab(:,1) z_lab(:,2) z_lab(:,3)], [z_fp_u(:,1) z_fp_u(:,2) z_fp_u(:,3)]);  %cross of  lab z ( 0 0 1) and  z_fp_u axes
n_lon_m = sqrt(n_lon(:,1).^2 + n_lon(:,2).^2 + n_lon(:,3).^2);
n_lon_u = n_lon./n_lon_m;


n_dot_x_lab = (n_lon_u(:,1).*x_lab(:,1) + n_lon_u(:,2).*x_lab(:,2) + n_lon_u(:,3).*x_lab(:,3));
n_dot_y_lab = (n_lon_u(:,1).*y_lab(:,1) + n_lon_u(:,2).*y_lab(:,2) + n_lon_u(:,3).*y_lab(:,3));

a_ab = atan2d(-n_dot_x_lab, n_dot_y_lab); % alpha ab  atan2d(y,x) i.e. four quadrature value of tan-1(y/x)


n_dot_x_fp_u = (n_lon_u(:,1).*x_fp_u(:,1) + n_lon_u(:,2).*x_fp_u(:,2) + n_lon_u(:,3).*x_fp_u(:,3));
n_dot_y_fp_u = (n_lon_u(:,1).*y_fp_u(:,1) + n_lon_u(:,2).*y_fp_u(:,2) + n_lon_u(:,3).*y_fp_u(:,3));
g_ab = atan2d(n_dot_x_fp_u, n_dot_y_fp_u); % alpha ab

cos_b_ab = (z_lab(:,1).*z_fp_u(:,1) + z_lab(:,2).*z_fp_u(:,2) + z_lab(:,3).*z_fp_u(:,3));
% aaa = (n_lon_u(:,1).*z_fp_u(:,1) + n_lon_u(:,2).*z_fp_u(:,2) + n_lon_u(:,3).*z_fp_u(:,3)); % to test cross product

%%  beta plot
nedges_x =-1:0.05:1;
nedges_y =0:bin_nv_an*3:180;
count_rec_f = histcounts2(cos_b_ab,theta_rec_n,nedges_x,nedges_y);  % this hitogram is exactly same as the histogram of previous  count_rec_f

count_rec_f_mod = max(count_rec_f,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;
close all
imagesc( nedges_x,nedges_y, ((count_rec_f_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis([-1 1 0 180])
xticks([-1 : 0.5: 1]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel('cos \beta_{CHBr-Br(2)}', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\theta_{CHBr-Br(2), Br(1)}/^{\circ}', 'FontWeight', 'normal','FontName', 'Arial');
% 
hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
  caxis([ 1 300]);

%% alpha plot
nedges_x = -180:bin_nv_an*5:180;
nedges_y =0:bin_nv_an*5:180;
count_rec_f = histcounts2(a_ab,theta_rec_n,nedges_x,nedges_y);  % this hitogram is exactly same as the histogram of previous  count_rec_f

count_rec_f_mod = max(count_rec_f,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;
close all
imagesc( nedges_x,nedges_y, ((count_rec_f_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis xy;
axis([-180 180 0 180])
xticks([-180:90: 180]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel(' \alpha_{CHBr-Br(2)}/^{\circ}', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\theta_{CHBr-Br(2), Br(1)}/^{\circ}', 'FontWeight', 'normal','FontName', 'Arial');
% 
hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);
%% gamma plot
nedges_x = -180:bin_nv_an*5:180;
nedges_y =0:bin_nv_an*5:180;
count_rec_f = histcounts2(g_ab,theta_rec_n,nedges_x,nedges_y);  % this hitogram is exactly same as the histogram of previous  count_rec_f

count_rec_f_mod = max(count_rec_f,1); 
% count_mod = count; 
% myColorMap = jet;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;
close all
imagesc( nedges_x,nedges_y, ((count_rec_f_mod')));
  
%colorbar('FontSize', 20,'Location','west');
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
% cb = colorbar; 
% set(cb,'position',[0.68 .5 .01 .4]) %[xposition yposition width height].

axis xy;
axis xy;
axis([-180 180 0 180])
xticks([-180:90: 180]);
yticks([0:60:180]);

set(gca,'FontSize',30)
% axis equal;
pbaspect([1 1 1]);

xlabel(' \gamma_{CHBr-Br(2)}/^{\circ}', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)}/^{\circ}', 'FontWeight', 'normal','FontName', 'Arial');
% 
hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);
%%
p23_m=sqrt(p23_x_fp.^2 + p23_y_fp.^2 + p23_z_fp.^2);
p23_x_fp_n = (p23_m).*sin(theta_rec_n.*pi/180);
p23_y_fp_n = -(p23_m).*cos(theta_rec_n.*pi/180);
p23_z_fp_n= p23_z_fp;

p23_lab = repmat(p23_x_fp_n,1,3).*x_fp_u + repmat(p23_y_fp_n,1,3).*y_fp_u + repmat(p23_z_fp_n,1,3).*z_fp_u;

p23_x_f = p23_lab(:,1);
p23_y_f = p23_lab(:,2);
p23_z_f = p23_lab(:,3);

%% total reconstructed sequential

pxs = px1_rec_f+px2_rec_f +px3_rec_f;
pys = py1_rec_f+py2_rec_f +py3_rec_f;
pzs = pz1_rec_f+pz2_rec_f +pz3_rec_f;
ke_23_f  = (p23_x_f.^2 + p23_y_f.^2 + p23_z_f.^2)./(2*mu_23)*hatoev;


%% 
close all
count_rec_f = histcounts2(ke_23_f,theta_rec_n,Xedges,Yedges);  % this hitogram is exactly same as the histogram of previous  count_rec_f

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

xlabel('KER_{CHBr-Br(2)}/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('\Theta_{CHBr-Br(2), Br(1)}', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',40)
set(gca,'colorscale','log');
caxis([ 1 200]);



%%
%     PAx= ( m(2)/m_123 )*pxs-(m(2)/(m(3)+m(2)))*(p23_1_x_n)  - (p23_x_f);
%     PAy= ( m(2)/m_123  )*pys-(m(2)/(m(3)+m(2)))*(p23_1_y_n) - (p23_y_f);
%     PAz= ( m(2)/m_123  )*pzs-(m(2)/(m(3)+m(2)))*(p23_1_z_n) - (p23_z_f);
% 
% 
%     PBx=( m(3)/m_123  )*pxs-(m(3)/(m(3)+m(2)))*(p23_1_x_n) + (p23_x_f);
%     PBy=( m(3)/m_123  )*pys-(m(3)/(m(3)+m(2)))*(p23_1_y_n) + (p23_y_f);
%     PBz=( m(3)/m_123  )*pzs-(m(3)/(m(3)+m(2)))*(p23_1_z_n) + (p23_z_f);
% 
%     PCx=(p23_1_x_n) ;
%     PCy=(p23_1_y_n)  ;
%     PCz=(p23_1_z_n)  ;

    
    PAx= -(m(2)/(m(3)+m(2)))*(p23_1_x_n)  - (p23_x_f);
    PAy= -(m(2)/(m(3)+m(2)))*(p23_1_y_n) - (p23_y_f);
    PAz= -(m(2)/(m(3)+m(2)))*(p23_1_z_n) - (p23_z_f);


    PBx= -(m(3)/(m(3)+m(2)))*(p23_1_x_n) + (p23_x_f);
    PBy= -(m(3)/(m(3)+m(2)))*(p23_1_y_n) + (p23_y_f);
    PBz= -(m(3)/(m(3)+m(2)))*(p23_1_z_n) + (p23_z_f);

    PCx=(p23_1_x_n) ;
    PCy=(p23_1_y_n)  ;
    PCz=(p23_1_z_n)  ;
    
    %%  Newton diagram
%    px2=PAx;py2=PAy;pz2=PAz; % because AB is the intermeduate i.e. 2nd and 3rd hit
%    px3=PBx;py3=PBy;pz3=PBz;
%    px1=PCx;py1=PCy;pz1=PCz; 
% % 
   px2=[PAx;px2_xx]; py2=[PAy;py2_xx]; pz2=[PAz;pz2_xx]; % because AB is the intermeduate i.e. 2nd and 3rd hit px_1 py_1 pz_1 from others
   px3=[PBx;px3_xx]; py3=[PBy;py3_xx]; pz3=[PBz;pz3_xx];
   px1=[PCx;px1_xx]; py1=[PCy;py1_xx]; pz1=[PCz;pz1_xx]; 
   
%    px2=[px2_xx]; py2=[py2_xx]; pz2=[pz2_xx]; % because AB is the intermeduate for newton plot from other set of analysis
%    px3=[px3_xx]; py3=[py3_xx]; pz3=[pz3_xx];
%    px1=[px1_xx]; py1=[py1_xx]; pz1=[pz1_xx]; 
%    
%% newton_diagram with respect to 1st hit without the boundary 


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


p2_m_np1=p2_m./p1_m;
p3_m_np1=p3_m./p1_m;

p2_m_np1_x = -(p2_m_np1).*cos(acute_theta_12.*pi/180); %momentum coordinate x
p2_m_np1_y =  (p2_m_np1).*sin(acute_theta_12.*pi/180); %momentum coordinate y

p3_m_np1_x = -(p3_m_np1).*cos(acute_theta_13.*pi/180);
p3_m_np1_y = -(p3_m_np1).*sin(acute_theta_13.*pi/180);

p23_x=[p2_m_np1_x;p3_m_np1_x];
p23_y=[p2_m_np1_y;p3_m_np1_y];

%making 2d histogram
% Xedges = min(p23_x(:,1)):(max(p23_x(:,1)) - min(p23_x(:,1)) )/5000:max(p23_x(:,1));
% Yedges = min(p23_y(:,1)):(max(p23_y(:,1)) - min(p23_y(:,1)) )/5000:max(p23_y(:,1));

Xedges = -20: bin_nw: 20;
Yedges = -20: bin_nw: 20;

count_rec_nt = histcounts2(p23_x(:,1),p23_y(:,1),Xedges,Yedges);


count_mod = max(count_rec_nt ,1); 
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
set(gca,'FontSize',30)
axis equal;
pbaspect([1 1 1]);

xlabel('rel. P_x', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('rel. P_y', 'FontWeight', 'normal','FontName', 'Arial');

hold on;
set(gca,'FontSize',50)
set(gca,'colorscale','log');

    
% subtratced newton's diagram
close all


  count_sub_nt=  count_rec_nt ;     % subtraction ------------------------------------
% count_sub_nt=count_allevt_nt -   count_rec_nt ;% 



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


% x_gate_1 = ([-0.76 -0.62 -.42  -0.76 ])';
% y_gate_1 = ([1.02 1.31 1.03  1.02])';
% 
% % [in_1,on_1] = inpolygon(p_f(:,13),p_f(:,15),x_gate,z_gate); %inside and on the edge
% 
% plot(x_gate_1,y_gate_1,'g', 'Linewidth', 4) % polygon
% 
% x_gate_2 = ([-.6 -0.38 -.21 -.6 ])';
% y_gate_2 = ([-1.03 -1.29 -1  -1.03])';
% 
% 
% plot(x_gate_2,y_gate_2,'g', 'Linewidth', 4) % polygon
% 
% 
% scatter([-0.489675374 -0.510307908 ],[0.871919914 -0.871939817],300,'x','g','LineWidth',3) %concerted
% p_wrt1_seq=dlmread('p_wrt1_seq.csv');
% scatter(p_wrt1_seq(:,1),p_wrt1_seq(:,2),'.','b','LineWidth',3)


%%  angular distribution
% ke and angular distribution
clc
    KE1 = (px1.*px1 + py1.*py1 + pz1.*pz1)/(2*frag_m(1))*hatoev;
    KE2 = (px2.*px2 + py2.*py2 + pz2.*pz2)/(2*frag_m(2))*hatoev;
    KE3 = (px3.*px3 + py3.*py3 + pz3.*pz3)/(2*frag_m(3))*hatoev;
    KER=KE1+KE2+KE3;

   %binsize_ke=0.05; %eV  
   binsize_ke=0.2; %eV   
   edges_ke=[0:binsize_ke:20];
   i_ke=1:length(edges_ke)-1;
   bincent_ke=[];
   bincent_ke(i_ke)=(edges_ke(i_ke)+edges_ke(i_ke+1))/2;
  
   % ke_raw=[k1 k2 k3 k];
   % dlmwrite('ke_raw.csv',ke_raw);
    
   [KE1_counts,KE1_edges]=histcounts(KE1,edges_ke);
   [KE2_counts,KE2_edges]=histcounts(KE2,edges_ke);
   [KE3_counts,KE3_edges]=histcounts(KE3,edges_ke);
   [KER_counts,KER_edges]=histcounts(KER,edges_ke);
  
   
   ke_hist=[bincent_ke' KE1_counts' KE2_counts' KE3_counts' KER_counts'];
%    dlmwrite('ke_hist.csv',ke_hist);
  
   figure 
   close all;
   plot(bincent_ke,KE1_counts,'b', bincent_ke,KE2_counts,'r', bincent_ke,KE3_counts,'g', bincent_ke,KER_counts,'k','LineWidth',2)   
   set(gca, 'XScale', 'linear','YScale', 'linear','FontSize',30 );
   grid on
   xlim([0 edges_ke(end)]);
%    pbaspect([1 1 1]);
   xlabel('KER /eV', 'FontWeight', 'normal','FontName', 'Arial');
   ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');
   legend({frag_m_z_str(1),frag_m_z_str(2),frag_m_z_str(3),'KER'}, 'FontSize', 20, 'FontWeight', 'normal')


    %sb
%     ke=[KE1_bin'  KE1_counts KE2_bin'  KE2_counts KE3_bin'  KE3_counts  KER_bin'  KER_counts];
%     dlmwrite('ke.csv',ke);
    %sb
    
%%     angular distribution for first hit
    
costheta1z= pz1./p1_m;
costheta2z= pz2./p2_m;
costheta3z= pz3./p3_m;

binsize_ke=0.2; %eV   
edges_ke=[0:binsize_ke:28];

i_ke=1:length(edges_ke)-1;
bincent_ke=[];
bincent_ke(i_ke)=(edges_ke(i_ke)+edges_ke(i_ke+1))/2;

binsize_costheta=0.02; %eV   
edges_costheta=[-1:binsize_costheta:1];

i_costheta=1:length(edges_costheta)-1;
bincent_costheta=[];
bincent_costheta(i_costheta)=(edges_costheta(i_costheta)+edges_costheta(i_costheta+1))/2;

close all
count_1z = histcounts2(KE1(:,1),costheta1z(:,1),edges_ke,edges_costheta);
count_mod_1z = max(count_1z,1); 

figure;
myColorMap=flipud(hot);
myColorMap(1,:) = 1;
imagesc( edges_ke,edges_costheta, (count_mod_1z'));  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 15 -1 1])
set(gca,'FontSize',30)
set(gca,'colorscale','log');
pbaspect([1 1 1]);
xlabel('KE1 (Br(1))/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('cos\theta\_zaxis', 'FontWeight', 'normal','FontName', 'Arial');
pbaspect([1 1 1]);

%%
close all
figure
aaaa=[bincent_costheta',sum(count_mod_1z)'];
plot(bincent_costheta,sum(count_mod_1z),'k','LineWidth',3);
xlabel('cos\theta','FontWeight', 'normal','FontName', 'Arial');
% ylabel('Br(1) counts', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('counts', 'FontWeight', 'normal','FontName', 'Arial');

set(gca, 'XScale', 'linear')
 %set(gca, 'YScale', 'log')
set(gca,'FontSize',40)
xlim([0 1]);
ylim([0 5000]);
%ylim([0 3000]);
ylim([0 700]);
xticks([0:0.2:1])
yticks([0 : 200: 700])

pbaspect([1 1 1]);

%% %%     angular distribution for second hit

close all
count_2z = histcounts2(KE2(:,1),costheta2z(:,1),edges_ke,edges_costheta);
count_mod_2z = max(count_2z,1); 

figure;
imagesc( edges_ke,edges_costheta, (count_mod_2z'));  
colorbar('FontSize', 20);
colormap(myColorMap);
colorbar 
axis xy;
axis([0 15 -1 1])
axis([0 15 -1 1])
set(gca,'FontSize',30)
set(gca,'colorscale','log');
pbaspect([1 1 1]);
xlabel('KE2 (Br(2))/eV', 'FontWeight', 'normal','FontName', 'Arial');
ylabel('cos\theta\_zaxis', 'FontWeight', 'normal','FontName', 'Arial');
pbaspect([1 1 1]);

%%
close all
figure

plot(bincent_costheta,sum(count_mod_2z),'r','LineWidth',2);
xlabel('cos\theta\_zaxis','FontWeight', 'normal','FontName', 'Arial');
ylabel('Br(2) counts', 'FontWeight', 'normal','FontName', 'Arial');
set(gca, 'XScale', 'linear')
% set(gca, 'YScale', 'log')
set(gca,'FontSize',30)
xlim([0 1]);
% ylim([0 max(sum(count_mod_2z))*1.1]);
ylim([0 3000]);
pbaspect([1 1 1]);

