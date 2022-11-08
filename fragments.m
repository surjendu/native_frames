function [frag_m,frag_m_z,frag_m_z_str]=fragments(frag)

mmns_to_au_conv=0.4571; %mm/ns to a.u. conversion
har_to_ev=  27.211396;
mpvsme=1822.888486424682;
H = 1.00782503223;
D = 2.01410177812;
He =  4.00260325413;
Li = 7.0160034366;
C = 12;
N = 14.00307400443;
O = 15.994914619;
F = 18.99840316273;
Ne =  19.9924401762;
S = 31.9720711744;
Cl = 34.968852682;
Ar = 39.9623831237;
Br79 = 78.9183376;
Br81 = 80.9162906;
I = 126.9044719;
C13 = 13.00335483521;

% m=(C+H+Br79*3)*mpvsme-q;


frag_m=[]; 
frag_m_z=[];
charge_z=[];
 frag_m_z_str=string.empty;
 for i=1:length(frag)
     frag_m_z_int_sum=0;
    frag_m_z_int_str=[];
     for j=1:(length(frag{i})-1)/2
          frag_m_z_int=frag{i}{2*j}*eval(frag{i}{2*j-1});
          frag_m_z_int_sum = frag_m_z_int_sum + frag_m_z_int;
          frag_m_z_int=[(frag{i}{2*j-1}),'_',num2str(frag{i}{2*j})];
          frag_m_z_int_str=[frag_m_z_int_str,frag_m_z_int];

     end
     charge_z=[charge_z,frag{i}{2*j+1}];
     frag_m=[frag_m, frag_m_z_int_sum];
     frag_m_z=[frag_m_z, frag_m_z_int_sum/frag{i}{2*j+1}];
     frag_m_z_str(i)=[frag_m_z_int_str '^' num2str(frag{i}{2*j+1}) '^' '+' ];
     
 end
 [charge_z;frag_m;frag_m_z]
 frag_m=frag_m*mpvsme-1;%q=1
 frag_m_z=frag_m_z*mpvsme-1;%q=1

end