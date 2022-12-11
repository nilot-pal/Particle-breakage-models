clear
clc
% % fid = fopen('Piecewise_linear_fitting_for_init_reb_PSD.txt','a+');
%initial_PSD = readmatrix("Rolls Royce project/Rolls_Royce_project/Particle_breakage/Stuttgart_2021/Datasets/Particle size distribution/CDF_bef_impact_Stainless_steel_target_speed_OP3.txt");
initial_PSD = readmatrix("Input_files/Particle size distribution/CDF_bef_impact_Stainless_Steel_target_speed_OP3.txt");
d_i = initial_PSD(:,1);
CDF_i = initial_PSD(:,2);
%PDF_i = diff([0 CDF_i]);
rebound_PSD = readmatrix("Input_files/Particle size distribution/CDF_aft_impact_Stainless_Steel_target_speed_OP3_angle_90_deg.txt");
%rebound_PSD = readmatrix("Rolls Royce project/Rolls_Royce_project/Particle_breakage/Stuttgart_2021/Datasets/Particle size distribution/CDF_aft_impact_Stainless_steel_target_speed_OP3_angle_90_deg.txt");
d_R = rebound_PSD(:,1);
CDF_r = rebound_PSD(:,2);
% % objective function = sum of squared errors
fun = @(x) sum((x(1)*x(2:17) + (1-x(1))*CDF_i(:) - CDF_r(:)).^2);
% Lower and upper bounds
lb = zeros(1,17)';
%lb(1) = 0.95;
ub = ones(1,17)';
% Linear inequality constraint
A = diag(ones(1,17)) + diag(-ones(1,16),1);
A(1,2) = 0;
A(end) = 0;
b = [1,zeros(1,16)]';
Aeq = zeros(17,17);
Aeq(end) = 1;
beq = [zeros(1,16),1]';
% Starting point
x0 = [0.5,linspace(0.05,0.95,16)]';
%x0 = [0.5,linspace(0,1,16)]';
% optimization
nonlcon = [];
options = optimoptions('fmincon','Algorithm','sqp');
[x,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
CDF_b = x(2:end);
CDF_r_dash = x(1)*x(2:17) + (1-x(1))*CDF_i(:);
plot(d_i,CDF_i,'o',d_R,CDF_b,d_R,CDF_r,'^',d_R,CDF_r_dash,'--','Linewidth',2);
legend('initial PSD','broken PSD','rebound PSD expt.','rebound PSD predicted');
xlabel('Particle diameter(\mum)');
ylabel('CDF');
% % disp(x(1));
n_i = 1e5;
%d_i = initial_PSD(:,1);
%d_R = rebound_PSD(:,1);
n = 200;
dq_i = linspace(d_i(1),d_i(end),n);
vq_i = interp1(d_i,CDF_i,dq_i);
% figure
% plot(d_i,CDF_i,'o',dq_i,vq_i,':.','LineWidth',2);
dq_r = linspace(d_R(1),d_R(end),n);
vq_r = interp1(d_R,CDF_r,dq_r);
% figure
% plot(d_R,CDF_r,'o',dq_r,vq_r,':.','LineWidth',2);
% dq_b = linspace(d_b(1),d_b(end),n);
% vq_b = interp1(d_b,CDF_b,dq_b);
% figure
% plot(d_b,CDF_b,'o',dq_b,vq_b,':.','LineWidth',2);
dd_i = (d_i(end)-d_i(1))/(n-1);
dd_r = (d_R(end)-d_R(1))/(n-1);
int_PSD_i = pi/6*(d_i(end)^3 - 3*sum(vq_i.*dq_i.^2*dd_i));
int_PSD_r = pi/6*(d_R(end)^3 - 3*sum(vq_r.*dq_r.^2*dd_r));
n_R = n_i*int_PSD_i/int_PSD_r;
B_n = x(1);
nR_ni = n_R/n_i;
S_i = 1 - (1-B_n)*(n_R/n_i);
disp(S_i);
disp(B_n)
%%
% fprintf( fid, '%f   %f   %f\n', B_n,nR_ni,S_i);
% fclose(fid);
%% pLOTTING S_i and n_R/n_i trends
fileID = fopen("Input_files/TiAl6V4_Piecewise_linear_fitting_for_init_reb_PSD - edited.txt","r");
formatSpec = ['%f' '%f' '%f' '%f' '%f'];
property = textscan(fileID,formatSpec);
v_i = property{1,1};
alpha = property{1,2};
v_ni = v_i.*sin(deg2rad(alpha));
vi_vni = v_i.*v_ni;
nR_ni = property{1,4};
S_i = property{1,5};
B_n = property{1,3};
% fun = @(params,v_ni) 1 - 1./(1+(v_ni/params(1)).^params(2));
% params_0 = [1,0.001];
% params = lsqcurvefit(fun,params_0,v_ni,S_i);
% S_i_analytic = 1 - 1./(1+(v_ni/params(1)).^params(2));
figure
plot(alpha(1:7),B_n(1:7),'-o',alpha(8:15),B_n(8:15),':s',alpha(16:23),B_n(16:23),'--^','LineWidth',2);
figure
plot(v_ni(1:7),B_n(1:7),'o',v_ni(8:15),B_n(8:15),'s',v_ni(16:23),B_n(16:23),'^','LineWidth',2);
figure
plot(v_ni(1:8),nR_ni(1:8),'bo',v_ni(9:16),nR_ni(9:16),'bs',v_ni(17:24),nR_ni(17:24),'r^','LineWidth',2);
figure
plot(vi_vni(1:8),nR_ni(1:8),'bo',vi_vni(9:16),nR_ni(9:16),'bs',vi_vni(17:24),nR_ni(17:24),'r^','LineWidth',2);
