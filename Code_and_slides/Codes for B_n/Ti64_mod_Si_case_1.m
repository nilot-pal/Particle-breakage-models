clear
clc

vmax = 220;
%y = load('TiAl6V4_Piecewise_linear_fitting_for_init_reb_PSD - edited.txt');
y = load('Input_files/SS_Piecewise_linear_fitting_for_init_reb_PSD-edited.txt');
op1 = y(1:7,:);
op2 = y(8:15,:);
op3 = y(16:23,:);
n_r_ni = zeros(23,1);
n_r_ni(1:7) = (1-op1(:,4))./(1-op1(:,3));
n_r_ni(8:15) = (1-op2(:,4))./(1-op2(:,3));
n_r_ni(16:23) = (1-op3(:,4))./(1-op3(:,3));

% gamma1_0 = 1;
% gamma2_0 = 1;
% gamma3_0 = 1;

% gamma1_0 = 1.0;
% gamma2_0 = 0.4;
% gamma3_0 = 0.1;

gamma3_0 = 0.1;
gamma2_0 = gamma3_0*1.7;
gamma1_0 = gamma3_0*4.35;

% gamma1 = gamma1_0.*sin(op1(:,2)/180*pi);
% gamma2 = gamma2_0.*sin(op2(:,2)/180*pi);
% gamma3 = gamma3_0.*sin(op3(:,2)/180*pi);

% gamma1 = gamma1_0;
% gamma2 = gamma2_0;
% gamma3 = gamma3_0;

gamma1 = gamma1_0./sin(op1(:,2)/180*pi);
gamma2 = gamma2_0./sin(op2(:,2)/180*pi);
gamma3 = gamma3_0./sin(op3(:,2)/180*pi);

nR_ni_op1 =(1-op1(:,4))./(1-op1(:,3)) .* gamma1 ;
nR_ni_op2 =(1-op2(:,4))./(1-op2(:,3)) .* gamma2 ;
nR_ni_op3 =(1-op3(:,4))./(1-op3(:,3)) .* gamma3 ;

op1(:,4) = 1-(1-op1(:,3)) .* nR_ni_op1;  %S_i
op2(:,4) = 1-(1-op2(:,3)) .* nR_ni_op2;
op3(:,4) = 1-(1-op3(:,3)) .* nR_ni_op3;

figure(1)
% plot(op1(:,2),op1(:,4),'linewidth',2,'bo-','markerfacecolor','b'); hold on
% plot(op2(:,2),op2(:,4),'linewidth',2,'rs--','markerfacecolor','w');
% plot(op3(:,2),op3(:,4),'linewidth',2,'k^:','markerfacecolor','g');
plot(op1(:,2),op1(:,4),'bo-.','linewidth',2); hold on
plot(op2(:,2),op2(:,4),'rs--','linewidth',2);
plot(op3(:,2),op3(:,4),'k^:','linewidth',2);
axis([0 100 0 1]);
xlabel('impact angle (degrees)');
ylabel('S_i (-)');
legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s, g_1 = %3.2f',op1(1,1),gamma1_0) ],...
       ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s, g_2 = %3.2f',op2(1,1),gamma2_0)],...
       ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s, g_3 = %3.2f',op3(1,1),gamma3_0)],'location','southeast');

% figure(2)
% subplot(121)
% plot(op1(:,1).^2.*sin(op1(:,2)/180*pi),nR_ni_op1,'bo','markerfacecolor','b'); hold on
% plot(op2(:,1).^2.*sin(op2(:,2)/180*pi),nR_ni_op2,'bs','markerfacecolor','w');
% plot(op3(:,1).^2.*sin(op3(:,2)/180*pi),nR_ni_op3,'r^','markerfacecolor','r');
% legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
%        ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
%        ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],'location','northwest');
% xlabel('v_{N,i}*v_i');
% ylabel('n_R/n_i');
% xlim([0 vmax^2]);
% subplot(122)
% plot(op1(:,1).*sin(op1(:,2)/180*pi),nR_ni_op1,'bo','markerfacecolor','b'); hold on
% plot(op2(:,1).*sin(op2(:,2)/180*pi),nR_ni_op2,'bs','markerfacecolor','w');
% plot(op3(:,1).*sin(op3(:,2)/180*pi),nR_ni_op3,'r^','markerfacecolor','r');
% legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
%        ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
%        ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],'location','northwest');
% xlabel('v_{N,i}');
% ylabel('n_R/n_i');
% xlim([0 vmax]);

% figure
% 
% dpm_op1   = (op1(:,3))./(nR_ni_op1.^-1-(1-op1(:,3))); % daughter per mother
% dpm_op2   = (op2(:,3))./(nR_ni_op2.^-1-(1-op2(:,3)));
% dpm_op3   = (op3(:,3))./(nR_ni_op3.^-1-(1-op3(:,3)));
% 
% subplot(121)
% plot(op1(:,1).^2.*sin(op1(:,2)/180*pi),dpm_op1,'bo','markerfacecolor','b'); hold on
% plot(op2(:,1).^2.*sin(op2(:,2)/180*pi),dpm_op2,'bs','markerfacecolor','w');
% plot(op3(:,1).^2.*sin(op3(:,2)/180*pi),dpm_op3,'r^','markerfacecolor','r');
% legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
%        ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
%        ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],'location','northwest');
% xlabel('v_{N,i}*v_i');
% ylabel('visible daughter per mother');
% xlim([0 vmax^2]);
% 
% subplot(122)
% plot(op1(:,1).*sin(op1(:,2)/180*pi),dpm_op1,'bo','markerfacecolor','b'); hold on
% plot(op2(:,1).*sin(op2(:,2)/180*pi),dpm_op2,'bs','markerfacecolor','w');
% plot(op3(:,1).*sin(op3(:,2)/180*pi),dpm_op3,'r^','markerfacecolor','r');
% legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
%        ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
%        ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],'location','northwest');
% xlabel('v_{N,i}');
% ylabel('visible daughter per mother');
% xlim([0 vmax]);



figure(2)

plot(op1(:,1).*sin(op1(:,2)/180*pi),op1(:,4),'bo','markerfacecolor','b'); hold on
plot(op2(:,1).*sin(op2(:,2)/180*pi),op2(:,4),'bs','markerfacecolor','w');
plot(op3(:,1).*sin(op3(:,2)/180*pi),op3(:,4),'b^','markerfacecolor','g');
vv = linspace(0,300,3000);
Dv_guess = 4;
v50_guess= 30;
a = y(:,1).*sin(y(:,2)/180*pi); % v_ni
b = [op1(:,4)          % modified S_i using gamma
     op2(:,4)
     op3(:,4)];

S_i_data = [a b];

%x0 = [v50_guess Dv_guess];
%[x, obj, info, iter, nf, lambda] = sqp(x0,@stuttgart_data_fitting); % minimize the objective function - sum of squared errors with initial value x0
% % Minimize fun using fmincon
A = [];
B = [];
S_i_fit = @(x) 1-1./(1+(S_i_data(:,1)/x(1)).^x(2));
diff = @(x) (S_i_fit(x)-S_i_data(:,2)).^2;
obj  = @(x) sum(diff(x));
x0 = [v50_guess Dv_guess];
x = fmincon(obj,x0,A,B);
S_i_fit = 1-1./(1+(vv/x(1)).^x(2));
plot(vv,S_i_fit,'r','linewidth',2);
legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
       ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
       ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],...
       'S_i = 1- 1/(1+(v_{ni}/v50)^{Dv})','location','southeast');
text(15,0.7,sprintf('v50 = %5.1f m/s, Dv = %3.1f',x(1),x(2)));
xlim([0 250]);
grid on
xlabel('normal incident speed v_{ni} (m/s)');
ylabel('S_i (-)');


% figure
% %pkg load statistics
% subplot(121)
% plot(op1(:,1).^2.*sin(op1(:,2)/180*pi),op1(:,4),'bo','markerfacecolor','b'); hold on
% plot(op2(:,1).^2.*sin(op2(:,2)/180*pi),op2(:,4),'bs','markerfacecolor','w');
% plot(op3(:,1).^2.*sin(op3(:,2)/180*pi),op3(:,4),'b^','markerfacecolor','g');
% vv = linspace(0,300^2,3000);
% mu = log(100^2*0.912);
% sigma = 1;
% median_val = exp(mu);
% variance = (exp(sigma^2)-1)*exp(2*mu+sigma^2);
% std_val  = sqrt(variance);
% cdf_lognorm = logncdf (vv,mu,sigma);
% plot(vv,cdf_lognorm,'r-');
% legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
%        ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
%        ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],...
%        'S_i = log_normal','location','southeast');
% text(4e4,0.4,sprintf('median = %5.1e (m/s)^2,  std. = %5.1e',median_val,std_val));
% 
% grid on
% xlabel('v_i \times n_{ni} (m/s)^2');
% ylabel('S_i (-)');
% subplot(122)
% title('lognorm distribution fitting, S_i (v_{Ni})');
% pdf_lognorm = lognpdf(vv,mu,sigma);
% semilogx(vv,pdf_lognorm);
% xlim([0 2e5]);
% xlabel('v_i \times n_{ni} (m/s)^2');
% ylabel('PDF of S_i (m/2)^{-2}');
% 
% figure
% %pkg load statistics
% subplot(121)
% plot(op1(:,1).*sin(op1(:,2)/180*pi),op1(:,4),'bo','markerfacecolor','b'); hold on
% plot(op2(:,1).*sin(op2(:,2)/180*pi),op2(:,4),'bs','markerfacecolor','w');
% plot(op3(:,1).*sin(op3(:,2)/180*pi),op3(:,4),'b^','markerfacecolor','g');
% vv = linspace(0,300,3000);
% mu = log(85);
% sigma = 0.4;
% median_val = exp(mu);
% variance = (exp(sigma^2)-1)*exp(2*mu+sigma^2);
% std_val  = sqrt(variance);
% cdf_lognorm = logncdf (vv,mu,sigma);
% plot(vv,cdf_lognorm,'r-');
% legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
%        ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
%        ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],...
%        'S_i = log_normal','location','southeast');
% %text(4e4,0.4,sprintf('median = %5.1e (m/s),  std. = %5.1e',median_val,std_val));
% 
% grid on
% xlabel('v_{ni} (m/s)');
% ylabel('S_i (-)');
% subplot(122)
% title('lognorm distribution fitting, S_i (v_{Ni})');
% pdf_lognorm = lognpdf(vv,mu,sigma);
% semilogx(vv,pdf_lognorm);
% xlim([8 300]);
% xlabel('v_{ni} (m/s)');
% ylabel('PDF of S_i (m/2)}');
% 
% 
% figure
% subplot(121)
% semilogy(op1(:,1).^2.*sin(op1(:,2)/180*pi),1-op1(:,4),'bo','markerfacecolor','b'); hold on
% semilogy(op2(:,1).^2.*sin(op2(:,2)/180*pi),1-op2(:,4),'bs','markerfacecolor','w');
% semilogy(op3(:,1).^2.*sin(op3(:,2)/180*pi),1-op3(:,4),'b^','markerfacecolor','g');
% legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
%        ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
%        ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],'location','southeast');
% 
% grid on
% xlabel('v_i \times n_{ni} (m/s)^2');
% ylabel('log(1-S_i) (-)');
% 
% subplot(122)
% title('exponential fitting, S_i = 1-exp(a*b)');
% semilogy(op1(:,1).*sin(op1(:,2)/180*pi),1-op1(:,4),'bo','markerfacecolor','b'); hold on
% semilogy(op2(:,1).*sin(op2(:,2)/180*pi),1-op2(:,4),'bs','markerfacecolor','w');
% semilogy(op3(:,1).*sin(op3(:,2)/180*pi),1-op3(:,4),'b^','markerfacecolor','g');
% legend(['Stuttgart-OP1 ' sprintf('v_i = %5.1f m/s',op1(1,1))],...
%        ['Stuttgart-OP2 ' sprintf('v_i = %5.1f m/s',op2(1,1))],...
%        ['Stuttgart-OP3 ' sprintf('v_i = %5.1f m/s',op3(1,1))],'location','southeast');
% grid on
% xlabel('normal incident speed v_{ni} (m/s)');
% ylabel('log(1-S_i) (-)');

% function obj = stuttgart_data_fitting (x,S_i_data)
%   data = S_i_data;
%    Dv = x(2);
%    v50= x(1);
%    S_i_fit = 1-1./(1+(data(:,1)/v50).^Dv);
%    diff = (S_i_fit-data(:,2)).^2;
%    obj  = sum(diff);
% end

