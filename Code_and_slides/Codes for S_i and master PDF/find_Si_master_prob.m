clear
clc
%% step 1: read in CDF of incident and rebound population and interpolate to fine, uniform slots
%raw_cdf_in      = load('incident_CDF_exp.txt');
%raw_cdf_out     = load('rebound_CDF_exp.txt');

%raw_cdf_in      = dlmread("Input_files/Particle size distribution/CDF_bef_impact_Stainless_steel_target_speed_OP3.txt",',',1,0);
%raw_cdf_out     = dlmread("Input_files/Particle size distribution/CDF_aft_impact_Stainless_steel_target_speed_OP3_angle_90_deg.txt",',',1,0);

%raw_cdf_in      = dlmread("Input_files/Particle size distribution/CDF_bef_impact_Ti6AlV4_target_speed_OP3.txt",',',1,0);
%raw_cdf_out     = dlmread("Input_files/Particle size distribution/CDF_aft_impact_Ti6AlV4_target_speed_OP3_angle_90_deg.txt",',',1,0);

raw_cdf_in = load("../../Particle_breakage/Levy_particle_breakage_quartz/Datasets/Particle_size_distribution/quartz_on_chromium_steel_PSD_bef_impact.txt");
cd Input_files/Levy_particle_breakage_quartz/Datasets/Particle_size_distribution;
cd ../../Particle_breakage/Levy_particle_breakage_quartz/Datasets/Particle_size_distribution
%raw_cdf_in = load('quartz_on_chromium_steel_PSD_bef_impact.txt');
%raw_cdf_out = load('quartz_on_chromium_steel_PSD_aft_impact_128_m_s.txt');
raw_cdf_in = load('quartz_on_silicate_glass_PSD_bef_impact.txt');
raw_cdf_out = load('quartz_on_silicate_glass_PSD_aft_impact_98_m_s.txt');

extra_row = [19,0];
k = 0;
raw_cdf_in = [raw_cdf_in(1:k,:); extra_row; raw_cdf_in(k+1:end,:)];
raw_cdf_out = [raw_cdf_out(1:k,:); extra_row; raw_cdf_out(k+1:end,:)];

raw_cdf_in(:,1) = (raw_cdf_in(:,1)+raw_cdf_out(:,1))/2; % make slot position the same
raw_cdf_out(:,1)= raw_cdf_in(:,1);                      % make slot position the same

nraw            = length(raw_cdf_in(:,1));              % number of slot edges originally

global cdf_in cdf_out nslot d1 d2_d1 nbin_master
% cdf_in              interpolated CDF of incident population
% cdf_out             interpolated CDF of rebound population
nslot  = 30;          % number of slots for binning of particles after refinning bin size
d1_min = raw_cdf_in(1,1);
d1_max = raw_cdf_in(nraw,1);
d1_slot_size = (d1_max-d1_min)/nslot;
d1     = zeros(nslot,1);    % nominal diameter in each slots
for i=1:nslot
  d1(i) = d1_slot_size*(i-0.5) + d1_min;
end

nbin_master = 16;           % mass-based master probability profile
d2_d1       = 0.5/nbin_master:1.0/nbin_master:1-0.5/nbin_master;
d2_d1       = d2_d1';       % nominal diameter of daughter particle relative to mother particle

cdf_in      = zeros(nslot+1,2);
cdf_out     = cdf_in;
cdf_in(:,1) = linspace(d1_min,d1_max,nslot+1);
cdf_out(:,1)= cdf_in(:,1);
for i=1:nslot+1
  cdf_in(i,2) = interp1(raw_cdf_in(:,1),raw_cdf_in(:,2),(i-1)*d1_slot_size + d1_min,'linear');
  cdf_out(i,2)= interp1(raw_cdf_out(:,1),raw_cdf_out(:,2),(i-1)*d1_slot_size + d1_min,'linear');
end
%figure(1)
%plot(raw_cdf_in(:,1),raw_cdf_in(:,2),'bs'); hold on;
%plot(cdf_in(:,1),cdf_in(:,2),'b+-'); hold on
%plot(raw_cdf_out(:,1),raw_cdf_out(:,2),'ro'); hold on;
%plot(cdf_out(:,1),cdf_out(:,2),'r+-');
%grid
%xlabel('particle diameter (\mum)');
%ylabel('truncated CDF');
%title('Linearly interpolated Stuttgart data, OP1, 40^o impact');
%legend('incident','rebound');
%rectangle('position',[0 0 cdf_in(1,1) 1],'EdgeColor','k','FaceColor',[0.86 0.86 0.86]);
%text(cdf_in(1,1)/3,0.4,'invisible');

%% Step 2: find the solution
Si    = 0.6;
guess = 0.4;
initial_guess = zeros(nbin_master+1,1);
initial_guess(nbin_master+1-4) = guess/4;
initial_guess(nbin_master+1-3) = guess;
initial_guess(nbin_master+1)   = 1-guess*1.25;

lowerBound = zeros(nbin_master+1,1);
upperBound = ones(nbin_master+1,1);
%lowerBound(1) = 0.99; % Varying the lower/upper bound to see how fitting may behave , sensitivity study
%upperBound(1) = 0.2;   % sensitivity study

[x, obj, info, iter, nf, lambda] = sqp (initial_guess, @measurement_prediction_CDF, ...
                                                       @unityConstraint, [], ...
                                                       lowerBound,upperBound,500);
%% step 3: show solution
Si       = x(1);
massProb = x(2:nbin_master+1);
d_grid   = cdf_in(:,1);       % edge of particle bin size

% reconstruct CDF from Si and mass-based master probability
n_in = 1e6;
N = zeros(nslot,nbin_master+1);
for k=1:nslot
  num_k = (cdf_in(k+1,2)-cdf_in(k,2))*n_in;
  N(k,nbin_master+1) = (1-Si)*num_k;
  for i=1:nbin_master
    N(k,i) = Si*num_k*massProb(i)*(1/d2_d1(i))^3;  % conversion from mass to number
  end
end
n_rebound = zeros(nslot,1);
cdf_pred  = zeros(size(cdf_out(:,2)));
for k=1:nslot  % each slot "k" of the REBOUND population
  low = d_grid(k);
  high= d_grid(k+1);

  n_rebound(k) = N(k,nbin_master+1);
  for m=k+1:nslot	% contribution from mother particles in slots with larger diameters
    for i=1:nbin_master
      dia = d2_d1(i)*d1(m);
        if(dia< high && dia > low)
          n_rebound(k) = n_rebound(k) + N(m,i);
        end
     end
  end

  cdf_pred(k+1) = cdf_pred(k) + n_rebound(k);
end
cdf_pred = cdf_pred/cdf_pred(nslot+1);
error = sum((cdf_pred-cdf_out(:,2)).^2);
display(Si)
display(sqrt(error/nslot))

% small probability is rounded to 0
for i=1:nbin_master
  if(abs(massProb(i)) <1e-5)
    massProb(i) = 0;
  end
end

figure(2)
subplot(121)
bar(d2_d1,massProb,1.0,'facecolor', [0.86 0.86 0.86], 'edgecolor', 'b'),
xlabel('d_{2}/d_{1}','FontSize',16);
ylabel('probability - mass based','FontSize',16);
title(['master curve, sum of probability = ' num2str(sum(massProb))],'FontSize',16);
set(gca,'FontSize',20)

subplot(122)
plot(d_grid,cdf_pred,'k-');
hold on
plot(d_grid,cdf_out(:,2),'rs');
xlabel('particle diameter (\mum)','FontSize',16);
ylabel('CDF','FontSize',16);
legend('reconstructed','measured','location','southeast','FontSize',16);
title(sprintf('S_i = %4.2f%% error = %4.2e',Si*100, sqrt(error/nslot)),'FontSize',16);
set(gca,'FontSize',20)

% show in terms of number-based probability
%num_prob = zeros(nbin_master,1);
%num_prob = d2_d1.^-3.*massProb;
%num_prob = num_prob/sum(num_prob);
%figure(3)
%bar(d2_d1,num_prob,1.0,'facecolor', [0.86 0.86 0.86], 'edgecolor', 'b'),
%loglog(d2_d1,num_prob,1.0,'bs-'),
%xlabel('d_{2}/d_{1}');
%ylabel('count - number');
%title(['master curve, sum of probability = ' num2str(sum(massProb))]);
