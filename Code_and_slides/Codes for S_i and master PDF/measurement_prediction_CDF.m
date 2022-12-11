function obj = measurement_prediction_CDF(Si_masterProb)

global cdf_in cdf_out nslot d1 d2_d1 nbin_master

N = zeros(nslot,nbin_master+1); % # of particles generated from each incident slot into different daughter sizes
Si= Si_masterProb(1);
massProb = Si_masterProb(2:nbin_master+1);

n_in = 1e6;         % number of incident particles, any number will do

% generate number of particles based on master probability
for k=1:nslot
  % number of incident particle in slot k
  num_k = (cdf_in(k+1,2)-cdf_in(k,2))*n_in;

  % intact particles
  N(k,nbin_master+1) = (1-Si)*num_k;
  for i=1:nbin_master
    N(k,i) = Si*num_k*massProb(i)*(1/d2_d1(i))^3;  % conversion from mass to number
  end
end

% generate CDF in experimental slots
d_grid    = cdf_in(:,1);       % edge of particle bin size
n_rebound = zeros(nslot,1);
cdf_pred  = zeros(size(cdf_out(:,2)));
for k=1:nslot  % each slot "k" of the REBOUND population
  low = d_grid(k);
  high= d_grid(k+1);

  n_rebound(k) = N(k,nbin_master+1); % contribution from unbroken mother particles
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
obj = sum((cdf_pred-cdf_out(:,2)).^2);

endfunction
