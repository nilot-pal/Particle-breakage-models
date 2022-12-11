function sum_of_prob = unityConstraint(Si_masterProb)

  global cdf_in cdf_out nslot d1 d2_d1 nbin_master

  sum_of_prob = sum(Si_masterProb(2:nbin_master+1)) - 1.0;  % sum of probability is 1.0


endfunction
