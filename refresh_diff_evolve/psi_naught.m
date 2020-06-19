function psi = psi_naught(bath_terms,beta,spin_state)

mode_counts = bath_terms(3,:);
omega_vec = bath_terms(1,:);
osc_count = length(omega_vec);
bath_count = prod(mode_counts);

%construct the bath states
psi = zeros(bath_count,1);
psi_last = zeros(bath_count,1);
counter_vec = zeros(1,osc_count);
for ii = 1:(bath_count-1)
     %increment counter by one and then iter. through terms
     counter_vec(1) = counter_vec(1) + 1;
    for jj = 1:1:osc_count
      if counter_vec(jj)== mode_counts(jj)
          counter_vec(jj) = 0;
          counter_vec(jj+1) = counter_vec(jj+1) +1;
          
      end 
    end
    counter_vec
    %prepare for direct product of osc1*(osc2*(...*oscn))
    tmp_state = zeros(mode_counts(end),1);
    tmp_state(counter_vec(end)+1) = 1;
    tmp_energy = bath_terms(1,end).*(counter_vec(end)+1/2);
    therm_prob = exp(-beta.*tmp_energy);
    tmp_state = sqrt(therm_prob).*tmp_state;
    psi = tmp_state;
    for jj = (osc_count-1):-1:1
        
        tmp_state = zeros(mode_counts(jj),1);
        tmp_state(counter_vec(jj)+1) = 1;
        tmp_energy = bath_terms(1,jj).*(counter_vec(jj)+1/2);
        therm_prob = exp(-beta.*tmp_energy);
        tmp_state = sqrt(therm_prob).*tmp_state;
        
        psi = kron(tmp_state,psi);
    end
    psi_last = psi_last + psi;

end

psi = 1/norm(psi_last).*psi_last;
psi = kron(spin_state,psi);

end