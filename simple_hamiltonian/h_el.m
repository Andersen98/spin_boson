function h_full = h_el(spin_b,osc_b_str,spin_k,osc_k_str,params)
%H_TERM calculates the <bra|H|ket> of the spin-boson Hamiltonian

osc_b = zeros(1,params.n_osc);
osc_k = zeros(1,params.n_osc);
for ii = 1:params.n_osc
    osc_b(ii) = str2int(osc_b_str(ii));
    osc_k(ii) = str2int(osc_k_str(ii));
end

%find coupling term
h_couple = 0;
if(spin_b == spin_k)
    for q = 1:params.n_osc
        if(osc_b(q)==(osc_k(q)-1))
            h_couple = h_couple + spin_k.*params.coupling(q).*sqrt(osc_k(q));
        elseif (osc_b(q) == (osc_k(q)+1)) 
            h_couple = h_couple + spin_k.*params.coupling(q).*sqrt(osc_k(q)+1);
        end
    end
    h_couple = spin_k * h_couple;
end
%find sho term
h_sho = 0;
if spin_b == spin_k
    
    for ii = 1:params.n_osc
        if strcmp(osc_b_str,osc_k_str)
        h_sho = h_sho + osc_b(ii)*params.omega;
        end
    end
end

%find Hopping term
h_hop = 0;
if spin_b ~= spin_k && strcmp(osc_b_str,osc_k_str)
    h_hop = params.delta;
end
h_full = h_couple + h_sho + h_hop;



end

