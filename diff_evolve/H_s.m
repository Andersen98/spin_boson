function dydt = H_s(t,y,n_max,bath_elements,epsilon)
%H_s Gives computes the output of H_s * y

%define variables
mode_count = n_max +1;
osc_count = length(bath_elements(1,:));
n_bath = mode_count^osc_count;
b_ones = ones(1,n_bath);
n_total = 2*n_bath;
%define initial derrivative
dydt = zeros(n_total,1);
%(\sigma x + epsilon \sigma_z)psi
    %=SIMGAX*PSI +SIGMAZ*PSI

%Sigma z*epsilon*
dydt(1:n_bath,1) = epsilon.*dydt(1:n_bath,1);
dydt((n_bath+1):n_total,1) = -1.*epsilon.*dydt((n_bath+1):n_total,1);

%sigma x just swaps the vector
tmpvec = zeros(n_total,1);
tmpvec(1:n_bath,1) = y(n_bath+1:n_total,1);
tmpvec(n_bath+1:n_total,1) = y(1:n_bath);

dydt = dydt + tmpvec;
        
end
end

