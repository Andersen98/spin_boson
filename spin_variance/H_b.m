function H_bath = H_b(n_max,bath_elements)
%outputs the oscillator bath hamiltonian (truncated of course)

%define variables
mode_count = n_max +1;
osc_count = length(bath_elements(1,:));
n_bath = mode_count^osc_count;
n_total = 2*n_bath;
H_bath = zeros(n_total,n_total);
my_ket = zeros(1,osc_count);
for ii = 1:n_bath
    dec = ii-1;
    %change the base into radix = mode_count
    for jj = 1:osc_count
        quo = floor(dec./mode_count);
        rem = mod(dec,mode_count);
        my_ket(jj) = rem;
        dec = quo;
    end
    %now my_ket is translated from index ii
    diagonal = sum(bath_elements(1,:).*my_ket) + sum(bath_elements(1,:)./2);
    H_bath(ii,ii) = diagonal;
    H_bath(ii+n_bath,ii+n_bath) = diagonal;
    
end
end

