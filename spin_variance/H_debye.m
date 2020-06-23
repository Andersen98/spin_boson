%Nb is the number of oscillator bath states
%N is the total number of state (Nb/2)
%input is a vector
function H_bc = H_debye(n_max, bath_terms)

%iterate over each ket state
mode_count = n_max +1;
osc_count = length(bath_terms(1,:));
n_bath = mode_count^osc_count;
n_total = 2*n_bath;
H_bc = zeros(n_total,n_total);
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
    
    %---------------raising operator start------------------------%
        % a^dag |n> = sqrt(1/(2w)) sqrt((n+1)%mode_count) |n+1>
    ket_prefactors =(bath_terms(2,:).*sqrt( mod(my_ket+1,mode_count))./sqrt(2.*bath_terms(1,:)));
    %new_kets = repelem(my_ket,osc_count,1) + diag(ones(1,osc_count));
    
    new_indices = ii + mode_count.^(0:(osc_count-1));
    for jj = 1:osc_count
        prefactor = ket_prefactors(jj);
        row = new_indices(jj);
        col = ii;
        if(prefactor ~= 0)
            H_bc(row,col) = prefactor;
            H_bc(row + n_bath, col + n_bath) = -prefactor;
        end
    end
    
    
    %--------------lowering operator start------------------------%
        % a |n> = sqrt(1/(2w)) sqrt((n)) |n-1>
    ket_prefactors = bath_terms(2,:)./sqrt((2.*bath_terms(1,:))).*sqrt(my_ket);
    new_indices = ii -mode_count.^(0:(osc_count-1));
    for jj = 1:osc_count
        prefactor = ket_prefactors(jj);
        row = new_indices(jj);
        col = ii;
        if(prefactor ~= 0)
            H_bc(row,col) = prefactor;
            H_bc(row + n_bath, col + n_bath) = -prefactor;
        end
    end
    
    
    

    
end

end