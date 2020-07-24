clear variables;
%energy is in units of Delta so we are left with the
%physics params - (k_b=h_bar=1, energy in units of delta)
params.epsilon = 0;
params.lambda = 2.5;        
params.T = .2; 

%simulation parameters
    %osc_count is the number of oscillators to be used in the bath
    %n_max is the maximal mode for each oscillator i,e, can range from
    %0,1,2,3...n
params.osc_count = 4;
params.mode_count = 2;

%Descretized Debey spectral density - need us to choose a
    %charactaristic frequency for the bath wc
    %core frequency range: min to max [w0,wf]
    %coupling strength (renormalizatin energy) lambda
params.wc= .25;
params.wf = 2*params.wc;
params.w0 = params.wf./(2.*params.osc_count);
params.lambda = 2.5;

%construct the bath terms (frequency and coupling
bath_elements = bath_terms(params.w0,params.wf,params.osc_count,params.wc,params.lambda,params.mode_count);
num_bath = prod(bath_elements(3,:));
initial_state = psi_naught(bath_elements,1/params.T);
 %H_sb(bath_elements,initial_state)
%H_s(bath_elements,initial_state,params.epsilon)
%H_b(bath_elements,initial_state)
 H_sb(bath_elements,initial_state)+ H_s(bath_elements,initial_state,params.epsilon) +H_b(bath_elements,initial_state)
 H_sb_w = @(psi) H_sb(bath_elements,psi);
 H_b_w = @(psi) H_b(bath_elements,psi);
 H_s_w = @(psi) H_s(bath_elements, psi, params.epsilon);

     
 FULL_H = @(t,y) -1i.*(anc_H(H_s_w,H_b_w,H_sb_w,y,num_bath));
[t,y] = ode15s(@(t,y) FULL_H(t,y),[0,12],initial_state);

pt = zeros(1,length(t));
for ii = 1:length(t)
    rho_spin = TrX(y(ii,:)',[1,3],[num_bath,2,num_bath]);
    pt(ii) = trace(rho_spin*[1,0;0,-1]);
end
%plot the solution
plot(t',pt,'-o');
title('Poulation difference Tr[\rho_s \sigma_z] solved with ODE 45')
xlabel('Time t');
ylabel('Population Difference');
