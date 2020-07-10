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
params.osc_count = 10;
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

initial_state = psi_naught(bath_elements,1/params.T,[1;0]);
FULL_H = @(t,y) -1i.*(H_s(bath_elements,y,params.epsilon) + H_b(bath_elements,y) + H_sb(bath_elements,y));  
[t,y] = ode15s(FULL_H,[0,15],initial_state);

pt = zeros(1,length(t));
for ii = 1:length(t)
    rho_spin = TrX(y(ii,:)',2,[2,prod(bath_elements(3,:))]);
    pt(ii) = trace(rho_spin*[1,0;0,-1]);
end
%plot the solution
plot(t',pt,'-o');
title('Poulation difference Tr[\rho_s \sigma_z] solved with ODE 45')
xlabel('Time t');
ylabel('Population Difference');

