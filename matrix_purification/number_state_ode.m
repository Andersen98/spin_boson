function [pt,t] = number_state_ode(odehandle,params,time_span)
%Given an ode handle and a set of parameters for spin boson system, as well
%as a time span,
%outputs solution to population over time span 


%construct the bath terms (frequency and coupling
bath_elements = bath_terms(params.w0,params.wf,params.osc_count,params.wc,params.lambda,params.mode_count);

initial_state = psi_naught(bath_elements,1/params.T,[1;0]);

FULL_diff = @(t,y) -1i.*(H_s(bath_elements,y,params.epsilon) + H_b(bath_elements,y) + H_sb(bath_elements,y));  

[t,y] = odehandle(FULL_diff,time_span,initial_state);


pt = zeros(1,length(t));
for ii = 1:length(t)
    rho_spin = TrX(y(ii,:)',2,[2,prod(bath_elements(3,:))]);
    pt(ii) = trace(rho_spin*[1,0;0,-1]);
end

end

