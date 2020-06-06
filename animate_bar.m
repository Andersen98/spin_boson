clear ket
clear xTickLabel
close all
ket = zeros(N,1);
ket(1:N/2) = 1;
ket(N/2+1) = 0;
ket = (1/norm(ket)).*ket;
h = figure;
axis tight manual %ensures getframe() returns consistent size
flag = 1;
filename = '3_osc_weak_coupling_spinUp_uniform.gif';
xTickLabel = cell(N/2,1);

%time setup
t0 = 0;
tf = 3000;
tcount = 400;
timevec =  linspace(t0,tf,tcount);

%scatter bin setup
spin_hist = zeros(2,tcount);

for ii = 1:(N/2)
    
    xTickLabel(ii,1) = {strcat('|',dec2base(ii-1,params.n_max +1,params.n_osc),'>')};
end
for n = 1:tcount
    
    %set current time
    t_c = timevec(n);
   
    %time evolve global state
    x = abs(expm(1i.*t_c.*H)*ket).^2;
    
    %create mask to select our up and down components
    spinDownMask = zeros(1,N);
    spinDownMask(1:N/2) = 1;
    spinDownMask = spinDownMask==1;
    
    %create 2-row figure 
    layout_ax = tiledlayout(2,1);
    
    %populate top plot
    ax1 = nexttile;
    
    %create grouped data
    Y = zeros(2,N/2);
    spindowntmp = x(~spinDownMask)';
    spinuptmp = x(spinDownMask)';
    Y(1,:) = spinuptmp;
    Y(2,:) = spindowntmp;
    xrng = 0:(N/2-1);
    hb1 = bar(ax1,xrng,Y);
    
    
    %populate bottom plot 
    %row 1 is spin up, row 2 is spin down
    spin_hist(1,n) = sum(x(spinDownMask));
    spin_hist(2,n) = sum(x(~spinDownMask));
    spin_hist(3,:) = abs(spin_hist(1,:)-spin_hist(2,:));
    ax2 = nexttile;
    up_s = plot(ax2,timevec(1:n),spin_hist(1,1:n));
    hold on
    down_s = plot(ax2,timevec(1:n),spin_hist(2,1:n));
    %diff_s = plot(ax2,timevec(1:n),spin_hist(3,1:n));
    hold off
    
    
    %apply formatting top plot
    set(hb1(1),'FaceColor','r');
    set(hb1(2),'FaceColor','b');
    set(hb1(1),'BarWidth',1)
    set(hb1(2),'BarWidth',1)
    title(ax1,'Time evolution of state magnitudes (3 Oscillators w/ range n=0,1,2)')
    xlabel(ax1,'States: red m=1/2 blue m=-1/2') 
    ylabel(ax1,'Magnitude of State Vector') 
    %set ticks to give osc state
    set(ax1,'XTickLabel',xTickLabel);
    set(ax1,'XTick',0:(N/2-1))
    ylim(ax1,[0,1])
    
    
    
    %apply formatting to bottom plot
    ylim(ax2,[0,max([spin_hist(3,:),1])]);
    xlim(ax2,[timevec(1),timevec(end)]);
    set(up_s(1),'Color','r');
    set(down_s(1),'Color','b');
    ylabel(ax2,'Probability of spin')
    xlabel(ax2,'Time (arbitrary)')
    
    
    %apply full tiled layout formatting
    layout_ax.Parent.Position=[350 421 1521 548];
    
    if(flag ==1)
        
       % set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        flag = 0;
    end
    drawnow
    %capture the plot as image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    %write to gif
    if n==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
    end
    [norm(spinuptmp), norm(spindowntmp)]
    
end
