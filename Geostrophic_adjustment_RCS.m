%==========================================================================
% Author: Ryan C. Scott, Scripps Institution of Oceanography, UCSD
% Email: rscot025@gmail.com
%==========================================================================
%The 1D Rossby geostrophic adjustment problem via numerical 
%solution of the linear shallow water equations:
%
%                           u_t - f*v   = -g*h_x
%                           v_t + f*u   = 0
%                           h_t + H*u_x = 0
%                        u(x,t=0) = v(x,t=0) = 0
%                       h(x,0) is specified below
%
%==========================================================================
%The solution to the non-dimensional form of these equations is computed
%using explicit forward-backward time-stepping on a staggered spatial grid 
%with periodic boundaries.
%==========================================================================
clear all;
format long g
%==========================================================================
%parameters of the problem
%==========================================================================
%g = 10 m/s^2
%f = 10^-4 1/s
%c = 10 m/s = sqrt(gH)
%H = 1 km 
%a = c/f        barotropic deformation radius
%L = a;         characteristic length scale = 100 km
%T = 1/f        rotational time scale = 10^4 s
%N = fUL/g      geostrophic height deviation scale
%U = gN/fL      geostrophic velocity scale
Bu = 1;         %Burger number (c.f., Vallis p. 200)
%==========================================================================
%grid set up
%==========================================================================
dx = 0.025;     
dt = 5e-3;
t = 0:dt:25;           
x = -25:dx:25;
N = length(t);
J = length(x);
sigma = dt/dx;
%==========================================================================
%initialize empty matrices to be filled by time stepping algorithm (no 
%ghost points yet)
%==========================================================================
u = zeros(N,2*J); 
v = zeros(N,2*J); 
h = zeros(N,2*J); 
%==========================================================================
%u,v defined at odd indicies, h even. The approach taken is to set columns 
%where u,v and h are not defined to NaN. u_0,h_0.5,u_1,h_1.5,... correspond
%to matrix indicies 1,2,3,4,...,2*J-1,2*J
%==========================================================================
for j = 2:2:length(u(1,:))
    u(:,j) = NaN;
    v(:,j) = NaN;
end
for j = 1:2:length(h(1,:))
    h(1,j) = NaN;
end

%==========================================================================
%                                  Set IC
%               different options (A-D) may be specified by user
%==========================================================================
% (A) box at origin of given height
%==========================================================================
%  for j = J-11:2:J+11
%     h(1,j) = 1;
%  end
%==========================================================================
% (B) unit step at origin
%==========================================================================
%  for j = 1:2:length(h(1,:))
%     if j < J
%     h(1,j+1)=1/2;
%     else h(1,j+1) = -1/2;
%     end
%  end
%==========================================================================
% (C) unit gaussian
%==========================================================================
%v-component
%    disp('IC: v(x,0)=cos(x)*exp(-x^2)') 
%    p = exp(-x.^2);
%    v(1,3:2:2*J+2) = p;
 
 disp('IC: h(x,0)=exp(-x^2)')
 %free surface 
 g = exp(-x.^2);

 %g = (1.*sech(x)).^2;
 h(1,4:2:2*J+2) = g;

%==========================================================================
% (D) sech^2(x) profile
%==========================================================================
% g = (1.*sech(x)).^2;
% p = exp(-x.^2);
% 
% h(1,4:2:2*J+2) = g;
% %v(1,3:2:2*J+2) = g;

%==========================================================================
%first two and last two points (4 total) ghost points for (u,v) and h
%append ghost points to grid
%u(n,1) and u(n,2*J+3)
%h(n,2) and h(n,2*J+4)
%==========================================================================
a0 = zeros(length(t),1);
b0 = NaN(length(t),1);
u = [a0 b0 u a0 b0];
v = [a0 b0 v a0 b0];
h = [b0 a0 h b0 a0];
 
h(1,2)= h(1,2*J+2);
h(1,2*J+4) = h(1,4);
u(1,1)= u(1,2*J+1);
u(1,2*J+3) = u(1,3);
v(1,1)= v(1,2*J+1);
v(1,2*J+3) = v(1,3); 
%==========================================================================
%implementation of forward-backward time-stepping on staggered spatial grid
%==========================================================================
 for n = 1:N-1
     for j = 3:2:(2*J+1)   
     u(n+1,j) = u(n,j)+dt.*v(n,j)-(dt./dx).*(h(n,j+1)-h(n,j-1));
     v(n+1,j) = v(n,j)-dt.*u(n+1,j);
     end
     for j = 4:2:(2*J+1)
     h(n+1,j) = h(n,j)-Bu.*(dt./dx).*(u(n+1,j+1)-u(n+1,j-1));
     end %periodic BCs
   h(n+1,2)= h(n+1,2*J+2);
   h(n+1,2*J+4) = h(n+1,4);
   u(n+1,1)= u(n+1,2*J+1);
   u(n+1,2*J+3) = u(n+1,3);
   v(n+1,1)= v(n+1,2*J+1);
   v(n+1,2*J+3) = v(n+1,3);
 end
%==========================================================================
%extract solution where defined, ignoring NaN
%==========================================================================
u0 = zeros(N,J);
v0 = zeros(N,J);
h0 = zeros(N,J);
for j = 1:J
    u0(:,j) = u(:,2*j+1);
    v0(:,j) = v(:,2*j+1);
    h0(:,j) = h(:,2*j+2);
end
%==========================================================================
%Figures  
%==========================================================================
%(1)
figure(1), clf
subplot(3,1,1)
imagesc(h0); colorbar
title('\eta(x,t)')
xtick = [1 501 1001 1501 2001];
set(gca,'XTick',xtick,'XTickLabel',{'-25','-12.5','0','12.5','25'})
xlabel('x/\lambda_R')
ytick= [1 1001 2001 3001 4001 5001 6001 7001 8001];
set(gca,'YTick',ytick,'YTickLabel',{'0','5','10','15','20','25','30','35','40'})
ylabel('ft')
subplot(3,1,2)
imagesc(v0); colorbar
title('v(x,t)')
xtick = [1 501 1001 1501 2001];
set(gca,'XTick',xtick,'XTickLabel',{'-25','-12.5','0','12.5','25'})
xlabel('x/\lambda_R')
ytick= [1 1001 2001 3001 4001 5001 6001 7001 8001];
set(gca,'YTick',ytick,'YTickLabel',{'0','5','10','15','20','25','30','35','40'})
ylabel('ft')
subplot(3,1,3)
imagesc(u0); colorbar
title('u(x,t)')
xtick = [1 501 1001 1501 2001];
set(gca,'XTick',xtick,'XTickLabel',{'-25','-12.5','0','12.5','25'})
xlabel('x/\lambda_R')
ytick= [1 1001 2001 3001 4001 5001 6001 7001 8001];
set(gca,'YTick',ytick,'YTickLabel',{'0','5','10','15','20','25','30','35','40'})
ylabel('ft')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',16,'fontName','Consolas')
%==========================================================================
% (2)
figure(2), clf
aa = ceil(0.25*(J+1));
bb = ceil(0.50*(J+1));
cc = ceil(0.75*(J+1));
dd = ceil((1/3)*(J+1));
ee = ceil((1/8)*(J+1));
subplot(3,1,1)
plot(x,h0(1,:),'k',x,h0(ee,:),'b',x,h0(aa,:),'m',x,h0(dd,:),'c',x,h0(bb,:)...
     ,'g',x,h0(cc,:),'r',x,h0(end,:),'k--')
xlim([-10 10])
box off
title('\eta(x,t)')
subplot(3,1,2)
plot(x,v0(1,:),'k',x,v0(ee,:),'b',x,v0(aa,:),'m',x,v0(dd,:),'c',x,v0(bb,:)...
    ,'g',x,v0(cc,:),'r',x,v0(end,:),'k--')
title('v(x,t)')
legend('IC','ft=1','ft=2','ft=3','ft=4','ft=5')
xlim([-10 10])
box off
subplot(3,1,3)
plot(x,u0(1,:),'k',x,u0(ee,:),'b',x,u0(aa,:),'m',x,u0(dd,:),'c',x,u0(bb,:)...
    ,'g',x,u0(cc,:),'r',x,u0(end,:),'k--')
title('u(x,t)')
xlabel('x/\lambda_R')
xlim([-10 10])
box off
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontName','Consolas')
%==========================================================================
% (3)
figure(3), clf
subplot(3,1,1)
plot(t,h0(:,0.5*(J+1)))
title('\eta(0,t)')
xlim([0 15])
box off
subplot(3,1,2)
plot(t,v0(:,0.5*(J+1)))
title('v(0,t)')
xlim([0 15])
box off
subplot(3,1,3)
plot(t,u0(:,0.5*(J+1)))
title('u(0,t)')
xlabel('ft')
xlim([0 15])
box off
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',14,'fontName','Consolas')
%==========================================================================
%(4)
figure(4), clf
for n=1:100:5001
    plot(x,u0(n,:),'r','LineWidth',1.00008)
    hold on
    plot(x,v0(n,:),'b','LineWidth',1.00008)
    plot(x,h0(n,:),'k','LineWidth',1.00008)
    hold off
    %title('rossby adjustment problem')
    text(-23,0.8,'ft=')
    text(-20,0.8,num2str(t(n)))
    text(-23,0.91,'h(x,0)=e^{-x^2}, u(x,0)=v(x,0)=0')
    xlim([x(1) x(end)])
    xlabel('x/\lambda_R')
    %ylim needs to be changed for each IC
    ylim([-1 1])
    ylabel('\eta (black), v (blue), and u (red)')
    box off
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14,'fontName',...
        'Consolas')
    %to plot every frame:
    %M(n) = getframe(gcf);
    %to plot every kth frame
    l=(n-1)/100+1;
    M(l) = getframe(gcf);
end
    movie(M)
    movie2avi(M,'movie.avi','quality',100,'fps',13)