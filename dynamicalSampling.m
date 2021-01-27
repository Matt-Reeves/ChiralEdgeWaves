clear
clc
close all

N = 1000;
R = 1;
Omega = N/R^2;
gamma0 = 0.005;
rho0 = N/pi/R^2;
Rt = @(tt) R*sqrt(1 + 2*pi*rho0*gamma0*tt);
kappa = ones(1,N);
theta = linspace(0,2*pi,100);
z0 = R*sqrt(rand(N,1)).*exp(1i*2*pi*rand(N,1));

figure(1)
clf
plot(z0,'.')
xlim([-1 1]*R*1.2)
ylim(xlim)
hold on
plot(exp(1i*theta),'--k')

%%
Nt = 1;
tf = 12/Omega/gamma0;
ss = @(tt) log(1+rho0*gamma0*2*pi*tt)/R^2/rho0/gamma0/2/pi;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
s = linspace(0,tf,Nt);

    [t,ztemp] = ode45(@(t,z) PointVortexPlane(t,z,kappa,gamma0,N,rho0,R),[0 s],z0,options);
    z = ztemp(end,:).';
    %[~,I] = max(abs(z));
    %z(I) = (abs(z(I))+0.1)*exp(1i*angle(z(I)));


%%
z0 = z;
Nt = 1000
s = linspace(0,0.1*tf,Nt);
z = zeros(Nt,N);
for ii = 1:Nt
    disp(num2str(ii))
[t,ztemp] = ode45(@(t,z) PointVortexPlane(t,z,kappa,0,N,rho0,R),[0 s(2)-s(1)],z0,options);
z0 = ztemp(end,:);
z(ii,:) = ztemp(end,:);
end
%%
% rad(kk) = max(abs(z(end,:)))
% E(kk) = Energy(z(end,:),N)/N^2
% M(kk) = mean(abs(z(end,:)).^2)

h = figure(3);
axis square
for jj = 1:1:Nt
    cla
    %plot(z(jj,:)/Rt(t(jj))*(1+rho0*gamma0*2*pi*t(jj))^(-1i*Omega/gamma0/2/pi/rho0),'.')
    plot(z(jj,:),'.r')
    hold on
    plot(exp(1i*theta),'-k')
    xlim([-1 1]*1.5)
    ylim(xlim)
    title(num2str(jj))
    drawnow
    F(jj) = getframe(gcf);
    E(jj) = Energy(z(jj,:),N);
end
%%
Nbins = 100;
dx = 1/Nbins
[rho,r] = hist(abs(z(:)),Nbins)

dr = r(2)-r(1);
rho = rho./r/dr/sum(rho);
figure
plot(r,rho,'.--')

% figure(10)
% hold on
% plot(Omega*t,E/N^2)

%%
% % create the video writer with 1 fps
%   writerObj = VideoWriter('vortexmatter.avi');
%   writerObj.FrameRate = 60;
%   % set the seconds per image
%     % open the video writer
%     open(writerObj);
%     % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

%end



    