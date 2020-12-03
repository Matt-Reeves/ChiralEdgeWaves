clear
clc
close all
Nvec = [2 5 10 20 50 100 200 500 1000 2000]
%for kk =1:length(Nvec)
savedat = struct;
for kk = 1:length(Nvec)
N = Nvec(kk);
R = 1;
Omega = N/R^2;
gamma0 = 0.01;
rho0 = N/pi/R^2;
Rt = @(tt) R*sqrt(1 + 2*pi*rho0*gamma0*tt);
kappa = ones(1,N);
theta = linspace(0,2*pi,100);

%for run = 1
%disp(num2str(run))
z0 = R*sqrt(rand(N,1)).*exp(1i*2*pi*rand(N,1));
%z0 = 0.25*(randn(N,1) + 1i*randn(N,1));

figure(1)
clf
plot(z0,'.')
xlim([-1 1]*R*1.2)
ylim(xlim)
hold on
plot(exp(1i*theta),'--k')

%%
Nt = 2000/8;
tf = 4/Omega/gamma0;
ss = @(tt) log(1+rho0*gamma0*2*pi*tt)/R^2/rho0/gamma0/2/pi;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
s = linspace(0,tf,Nt);
%s = ss(t);
[t,z] = ode45(@(t,z) PointVortexPlane(t,z,kappa,gamma0,N,rho0,R),s,z0,options);
eval(['savedat.N' num2str(N) ' = z(end,:)'])
%z_save(:,kk) = z(end,:);
%u_save(:,kk) = getVelocity(z_save(:,kk),kappa,N);
%%

% zbin = 0:0.005:1;
% for jj = 1:length(zbin)-1
% ind = z_save(:) > zbin(jj) & z_save(:) < zbin(jj+1);
% U(jj) = mean(abs(u_save(ind)));
% U2(jj) = mean(angle(u_save(ind)));
% end

% rad(kk) = max(abs(z(end,:)))
% E(kk) = Energy(z(end,:),N)/N^2
% M(kk) = mean(abs(z(end,:)).^2)
figure(3)
plot(z(end,:),'.r')
drawnow
disp(num2str(kk))
end
% h = figure(3);
% axis square
% for jj = 1:1:Nt
%     cla
%     %plot(z(jj,:)/Rt(t(jj))*(1+rho0*gamma0*2*pi*t(jj))^(-1i*Omega/gamma0/2/pi/rho0),'.')
%     plot(z(jj,:),'.r')
%     hold on
%     plot(exp(1i*theta),'-k')
%     xlim([-1 1]*1.5)
%     ylim(xlim)
%     title(num2str(jj))
%     drawnow
%     F(jj) = getframe(gcf);
%     E(jj) = Energy(z(jj,:),N);
% end

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



    