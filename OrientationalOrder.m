%clear 
close all
clc

%load 0.5tf.mat;
%load strong_edgewave.mat;
%load crystallization.mat;

for jj = 1:1
    for zz = 1:1000
  
%z0 = z_save(zz,:,jj);
z0 = z(zz,:);

[~,ind] = min(abs(z0));
zi = z0(ind);
K = 6;

temp = abs(z0-zi);
temp(ind) = NaN;
IDX = zeros(1,6);
for kk = 1:K
[~,ind2] = nanmin(temp);
IDX(kk) = ind2;
temp(ind2) = NaN;
end

[~,IDX2] = sort(angle(z0(IDX)));
IDX = IDX(IDX2);

% 
figure(1)
 clf
 plot(z0,'.')
 hold on
 plot(z0(ind),'ok');
 hold on
 %plot(z0(IDX),'sqr')
 plot(z0(IDX(1)),'pm')
 plot(z0(IDX(2)),'vg')
 plot(z0(IDX(3)),'^b')
 plot(z0(IDX(4)),'^k')
 plot(z0(IDX(5)),'^r')
 plot(z0(IDX(6)),'^y')
 
 angles = wrapTo2Pi(angle(z0([IDX(2:end) IDX(1)])) - angle(z0(IDX)));

psi6(zz,jj) = 1/6*sum(exp(1i*6*angles));

    end
    disp(num2str(jj))
end

figure
hist(abs(psi6),40)