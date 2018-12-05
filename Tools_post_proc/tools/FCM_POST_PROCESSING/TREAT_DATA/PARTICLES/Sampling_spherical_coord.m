clear all
close all


L = 2*pi;
a = 0.081;
Ntheta = 63;
Nphi = 31;
Nr = 60;

r_reg = linspace(2*a,8*a,Nr);
phi_reg = linspace(0.1,2*pi,Nphi);
theta_reg = linspace(0.1,pi-0.1,Ntheta);

dr_reg = r_reg(2) - r_reg(1);
dphi_reg = phi_reg(2) - phi_reg(1);
dtheta_reg = theta_reg(2) - theta_reg(1);



Vol_reg = zeros(Nr,Nphi,Ntheta);
Vol_shell = zeros(1,Nr);
Vol_ring = zeros(Nr,Ntheta);
Vol_ring_ana = zeros(Nr,Ntheta);
for i = 1:Nr
    Vol_shell(i) = 0;    
        for k = 1:Ntheta            
            Vol_ring(i,k) = Vol_ring(i,k) + 2*pi*r_reg(i)^2*sin(theta_reg(k))*dr_reg*dtheta_reg;
            Vol_ring_ana(i,k) = 2*pi*1/3*((r_reg(i)+dr_reg)^3-r_reg(i)^3)*(-cos(theta_reg(k)+dtheta_reg)+cos(theta_reg(k)));
           for j = 1:Nphi            
            Vol_shell(i) = Vol_shell(i) + r_reg(i)^2*sin(theta_reg(k))*dr_reg*dtheta_reg*dphi_reg;
            dtheta_irreg(i,k) =  1/(sin(theta_reg(k))*r_reg(i)^2);
            Vol_reg(i,j,k) = r_reg(i)^2*sin(theta_reg(k))*dr_reg*dtheta_reg*dphi_reg;
            Vol_irreg(i,j,k) = r_reg(i)^2*sin(theta_reg(k))*dr_reg*dtheta_irreg(i,k)*dphi_reg;
        end
    end
end

figure
pcolor(dtheta_irreg)
colorbar
max(max(dtheta_irreg))
min(min(dtheta_irreg))

figure
pcolor(squeeze(Vol_reg(floor(Nr/2),:,:)))
colorbar
max(max(max(Vol_reg)))
min(min(min(Vol_reg)))

figure
pcolor(squeeze(Vol_irreg(floor(Nr/2),:,:)))
colorbar

max(max(max(Vol_irreg)))
min(min(min(Vol_irreg)))


figure
hold on
plot(r_reg,Vol_shell)
plot(r_reg,4/3*pi*((r_reg+dr_reg).^3-r_reg.^3),'+k')

figure
hold on
plot(r_reg,Vol_shell)
plot(r_reg,4/3*pi*((r_reg+dr_reg).^3-r_reg.^3),'+k')
plot(r_reg,4*pi*(r_reg.^2)*dr_reg,'ok')

figure
pcolor(r_reg,theta_reg,Vol_ring')
colorbar
figure
pcolor(r_reg,theta_reg,Vol_ring_ana')
colorbar