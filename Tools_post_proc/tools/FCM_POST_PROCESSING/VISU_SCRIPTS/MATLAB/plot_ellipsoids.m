%  plot_ellipsoids.m
clear all
close all

radii_data = importdata('FCM_PART_RADII.dat');
centers_data = importdata('FCM_PART_POS.dat');
quats_data = importdata('FCM_PART_ORIENT.dat');


Np = radii_data.data(1)

for i = 1:Np
    radii(i,1:3)=radii_data.data((i-1)*3+2:i*3+1)';
    centers(i,1:3)=centers_data.data((i-1)*3+2:i*3+1)';
    quats(i,1:4)=quats_data.data((i-1)*4+2:i*4+1)';
end


LXMAX = 2*pi;
LYMAX = LXMAX;
LZMAX = LXMAX;
 


colors = ones(21);


close all
figure(1)
ifoff = 0;


for i = 1:Np

        
        solid_volume(i) = 4/3*prod(radii(i,1:3))*pi;
        
        q = quats(i,1:4);
        U = (2*q(1)^2 - 1)*eye(3) + 2*(q(2:4)'*q(2:4) + q(1)*skew(q(2:4)')) ; 

        
        [x ,y, z] = ellipsoid(0,0,0,radii(i,1), radii(i,2), radii(i,3));
        
        a = kron(U(:,1),x);
        b = kron(U(:,2),y);
        c = kron(U(:,3),z);

        
        data = a+b+c;
        n = size(data,2);
        
        
        x = data(1:n,:)+centers(i,1); 
        y = data(n+1:2*n,:)+centers(i,2);
        z = data(2*n+1:end,:)+centers(i,3);

        surf(x,y,z,i*colors)
       
   
        if ifoff==0; hold on; ifoff=1; end
        
        plot3( centers(i,1), centers(i,2), centers(i,3),'+k' )
        text(centers(i,1), centers(i,2), centers(i,3), num2str(i))

        pause
        
%         (ir-1)/7
        
%    U
%    pause
 %       if mod(ir-1,15)==0; pause; hold off; ifoff=0; end
end



xlabel('x')
ylabel('y')
zlabel('z')

alpha_p = sum(solid_volume)/(LXMAX*LYMAX*LZMAX)

axis([0 LXMAX 0 LYMAX 0 LZMAX])
    
shading interp
camlight
lighting gouraud
 alpha(0.3)
 view(3)
 
 title(['N_p = ' num2str(Np) ', \alpha_p = ' num2str(alpha_p) ', radii = ' num2str(radii(i,:)) ])