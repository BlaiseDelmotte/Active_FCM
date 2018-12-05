clear all
close all

plot_3D = 1;

posdatfiles = dir(fullfile('./', 'FCM_PART_POS_t*.dat'))
pswimdatfiles = dir(fullfile('./', 'FCM_PART_PSWIM_t*.dat'))
Nfiles = size(posdatfiles,1);

f1 = importdata(posdatfiles(1).name);
f2 = importdata(pswimdatfiles(1).name);

rad = 0.647700951441362;
LX = 2*pi;
LY=LX;
LZ=LY;

offset = 1;
Npart = f1.data(offset);

fraction = Npart*4/3*pi*rad^3/(LX*LY*LZ)

XP = zeros(Nfiles,Npart);
YP = zeros(Nfiles,Npart);
ZP = zeros(Nfiles,Npart);

P1x = zeros(Nfiles,Npart);
P1y = zeros(Nfiles,Npart);
P1z = zeros(Nfiles,Npart);

dist_surf = zeros(Npart,Npart);
min_dist_surf = zeros(Nfiles,1);

for ind_files = 1:Nfiles
    f1 = importdata(posdatfiles(ind_files).name);
    f2 = importdata(pswimdatfiles(ind_files).name);
    
    for ind_p = 1: Npart
        XP(ind_files,ind_p) = f1.data(offset + (ind_p-1)*3 + 1);
        YP(ind_files,ind_p) = f1.data(offset + (ind_p-1)*3 + 2);
        ZP(ind_files,ind_p) = f1.data(offset + (ind_p-1)*3 + 3);  
        
        P1x(ind_files,ind_p) = f2.data(offset + (ind_p-1)*3 + 1);
        P1y(ind_files,ind_p) = f2.data(offset + (ind_p-1)*3 + 2);
        P1z(ind_files,ind_p) = f2.data(offset + (ind_p-1)*3 + 3);
        
    end

    
    for ind_p1 = 1:Npart
        for ind_p2 = ind_p1+1:Npart
         dist_x = XP(ind_files,ind_p1) - XP(ind_files,ind_p2) ;
         dist_x = dist_x - LX*fix(dist_x/(LX/2));
         
         dist_y = YP(ind_files,ind_p1) - YP(ind_files,ind_p2);
         dist_y = dist_y - LY*fix(dist_y/(LY/2));
         
         dist_z = ZP(ind_files,ind_p1) - ZP(ind_files,ind_p2);
         dist_z = dist_z - LZ*fix(dist_z/(LZ/2));
         
         dist = sqrt(dist_x^2 + dist_y^2 + dist_z^2);
         dist_surf(ind_p1,ind_p2) = dist - 2*rad;
         
        end        
    end
    
   if isempty(dist_surf(abs(dist_surf)>0))==0
    min_dist_surf(ind_files) = min(min(dist_surf(abs(dist_surf)>0)));
   end
    
    
end



figure
plot(1:Nfiles,min_dist_surf/rad,'k')
xlabel('Dump number')
ylabel('Min dist surf/a')

figure
hold on
color = {'k','b','r','m'};
for ind_p = 1: Npart
    plot3(XP(1,ind_p)/rad,YP(1,ind_p)/rad,ZP(1,ind_p)/rad,['+' color{mod(ind_p,length(color))+1}])
    plot3(XP(Nfiles,ind_p)/rad,YP(Nfiles,ind_p)/rad,ZP(Nfiles,ind_p)/rad,['o' color{mod(ind_p,length(color))+1}])
    plot3(XP(:,ind_p)/rad,YP(:,ind_p)/rad,ZP(:,ind_p)/rad,color{mod(ind_p,length(color))+1})
end
xlabel('x/a')
ylabel('y/a')
zlabel('z/a')
axis equal



if plot_3D==1
    k=0
    [Xs,Ys,Zs]=sphere(20);
    jump_mov=1;
    figure
    
    for n=1:jump_mov:Nfiles

         k=k+1;

        for i=1:Npart
            

             s=surf(XP(n,i)/rad-Xs,YP(n,i)/rad-Ys,ZP(n,i)/rad-Zs,0.5*ones(size(Xs)));
             if mod(i,5)==1
                set(s,'FaceColor','green','EdgeColor','none');
             elseif mod(i,5)==2
                set(s,'FaceColor','red','EdgeColor','none');
             elseif mod(i,5)==3
                set(s,'FaceColor','yellow','EdgeColor','none');
             elseif mod(i,5)==4
                set(s,'FaceColor','magenta','EdgeColor','none');      
             elseif mod(i,5)==0
                set(s,'FaceColor','cyan','EdgeColor','none');      
             end

             hold on
        end

        quiver3(XP(n,:)/rad,YP(n,:)/rad,ZP(n,:)/rad,P1x(n,:),P1y(n,:),P1z(n,:),1,'LineWidth',2,'MaxHeadSize',0.1)
 
         axis([0 LX/rad 0 LY/rad 0 LZ/rad])
         
       view(3); 
       camlight
       lighting gouraud

        grid on
        box on

          M(k)=getframe;

%         pause
        clf
    end
    
end



save_pos  = input('Save ?')
if save_pos==1
    save('2_sq_positions_Maury.mat')
end