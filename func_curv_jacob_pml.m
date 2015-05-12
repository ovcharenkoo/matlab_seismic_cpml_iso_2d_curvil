%Constructs curvilinear mesh with regular regions for pmls and its Cartesian analogus. Calculates
%Jacobian 
% J=[dksi_dx dksi_dy;
%    deta_dx deta_dy];
% Ji=[dx_dksi dx_deta;
%     dy_dksi dy_deta];

% Input arguments:
% nx - number of nx grid points
% ny - number of ny grid points
% xmin, xmax - min and max values over OX
% ymin, ymax - min and max values over OY
% argument - string, for sin(argument), as function of x, that depend on xmin and xmax
% dxx - x spacing for Cartesian grid
% dyy - y spacing for Cartesian grid

%Output:
%ksi - Cartesian x
%eta - Cartesian y

%xx - curvil x
%yy - curvil y

%To test:
% [x,y, ksi,eta,J] = func_curv_jacob_pml(30,30,5,0,100, 0, 100,'-(2*pi*x/max(x)+0.25*pi)',10,10, 0.1,true);
%[curvil x, curvil y, cartes x, cartes y, J of curvil ---> cartes]= ...

function [xx,yy,ksi,eta,J] = func_curv_jacob_pml(nx,ny,npml,xmin,xmax,ymin,ymax,argument,dxx,dyy,curvature,showornot)
%Step of grid
dxxt=xmax/nx;
dyyt=ymax/ny;

% npml=npml+1;

%Actual max yvalue, without pmls
ymaxr=ymax-2*dyyt*npml;
ymidr=(ymaxr+ymin)/2;

xmaxr=xmax-2*npml*dxxt;
nxr=nx-2*npml;
nyr=ny-2*npml;

x=linspace(xmin,xmaxr,nxr+1);

phi=inline(argument,'x');

%curvature=0.1;
%Bottom part
y=curvature*ymaxr*sin(phi(x));
y=ymidr+y;
%Top part, shifted by pi to be in phase
y1=-ymidr-curvature*ymaxr*sin(phi(x)-pi);

%Phys. domain
xx=zeros(nxr+1, nyr+1);
yy=zeros(nxr+1, nyr+1);

% %Comp. domain
ksi=zeros(nx+1,ny+1);
eta=zeros(nx+1,ny+1);
% dxx=1.0;
% dyy=1.0;

for i=1:nx+1
    for j=1:ny+1
        ksi(i,j)=dxx*(i-1);
        eta(i,j)=dyy*(j-1);
    end
end


%Calculate x-ksi in phys domain
for i=1:nxr+1
        xx(i,:)=x(1,i);
end

%cksi=x;  %x in real space
nyro=round(nyr/2);
for i=1:nyr+1
    if i<=nyro
        ceta=(i-1)*y/nyro; %y in real space
        yy(:,i)=ceta(1,:);
    else
        ceta=ymaxr+(nyr-i+1)*y1/(nyr-nyro); %y in real space
        yy(:,i)=ceta(1,:);
    end
end


%for replacing of middle part of cartesian by curvil
xxt=zeros(nx+1,ny+1);
yyt=zeros(nx+1,ny+1);

for i=1:nx+1
    for j=1:ny+1
        xxt(i,j)=dxxt*(i-1);
        yyt(i,j)=dyyt*(j-1);
    end
end

xx=(xx+(dxxt*npml));
yy=(yy+(dyyt*npml));
xxt((npml+1):(npml+size(xx,1)),(npml+1):(npml+size(xx,2)))=xx;
yyt((npml+1):(npml+size(xx,1)),(npml+1):(npml+size(xx,2)))=yy;

xx=xxt;
yy=yyt;
%-------------------------------------------------------------------------
if showornot 
    %Plot Cartesian grid
    subplot(1,2,2)
    title('Cartesian grid')
    %plot grid lines
    for i=2:nx+1
        for j=2:ny+1
            line([ksi(i-1,j-1) ksi(i-1,j-1)],[eta(i-1,j-1) eta(i-1,j)]); hold on; %vertical
            line([ksi(i-1,j-1) ksi(i,j-1)],[eta(i-1,j-1) eta(i-1,j-1)]); hold on; %horizontal
        end
    end
    line([ksi(1,1) ksi(end,1)],[max(eta(1,:)) max(eta(1,:))]);
    line([ksi(end,1) ksi(end,1)],[eta(1,1) eta(1,end)]);

    %----------------------------------------------------------
    %Plot curvilinear coordinates
    subplot(1,2,1)
    for i=1:ny+1
        plot(xx(:,i),yy(:,i)); hold on; %horizontal
    end
    for i=1:nx+1
        line([xx(i,1),xx(i,1)],[0,yy(i,end)]); hold on; %vertical
    end
    title('Curvilinear grid')
% 
%     for i=1:nx+1
%         for j=1:ny+1
%             plot(ksi(i,j),eta(i,j),'*'); drawnow; hold on;
%         end
%     end
end
%--------------------------------------------------------------------------
%Derivatives and Jacobian
 J=cell(nx,ny);
%  Ji=cell(nx,ny);
% Ji=zeros(nx,ny);
% J=zeros(nx,ny);
for i=2:nx
     for j=2:ny
        dksi_dx=(ksi(i+1,j)-ksi(i-1,j))/(xx(i+1,j)-xx(i-1,j));
        dksi_dy=(ksi(i,j+1)-ksi(i,j-1))/(yy(i,j+1)-yy(i,j-1));
        deta_dx=(eta(i+1,j)-eta(i-1,j))/(xx(i+1,j)-xx(i-1,j));
        deta_dy=(eta(i,j+1)-eta(i,j-1))/(yy(i,j+1)-yy(i,j-1));
        J{i,j}=[dksi_dy deta_dy; dksi_dx deta_dx];
%         dksi_dx=(ksi(i+1,j)-ksi(i-1,j))/(xx(i+1,j)-xx(i-1,j));
%         dksi_dy=(ksi(i,j+1)-ksi(i,j-1))/(yy(i,j+1)-yy(i,j-1));
%         deta_dx=(eta(i+1,j)-eta(i-1,j))/(xx(i+1,j)-xx(i-1,j));
%         deta_dy=(eta(i,j+1)-eta(i,j-1))/(yy(i,j+1)-yy(i,j-1));
%         J{i-1,j-1}=[dksi_dy deta_dy; dksi_dx deta_dx];
        
%         dx_dksi=(xx(i+1,j)-xx(i-1,j))/(ksi(i+1,j)-ksi(i-1,j));
%         dx_deta=(xx(i,j+1)-xx(i,j-1))/(eta(i,j+1)-eta(i,j-1));
%         dy_dksi=(yy(i+1,j)-yy(i-1,j))/(ksi(i+1,j)-ksi(i-1,j));
%         dy_deta=(yy(i,j+1)-yy(i,j-1))/(eta(i,j+1)-eta(i,j-1));
        
%          J{i-1,j-1}=[dx_deta dy_deta; dx_dksi dy_dksi];
%         Ji{i-1,j-1}=[dksi_dx deta_dx; dksi_dy deta_dy];

%         J(i-1,j-1) =dksi_dx*deta_dy-deta_dx*dksi_dy;
%         Ji(i-1,j-1)=dx_dksi*dy_deta-dx_deta*dy_dksi;
        
    end
end

fprintf('Jacobian cell-array %d x %d created\n',size(J,1),size(J,2));
disp('Finished');
end




