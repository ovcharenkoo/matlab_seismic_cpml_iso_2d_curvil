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
%ksi - curvilinear x
%eta - curvilinear y

%xx - cartesian x
%yy - cartesian y

%To test:
%[ksi,eta, gr_x,gr_y,J, Ji] = func_curv_jacob_pml(30,30,5,0,100, 0, 100,'-(2*pi*x/max(x)+0.25*pi)',10,10,true)

function [ksi,eta,xx,yy,J,Ji] = func_curv_jacob_pml(nx,ny,npml,xmin,xmax,ymin,ymax,argument,dxx,dyy,showornot)

dxxt=xmax/nx;
dyyt=ymax/ny;
ymaxr=ymax-2*dyyt*npml;
ymidr=(ymaxr+ymin)/2;

xmaxr=xmax-2*npml*dxxt;
nxr=nx-2*npml;
nyr=ny-2*npml;

x=linspace(xmin,xmaxr,nxr+1);
%argument of sin
%phi=2*pi*x/max(x);
%phi=-(1.25*pi*x/max(x)+0.25*pi);
phi=inline(argument,'x');

%Bottom part
y=0.2*ymaxr*sin(phi(x));
y=ymidr+y;
%Top part, shifted by pi to be in phase
y1=-ymidr-0.2*ymaxr*sin(phi(x)-pi);

%Phys. domain
ksi=zeros(nxr+1, nyr+1);
eta=zeros(nxr+1, nyr+1);

% %Comp. domain
xx=zeros(nx+1,ny+1);
yy=zeros(nx+1,ny+1);
% dxx=1.0;
% dyy=1.0;

for i=1:nx+1
    for j=1:ny+1
        xx(i,j)=dxx*(i-1);
        yy(i,j)=dyy*(j-1);
    end
end


%Calculate x-ksi in phys domain
for i=1:nxr+1
        ksi(i,:)=x(1,i);
end

%cksi=x;  %x in real space
nyro=round(nyr/2);
for i=1:nyr+1
    if i<=nyro
        ceta=(i-1)*y/nyro; %y in real space
        eta(:,i)=ceta(1,:);
    else
        ceta=ymaxr+(nyr-i+1)*y1/(nyr-nyro); %y in real space
        eta(:,i)=ceta(1,:);
    end
end

xxt=zeros(nx+1,ny+1);
yyt=zeros(nx+1,ny+1);

for i=1:nx+1
    for j=1:ny+1
        xxt(i,j)=dxxt*(i-1);
        yyt(i,j)=dyyt*(j-1);
    end
end

ksi=(ksi+(dxxt*npml));
eta=(eta+(dyyt*npml));
xxt((npml+1):(npml+size(ksi,1)),(npml+1):(npml+size(ksi,2)))=ksi;
yyt((npml+1):(npml+size(ksi,1)),(npml+1):(npml+size(ksi,2)))=eta;

ksi=xxt;
eta=yyt;
%-------------------------------------------------------------------------
if showornot 
    %Plot Cartesian grid
    subplot(1,2,2)
    title('Cartesian grid')
    %plot grid lines
    for i=2:nx+1
        for j=2:ny+1
            line([xx(i-1,j-1) xx(i-1,j-1)],[yy(i-1,j-1) yy(i-1,j)]); hold on; %vertical
            line([xx(i-1,j-1) xx(i,j-1)],[yy(i-1,j-1) yy(i-1,j-1)]); hold on; %horizontal
        end
    end
    line([xx(1,1) xx(end,1)],[max(yy(1,:)) max(yy(1,:))]);
    line([xx(end,1) xx(end,1)],[yy(1,1) yy(1,end)]);

    %----------------------------------------------------------
    %Plot curvilinear coordinates
    subplot(1,2,1)
    for i=1:ny+1
        plot(ksi(:,i),eta(:,i)); hold on; %horizontal
    end
    for i=1:nx+1
        line([ksi(i,1),ksi(i,1)],[0,eta(i,end)]); hold on; %vertical
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
 Ji=cell(nx,ny);
% Ji=zeros(nx,ny);
% J=zeros(nx,ny);
for i=2:nx+1
    for j=2:ny+1
        dksi_dx=(ksi(i,j)-ksi(i-1,j))/(xx(i,j)-xx(i-1,j));
        dksi_dy=(ksi(i,j)-ksi(i,j-1))/(yy(i,j)-yy(i,j-1));
        deta_dx=(eta(i,j)-eta(i-1,j))/(xx(i,j)-xx(i-1,j));
        deta_dy=(eta(i,j)-eta(i,j-1))/(yy(i,j)-yy(i,j-1));
         %J{i-1,j-1}=[dksi_dx dksi_dy; deta_dx deta_dy];
        
        dx_dksi=(xx(i,j)-xx(i-1,j))/(ksi(i,j)-ksi(i-1,j));
        dy_dksi=(yy(i,j)-yy(i,j-1))/(ksi(i,j)-ksi(i,j-1));
        dx_deta=(xx(i,j)-xx(i-1,j))/(eta(i,j)-eta(i-1,j));
        dy_deta=(yy(i,j)-yy(i,j-1))/(eta(i,j)-eta(i,j-1));
        
        J{i-1,j-1}=[dx_dksi dy_dksi; dx_deta dy_deta];
        Ji{i-1,j-1}=[dksi_dx deta_dx; dksi_dy deta_dy];

%         J(i-1,j-1) =dksi_dx*deta_dy-deta_dx*dksi_dy;
%         Ji(i-1,j-1)=dx_dksi*dy_deta-dx_deta*dy_dksi;
        
    end
end
fprintf('Jacobian cell-array %d x %d created\n',size(J,1),size(J,2));
disp('Finished');
end




