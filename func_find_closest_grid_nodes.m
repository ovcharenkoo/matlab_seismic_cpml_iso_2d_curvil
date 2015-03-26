%Function as output gives arrays of markers of nearby grid points that then
%can beb visualized by: 
% for i=1:nx+1
%     for j=1:ny+1
%         if markers(i,j)==1
%             scatter(gr_x(i,j),gr_y(i,j),'r','filled'); drawnow; hold on;
%         end
%     end
% end
% + discretized curve and it's normals
% To plot normals use:
%%line([x1n x2n],[y1n y2n]);
%line([cxmn cxmn-1000*cnorm(1)],[cymn cymn-0.1*max(y_topo)*cnorm(2)],'Color','m'); drawnow; hold on;

function [markers, xt_dis, yt_dis, nvec, xmn, ymn] = func_find_closest_grid_nodes(nx,ny,nc,gr_x,gr_y ,x_topo, y_topo)
fprintf('Started looking for markers, x,y-t_dis, nvec, xmn,ymn with\n');
fprintf('gridsize NX =%d, NY=%d\n', nx, ny);
topo_szx=size(x_topo,2);
fprintf('topo_szx = %.2f - number of topography points\n', topo_szx);
tgrx=topo_szx/nx;
fprintf('tgrx = %.2f - points in one dx spacing |___|\n', tgrx);
nct=round(tgrx/nc);
fprintf('Nc=%d - fine division of dx to nc small samples\n nct = %d - size in points of small topography sample\n',nc, nct);

xt_dis=[]; %x of discretized topography
yt_dis=[]; %y of discretized topography

% nu0x=[];
% nu0y=[];
% 
% nu1x=[];
% nu1y=[];

markers=zeros(nx+1,ny+1);
for i=1:nx
    x_trial=(1+(i-1)*tgrx):i*tgrx; %vector of curve x values
    for j=1:ny
        y_trial=y_topo(x_trial);
        if ~isempty(find(y_trial<gr_y(i,j+1), 1)) && ~isempty(find(y_trial>gr_y(i,j), 1))
            %fprintf('Involved cell ij point: %d %d\n',i,j);
            %set involved grid nodes
            markers(i,j)=1;
            markers(i,j+1)=1;
            markers(i+1,j)=1;
            markers(i+1,j+1)=1;
            %discretize curve inside cells
            yvert=linspace(gr_y(1,j),gr_y(1,j+1),size(x_trial,2));
            for k=1:nc
                xvert=x_topo(1,k*nct+(i-1)*nc*nct)*ones(1,size(x_trial,2));
                %plot(xvert,yvert); drawnow; hold on; %pause(0.02);
                [x0,y0]=curveintersect(xvert,yvert,x_topo, y_topo);
                xt_dis=[xt_dis; x0];
                yt_dis=[yt_dis; y0];
                %scatter(x0,y0,'g','filled'); drawnow; hold on;
            end              

        end  
    end                
end
fprintf('Total number of involved near-boundary grid points = %d  %.2f%%\n', nnz(markers), nnz(markers)*100/(nx*ny));
fprintf('All subarrays created\n');
%Sort xt_dis  and ut_dis vector to get rid of loops on curve
 [xt_dis_sorted, xt_sortindex] = sort(xt_dis);
  yt_dis_sorted = yt_dis(xt_sortindex);
  
 xt_dis=xt_dis_sorted;
 yt_dis=yt_dis_sorted;
 
nvec=[];
xmn=[];
ymn=[];

%subplot(2,2,4);%subplot(2,2,3);
%Find normal vectors
fprintf('Calculating of normals started\n');
for i=2:max(size(xt_dis))
    y1n=yt_dis(i-1); %from discr topography take first point of vector
    y2n=yt_dis(i);
    x1n=xt_dis(i-1);
    x2n=xt_dis(i);

    cxmn=(x1n+x2n)/2; %find x center of current pattern
    cymn=(y1n+y2n)/2; %find y center of current pattern
    
    xmn=[xmn; cxmn]; %x array of middles of the vectors
    ymn=[ymn; cymn]; %y 
    
    x0n=x2n-x1n;
    y0n=y2n-y1n;
    
    cnorm=[y0n, -x0n]/sqrt(x0n^2+y0n^2);
    nvec=[nvec; cnorm];

    %scatter(xmn,ymn,'o','filled');
    %line([x1n x2n],[y1n y2n]);
    %line([cxmn cxmn-1000*cnorm(1)],[cymn cymn-0.1*max(y_topo)*cnorm(2)],'Color','m'); drawnow; hold on;
end
fprintf('Calculating of normals finished\n');
end