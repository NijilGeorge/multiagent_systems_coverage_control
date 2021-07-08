clear;
clc;
close all;
global X
global Y
global M
global cost_list
global sensing %changed inside cost_cell
global binary_obst
global xD_list
global yD_list
global Target_coordinates
global Current_coordinates
global init_heading
global init_heading_list
global a_mat
global clambda_mat
global lambda_mat
global kappa
global phi
global alpha_rate
global a_min
global Gamma_mat
global gamma
global avg_para_err
global a_true
sensing_list=zeros(size(M));
showPlot = true;
numIterations  = 100;
xrange = 10;  %region size
yrange = 10;
n = 15; %number of robots  (changing the number of robots is interesting)
[X,Y] = meshgrid(0:.1:max(xrange,yrange));
% M =gauss_multiple(X, Y, [5 5],[3 3;9 9],[3 0]);
% M1 =gauss_multiple(X, Y, [2 5],[3 3;9 9],[2 0]);
% M=M-M1;
M=ones(size(X));
Target_coordinates=[1 1];
Current_coordinates=[0 0];
init_heading=0;
init_heading_list=zeros(n,1);
sim('single_integrator.slx');
cost_list=zeros(numIterations,1);
% Generate and Place  n stationary robots
Px = 0.01*mod(1:n,ceil(sqrt(n)))'*xrange; %start the robots in a small grid
Py = 0.01*floor((1:n)/sqrt(n))'*yrange;
% rand('seed', 2);
% xy=rand(n,2);
% Px=round(xy(:,1)*xrange,1);
% Py=round(xy(:,2)*yrange,1);
% crs = [ 0, 0;    
%     0, yrange;
%     xrange, yrange;
%     xrange, 0];

% crs = [ 0, 0;    
%     0, yrange;
%     1/3*xrange, yrange;  % a world with a narrow passage
%     1/3*xrange, 3/6*yrange;
%     2/3*xrange, 3/6*yrange;
%     2/3*xrange, yrange;
%     xrange, yrange;
%     xrange, 0];

% crs = [ 0, 0;    
%     0, yrange;
%     1/3*xrange, yrange;  % a world with a narrow passage
%     1/3*xrange, 3/6*yrange;
%     2/3*xrange, 3/6*yrange;
%     xrange, 3/6*yrange;
%     xrange, 0];

% crs = [ 0, 0;    
%     0, yrange;
%     0, 1/3*yrange;  % a world with a narrow passage
%     2/3*xrange, 1/3*yrange;
%     2/3*xrange, 0.5*yrange;
%     1/3*xrange, 0.5*yrange;
%     1/3*xrange, yrange;
%     xrange, yrange;
%     xrange, 0];
crs = [ 0, 0;    
    0, yrange;
    0, 1/3*yrange;  % a world with a narrow passage
    1/3*xrange, 1/3*yrange;
    1/3*xrange, 0.5*yrange;
    0.25*xrange, 0.5*yrange;
    0.25*xrange, yrange;
    xrange, yrange;
    xrange, 0];

% surf(X,Y,ones(size(valid_reg)).*valid_reg)
valid_reg=valid_region(crs);
M=M.*valid_reg;
% M=M/sum(M,'all');
m = 9;
avg_para_err=zeros(numIterations,1);
clambda=zeros(m+1,m+1);
lambda=zeros(m+1,1);
alpha_rate=1;
Gamma_mat=eye(m+1)*1e-2;
gamma=20%2e7;%20
sqrtm=sqrt(m);
max_dim=max(xrange,yrange);
%----------Defining a_true---------------
a_min=0.1
a_true=ones(m+1,1)*a_min;
a_true(1)=1000;
a_true(m)=1000;
% a_true(m+1)=1000;
%----------Defining initial values for a_is,clambda_i,lambda_i---------------
a_mat=ones(m+1,n)*a_min;
clambda_mat=zeros(m+1,m+1,n);
lambda_mat=zeros(m+1,n);
%-----------Defining Kappa---------------
kappa=zeros([size(M) m]);
kappa=cat(3,kappa,ones(size(kappa(:,:,1)))/sum(ones(size(kappa(:,:,1))),'all'));%adding one more basis vector(constant through out) /sum(ones(size(kappa(:,:,1))),'all')

[xx,yy]=meshgrid([max_dim/(2*sqrtm):max_dim/sqrtm:(2*sqrtm*max_dim-max_dim)/(2*sqrtm)]); 
    %mesh for centers of basis functions
xx=reshape(xx,[1 m]);
yy=reshape(yy,[1 m]);
M=zeros(size(X));
for j=1:m
    kappa(:,:,j)=gauss_multiple(X, Y, [1.8],[xx(j) yy(j)],[1]);%x,y,sigma,centers,weights
    %kappa(:,:,j)=kappa(:,:,j)/sum(kappa(:,:,j),'all');
    M=M+kappa(:,:,j)*a_true(j);
end
%----------Defining density function for visualization--------------
M=M+kappa(:,:,end)*a_true(m+1);
phi_i=sum(kappa,3)*a_min;
phi=ones([size(X),n]).*phi_i;
%comment out these
M=ones(size(X));
M=M.*valid_reg;
M=M/sum(M,'all');
% M=M/sum(M,'all');
%----------
binary_obst=ones(size(valid_reg))-valid_reg.*ones(size(valid_reg));
binary_obst(:,1)=ones(size(binary_obst(:,1)));
binary_obst(:,end)=ones(size(binary_obst(:,end)));
binary_obst(1,:)=ones(size(binary_obst(1,:)));
binary_obst(end,:)=ones(size(binary_obst(end,:)));
binary_obst=bwdist(binary_obst,'euclidean');
% binary_obst=binary_obst/sum(binary_obst,'all');
for i = 1:numel(Px)  
    while ~inpolygon(Px(i),Py(i),crs(:,1),crs(:,2))% ensure robots are inside the boundary
        Px(i) = rand(1,1)*xrange; 
        Py(i) = rand(1,1)*yrange;
    end
end


[Px, Py]=lloydsAlgorithm(Px,Py, crs, numIterations, showPlot);
% plot(cost_list)
writematrix(a_mat,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/a_mat.csv')


function [Px, Py] = lloydsAlgorithm(Px,Py, crs, numIterations, showPlot)
% LLOYDSALGORITHM runs Lloyd's algorithm on the particles at xy positions 
% (Px,Py) within the boundary polygon crs for numIterations iterations
% showPlot = true will display the results graphically.  
% Run with no input to see example.  To initialize a square with 50 robots 
% in left middle, run:
%lloydsAlgorithm(0.01*rand(50,1),zeros(50,1)+1/2, [0,0;0,1;1,1;1,0], 200, true)

close all
format compact
global cost_list
global sensing %changed inside cost_cell
global M
global X
global Y
global yD_list
global xD_list
global Target_coordinates
global Current_coordinates
global init_heading
global init_heading_list
global phi
global avg_para_err
% initialize random generator in repeatable fashion

if nargin < 1   % demo mode
    showPlot = true;
    numIterations  = 200;
    xrange = 10;  %region size
    yrange = 10;
    n = 20; %number of robots  (changing the number of robots is interesting)

% Generate and Place  n stationary robots
    Px = 0.01*mod(1:n,ceil(sqrt(n)))'*xrange; %start the robots in a small grid
    Py = 0.01*floor((1:n)/sqrt(n))'*yrange;
    
%     Px = 0.1*rand(n,1)*xrange; % place n  robots randomly
%     Py = 0.1*rand(n,1)*yrange;
%     
%     crs = [ 0, 0;    
%         0, yrange;
%         1/3*xrange, yrange;  % a world with a narrow passage
%         1/3*xrange, 1/4*yrange;
%         2/3*xrange, 1/4*yrange;
%         2/3*xrange, yrange;
%         xrange, yrange;
%         xrange, 0];
    crs = [ 0, 0;    
        0, yrange;
        xrange, yrange;
        xrange, 0];
    
    for i = 1:numel(Px)  
        while ~inpolygon(Px(i),Py(i),crs(:,1),crs(:,2))% ensure robots are inside the boundary
            Px(i) = rand(1,1)*xrange; 
            Py(i) = rand(1,1)*yrange;
        end
    end
else
    xrange = max(crs(:,1));
    yrange = max(crs(:,2));
    n = numel(Px); %number of robots  
end
Prev_Px=Px;
Prev_Py=Py;

%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showPlot
    verCellHandle = zeros(n,1);
    cellColors = cool(n);
    for i = 1:numel(Px) % color according to
        verCellHandle(i)  = patch(Px(i),Py(i),cellColors(i,:)); % use color i  -- no robot assigned yet
        hold on
    end
    pathHandle = zeros(n,1);    
    numHandle = zeros(n,1);    
    for i = 1:numel(Px) % color according to
       pathHandle(i)  = plot(Px(i),Py(i),'-','color',cellColors(i,:)*.8);
       numHandle(i) = text(Px(i),Py(i),num2str(i));
    end
    goalHandle = plot(Px,Py,'+','linewidth',2);
    currHandle = plot(Px,Py,'o','linewidth',2);
    titleHandle = title(['o = Robots, + = Goals, Iteration ', num2str(0)]);
end
%%%%%%%%%%%%%%%%%%%%%%%% END VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteratively Apply LLYOD's Algorithm
for counter = 1:numIterations
    sensing=zeros(size(M));
    avg_para_err(counter)=find_avg_parameter_err(numel(Px));
    %[v,c]=VoronoiLimit(Px,Py, crs, false);
    [v,c]=VoronoiBounded(Px,Py, crs);
    
    if showPlot
        set(currHandle,'XData',Px,'YData',Py);%plot current position
        for i = 1:numel(Px) % color according to
            
            Current_coordinates=[Prev_Px(i); Prev_Py(i)];
            Target_coordinates=[Px(i);Py(i)];
            init_heading=init_heading_list(i);
            if  Current_coordinates ~= Target_coordinates
%             if 0 
                if (i<=numel(Px)/3)
%                     sim('single_integrator.slx');
                    sim('unicycle_polar_on_spot_simulink.slx');
                    path_from_model=ans.path_coordinates.data;
                    b=ans.heading.data;
                    init_heading_list(i)=b(end);
%                     sim('single_integrator.slx');
                    update_lambda(path_from_model,i);
                elseif(i<=2*numel(Px)/3)
                    sim('double_integrator_pd.slx');
%                     sim('unicycle_simulink_const_input.slx');
                    path_from_model=ans.path_coordinates.data;
%                     b=ans.heading.data;
%                     init_heading_list(i)=b(end); 
                    update_lambda(path_from_model,i);
                else
                    sim('single_integrator.slx');
%                     sim('unicycle_polar_on_spot_simulink.slx');
                    path_from_model=ans.path_coordinates.data;
%                     b=ans.heading.data;
%                     init_heading_list(i)=b(end); 
                    update_lambda(path_from_model,i);
                end
                
%             wait('unicycle_polar_on_spot_simulink.slx','stateOfJob',timeout);
                xD = [get(pathHandle(i),'XData'),path_from_model(:,1)'];
                yD = [get(pathHandle(i),'YData'),path_from_model(:,2)'];
                xD_list=xD;
                yD_list=yD;
                set(pathHandle(i),'XData',xD,'YData',yD);%plot path position
            else
                xD = [get(pathHandle(i),'XData'),Px(i)];
                yD = [get(pathHandle(i),'YData'),Py(i)];
                set(pathHandle(i),'XData',xD,'YData',yD);
            end
           set(numHandle(i),'Position',[ Px(i),Py(i)]);
        end 
    end
    costfn=0;
    for i = 1:numel(c) %calculate the centroid of each cell 
        phi(:,:,i)=find_phi_i(i);
        [cx,cy] = WtPolyCentroid(v(c{i},1),v(c{i},2),counter,i);
        costfn=costfn+cost_cell(v(c{i},1),v(c{i},2),cx,cy);
        cx = min(xrange,max(0, cx));
        cy = min(yrange,max(0, cy));
        update_a(v(c{i},1),v(c{i},2),i);
        if ~isnan(cx) && inpolygon(cx,cy,crs(:,1),crs(:,2))
            Prev_Px(i)=Px(i);
            Prev_Py(i)=Py(i);
            Px(i) = cx;  %don't update if goal is outside the polygon
            Py(i) = cy;
        end
    end
    cost_list(counter)=costfn;
    if showPlot
        
        f1=figure(1);
        for i = 1:numel(c) % update Voronoi cells
            set(verCellHandle(i), 'XData',v(c{i},1),'YData',v(c{i},2));
        end

        set(titleHandle,'string',['o = Robots, + = Goals, Iteration ', num2str(counter,'%3d'),' of ',num2str(numIterations,'%3d')]);
        set(goalHandle,'XData',Px,'YData',Py);%plot goal position
        
        
        drawnow
        filename1 = '/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/path.gif';
        savegif(f1,filename1,counter);
        if counter==1
            saveas(f1,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/path_initial.jpg');
        end
        f2=figure(2);
        plot(cost_list);
        text(counter,cost_list(counter),'.');
        text(counter,cost_list(counter),num2str(cost_list(counter)));
        xlabel('Iteration');
        ylabel('Cost Function');
        title('Cost Function v/s Iteration');
        filename3 = '/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/costvsiter.gif';
        savegif(f2,filename3,counter);
        
        f3=figure(3);
        surf(X,Y,sensing,'EdgeColor', [0 0 0]);
        view(2);
        xlabel('x');
        ylabel('y');
        zlabel('Cost function');
        title(['Cost Function Iteration ', num2str(counter,'%3d'),' of ',num2str(numIterations,'%3d')]);
        colorbar
        filename2 = '/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/cost_function.gif';
        savegif(f3,filename2,counter);
         if counter==1
            saveas(f3,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/cost_initial.jpg');
        end
%         if mod(counter,50) ==0
%             pause
%             %pause(0.1)
%         end
        f4=figure(4);
        surf(X,Y,M);
        xlabel('x');
        ylabel('y');
        zlabel('Density function');
        title('Density Function');
        colorbar;
    end
    %find cost function value here and store\
    
end
saveas(f1,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/path_final.jpg');
saveas(f2,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/cost_vs_iter.jpg');
saveas(f4,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/density_3d.jpg');
f3=figure(3);
surf(X,Y,sensing,'EdgeColor', [0 0 0]);
view(2);
xlabel('x');
ylabel('y');
zlabel('Cost function');
title(['Cost Function Iteration ', num2str(counter,'%3d'),' of ',num2str(numIterations,'%3d')]);
colorbar
saveas(f3,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/cost_2d.jpg');
f5=figure(5);
surf(X,Y,M);
view(2);
xlabel('x');
ylabel('y');
zlabel('Density function');
title('Density Function');
colorbar;
saveas(f5,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/density_2d.jpg');
writematrix(cost_list,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/cost_vs_iter.csv')
f7=figure(7);
plot(avg_para_err);
xlabel('Iteration');
ylabel('Average Parameter Error');
title('Average Parameter Error');
saveas(f7,'/home/nijil/Documents/Nijil/MTP/Report/Simultaion Figures/adaptive_heterogeneous/avg_para_err.jpg');
end

function cost_i=cost_cell(Xa,Ya,cx,cy) 
    %cost incurred by the agent at position (cx,cy) with voronoi region
    %with vertices given by Xa,Ya
    global X;
    global Y;
    global M;
    global sensing;
    Xn=X(:);
    Yn=Y(:);
    Xa = [Xa(2:end);Xa(1)];
    Ya = [Ya(2:end);Ya(1)];
    in=inpolygon(Xn,Yn,Xa,Ya); %getting points that are inside the given polygon
    in=reshape(in,size(M));
    costx=(X-ones(size(X)).*cx);
    costy=(Y-ones(size(Y)).*cy);
    cost_i=sum(((costx.*costx+costy.*costy).*M).*in,'all')/2;
    sensing=sensing+((costx.*costx+costy.*costy).*M).*in;
end

function [Cx,Cy] = WtPolyCentroid(Xa,Ya,iter,i)
    global X;
    global Y;
    global M;
    global binary_obst;
    global phi
    M1=M+0*(exp(-0.5*iter))*binary_obst;
%     M1=phi(:,:,i)+0*(exp(-0.5*iter))*binary_obst;
    Xn=X(:);
    Yn=Y(:);
    Xa = [Xa(2:end);Xa(1)];
    Ya = [Ya(2:end);Ya(1)];
    in=inpolygon(Xn,Yn,Xa,Ya); %getting points that are inside the given polygon
    in=reshape(in,size(X));
    Cx=sum(M1.*X.*in,'all')/sum(M1.*in,'all');
    Cy=sum(M1.*Y.*in,'all')/sum(M1.*in,'all');
end

function val = gauss_multiple(x, y, sigma, centers,w)
%x,y=mesh sigma center w=>weight
    nc=size(centers);
    Cx=centers(:,1);
    Cy=centers(:,2);
    val=zeros(size(x));
    for ii=1:nc(1)
        xc = Cx(ii);
        yc = Cy(ii);
        exponent = ((x-xc).^2 + (y-yc).^2)/(2*sigma(ii)*sigma(ii));
        val       =val+w(ii)*(exp(-exponent))/(sqrt(2*pi)*sigma(ii));
    end
end

function [V,C]=VoronoiBounded(x,y, crs)
    % VORONOIBOUNDED computes the Voronoi cells about the points (x,y) inside
    % the bounding box (a polygon) crs.  If crs is not supplied, an
    % axis-aligned box containing (x,y) is used.

    bnd=[min(x) max(x) min(y) max(y)]; %data bounds
    if nargin < 3
        crs=double([bnd(1) bnd(4);bnd(2) bnd(4);bnd(2) bnd(3);bnd(1) bnd(3);bnd(1) bnd(4)]);
    end

    rgx = max(crs(:,1))-min(crs(:,1));
    rgy = max(crs(:,2))-min(crs(:,2));
    rg = max(rgx,rgy);
    midx = (max(crs(:,1))+min(crs(:,1)))/2;
    midy = (max(crs(:,2))+min(crs(:,2)))/2;

    % add 4 additional edges
    xA = [x; midx + [0;0;-5*rg;+5*rg]];
    yA = [y; midy + [-5*rg;+5*rg;0;0]];

    [vi,ci]=voronoin([xA,yA]);

    % remove the last 4 cells
    C = ci(1:end-4);
    V = vi;
    % use Polybool to crop the cells
    %Polybool for restriction of polygons to domain.

    for ij=1:length(C)
            % thanks to http://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
            % first convert the contour coordinate to clockwise order:
            [X2, Y2] = poly2cw(V(C{ij},1),V(C{ij},2));
            [xb, yb] = polybool('intersection',crs(:,1),crs(:,2),X2,Y2);
            ix=nan(1,length(xb));
            for il=1:length(xb)
                if any(V(:,1)==xb(il)) && any(V(:,2)==yb(il))
                    ix1=find(V(:,1)==xb(il));
                    ix2=find(V(:,2)==yb(il));
                    for ib=1:length(ix1)
                        if any(ix1(ib)==ix2)
                            ix(il)=ix1(ib);
                        end
                    end
                    if isnan(ix(il))==1
                        lv=length(V);
                        V(lv+1,1)=xb(il);
                        V(lv+1,2)=yb(il);
                        ix(il)=lv+1;
                    end
                else
                    lv=length(V);
                    V(lv+1,1)=xb(il);
                    V(lv+1,2)=yb(il);
                    ix(il)=lv+1;
                end
            end
            C{ij}=ix;
    end
end

function returning=savegif(figure,filename,n)
    h = figure;
    axis tight manual; % this ensures that getframe() returns a consistent size
    %filename = '/home/nijil/Documents/Nijil/Matlab_simulations/testAnimated.gif';
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.5); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.5); 
    end
    returning=1;
end

function vr=valid_region(crs)
    global X;
    global Y;
    Xn=X(:);
    Yn=Y(:);
    c=[crs;crs(1,:)];
    vr=inpolygon(Xn,Yn,c(:,1),c(:,2));
    vr=reshape(vr,size(X));
end

function phi_i=find_phi_i(iter)
    global kappa;
    global a_mat
    phi_i=zeros(size(kappa(:,:,1)));
    for jj=1:size(a_mat(:,iter))
        phi_i=phi_i+kappa(:,:,jj)*a_mat(jj,iter);
    end
end

function [clambda_i,lambda_i]=update_lambda(path,iter)
%updates the clambda and lambda for ith agent, also updates in the global
%matrices
global X;
global Y;
global kappa;
global clambda_mat;
global lambda_mat;
global alpha_rate;
global M;
clambda_i=clambda_mat(:,:,iter);
lambda_i=lambda_mat(:,iter);
path_length=size(path(:,1));
path_length=path_length(1);
    for j=1:path_length
        xcord=find(round(X(1,:),1)-round(path(j,1),1)==0);
%         isempty(xcord)
        ycord=find(round(Y(:,1),1)-round(path(j,2),1)==0);
%         isempty(ycord)
        if not(or(isempty(xcord),isempty(ycord)))
            clambda_i=(1-alpha_rate)*clambda_i+squeeze(kappa(xcord,ycord,:))*squeeze(kappa(xcord,ycord,:))';
            lambda_i=(1-alpha_rate)*lambda_i+squeeze(kappa(xcord,ycord,:))*M(xcord,ycord);
        end
    %     templist(j)=lambda;
    end
clambda_mat(:,:,iter)=clambda_i;
lambda_mat(:,iter)=lambda_i;
end

function a_i=update_a(Xa,Ya,iter)
    global X;
    global Y;
    global Current_coordinates;
    global Target_coordinates;
    global kappa;
    global a_mat;
    global a_min;
    global Gamma_mat;
    global clambda_mat;
    global lambda_mat;
    global gamma;
    a_i=a_mat(:,iter);
    Xn=X(:);
    Yn=Y(:);
    Xa = [Xa(2:end);Xa(1)];
    Ya = [Ya(2:end);Ya(1)];
    in=inpolygon(Xn,Yn,Xa,Ya); %getting points that are inside the given polygon
    in=reshape(in,size(X));
    F=[sum(kappa.*(round(X,1)-round(Current_coordinates(1),1)).*in,[1 2]),sum(kappa.*(round(Y,1)-round(Current_coordinates(2),1)).*in,[1 2])];
    F=squeeze(F)';
    b=-F*(Target_coordinates-Current_coordinates)-gamma*(clambda_mat(:,:,iter)*a_i-lambda_mat(:,iter));
    for jj=1:size(b)
        if not(or((a_i(jj)>a_min),and(a_i(jj)==a_min,b(jj)>=0)))
            b(jj)=0;
        end
    end
    a_i=a_i+Gamma_mat*b
    a_mat(:,iter)=a_i;
end

function avg_err=find_avg_parameter_err(n)
global a_true
global a_mat
avg_err=sum(vecnorm(a_true-a_mat))/n;
end