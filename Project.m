clc
clear
close all
warning off
data = readtable("BBD_2021_VERTICAL.csv");


%% true time in juliandate
time = table2array(data(1,8:end));
time = time';
year = round(time/10000);
month = round((time - year*10000)/100);
day = time - year*10000 - month*100;

w = sum(numel(find(year == 2017)))/2;


Point_data = [];


for ii = 1:size(data,1)-1
    t = juliandate([year month day]);
    p = data(ii+1,:);
    xyz = table2array(p(1,4:6));
    y = table2array(p(1,8:end));
    k = find(isnan(y) == 0);
    y = y(k);
    t = t(k);
    t = t-min(t);
    n = numel(y);
    y = y';

            figure()
            plot(t,y,'.')
    Point_data(ii).xyz = xyz;
    Point_data(ii).obs = y;
    Point_data(ii).time = t;



    %% eliminate Linear Part

    A = [ones(n,1) t];
    x = (A'*A)^-1*(A'*y);
    a = x(2);
    b1 = x(1);
    y_linear = a*t + b1;
    
            hold on
            plot(t,y_linear,'.')

    yn = y - y_linear;
    Point_data(ii).linearX = x;
    Point_data(ii).DelLinearY = yn;
    yn_2 = y;

    %% LSSA
    T = 0.1:0.1:w/2;
    Q = 0.25*eye(n);
    W = Q^-1;
    df = n - 2;
    repetition = 1;
    yn_1 = yn;
    %
        [spectrum, chats, CritVal, g, xhat1, COVcoeff, nc]=...
                         LSSA(t, y, W, T, [0,0,0,0,0], [], 0.01);


    while 1
        for i = 1:length(T)

            A = [cos(2*pi*t/T(i))  sin(2*pi*t/T(i))];
            xhat = (A'*W*A)^-1*A'*W*yn;
            S(i) = (yn'*W*A*xhat)/(yn'*W*yn);
        end

                        figure;	plot(T,S,'r');hold on;grid
                        title(texlabel('Spectrum computation'),'FontSize',14,'FontName','Arial','FontWeight','bold');
                        xlabel('T (day)','FontSize',12,'FontName','Arial');
                        hold on
%         statistical testing
        if max(S)<=(1 + 0.5*df*finv(0.01,df,2) )^-1
            disp('H0 Accept')
            break
        else
            disp('H0 reject')
            j = find(S==max(S));
            j = j(1);
            TT(repetition) = T(j);
            A = [cos(2*pi*t/T(j))  sin(2*pi*t/T(j))];
            xhat = (A'*W*A)^-1*A'*W*yn;
            figure
                        plot(T(j),S(j),'*b')
            yn = yn - A*xhat;
            yn_2 = yn_2 - A*xhat;
            repetition = repetition + 1;
        end

    end

    Point_data(ii).w = TT;


    %
    A = zeros(n,length(TT)*2);
    for i = 1:length(TT)
        A(:,2*i-1:2*i) = [cos(2*pi*t/TT(i))  sin(2*pi*t/TT(i))];
    end
    xhat = inv(A'*W*A)*A'*W*yn_1;
    Qxhat =  inv(A'*W*A);
    sigma_xhat = sqrt(diag(Qxhat));

    Point_data(ii).furie = xhat;
    figure;plot(t,yn_1-A*xhat,'.r');hold on;grid;plot(t,yn_1-A*xhat,'b')
    title(texlabel('e_h_a_t & y_h_a_t'),'FontSize',14,'FontName','Arial','FontWeight','bold');
    xlabel('t (year)','FontSize',12,'FontName','Arial');
    e = yn_1-A*xhat;

            figure
            plot(t,e,'.y');
            close all


    A = [ones(n,1) t];
    x2 = (A'*A)^-1*(A'*yn_2);
    a = x2(2);
    b1 = x2(1);
    y_linear2 = a*t + b1;
    coordinate(ii,:) = [xyz x' x2'];




end
%%
x = coordinate(:,1);
y = coordinate(:,2);
z = coordinate(:,3);
figure

b1 = coordinate(:,end)*365;
% triangulate and plot
for i = 1:0.5:20
    z1 = z + b1*i;
    tri = delaunay(x, y);
    trisurf(tri, x, y, z1);
    view(90-285,50 )
    colorbar
    zlim([-100 300])
%     clim([-100 125])
    % optional, could help make the plot look nicer
    shading interp
    pause(0.1)

end
%%

figure()
subplot(2,2,[1,2])
scatter(x,y,70,b1,'filled')
colorbar
colormap turbo
title('Estimated')
xlabel('Easting')
ylabel('Northing')
grid on
axis equal

b2 = table2array(data(2:end,2));
subplot(2,2,3)
scatter(x,y,70,b2,'filled')
colorbar
colormap turbo
title('Real')
xlabel('Easting')
ylabel('Northing')
grid on
axis equal

subplot(2,2,4)
scatter(x,y,70,b2-b1,'filled')
colorbar
colormap turbo
title('Estimated - Real')
xlabel('Easting')
ylabel('Northing')
grid on
axis equal

%%
figure('Renderer' ,'painters','Position',[300 100 1500 500])

subplot(121)
z1 = + b1(:,end)*0;
tri = delaunay(x, y);
trisurf(tri, x, y, z1);
view(45,50 )
colorbar
zlim([-600 0])
clim([-370 0])
% optional, could help make the plot look nicer
shading interp
title('First year of zero heght hypothese')
xlabel('Easting')
ylabel('Northing')
zlabel('Zenith')
subplot(122)
% triangulate and plot

z1 = + b1(:,end)*40;
tri = delaunay(x, y);
trisurf(tri, x, y, z1);
view(45,50 )
colorbar
zlim([-600 0])
clim([-370 0])
% optional, could help make the plot look nicer
shading interp

title('40 years later of zero heght hypothese')
xlabel('Easting')
ylabel('Northing')
zlabel('Zenith')


saveas(gcf,'40.png')