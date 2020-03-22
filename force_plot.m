filename = 'city0_Force.dat';
f = fopen(filename);
Force_0 = fscanf(f,'%f');
fclose(f);
total = max(size(Force_0));
Force_0 = reshape(Force_0,[8,total/8]);
Force_0 = Force_0';

filename = 'city1_Force.dat';
f = fopen(filename);
Force_1 = fscanf(f,'%f');
fclose(f);
total = max(size(Force_1));
Force_1 = reshape(Force_1,[8,total/8]);
Force_1 = Force_1';

filename = 'city2_Force.dat';
f = fopen(filename);
Force_2 = fscanf(f,'%f');
fclose(f);
total = max(size(Force_2));
Force_2 = reshape(Force_2,[8,total/8]);
Force_2 = Force_2';

for i = 1:4
    figure;  % Fx
    Nt = max(size(Force_0));
    t = linspace(0,2,Nt);
    plot(t, Force_0(:,2*i-1),'LineWidth',3,'DisplayName','city0');
    hold on;
    Nt = max(size(Force_1));
    t = linspace(0,2,Nt);    
    plot(t, Force_1(:,2*i-1),'LineWidth',3,'DisplayName','city1');
    Nt = max(size(Force_2));
    t = linspace(0,2,Nt);    
    plot(t, Force_2(:,2*i-1),'LineWidth',3,'DisplayName','city2');
    legend;
    xlabel('T');
    ylabel('F_x (N)');
    title(['B' num2str(i)]);
    ax = gca;
    ax.FontSize = 20;
    
    figure;  % Fy
    Nt = max(size(Force_0));
    t = linspace(0,2,Nt);
    plot(t, Force_0(:,2*i),'LineWidth',3,'DisplayName','city0');
    hold on;
    Nt = max(size(Force_1));
    t = linspace(0,2,Nt);    
    plot(t, Force_1(:,2*i),'LineWidth',3,'DisplayName','city1');
    Nt = max(size(Force_2));
    t = linspace(0,2,Nt);    
    plot(t, Force_2(:,2*i),'LineWidth',3,'DisplayName','city2');
    legend;
    xlabel('T');
    ylabel('F_y (N)');
    title(['B' num2str(i)]);
    ax = gca;
    ax.FontSize = 20;    
end


