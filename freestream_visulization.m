clear;

%% Full State boundary
filename = 'fullstate_test_Res.dat';
f = fopen(filename);
Residual = fscanf(f,'%f');
fclose(f);

figure;
semilogy(abs(Residual),'LineWidth',3);
xlabel('Number of iteration');
ylabel('Residual');
title('Full State boundary condition');
ax = gca;
ax.FontSize = 20;

filename = 'fullstate_test_Force.dat';
f = fopen(filename);
Force = fscanf(f,'%f');
fclose(f);
total = max(size(Force));
Force = reshape(Force,[8,total/8]);
Force = Force';

figure;
for i = 1:4
    semilogy(abs(Force(:,2*i-1)),'LineWidth',3,'DisplayName',['B' num2str(i)]);
    hold on;
end
legend;
xlabel('Number of iteration');
ylabel('F_x (N)');
title('Full State boundary condition');
ax = gca;
ax.FontSize = 20;

figure;
for i = 1:4
    semilogy(abs(Force(:,2*i)),'LineWidth',3,'DisplayName',['B' num2str(i)]);
    hold on;
end
legend;
xlabel('Number of iteration');
ylabel('F_y (N)');
title('Full State boundary condition');
ax = gca;
ax.FontSize = 20;


%% Wall boundary

filename = 'wall_test_Res.dat';
f = fopen(filename);
Residual = fscanf(f,'%f');
fclose(f);

figure;
semilogy(abs(Residual),'LineWidth',3);
xlabel('Number of iteration');
ylabel('Residual');
title('Wall boundary condition');
ax = gca;
ax.FontSize = 20;

filename = 'wall_test_Force.dat';
f = fopen(filename);
Force = fscanf(f,'%f');
fclose(f);
total = max(size(Force));
Force = reshape(Force,[8,total/8]);
Force = Force';

figure;
for i = 1:4
    semilogy(abs(Force(:,2*i-1)),'LineWidth',3,'DisplayName',['B' num2str(i)]);
    hold on;
end
legend;
xlabel('Number of iteration');
ylabel('F_x (N)');
title('Wall boundary condition');
ax = gca;
ax.FontSize = 20;

figure;
for i = 1:4
    semilogy(abs(Force(:,2*i)),'LineWidth',3,'DisplayName',['B' num2str(i)]);
    hold on;
end
legend;
xlabel('Number of iteration');
ylabel('F_y (N)');
title('Wall boundary condition');
ax = gca;
ax.FontSize = 20;

