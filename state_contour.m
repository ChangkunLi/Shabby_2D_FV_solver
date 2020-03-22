clear;

screen = get(0,'ScreenSize');
W_g = screen(3); H_g = screen(4);

for i = 0:0
    mesh = read_gri(['city' num2str(i) '.gri']);
    
    Elem = mesh.Elem;
    Node = mesh.Node;
    X = zeros(3,mesh.nElem);
    Y = zeros(3,mesh.nElem);
    for ii = 1:mesh.nElem
        for jj = 1:3
            nj = Elem(ii,jj);
            X(jj,ii) = Node(nj,1);
            Y(jj,ii) = Node(nj,2);
        end
    end
    
    for j = 0:5
        filename = ['city' num2str(i) '_State_' num2str(j) '.dat'];
        f = fopen(filename);
        U = fscanf(f,'%f');
        fclose(f);
        total = max(size(U));
        U = reshape(U,[total/3,3]);
        U = U';     
        
        h = U(1,:);
        
        w_g = 0.5*W_g;              % width of the window
        h_g = 0.5*H_g;              % hight of the window

        figure('Color',[1 1 1],'Position',[0,0,w_g,h_g]);
        patch(X,Y,h);    
        axis equal;
        axis tight;
        xlabel('X');
        ylabel('Y');
        title(['city' num2str(i) ', t = ' num2str(j*0.05)]);
        ax = gca;
        ax.FontSize = 20;
        colorbar;
        colormap jet;
        if j == 0
            caxis([0.7 1.3]);
        else
            caxis([0.9 1.1]);
        end
    end
end