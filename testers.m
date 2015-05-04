%%
%Check J and Jd
close all;
clc;
for i=2:NX
    for j=2:NY
        subplot(1,2,1);
        pcolor(J{i,j});
        J{i,j}
        colorbar();
        title(['J. i=' num2str(i) ' j=' num2str(j)]);
        subplot(1,2,2);
        pcolor(Jd{i,j});
        Jd{i,j}
        colorbar();
        title(['Jd. i=' num2str(i) ' j=' num2str(j)]);
        input('Next?');
    end
end
%%
%Visualy compare determinant destribution of J and Jd
close all;
clc;

for i=2:NX-1
    for j=2:NY-1
       Jm(i,j)=det(J{i,j});
       Jdm(i,j)=det(J2{i,j});
    end
end
subplot(1,2,1)
imagesc(Jm');
colorbar();
title('det J');
subplot(1,2,2)
imagesc(Jdm');
colorbar();
title('det Jd');
%%