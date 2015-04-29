function oper=create_oper_for_derivatives_from_taylor(dx,dz)
ctr=0;

% dx=1.d0;
% dz=1.d0;

for i=1:-1:-1
    for j=1:-1:-1
        r=i;
        s=j;
        ctr=ctr+1;
        A=[1.d0 r*dx s*dz (r*dx)^2.d0/2.d0 (s*dz)^2.d0/2.d0 r*s*dx*dz; ...
            0 1.d0 0 r*dx 0 s*dz; ...
            0 0 1.d0 0 s*dz r*dx; ...
            0 0 0 1.d0 0 0; ...
            0 0 0 0 1.d0 0; ...
            0 0 0 0 0 1.d0];
        fprintf('i=%d j=%d\n',i,j);
        biOp(ctr,:)=A(1,:);
    end
end
biOp
oper=svdinv(biOp)
% 
% subplot(1,2,1);
% pcolor(biOp);
% colorbar();
% title('biOp');
% 
% subplot(1,2,2);
% pcolor(oper);
% colorbar();
% title('oper');
end