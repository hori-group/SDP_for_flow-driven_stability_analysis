%
% Fig. 2A
%
close all;
clear;

a = 0.01:0.002:0.25;
b = 0:0.002:0.07;
v_norm = 0;
d = 6;

res = zeros(length(a),length(b));
sol1 = zeros(length(a),length(b));

filename_head = sprintf('Figure2A');
filename_mat = strcat(filename_head,'.mat');
filename_fig = strcat(filename_head,'.fig');
filename_png = strcat(filename_head,'.png');

zeta = sym('zeta','real');
syms a11 a12 a21 a22 d1 d2 v1 v2 real;

A = [a11 a12; a21 a22];
D = [d1 0; 0 d2];
V = [v1 0; 0 v2];

s = sym('s');

p = det(1i*s*eye(2) - (A-zeta^2*D + 1i*zeta*V));

p_coeff = coeffs(p,s);

p0 = real(p_coeff(1));
q0 = imag(p_coeff(1));
p1 = real(p_coeff(2));
q1 = imag(p_coeff(2));
p2 = real(p_coeff(3));
q2 = imag(p_coeff(3));

D1 = det([q2 q1; p2 p1]);
D2 = det([q2 q1 q0 0; p2 p1 p0 0; 0 q2 q1 q0; 0 p2 p1 p0]);


ops = sdpsettings('verbose',0,'solver','sedumi','dualize',1);


for i=1:length(a)
    for j=1:length(b)
         z = 1 - 4*(a(i)+b(j))^2 / a(i);

        if(z<0 || a(i) == 0)
            sol1(i,j) = -1;
            continue;
        end

        c1 = (1-sqrt(z))/2;
        c2 = a(i)/(2*(a(i)+b(j)))*(1+sqrt(z));

        par = [-a(i)-c2^2 -2*c1*c2 c2^2 -(a(i)+b(j))+2*c1*c2 d 1 v_norm*d v_norm];
        
        dd1 = subs(D1,[a11 a12 a21 a22 d1 d2 v1 v2],par);
        dd2 = subs(D2,[a11 a12 a21 a22 d1 d2 v1 v2],par);
        
        mu1 = coeffs(dd1,'All'); mu2 = coeffs(dd2,'All');
        M1 = poly2mat(dd1);
        M2 = poly2mat(dd2);

        ell1 = ceil((length(coeffs(dd1,zeta,'All'))-1)/2);
        ell2 = ceil((length(coeffs(dd2,zeta,'All'))-1)/2);
        N1 = sdpvar(ell1+1,ell1+1,'symmetric');
        N2 = sdpvar(ell2+1,ell2+1,'symmetric');

        F = [M1 + N1 >=0, M2 + N2 >= 0, eq_constraint(N1), eq_constraint(N2)];
        sol_tmp = optimize(F,[],ops);
              
        sol1(i,j)=sol_tmp.problem;
        if(sol1(i,j) == 4)
            disp(['numerical error ','a=',num2str(a(i)),' b(j)=',num2str(b(j))]);
        end            
    end
    disp(a(i));
end

save(filename_mat,'a','b','v_norm','d','sol1');

hfig = figure('visible','off');
map = [0 0 1; 0 1 0; 1 0 0];
colormap(map)
imagesc(b,a, sol1);
xlim([0 b(length(b))]);
ylim([a(1) a(length(a))]);

axis square
set( gca, 'FontName','Times','FontSize',16 ); 
set( gca, 'YDir', 'normal')
xlabel('b', 'FontName','Times','FontSize',24);
ylabel('a', 'FontName','Times','FontSize',24);

saveas(hfig,filename_fig,'fig');
saveas(hfig,filename_png,'png');

function M = poly2mat(p)
    mu1 = coeffs(p,'All');
    mu1 = fliplr(mu1); 
    nmonomials = ceil((length(mu1)-1)/2)+1;
    M = zeros(nmonomials);
    for k=1:length(mu1)
        if(mod(k,2)) %zeta^2n
            M((k+1)/2,(k+1)/2) = mu1(k);
        else % zeta^(2n+1)
            M(k/2,k/2+1)=mu1(k)/2;
            M(k/2+1,k/2)=M(k/2,k/2+1);
        end
    end
end

function F = eq_constraint(M)
    %M: matrix
    F = [];
    nsize = length(M);
    for ell=2:2*nsize
        antidiag = {};
        if(ell <= nsize+1)
            for j=1:ell-1
            antidiag{end+1} = M(j,ell-j);
            end
        else
            for j=ell-nsize:nsize
            antidiag{end+1} = M(j,ell-j);
            end
        end
        F = [F sum([antidiag{:}]) == 0];
%       sdisplay(sum([antidiag{:}]));        
       
    end
end


