%
% Fig. 3
%
close all;
clear;

a = 0.06;
b = 0.04;
v_norm = 0.3162;
d = 6;

Xstart = 0.10;
Xend = 0.25;

res = zeros(length(a),length(b));
sol1 = zeros(length(a),length(b));
sol2 = zeros(length(a),length(b));
sol3 = zeros(length(a),length(b));

filename_head = sprintf('Figure4');
filename_mat = strcat(filename_head,'.mat');
filename_fig = strcat(filename_head,'.fig');

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
            sol1(i,j)= -1;
            sol2(i,j)= -1;
            sol3(i,j)= -1;
            res(i,j) = -1;
            continue;
        end

        c1 = (1-sqrt(z))/2;
        c2 = a(i)/(2*(a(i)+b(j)))*(1+sqrt(z));
        
        par = [-a(i)-c2^2 -2*c1*c2 c2^2 -(a(i)+b(j))+2*c1*c2 d 1 v_norm*d v_norm];

        dd1 = subs(D1,[a11 a12 a21 a22 d1 d2 v1 v2],par);
        dd2 = subs(D2,[a11 a12 a21 a22 d1 d2 v1 v2],par);

        %%
        %% solve SOS for [0,Xstart]
        %%
        tilde_d1 = subs(dd1,zeta,((Xstart-0)*zeta+(Xstart+0))/2);        
        tilde_d2 = subs(dd2,zeta,((Xstart-0)*zeta+(Xstart+0))/2);
        
        M1 = poly2mat(tilde_d1);
        M2 = poly2mat(tilde_d2);

        ell1 = ceil((length(coeffs(tilde_d1,zeta,'All'))-1)/2);
        ell2 = ceil((length(coeffs(tilde_d2,zeta,'All'))-1)/2);
        K1 = sdpvar(ell1+1); L1 = sdpvar(ell1);
        K2 = sdpvar(ell2+1); L2 = sdpvar(ell2);

        constraint = [];
        for m = 0:length(coeffs(tilde_d1,zeta,'All'))-1 %m: Delta_1(s) = \sum_{m} delta_m zeta^m
            %sdisplay(sum([-antidiag_sum(M1,m+2),antidiag_sum(K1, m+2),antidiag_sum(L1,m+2), -antidiag_sum(L1,m)]));
            constraint = [constraint, sum([-antidiag_sum(M1,m+2),antidiag_sum(K1, m+2),antidiag_sum(L1,m+2), -antidiag_sum(L1,m)])==0];
        end
        for m = 0:length(coeffs(tilde_d2,zeta,'All'))-1 %m: Delta_2(s) = \sum_{m} delta_m zeta^m
%            sdisplay(sum([-antidiag_sum(M2,m+2),antidiag_sum(K2, m+2),antidiag_sum(L2,m+2), -antidiag_sum(L2,m)]));
            constraint = [constraint, sum([-antidiag_sum(M2,m+2),antidiag_sum(K2, m+2),antidiag_sum(L2,m+2), -antidiag_sum(L2,m)])==0];
        end

        sol=optimize([constraint, K1 >= 0, L1 >= 0, K2 >= 0, L2 >= 0],[],ops);
        sol1(i,j)=sol.problem;

        %% validation
        if 0
           sdpvar x1;
           z1 = monolist(x1, length(K2)-1);
           z2 = monolist(x1, length(L2)-1);
           sdisplay(z1'*value(K2)*z1 + (1 - x1^2)*z2'*value(L2)*z2);
           coefficients(z1'*value(K2)*z1 + (1 - x1^2)*z2'*value(L2)*z2)
           double(coeffs(tilde_d2))
        end
        %value(K2)*10^4
        %value(L2)*10^4
        
        %%
        %% solve SOS for [Xstart,Xend]
        %%
        tilde_d1 = subs(dd1,zeta,zeta+Xstart);        
        tilde_d2 = subs(dd2,zeta,zeta+Xstart);

        M1 = poly2mat(tilde_d1);
        M2 = poly2mat(tilde_d2);

        ell1 = ceil((length(coeffs(tilde_d1,zeta,'All'))-1)/2);
        ell2 = ceil((length(coeffs(tilde_d2,zeta,'All'))-1)/2);
        K1 = sdpvar(ell1+1); L1 = sdpvar(ell1);
        K2 = sdpvar(ell2+1); L2 = sdpvar(ell2);

        constraint = [];
        for m = 0:length(coeffs(tilde_d1,zeta,'All'))-1 %m: Delta_1(s) = \sum_{m} delta_m zeta^m
            %sdisplay(sum([-antidiag_sum(M1,m+2),antidiag_sum(K1, m+2),antidiag_sum(L1,m+2), -antidiag_sum(L1,m)]));
            constraint = [constraint, sum([-antidiag_sum(M1,m+2),antidiag_sum(K1, m+2),antidiag_sum(L1,m+2), -antidiag_sum(L1,m)])==0];
        end
        for m = 0:length(coeffs(tilde_d2,zeta,'All'))-1 %m: Delta_2(s) = \sum_{m} delta_m zeta^m
%            sdisplay(sum([-antidiag_sum(M2,m+2),antidiag_sum(K2, m+2),antidiag_sum(L2,m+2), -antidiag_sum(L2,m)]));
            constraint = [constraint, sum([-antidiag_sum(M2,m+2),antidiag_sum(K2, m+2),antidiag_sum(L2,m+2), -antidiag_sum(L2,m)])==0];
        end

        sol=optimize([constraint, K1 >= 0, L1 >= 0, K2 >= 0, L2 >= 0],[],ops);
        sol2(i,j)=sol.problem;

        %%
        %% solve SOS for [Xend,infty)
        %%
        tilde_d1 = subs(dd1,zeta,zeta+Xend);        
        tilde_d2 = subs(dd2,zeta,zeta+Xend);

        M1 = poly2mat(tilde_d1);
        M2 = poly2mat(tilde_d2);

        ell1 = ceil((length(coeffs(tilde_d1,zeta,'All'))-1)/2);
        ell2 = ceil((length(coeffs(tilde_d2,zeta,'All'))-1)/2);
        K1 = sdpvar(ell1+1); 
        K2 = sdpvar(ell2+1); 

        if(mod(length(coeffs(tilde_d1,zeta,'All'))-1, 2) == 1) %deg(Delta_i(zeta)) == odd
            L1 = sdpvar(ell1+1);
        else
            L1 = sdpvar(ell1);
        end
        if(mod(length(coeffs(tilde_d2,zeta,'All'))-1, 2) == 1) %deg(Delta_i(zeta)) == odd
            L2 = sdpvar(ell2+1);
        else
            L2 = sdpvar(ell2);
        end        

        constraint = [];
        for m = 0:length(coeffs(tilde_d1,zeta,'All'))-1 %m: Delta_1(s) = \sum_{m} delta_m zeta^m
            %sdisplay(sum([-antidiag_sum(M1,m+2), antidiag_sum(K1, m+2), antidiag_sum(L1,m+1)]));
            constraint = [constraint, sum([-antidiag_sum(M1,m+2), antidiag_sum(K1, m+2), antidiag_sum(L1,m+1)])==0];
        end
        for m = 0:length(coeffs(tilde_d2,zeta,'All'))-1 %m: Delta_2(s) = \sum_{m} delta_m zeta^m
            %sdisplay(sum([-antidiag_sum(M2,m+2), antidiag_sum(K2, m+2), antidiag_sum(L2,m+1)]));
            constraint = [constraint, sum([-antidiag_sum(M2,m+2), antidiag_sum(K2, m+2), antidiag_sum(L2,m+1)])==0];
        end

        sol=optimize([constraint, K1 >= 0, L1 >= 0, K2 >= 0, L2 >= 0],[],ops);
        sol3(i,j)=sol.problem;

        %% validation
        if 0
           sdpvar x1;
            z1 = monolist(x1, length(K2)-1);
            z2 = monolist(x1, length(L2)-1);
            sdisplay(z1'*value(K2)*z1 + x1*z2'*value(L2)*z2);
            coefficients(z1'*value(K2)*z1 + x1*z2'*value(L2)*z2)
            double(coeffs(tilde_d2))
        end        
%        value(K2)
%        value(L2)
        %%
        %% result
        %%
        if sol1(i,j) == 0 && sol2(i,j) == 2 && sol3(i,j) == 0
            res(i,j) = 1;
        else
            res(i,j)=0;
        end
    end
    disp(a(i));
end


ezplot(dd2,[0,0.25])
ylim([-1e-3 4e-3])
hold on;

h=line([0 0.25],[0 0]);
h.LineStyle = '--';
h.Color='k';

axis square
title('');
ylabel('\Delta_2(\zeta)');

savefig('Figure3.fig');

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

function H = antidiag_sum(M,ell)
%M: matrix
nsize = length(M);
G= {};
if(ell <= nsize+1)
    for j=1:ell-1
        G{end+1} = M(j,ell-j);
    end
else
    for j=ell-nsize:nsize
        G{end+1} = M(j,ell-j);
    end
end
H = sum([G{:}]);
end


