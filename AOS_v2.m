
function u = AOS_v2(u,n,g,tau,mu,f,beta,dNu,b,alpha,varsigma)
%%%%%%%%% CMS Spencer, Chen AOS2 %%%%%%%%%%

global theta_tai

h = 1; m = n; N = m*n;

if sum(u(:))==0, fprintf('Broken All Zero'); return; end

btaylor = zeros(n,n);
for i = 1:n
    for j = 1:n
        if (u(i,j) <= varsigma && u(i,j) >= -varsigma) || (u(i,j) <= (1+varsigma) && u(i,j) >= (1-varsigma))
            btaylor(i,j) = b;
        else
            %             keyboard;
        end
    end
end
% btaylor = b*ones(n,n);

ON = speye(N,N);

bvecx = reshape(btaylor,N,1);
taybx = spdiags(bvecx,0,N,N);
bvecy = reshape(btaylor',N,1);
tayby = spdiags(bvecy,0,N,N);

Bx = inv(ON + tau.*alpha.*taybx);
By = inv(ON + tau.*alpha.*tayby);

R0 = reshape((-tau*(f) - tau*dNu),N,1);
R1 = Bx*R0;
R = reshape(R1,n,n);

ux = diff(u);
uy = diff(u')';
ux = [ux;zeros(1,m)];
uy = [uy zeros(n,1)];

ucx = u(3:n,:)-u(1:n-2,:);
ucx = [u(2,:)-u(1,:);ucx;u(n,:)-u(n-1,:)];
ucy = u(:,3:m)-u(:,1:m-2);
ucy = [u(:,2)-u(:,1) ucy u(:,m)-u(:,m-1)];
d1 = 1./sqrt((ux/h).^2+(ucy/(2*h)).^2+beta);
d2 = 1./sqrt((ucx/(2*h)).^2+(uy/h).^2+beta);

d11 = g.*d1;
a1 = zeros(N,1); a2 = zeros(N,1); a3 = zeros(N,1);
for j = 1:m
    j0 = (j-1)*n;
    d_1 = d11(1:n-1,j);
    a1(1+j0:n-1+j0) = +d_1;
    a3(2+j0:n+j0) = +d_1;
    a2(1+j0:n-1+j0) = -d_1;
    a2(2+j0:n+j0) = a2(2+j0:n+j0) - d_1;
    %     a2(1+j0:n+j0) = a2(1+j0:n+j0) - lambda_TV*w(:,j);
end
Ax = spdiags([a1,a2,a3],[-1,0,1],N,N);
A1 = speye(N,N)*(1-tau*theta_tai) - Bx*(2*tau*mu*Ax);

u1 = A1\(reshape((u+R),N,1));
u1 = reshape(u1,n,m);

d22 = g.*d2;
b1 = zeros(N,1); b2 = zeros(N,1); b3 = zeros(N,1);
for i = 1:n
    i0 = (i-1)*m;
    d_2 = d22(i,1:m-1)';
    b1(1+i0:m-1+i0) = +d_2;
    b3(2+i0:m+i0) = +d_2;
    b2(1+i0:m-1+i0) = -d_2;
    b2(2+i0:m+i0) = b2(2+i0:m+i0) - d_2;
    %     b2(1+i0:m+i0) = b2(1+i0:m+i0) - lambda_TV*w(i,:)';
end
Ay = spdiags([b1,b2,b3],[-1,0,1],N,N);
A2 = speye(N,N).*(1-tau*theta_tai) - By*(2*tau*mu*Ay);
u2 = A2\(reshape((u+R)',N,1));
u2 = reshape(u2,m,n);
u2 = u2';

u = (u1+u2)/2;