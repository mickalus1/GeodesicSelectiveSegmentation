function [TC,Iters_Final,Res_Final,u,cols,rows,u_data,R_min,R_max,C_min,C_max] = ConvexSeg_Run(output,z,c1,c2,lambda,tau,Iters,utol,theta)

global sigma1 CV_type lambda3 eps2 w lambda_TV theta_tai

t0 = cputime; [m,n] = size(z);

%% FLAG KEY
% 0 - Spencer-Chen
% 1 - Geodesic

for flag = [1]
    
    tic
    
    %% Geodesic Fitting Term
    
    u_data = struct();
    
    if (theta == 0) || (lambda ==0)
        SSF = zeros(size(z));
        u=z;
        
    else
        
        if flag == 0
            figure;imagesc(z);colormap(gray) % Me
            [Mask,cols,rows] = roipoly;
            close
            
            Mask = imdilate(Mask,strel('disk',3));
            SSF = theta.*Pd(Mask>0.5); % fitting terms
        elseif (flag ==1)
            
            [Pd1,cols,rows,R_min,R_max,C_min,C_max,z] = Geo_MarkerAndAntimarker( z );
            [n,m] = size(z);
            C = convhull(rows,cols);
            Mask = roipoly(Pd1,cols(C),rows(C));
            SSF = theta.*Pd1/max(Pd1(:));
            
        end
        
        u = Mask;
        
        if (size(cols,1)*size(cols,2) == 1) &&...
                (size(rows,1)*size(rows,2) == 1)
            Mask(rows:rows+1,cols:cols+1) = 1;
        end
        
    end
    
    %%%% Additional parameters
    mu = 1; % regularisation term
    varsigma = 1e-2; % parameter in penalty function
    as = 2; % multiply minimum alpha value by as
    beta = eps2; % parameter in curvature
    b = 161.7127690; % from Taylor expansion of fpen fn.
    
    z_sm = imgaussfilt(z,max(1,sqrt(100*sigma1)));
    [gx,gy] = gradient(z_sm);
    nab_z = sqrt(gx.^2+gy.^2);
    beta1 = 10;
    g = 1./ (1+ beta1*nab_z.^2);
    
    %%%% For recording progression of residual
    res = [];
    %%%% Calculate fitting term and set alpha
    
    c1 = sum(z(:).*(Mask(:)>0.5))/sum(Mask(:)>0.5);
    c2 = sum(z(:).*(Mask(:)<0.5))/sum(Mask(:)<0.5);
    
    if CV_type ==0 % Original Chan-Vese
        
        %             f1 = (((z-c1)/(c1+eps)).^2-lambda3 *((z-c2)/(c2+eps)).^2) + SSF;
        f1 = (z-c1).^2-lambda3 *(z-c2).^2 + SSF;
        
    elseif (CV_type==1) % Roberts-Spencer
        
        N = 3;
        K = [0,multithresh(z(:),N-1),1];
        
        % LOWER THRESHOLD
        L_vect = max(c1 - K,0);
        L_vect2 = L_vect;
        TF_vect = 1:size(L_vect,2);
        TF_vect(L_vect==0) = [];
        L_vect2(L_vect==0) = [];
        TF2 = find(L_vect2==min(L_vect2));
        L = K(TF_vect(TF2(1)));
        
        % UPPER THRESHOLD
        H_vect = max(K - c1,0);
        H_vect2 = H_vect;
        TF_vect = 1:size(L_vect,2);
        TF_vect(H_vect==0) = [];
        H_vect2(H_vect==0) = [];
        if isempty(H_vect2)
            H = K(end);
        else
            TF2 = find(H_vect2==min(H_vect2));
            H = K(TF_vect(TF2(1)));
        end
        
        gamma1 = c1 - L;
        gamma2 = H - c1;
        
        TF3 = (z>=c1-gamma1) & (z<=c1);
        TF4 = (z<=c1+gamma2) & (z>c1);
        
        f3 = ( 1+((z-c1)/gamma1) ).*TF3 +...
            ( 1-((z-c1)/gamma2) ).*TF4 ;
        
        f1 = ((z-c1).^2) - lambda3 * f3 + SSF ;
        
    end
    
    fprintf('\n outer = 1 -- c1 = %.3f, c2 = %.3f, do %d Iterations: ____',c1,c2,Iters);
    for l = 1:Iters
        
        if rem(l,50) ==0
            tau = max(1e-3,tau*0.9);
        end
        
        if rem(l,2)==0
            
            if (CV_type == 1)
                
                % LOWER THRESHOLD
                L_vect = max(c1 - K,0);
                L_vect2 = L_vect;
                TF_vect = 1:size(L_vect,2);
                TF_vect(L_vect==0) = [];
                L_vect2(L_vect==0) = [];
                TF2 = find(L_vect2==min(L_vect2));
                L = K(TF_vect(TF2(1)));
                
                % UPPER THRESHOLD
                H_vect = max(K - c1,0);
                H_vect2 = H_vect;
                TF_vect = 1:size(L_vect,2);
                TF_vect(H_vect==0) = [];
                H_vect2(H_vect==0) = [];
                
                if isempty(H_vect2)
                    H = K(end);
                else
                    TF2 = find(H_vect2==min(H_vect2));
                    H = K(TF_vect(TF2(1)));
                end
                
                gamma1 = c1 - L;
                gamma2 = H - c1;
                
                TF3 = (z>=c1-gamma1) & (z<=c1);
                TF4 = (z<=c1+gamma2) & (z>c1);
                
                f3 = ( 1+((z-c1)/gamma1) ).*TF3 +...
                    ( 1-((z-c1)/gamma2) ).*TF4 ;
                
                f1 = ((z-c1).^2) - lambda3 * f3 + SSF ;
                
            end
        end
        
        if rem(l,50)==0
            figure(1);
            subplot(121);imagesc(z);colormap(gray);hold on; contour(u,[0.5 0.5],'r','LineWidth',2);drawnow
            subplot(122);surf(u);drawnow
        end
        
        oldu = u;
        
        f0 = lambda*f1/max(f1(:));
        f11 = f0(:);
        A = (1/2)*(norm(f11,Inf)); % minimum for alpha from Chan paper
        alpha = as*A;
        
        %%%% Calculate penalty term
        N1 = sqrt((2*u-1).^2+varsigma)-1;
        Hnu = 1/2+1/pi*atan(N1./varsigma);
        dN = ((4*u-2)./(N1+1)).*((varsigma*N1)./(pi*(varsigma^2+N1.^2))+Hnu);
        dNu = alpha*dN;
        dNu = double(dNu);
        
        %%%% Perform AOS iteration
        %%%% v2 is from CMS 15 paper
        u = double(u);
        u = AOS_v2(u,n,g,tau,mu,f0,beta,dNu,b,alpha,varsigma);
        
        %%%% Check residual and store
        R = Residual(u.*(g>0.95),oldu.*(g>0.95)); res = [res,R];
        
        if rem(l,50)==0
            fprintf('It:%d, Res:%0.4e \n',l,R);
        end
        
        %%%% Stopping criterion for u
        if R < utol, break; end
        
        u_data(l).h = u;
        
    end % l loop
    fprintf(' == Done (%d iters) \n',l);
    toc
    
    Res_Final = res;
    Iters_Final = l;
    TC = 1;
    
    t0 = cputime-t0;
    
end

function R = Residual(new,old)

[n,m] = size(new);
new = reshape(new,n*m,1);
old = double(reshape(old,n*m,1));

Res = new - old;
R = norm(Res)/norm(old);

function z = normalize0a(z,a)
% Normalize to the range of [0,1]
fmin  = min(z(:));
fmax  = max(z(:));
z = (z-fmin)/(fmax-fmin);
z = a*z;

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