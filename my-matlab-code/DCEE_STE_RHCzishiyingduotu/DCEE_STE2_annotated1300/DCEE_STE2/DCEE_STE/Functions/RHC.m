function [var, npos] = RHC(HL,count,moveDist,origin_pos,pos,theta,Wpnorm,D,m,s,M,MM,i,P_k_store,N,PF_Memory,domain,xmin,xmax,ymin,ymax,zmin,zmax)

if HL == 0
    return
end

ynew = [moveDist,moveDist,0,-moveDist,-moveDist,-moveDist,0,moveDist];
xnew = [0,moveDist,moveDist,moveDist,0,-moveDist,-moveDist,-moveDist];
znew = zeros(1,8);

% Entrotaxis reward
indx = resampleStratified(Wpnorm,M);
d.x = theta.x(indx);
d.y = theta.y(indx);
d.z = theta.z(indx);
d.Q = theta.Q(indx);
d.u = theta.u(indx);
d.phi = theta.phi(indx);
d.ci = theta.ci(indx);
d.cii = theta.cii(indx);
%d.duration = 0;%theta.duration(indx);

Xneighbour = zeros(1,8);
Yneighbour = zeros(1,8);
Zneighbour = zeros(1,8);
var = zeros(1,8);
dist_theta = zeros(1,8);
theta_RMSE = zeros(1,8);

for j = 1:8
    npos.x_matrix = pos.x_matrix+xnew(j);
    npos.y_matrix = pos.y_matrix+ynew(j);
    npos.z_matrix = pos.z_matrix+znew(j);

    if npos.x_matrix<xmin || npos.x_matrix>xmax || ...
       npos.y_matrix<ymin || npos.y_matrix>ymax || ...
       npos.z_matrix<zmin || npos.z_matrix>zmax
        % invalid position
        var(j)=NaN;
        dist_theta(j)=NaN;
        theta_RMSE(j)=NaN;
        continue
    end

    % calculate the minimum distance for the position
    dist = count*moveDist;
    % calculate the distance between the position and the origin position
    distance = sqrt((npos.x_matrix-origin_pos.x_matrix).^2+(npos.y_matrix-origin_pos.y_matrix).^2+(npos.z_matrix-origin_pos.z_matrix).^2);

    if distance < dist
        % invalid distance
        var(j)=NaN;
        dist_theta(j)=NaN;
        theta_RMSE(j)=NaN;
        continue
    end

    % save the coordinates of the position
    Xneighbour(j) = npos.x_matrix;
    Yneighbour(j) = npos.y_matrix;
    Zneighbour(j) = npos.z_matrix;

    if HL-1 > 0
        % still need to calculate the next position

        %predict measurement
        ni = i + 1;
        Dsim=PredictMeasurement(s,npos,m.thresh);
        D(ni)=Dsim;
    
        %update particle filter
        [theta_updated,Wpnorm_updated]=UpdatePFPlume(D,theta,Wpnorm,npos,P_k_store,m,N,PF_Memory,domain);
        
        % calculate the next position recursively
        [nvar, ~] = RHC(HL-1,count+1,moveDist,origin_pos,npos,theta_updated,Wpnorm_updated,D,m,s,M,MM,i,P_k_store,N,PF_Memory,domain,xmin,xmax,ymin,ymax,zmin,zmax);
       
        % get the smallest variance
        next_var = min(nvar);
    else
        % complete the calculation of the next position
        next_var = NaN;
    end

    % calculate the variance of the position
    pC = RadioactiveDispersionModel(theta,npos);

    desC = RadioactiveDispersionModel(d,npos);
    designd = repmat(desC,1,MM);
    designd = designd+(1*designd.*randn(M,MM));
    designd(rand<0.3)=0;
    designd(designd<m.thresh)=0;

    for k = 1:M
        for kk = 1:MM
            dC = designd(k,kk);

            zWp= Likelihood_Like_Yee(pC, dC, Wpnorm, m);
            zWpnorm = zWp./sum(zWp);

            theta_mean_xy=N*[mean(theta.x.*zWpnorm),mean(theta.y.*zWpnorm)];
            theta_RMSE(j)=theta_RMSE(j)+sum(zWpnorm.*vecnorm([theta.x,theta.y]'-theta_mean_xy')'.^2)/(M*MM);              
            dist_theta(j)=dist_theta(j)+norm(theta_mean_xy-[Xneighbour(j),Yneighbour(j)])/(M*MM);
        end
    end
    var_j = dist_theta(j)+theta_RMSE(j);
    
    % get the minimum value of the variance between the position and the next position
    var(j) = min(next_var, var_j);
end

% get the coordinates of the next location
[~,ind] = min(var);
npos.x_matrix = Xneighbour(ind);
npos.y_matrix = Yneighbour(ind);
npos.z_matrix = Zneighbour(ind);

end

