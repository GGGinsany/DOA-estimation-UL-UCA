function err = cal_theta_match_error_ul(idea_theta, esti_theta)
L = length( idea_theta );

if isempty( esti_theta )|| (length(idea_theta)~=length(esti_theta))
    err = inf;
    return;
end
 
err=sum((sort(idea_theta)-sort(esti_theta)).^2);
err = sqrt( err./L );
end
