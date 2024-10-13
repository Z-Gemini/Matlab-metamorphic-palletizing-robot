function [t,v] = compute_time_velocity(x1,x2)
    tc = 0.2;
    vc = deg2rad(3);
    if(abs(x2-x1))/vc > tc
       t=abs(x2-x1)/vc;
       v=vc * (abs(x2-x1))/((x2-x1));
    else
        t=tc;
        v=(x2-x1)/tc;
    end
end
