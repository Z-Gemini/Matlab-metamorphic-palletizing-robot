function y = coordinate(x) 
    ab=0.6;
    ae=0.196;
    bd=0.455;
    de=0.35;
    
    xb=ab*cos(x);
    yb=ab*sin(x);
    xe=0;
    ye=-ae;
    
    for i=1:length(x)
        azi(i)=azimuth(xb(i),yb(i),xe,ye);
        be(i)=sqrt((xb(i)-xe).^2+(yb(i)-ye).^2);
        dbe(i)=acos((bd.^2+be(i).^2-de.^2)/(2*bd*be(i)));
        xd(i)=xb(i)+bd*cos(azi(i)-dbe(i));
        yd(i)=yb(i)+bd*sin(azi(i)-dbe(i));
        y(i)=azimuth(xe,ye,xd(i),yd(i));
    end
end