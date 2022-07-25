function y=G(eco,x)
F=x(1);
M1=x(2);
M2=x(3);
fy=x(4)*1000;

if length(x)==4
    h=eco(1);
    b=eco(2);
else
    h=x(5);
    b=x(6);
end

y=(1-(4*M1/(b*h*h*fy))-(4*M2/(b*b*h*fy))-(F*F/(b*h*fy)^2)); 
end

