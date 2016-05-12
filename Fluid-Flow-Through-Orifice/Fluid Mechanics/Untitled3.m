d = 0.01; 
D = 0.109;
k1 = 164; 
k2 = 1; 
g = 9.81; 
p = 998; 
e = 1.5*10^-7; 
l = .127;
u = 1.003*10^-3; 
h = .15; 
H = .044;
ts = 0.0001; 
rows = (h-H)/ts;
value = zeroes(rows, 4); 
solution = (3,1); 


%laminar 
a = k1*u/(p*d) + k2 + 1 + 6.4*u*l/(p*d^2);

%turbleunt 
b = (k1*u/(p*d) + k2)/2; 
k = ((e/d)/3.7)^1.11; 


% for l = lStart:.001:lEnd    
% 
%    for j = 1:4
%      Re(i,j) = (p*l*v(1,j))/u;
%    end 
%    i = i + 1; 
% end

%starting height 
value(1,1) = h; 


for t=1:rows  
c = -2*g*(value(t,1) - H); 

%start off with turbleunt assumption 
%get value of v(t) 
f = v*b + v*(l/(2*d))*(log(  ((e/d)/3.7)^1.11 + 6.9*u/(v*p*d))^2) + (v^2)/2 + c; 
solution = fsolve(f, [0,0]); 

%fuck it ttititititititi

end


   
display(value)

