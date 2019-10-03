function ln=leg(x,n)

ln=0.;
if n~= 0 && n~= 1
for l=0:n
a=factorial(n)./(factorial(l)*factorial(n-l));
b=factorial(n)./(factorial(n-l)*factorial(l));
c=(x-1).^l;
d=(x+1).^(n-l);
ln=ln+((1./(2^n)).*a.*b.*c.*d);
end
elseif n== 1
ln=x;
else
ln=1.;
end
