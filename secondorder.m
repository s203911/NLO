function adot=secondorder(z,a,K,g,dk)

adot=zeros(3,1);

% these are the three equations from chapter 6 in Tideman and Rotwitt 

adot(1)=1i*K*g*a(3)*conj(a(2))*exp(1i*dk*z);
adot(2)=1i*K*g*a(3)*conj(a(1))*exp(1i*dk*z);
adot(3)=1i*K*g*a(1)*a(2)*exp(-1i*dk*z);



end