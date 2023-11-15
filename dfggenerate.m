function adot=dfggenerate(K,g,a)

adot=zeros(3,1);

% these are the three equations from chapter 6 in Tideman and Rotwitt 
% page 150

adot(1) = 1i*K*g*a(3)
adot(2) = 1i*K*g*a(3)
adot(3) = 0

end