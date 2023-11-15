function result = neo(lambda, T)
    
    % sellmeier oefficients %
    % Choose data from Gayer paper;
    % 1 = ne for 5% MgO doped CLN
    % 2 = no for 5% MgO doped CLN
    % 3 = ne for 1% MgO doped SLN
    filename = 'sellmeierLnB.csv'; coef = csvread(filename);
    dataCol = 1;
    
    a1 = coef(dataCol, 1);
    a2 = coef(dataCol, 2);
    a3 = coef(dataCol, 3);
    a4 = coef(dataCol, 4);
    a5 = coef(dataCol, 5);
    a6 = coef(dataCol, 6);
    b1 = coef(dataCol, 7);
    b2 = coef(dataCol, 8);
    b3 = coef(dataCol, 9);
    b4 = coef(dataCol, 10);
    
    
    result = sqrt(a1 + b1*f(T) + (a2 + b2*f(T)) / (lambda^2-(a3+b3*f(T))^2) + (a4+b4*f(T)) / (lambda^2-a5^2) - a6*lambda^2);

end