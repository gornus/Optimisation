x0 = [5,5];
A =[]; B=[]; Aeq=[]; Beq=[]; lb=[0,0]; ub=[13.2,13.2];
[x,fval]=fmincon('obj',x0,A,B,Aeq,Beq,lb,ub,'nonlcon')
