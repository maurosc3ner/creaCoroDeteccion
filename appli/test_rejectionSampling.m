rand('seed',12345);
x = -10:.1:10;
% CREATE A "COMPLEX DISTRIBUTION" f(x) AS A MIXTURE OF TWO NORMAL
% DISTRIBUTIONS
f = inline('normpdf(x,3,2) + normpdf(x,-5,1)','x');
t = plot(x,f(x),'b','linewidth',2); hold on;
 
% PROPOSAL IS A CENTERED NORMAL DISTRIBUTION
q = inline('normpdf(x,0,4)','x');
 
% DETERMINE SCALING CONSTANT
c = max(f(x)./q(x))
 
%PLOT SCALED PROPOSAL/ENVELOP DISTRIBUTION
p = plot(x,c*q(x),'k--');
 
% DRAW A SAMPLE FROM q(x);
qx = normrnd(0,4);
fx = f(qx);
 
% PLOT THE RATIO OF f(q(x)) to cq(x)
a = plot([qx,qx],[0 fx],'g','Linewidth',2);
r = plot([qx,qx],[fx,c*q(qx)],'r','Linewidth',2);
legend([t,p,a,r],{'Target','Proposal','Accept','Reject'});
xlabel('x');