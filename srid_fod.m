% Step response invariant discretization of fractional order integrators
% 
% srid_fod function is prepared to compute a discrete-time finite dimensional
% (z) transfer function to approximate a continuous-time fractional order 
% integrator/differentiator function s^r, where "s" is the Laplace transform variable, and "r" is a
% real number in the range of (-1,1). s^r is called a fractional order
% differentiator if 0 < r < 1 and a fractional order integrator if -1 < r < 0.
%
% The proposed approximation keeps the step response "invariant"
%
% IN: 
%       r: the fractional order
%       Ts: the sampling period
%       norder: the finite order of the approximate z-transfer function 
%               (the orders of denominator and numerator z-polynomial are the same)
% OUT: 
%       sr: returns the LTI object that approximates the s^r in the sense
%       of step response. 
% TEST CODE
% dfod=srid_fod(-.5,.01,5);figure;pzmap(dfod)
%
% Reference: YangQuan Chen. "Impulse-invariant and step-invariant
% discretization of fractional order integrators and differentiators".
% August 2008. CSOIS AFC (Applied Fractional Calculus) Seminar.
% http://fractionalcalculus.googlepages.com/
% --------------------------------------------------------------------
% YangQuan Chen, Ph.D, Associate Professor and Graduate Coordinator
% Department of Electrical and Computer Engineering,
% Director, Center for Self-Organizing and Intelligent Systems (CSOIS)
% Utah State University, 4120 Old Main Hill, Logan, UT 84322-4120, USA
% E: yqchen@ece.usu.edu or yqchen@ieee.org, T/F: 1(435)797-0148/3054; 
% W: http://www.csois.usu.edu or http://yangquan.chen.googlepages.com 
% --------------------------------------------------------------------
%
% 9/6/2009 
% Only supports when r in (-1,0). That is fractional order integrator
% To get fractional order differentiator, use 1/sr.
%
% See also irid_fod.m at
% http://www.mathworks.com/matlabcentral/files/21342/irid_fod.m

function [sr]=srid_fod(r,Ts,norder)
if nargin<3; norder=5; end
if Ts < 0 , sprintf('%s','Sampling period has to be positive'),     return, end
if r>=0 | r<= -1, sprintf('%s','The fractional order should be in (-1,0)'), return, end
if norder<2, sprintf('%s','The order of the approximate transfer function has to be greater than 1'), return, end
% 
L=200; %number of points of the step response function h(n)
Taxis=[0:L-1]*Ts;r0=r;r=abs(r);n=0:L-1;h=[(Ts^r)*(n.^(r))/gamma(r)/r]; 
[b,a] = stmcb(h,ones(size(h)),norder,norder,100);sr=tf(b,a,Ts);

% Note that the generated "sr" LTI object might be nonminimum phase!
% although a good fitting is obtained

if 1  % change this to 0 if you do not want to see plots
% approximated h()
wmax0=2*pi/Ts/2; % rad./sec. Nyquist frequency
hhat=step(sr,Taxis); 
figure;plot(Taxis,hhat,'r');hold on;plot(Taxis,h,'ok')
xlabel('time');ylabel('step response'); legend(['approximated  for 1/s^{',num2str(abs(r)),'}'],'true')
figure;
wmax=floor(1+ log10(wmax0) ); wmin=wmax-5;
w=logspace(wmin,wmax,1000);
srfr=(j*w).^(-r); 
subplot(2,1,1)
semilogx(w,20*log10(abs(srfr)),'r');grid on
hold on;
srfrhat=freqresp(sr,w);semilogx(w,20*log10(abs(reshape(srfrhat, 1000, 1))),'k');grid on
xlabel('frequency in Hz');ylabel('dB');
legend('true mag. Bode','approximated mag. Bode')
subplot(2,1,2)
semilogx(w,(180/pi) * (angle(srfr)),'r');grid on;hold on
semilogx(w,(180/pi) * (angle(reshape(srfrhat, 1000, 1))),'k');grid on
xlabel('frequency in Hz');ylabel('degree');
legend('true phase Bode','approximated Phase Bode')
figure;pzmap(sr)
end % if 1
% get stable, minimum phase approximation.
[zz,pp,kk]=zpkdata(sr,'v');
for i=1:norder;
    if abs(zz(i)) > 1
        kk=kk*(-zz(i)); zz(i)=1/zz(i); 
        sprintf('%s','nonminimum phase approximation - forced minimum phase!!'), 
    end
    if abs(pp(i)) > 1
        kk=kk/(-pp(i)); pp(i)=1/pp(i); 
        sprintf('%s','unstable approximation - forced stable!!'), 
    end
end
sr1=zpk(zz,pp,kk,Ts);

if 1  % change this to 0 if you do not want to see plots
% approximated h()
wmax0=2*pi/Ts/2; % rad./sec. Nyquist frequency
hhat=step(sr1,Taxis); 
figure;plot(Taxis,hhat,'r');hold on;plot(Taxis,h,'ok')
xlabel('time');ylabel('step response'); legend(['approximated for 1/s^{',num2str(abs(r)),'}'],'true')
figure;
wmax=floor(1+ log10(wmax0) ); wmin=wmax-5;
w=logspace(wmin,wmax,1000);
srfr=(j*w).^(-r); 
subplot(2,1,1)
semilogx(w,20*log10(abs(srfr)),'r');grid on
hold on;
srfrhat=freqresp(sr1,w);semilogx(w,20*log10(abs(reshape(srfrhat, 1000, 1))),'k');grid on
xlabel('frequency in Hz');ylabel('dB');
legend('true mag. Bode','approximated mag. Bode')
subplot(2,1,2)
semilogx(w,(180/pi) * (angle(srfr)),'r');grid on;hold on
semilogx(w,(180/pi) * (angle(reshape(srfrhat, 1000, 1))),'k');grid on
xlabel('frequency in Hz');ylabel('degree');
legend('true phase Bode','approximated Phase Bode')
figure;pzmap(sr1)
end % if 1

 