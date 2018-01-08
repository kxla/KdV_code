function vkdv0
% programming: Alex Sheremet, 2017
% variable-coefficient KdV integrator, ver13.
% details in: 3. Sheremet, A., U. Gravois, and V. Shrira (2016), 
% Meteotsunami observations on the Louisiana shelf: lone solitons and soliton packs, 
% Nat Hazards DOI 10.1007/s11069-016-2446-2.

clear all
global problemData

	setXTGrids;			% define space-time grids (not scaled)
	setScales;			% define scaling parameters
	setScaledConfGrids;		% set configuration grids (scaled)
	setScaledPerturbation;		% define initial perturbation (scaled)
	plotProblemData;		% plot problem data
	ssSolver;			% integrate vKdV
	plotSolution;		% plot solution 
end
	
% -------------------------various functions-----------------------------------	
	
function setXTGrids 
global problemData
% flat bottom @ 20 m depth
	xb=[0,50000]; 		% begin/end along-ray coordinate (m)
	dx=20; 			% spatial increment (m)
	hb=[20,20]; 		% begin/end along-ray depth (m)
	db=[1,1]; 			% begin/end ray separation
	m=12; N=2^m; 		% Fourier transform power; 	
	TMax=100*60; 		% duration of time window (s)
% physical grid (in physical units)
	problemData.x=0:dx:xb(end); 
	problemData.h=interp1(xb,hb,problemData.x); 
	problemData.d=interp1(xb,db,problemData.x);
	problemData.t=linspace(0,TMax,N)'; 
end	
	
	
function setScales 
global problemData
	TScale=60; 			% time scale in sec
	AScale=1; 			% height scale im m
	NU=0; 			% kinematic viscosity
% save scales
	problemData.pars.A=AScale; 
	problemData.pars.T=TScale; 
	problemData.pars.NU=NU;
end
	
	
function setScaledConfGrids 	
global problemData
	[S,tShift]=xt2STGrid;	
	problemData.scl.S=S; 
	problemData.scl.tShift=tShift; 
end


function setScaledPerturbation 
global problemData
	t=problemData.t/problemData.pars.T; TMax=t(end);
	[U,~]=URNumber(problemData.h(1),1,problemData.pars.A,problemData.pars.T); 
% 2-soliton collision (standard KdV test)
	a1=1.2; 			% scaled amplitude of large soliton
	t1=0.6; 			% center time soliton1 on scaled time axis
	a2=0.3; 			% scaled amplitude of small soliton 
	t2=0.59; 			% center time soliton2 on scaled time axis 
	t1=t1*TMax; t2=t2*TMax;
	th=0.5*sqrt(U)*sqrt(a1/3)*(t-t1);  q1=a1*(sech(th)).^2;
	th=0.5*sqrt(U)*sqrt(a2/3)*(t-t2);  q2=a2*(sech(th)).^2;
	Q=q1+q2;
	problemData.Q=Q(:);
	problemData.QType='Soliton collision';
end
	

function plotProblemData 
global problemData
	ftag='ray configuration'; figWin(ftag);
	hax=sameAxSubV(2,1,ftag);
	hx=hax(1); axes(hx);
	plot(problemData.scl.S,-problemData.h);
	axis tight; set(hx,'xaxisloc','top');
	ylabel('Depth(m)'); xlabel('Scaled along-ray distance');
% plot initial perturbation
	hx=hax(2); axes(hx)
	plot(problemData.t/problemData.pars.T,problemData.Q); axis tight;
	xlabel('Scaled time'); ylabel('Initial perturbation'); 		
end


function plotSolution 
global problemData
	ftag='solution'; figWin(ftag);
	tLim=[55,65]*problemData.pars.T;
	kt=find(problemData.t>=tLim(1) & problemData.t<=tLim(2));
	pcolor(problemData.x/1000,problemData.t(kt)/60,problemData.Q(kt,:)); 
	shading interp;
	xlabel('Along-ray distance (km)');
	ylabel('Relative arrival time (min)')
	title(problemData.QType,'fontsi',12,'fontw','n');
end
	

function ssSolver	
% split step Fourier solver for the variable KdV Solver
global problemData OM
	t1=cputime;
% load problem parameters
	S=problemData.scl.S; ns=length(S); 
	dt=(problemData.t(2)-problemData.t(1))/problemData.pars.T;
	N=length(problemData.t);
	U=[]; R=[];
	for k=1:length(problemData.x), 
		[u,r]=URNumber(problemData.h(k),problemData.d(k),problemData.pars.A,problemData.pars.T); 
		U(k)=u; R(k)=r;
	end
	Nu=problemData.pars.NU;
	h=problemData.h; c=sqrt(9.81*h);
% caution - funny omega grid, handle with care
	dOM=2*pi/(N*dt); OM=[0:dOM:N/2*dOM,-(N/2-1)*dOM:dOM:-dOM]; OM=OM(:);	
	NU=ones(size(OM))*Nu; 

	opts=odeset('RelTol',1e-10);
%	fprintf('computing time evolution:%8d/%8d',0,ns);

	outQ=problemData.Q(:,1); 
	Q=fft(problemData.Q(:,1)); 
	for ks=2:ns
		if ks/50==floor(ks/50),
			fprintf('s: %g/%g h=%g (%d/%d)\n',S(ks),S(end),h(ks),ks,ns);
		end
		ds=S(ks)-S(ks-1);
		K=-i*OM.^3/U(ks)-0.75*NU/c(ks)/h(ks)/R(ks);
		Q1=Q.*exp(0.5*K*ds);
		q1=[real(Q1);imag(Q1)];
		q1=RK4(0,q1,ds); % nonlinear term handled by Runge-Kutta 4	
		Q1=q1(1:end/2)+i*q1(end/2+1:end);
		Q1=Q1.*exp(0.5*K*ds);
		Q=Q1;
		outQ=[outQ,real(ifft(Q(:)))];  
	end 
	fprintf('\n');
	problemData.Q=outQ;
	t2=cputime;
	fprintf('odeSolver time: %g min\n',(t2-t1)/60);
end


function [S,tShift]=xt2STGrid    
% xt  is the physical grid; ST is the configuration (integration-space) grid
global problemData

% compute S	
	dx=diff(problemData.x);
	c=sqrt(9.81*problemData.h);
% compute S	
	F=1.5./problemData.h./c./sqrt(problemData.d);
	F=0.5*(F(1:end-1)+F(2:end));
	S=[0,cumsum(F.*dx)]; S=S*problemData.pars.A/problemData.pars.T;
% compute local time:
	F=1./c;
	F=0.5*(F(1:end-1)+F(2:end));
	tShift=[0,cumsum(F.*dx)]; 
end


function dq=dQNL(x,q)	
global OM
	Q=q(1:end/2)+i*q(end/2+1:end);
	dQ=0.5*i*OM.*fft(real(ifft(Q)).^2);
	dq=[real(dQ);imag(dQ)];
end


function q=RK4(x,q,dx)	
% Runge-Kutta 4th order
 	k1 = dQNL(x,q);  q1 = q+k1*dx/2; 
  	k2 = dQNL(x,q1); q2 = q+k2*dx/2;      
  	k3 = dQNL(x,q2); q3 = q+k3*dx;       
  	k4 = dQNL(x,q3);     
  	q=q+(k1+2*k2+2*k3+k4)*dx/6;
end

	
function [U,R]=URNumber(h,d,A,T)	
	c=sqrt(9.81*h); d2=sqrt(d);
	U=9*9.81*A*T^2./(h.^2.*d2);
	R=1.5*A^2./(h.*c.*d2*T);	
end


function figWin(ftag)      
      if isempty(findobj('tag',ftag)), h=figure; set(h,'tag',ftag); 
      else, figure(findobj('tag',ftag)); clf;     
      end
end      
   
   
function hax=sameAxSubV(N,M,ftag,varargin)
      if ~isempty(varargin), kRange=varargin{1}; 
      else, kRange=1:N*M;
      end
% This sets the subplot axes
      Y0=0.92;  Dy=0.8/N; dy=0.93*Dy;
      X0=0.12; Dx=0.8/M; dx=0.93*Dx; 
      hax=[]; k=0;
      for iN=1:N, 
      for iM=1:M
		k=k+1; 
		if ismember(k,kRange),
			pos=[X0+(iM-1)*Dx,Y0-iN*Dy,dx,dy];
			hx=axes('par',findobj('tag',ftag),'pos',pos);		     
			hax=[hax,hx];
		end
	end 
	end      
end
