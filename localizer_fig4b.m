clear
kappa = 1.25e-7;

a=2; % structural lattice constant in nm

% # of holes along x,y-axis 
N_xH=24; 
N_yH=28;

% keep N_x,N_y>=4
N_y=80*N_yH/a;
N_x=round(N_y/N_yH*N_xH*sqrt(3)/2)*2; 
rounding_error=abs(N_y/N_yH*N_xH*sqrt(3)-round(N_y/N_yH*N_xH*sqrt(3)/2)*2);

nx_A0=N_x/2;

% supercell size
Lx=N_x*a; 
L=N_y*a;
           
shift_A=round(L/N_yH*(2/sqrt(3)-sqrt(3)/2))/a;

% Heaviside broadening --> hole radius
eta_H=a*N_y/N_yH/4;

%Hamiltonian parameters
gamma=1; %Unit of energy
t=-gamma; %hopping inside the sample
E0=0*gamma; %on-site potential

% units:
me=9.10938356e-31;
h=6.62607015e-34;
hbar=1.05457e-34;
c=299792458;
nm=1e-9;
cm=1e-2;
eV=1.602176634e-19; 
aa=a*nm;
a_in_nm=aa/nm;
meff=0.03; % 2DEG effective mass (GaAs)
t_eV=hbar^2/(meff*me*2*aa^2)/eV;
%%%%%%%%%%%%%%%%%%%
Vhole_eV=5; % infinite potential strength on holes
Vhole=Vhole_eV/t_eV;
%%%%%%%%%%%%%% add holes, checkers board %%%%%%
dx=a;
dy=dx; % underlying TB grid is square

V_H_onTB=zeros(N_x*N_y,1);
%%%%%%%%%%%%%%%%% add holes potential pilars, use rectangular graphene uc
tau1=[-L/sqrt(3) 0]/N_yH;
tau2=[L/sqrt(3)/2  -L/2]/N_yH;
tau3=[L/sqrt(3)/2 L/2]/N_yH;

dxyH=[Lx/N_xH L/N_yH];

xyA=[L/N_yH/sqrt(3)/2-shift_A*a L/N_yH*3/4-a*1]; 

xyH=zeros(N_xH*N_yH,2);
nH=0;
for ixH=1:N_xH
    for iyH=1:N_yH
nH=nH+1;
xyH(nH,:)=xyA+[ixH-1 iyH-1].*dxyH;
    end
end

dxy=[dx dy];

xyTB=zeros(N_x*N_y,2);
nTB=0;
for ix=1:N_x
    for iy=1:N_y
nTB=nTB+1;
        xyTB(nTB,:)=dxy/2+[ix-1 iy-1].*dxy;
    end
end

xTB=xyTB(:,1);
yTB=xyTB(:,2);

nH=0;
for ixH=1:N_xH
    for iyH=1:N_yH
        nH=nH+1;
        for iTB=1:N_x*N_y
        d_TBtoH_A(iTB,1)=norm(xyH(nH,:)-xyTB(iTB,:));
        d_TBtoH_A2(iTB,1)=norm(xyH(nH,:)+[Lx 0]-xyTB(iTB,:));
        d_TBtoH_A3(iTB,1)=norm(xyH(nH,:)+[0 -L]-xyTB(iTB,:));
        d_TBtoH_C(iTB,1)=norm(xyH(nH,:)+tau2-tau1-xyTB(iTB,:));
        d_TBtoH_C2(iTB,1)=norm(xyH(nH,:)+tau2-tau1+[0 L]-xyTB(iTB,:));
        d_TBtoH_C3(iTB,1)=norm(xyH(nH,:)+tau2-tau1+[Lx 0]-xyTB(iTB,:));
        end
        V_H_onTB=V_H_onTB+Vhole*heaviside(eta_H-d_TBtoH_A)+Vhole* ...
                 heaviside(eta_H-d_TBtoH_A2)+Vhole*heaviside(eta_H-d_TBtoH_A3);
     V_H_onTB=V_H_onTB+Vhole*heaviside(eta_H-d_TBtoH_C)+Vhole*heaviside(eta_H-d_TBtoH_C2)+Vhole*heaviside(eta_H-d_TBtoH_C3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

midX=(min(xTB)+max(xTB))/2;
midY=(min(yTB)+max(yTB))/2;

Hsize=N_x*N_y; %Size of the full Hamiltonian matrix
H=0*sparse(speye(Hsize,Hsize)); %Fill the full Hamiltonian matrix with zeros
  
H_ij=sparse(t*diag(ones(1,N_y),0)); 

for n=2:N_x
    H(1+(n-1)*N_y:n*N_y,1+(n-2)*N_y:(n-1)*N_y)=H_ij; 
    H(1+(n-2)*N_y:(n-1)*N_y,1+(n-1)*N_y:n*N_y)=H_ij'; 
end

%%%%%%%%%%%%%%
fact_obc=13.3333; % insures B[T]=0.046 for the particular flake size
p=1;

B=2*pi*p/N_x/a/fact_obc; 
                         %units_eV
B_in_Tsla=h/(eV)/aa^2*p/N_x/fact_obc

B_T=B_in_Tsla;
for n=1:N_x
    H(1+(n-1)*N_y:n*N_y,1+(n-1)*N_y:n*N_y)=sparse(E0*diag(ones(1,N_y))+t*diag(ones(1,N_y-1),1)*exp(1i*(n-0.5-nx_A0)*B*a)+t*diag(ones(1,N_y-1),-1)*exp(-1i*(n-0.5-nx_A0)*B*a)); 
end 

H=H+diag(sparse(V_H_onTB));
%%%%%%%%

X_op=diag(sparse(xTB-midX));
Y_op=diag(sparse(yTB-midY));

Hsparse = sparse(H);
Xsparse = sparse(X_op);
Ysparse = sparse(Y_op);
II = speye(size(Hsparse));
%%%%%%%%%%%%%%%%%%% % spectral localizer
% particular choice of xy near the middle  of flake
ix=N_x/2+1;
xx=a/2+(ix-1)*a-midX;
iy=N_y/2+1-22;
yy=a/2+(iy-1)*a-midY;

sigma_x=sparse([0  1
         1  0]);
sigma_y=sparse([0 -1i
         1i 0]);
sigma_z=sparse([1  0
         0 -1]);

HH = Hsparse;
XX = Xsparse;
YY = Ysparse;

energyVec = linspace(4,8,161);
energyVecAct = -4.0+energyVec/(t_eV*1000);

indVec = zeros(length(energyVec),1);
Neigs=10; % needs to be large enough that matlab can perform eig
          % without NaN - due to high degeneracy of LLs
gapVec = zeros(length(energyVec),Neigs);

for ie = 1:length(energyVec)
            [ie length(energyVec)]
        LL = sparse(kappa*kron((XX-xx*speye(size(XX,1))),sigma_x) + kappa*kron((YY-yy*speye(size(YY,1))),sigma_y) + kron((HH-energyVecAct(ie)*speye(size(HH,1))),sigma_z));

    indVec(ie) = (1/2)*signatureComplex(LL);
    gapVec(ie,:) = eigs(LL,Neigs,'sm');

end
loc_gap=squeeze(gapVec(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot localizer results as in fig 4b

set(0,'defaultAxesFontSize',25)
figure
plot(energyVec,-indVec,'r-','LineWidth',2)
hold on
plot(energyVec,1E5*abs(loc_gap),'b-','LineWidth',2)
xlabel 'E [meV]'
ylabel 'Localizer'
xlim([4 8])
title(['V_h=\infty, sc24x28, \kappa=1.25E-7, a=2nm, B=',num2str(round(B_T*1000)/1000),'T'])
grid on
legend('Index','|Gap|x1E5')
ylim([-10 15])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functios needed for signature
function sig = signatureComplex(X)

%
% Matlab does not have a sparse  LDLT factorization
% in the comples case.  A hint from mathexchage was used:
% embed M_n(C)  in  M_{2n}(R) and use 1/2 of the index 
% computed there.
%
% This actually works for a dense input.
%

if ~ishermitian(X)
	error('Input is not hermitian');
end

X2 = [real(X), imag(X); -imag(X), real(X)];

sig = (signatureReal(X2))/2;

end

function sig = signatureReal(X)

% Uses  Sylvester's law of inertia
% and the sparse LDL^adjoint factorization
% that comes with Matlab.
% It is assumed that X is real symmetric.
% If there is not reasonable 
% gap in the spectrum at 0 the answer returned
% is basically random.
%
% This actually works for a dense input.

if ~isreal(X)
	error('Input is not stored using real numbers');
end
if ~issymmetric(X)
	error('Input is not hermitian/symmetric');
end

n = length(X);
[L,D,p]=ldl(X,'vector');
sig = 0;
j=1;
while j <= n
	if j<n && D(j,j+1) ~= 0
		sig = sig + sum( eig( [D(j,j),D(j,j+1);D(j+1,j),D(j+1,j+1)] ) > 0 );
		j = j+2;
	else
		sig = sig + (D(j,j) > 0);
		j = j+1;
	end
end
sig = sig - (n - sig);

end
