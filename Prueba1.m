

%DEFINICION MATERIAL.
E = 70e9; nu=0.3;
G = E/(2*(1+nu));
L = 0.5; LadoCuadrado = 1e-2; A = LadoCuadrado^2; I = LadoCuadrado^4/12;
Iy = I; Iz = I; 
J=0.141*LadoCuadrado^4;
m = 7;

%DEFINICION GEOMETRICA
gdln = 6; 

ne1 = 5; 
ne2 = 2*ne1;
ne = ne1+ne2; nn = ne+1; nP = ne1+1; 
Px = -0.7071; Py = -0.7071 ; Pz = 0; 
L1 = 0.5*L; 
L2 = L;
cx_n = zeros(nn,1); cy_n = cx_n; cz_n = cx_n;
for n=0:(ne1+1);
   cx_n(n+1) = n*(L1/ne1);
   cy_n(n+1) = L2;
   cz_n(n+1) = L; 
end
for n=0:(nn-(ne1+1));
   cy_n(ne1+1+n) = L2-n*(L2/ne2);
   cx_n(ne1+1+n) = L1;
   cz_n(ne1+1+n) = L;
end
c_n(:,1) = cx_n(:);c_n(:,2) = cy_n(:);c_n(:,3) = cz_n(:);

%% 

gdl = nn*gdln; 
elementos = zeros(ne,2); 
gdl_e = zeros(ne,2*gdln); 
L_e =  zeros(ne,1); 
T_e =  cell(ne,1);
for e = 1:ne 
    elementos(e,:) = [e e+1];
    index = elementos(e,:); 
    gdl_e(e,:)=[index(1)*gdln-5 index(1)*gdln-4 index(1)*gdln-3 ...         
        index(1)*gdln-2 index(1)*gdln-1 index(1)*gdln...
        index(2)*gdln-5 index(2)*gdln-4 index(2)*gdln-3 ...
        index(2)*gdln-2 index(2)*gdln-1 index(2)*gdln];
    x1 = cx_n(index(1)); x2 = cx_n(index(2));
    y1 = cy_n(index(1)); y2 = cy_n(index(2));
    z1 = cz_n(index(1)); z2 = cz_n(index(2));
    dx = x2-x1; dy = y2-y1; dz = z2-z1;
    Le = sqrt(dx^2+dy^2+dz^2); L_e(e) = Le;

   
    if x1 == x2 && y1 == y2
        if z2 > z1
          
            T1 = [0 0 1 ; 0 1 0 ; -1 0 0];
        else
           
            T1 = [0 0 -1 ; 0 1 0 ; 1 0 0];
        end
    else
        CXx = (x2-x1)/Le; CYx = (y2-y1)/Le; CZx = (z2-z1)/Le;
        LD = sqrt(CXx*CXx + CYx*CYx);
        CXy = -CYx/LD; CYy = CXx/LD; CZy = 0; CXz = -CXx*CZx/LD;
        CYz = -CYx*CZx/LD; CZz = LD;
        T1 = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];        
    end
    T0 = [0 0 0; 0 0 0;0 0 0];    
    T = [T1 T0 T0 T0;T0 T1 T0 T0;T0 T0 T1 T0; T0 T0 T0 T1];
    T_e{e}= T;
end
%%%%%%% MATRIZ de RIGIDEZ
K0 = zeros(gdl); 
Ke = zeros(2*gdln,2*gdln,ne); 
K0_e =  cell(ne,1);
for e = 1:ne
    gdle = gdl_e(e,:);
    Le = L_e(e);
    T = T_e{e};
    k1 = E*A/Le;
    k2 = 12*E*Iz/Le^3; k3 = 6*E*Iz/Le^2; k4 = 4*E*Iz/Le; k5 = 2*E*Iz/Le;
    k6 = 12*E*Iy/Le^3; k7 = 6*E*Iy/Le^2; k8 = 4*E*Iy/Le; k9 = 2*E*Iy/Le;
    k10 = G*J/Le;
    a =[k1 0 0; 0 k2 0; 0 0 k6]; b =[ 0 0 0;0 0 k3; 0 -k7 0];
    c =[k10 0 0;0 k8 0; 0 0 k4]; d =[-k10 0 0;0 k9 0;0 0 k5];
    k0_viga = [a b -a b;b' c (-b)' d;(-a)' -b a -b;b' d' (-b)' c];
    Ke(:,:,e) = Ke(:,:,e)+k0_viga;
    k0 = T'*k0_viga*T; 
    K0_e{e}=round(k0);
    K0(gdle,gdle) = K0(gdle,gdle)+k0; 
end

Kcomparar = K0;
K0=round(K0);
U0 = zeros(gdl,1);
ubc = [1 2 3 4 5 6 gdl-5 gdl-4 gdl-3 gdl-2 gdl-1 gdl]; 
uL = setdiff((1:length(U0))',ubc); 
ux = (1:6:gdl); uy = (2:6:gdl); uz = (3:6:gdl); 
uyz = (4:6:gdl); uxz = (5:6:gdl); uxy = (6:6:gdl); 

F = zeros(gdl,1); 
F(nP*gdln-5) = Px; F(nP*gdln-4) = Py; F(nP*gdln-3) = Pz; 
F_L = F(uL);
KLL_0 = K0(uL,uL); 
Kbb_0 = K0(ubc,ubc);
U0_L = KLL_0\F_L;
U0(uL) = U0_L; 
F = K0*U0; 
U0bc = U0(ubc); U0x = U0(ux); U0y = U0(uy); U0z = U0(uz);U0xy = U0(uxy);
U_e = zeros(ne,gdln*2); 
F_e = zeros(ne,gdln*2);
St = zeros(ne,gdln); 
for e=1:ne 
    gdle = gdl_e(e,:); Le = L_e(e); T = T_e{e};
    U_e(e,:) = U0(gdle)'; 
    U_e(e,:) = (T*U_e(e,:)')'; 
    F_e(e,:) = (Ke(:,:,e)*U_e(e,:)')';
    St(e,1) = F_e(e,7); 

    St(e,2) = F_e(e,8);
    St(e,3) = F_e(e,9); 
    St(e,4) = F_e(e,10); 
    St(e,5) = F_e(e,11)+F_e(e,2)*Le; 
    St(e,6) = F_e(e,12)+F_e(e,3)*Le; 

end


%%%%%%%% MATRIZ de RIGIDEZ GEOMETRICA
KG = zeros(gdl); 
KGe = zeros(2*gdln,2*gdln,ne); 
for e=1:ne
    gdle = gdl_e(e,:); 
    Le = L_e(e); T = T_e{e}; 
    Fx = St(e,1); Vy = St(e,2); Vz = St(e,3); Mx = St(e,4);
     My1 = F_e(e,5)-F_e(e,8)*Le; My2 =F_e(e,11)-F_e(e,2)*Le;
    Mz1 = F_e(e,6)-F_e(e,9)*Le; Mz2 = F_e(e,12)-F_e(e,3)*Le; 
    Dg1 = [0 0 0 0 0 0;...
        0 6*Fx/(5*Le) 0 -My1/Le 0 Fx/10;...
        0 0 6*Fx/(5*Le) -Mz1/Le -Fx/10 0;...
        0 -My1/Le -Mz1/Le Fx*J/(Le*A) -Vy*Le/6 -Vz*Le/6;...
        0 0 -Fx/10 -Vy*Le/6 Fx*2*Le/15 0;...
        0 Fx/10 0 -Vz*Le/6 0 Fx*2*Le/15];
    Dg2 = [0 0 0 0 0 0;...
        0 -6*Fx/(5*Le) 0 -My2/Le 0 Fx/10;...
        0 0 -6*Fx/(5*Le) -Mz2/Le -Fx/10 0;...
        0 My1/Le Mz1/Le -Fx*J/(Le*A) Vy*Le/6 Vz*Le/6;...
        0 0 Fx/10 Vy*Le/6 -Fx*Le/30 0;...
        0 -Fx/10 0 Vz*Le/6 0 -Fx*Le/30];
    Dg3 = [0 0 0 0 0 0;...
        0 -6*Fx/(5*Le) 0 My1/Le 0 -Fx/10;...
        0 0 -6*Fx/(5*Le) Mz1/Le Fx/10 0;...
        0 -My2/Le -Mz2/Le -Fx*J/(Le*A) Vy*Le/6 Vz*Le/6;...
        0 0 -Fx/10 Vy*Le/6 -Fx*Le/30 0;...
        0 Fx/10 0 Vz*Le/6 0 -Fx*Le/30];
    Dg4 = [0 0 0 0 0 0;...
        0 6*Fx/(5*Le) 0 My2/Le 0 -Fx/10;...
        0 0 6*Fx/(5*Le) Mz2/Le Fx/10 0;...
        0 -My2/Le Mz2/Le Fx*J/(Le*A) -Vy*Le/6 -Vz*Le/6;...
        0 0 Fx/10 -Vy*Le/6 Fx*2*Le/15 0;...
        0 -Fx/10 0 -Vz*Le/6 0 Fx*2*Le/15];
    kg_viga = [Dg1 Dg2;Dg3 Dg4];
    kg = T'*kg_viga*T; 
    KG(gdle,gdle) = KG(gdle,gdle)+kg;
    KGe(:,:,e) = KGe(:,:,e)+kg_viga; 
end

%%%%%% COMPUTACION de AUTOVALORES
%Matriz de rigidez reducida por las condiciones de contorno
KLL_G = -KG(uL,uL); KLL_0 = K0(uL,uL);
[eigenvec,eigenval] = eigs(KLL_0,KLL_G,m,'sm'); %autovales y autovectores
[Pc iPc] = sort(diag(eigenval));
UG = zeros(gdl,m); UG(uL,:) = eigenvec(:,:); UGbc = UG(ubc,:);
UGX = UG(1:6:(nn*gdln),:);
UGY = UG(2:6:(nn*gdln),:);
UGZ = UG(3:6:(nn*gdln),:);

% Representación de los modos de PANDEO
for n = 1:m;
    %desplazamientos de apoyo
    autovalor = Pc(n);
    UGbc = UG(ubc,iPc(n));
    UGx = UGX(:,iPc(n));
    UGy = UGY(:,iPc(n));
    UGz = UGZ(:,iPc(n));
    
end


%% Matrices de rigidez del pórtico (sin reducir)

K0_quest = round(K0,0);
KG_quest = round(KG,0);

%% Variables de salida
% Cargas críticas de los dos primeros modos de pandeo con 3 decimales (comando "round" de matlab)
Pc1_quest = Pc(1,1);
Pc2_quest = Pc(2,1);
