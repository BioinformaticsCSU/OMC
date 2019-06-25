clear all
addpath('Datasets');
%% 1. Load Datesets
load Fdataset 
load Fdataset_Matrix_DrProA_DrugBank
load Fdataset_Matrix_DisProA_CTDomim             %The gold standard dataset
Wrr = drug;
Wdd = disease;
Wdr = didr;
Wrd = Wdr';
Wdit = Matrix_dis_pro_associ;
Wdrt = Matrix_dr_pro_associ;
Wtdi = Wdit';
Wtdr = Wdrt';
%% 2. OMC3 algorithm
alpha = 1;
beta = 10;
K = 10;
tol1 = 2*1e-3;
tol2 = 1*1e-5;
maxiter = 300;
[dn,dr] = size(Wdr);
P_TMat = Wdr;

row_no = find(sum(P_TMat, 2) == 0);
if isempty(row_no) == 0
    P_TMat_new1 = KNN_diseaseS(P_TMat, Wdd, K);          %KNN Preprocessing
    P_TMat_new = P_TMat_new1 + P_TMat;
else
    P_TMat_new =P_TMat;
end
T1 = [Wrr; Wtdr; P_TMat_new];
[t1, t2] = size(T1);
trIndex1 = double(T1 ~= 0);
[W1, iter1] = BNNR(alpha, beta, T1, trIndex1, tol1, tol2, maxiter, 0, 1);
M_ResultMat1 = W1((t1-dn+1):t1, 1:dr);

col_no = find(sum(P_TMat, 1) == 0);
if isempty(col_no) == 0
    P_TMat_new2 = KNN_drugS(P_TMat, Wrr, K);             %KNN Preprocessing
    P_TMat_new = P_TMat_new2 + P_TMat;
else
    P_TMat_new = P_TMat;
end
T2 = [P_TMat_new, Wdit, Wdd];
[t_1, t_2] = size(T2);
trIndex2 = double(T2 ~= 0);
[W2, iter2] = BNNR(alpha, beta, T2, trIndex2, tol1, tol2, maxiter, 0, 1);
M_ResultMat2 = W2(1:dn, 1:dr);

M_recovery = (M_ResultMat1 + M_ResultMat2) / 2;     

