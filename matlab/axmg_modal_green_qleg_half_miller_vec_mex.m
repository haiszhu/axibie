function [qm, qmd, qmdd, pm, pmd, pmdd] = axmg_modal_green_qleg_half_miller_vec_mex(nc, t, M, qm, qmd, qmdd, pm, pmd, pmdd)
M=double(M); Mp1=M+1; nt=size(t,2); ns=size(t,3);
t=double(t(:)); nc=numel(t);
qm=zeros(Mp1,nc); qmd=zeros(Mp1,nc); qmdd=zeros(Mp1,nc);
pm=zeros(Mp1,nc); pmd=zeros(Mp1,nc); pmdd=zeros(Mp1,nc);
mex_id_ = 'axmg_modal_green_qleg_half_miller_vec_r64(c i int64_t[x], c i double[x], c i int64_t[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx])';
[qm, qmd, qmdd, pm, pmd, pmdd] = AxiStokes3D_mex(mex_id_, nc, t, M, qm, qmd, qmdd, pm, pmd, pmdd, 1, nc, 1, Mp1, nc, Mp1, nc, Mp1, nc, Mp1, nc, Mp1, nc, Mp1, nc);
qm=reshape(qm,Mp1,nt,ns); qmd=reshape(qmd,Mp1,nt,ns); qmdd=reshape(qmdd,Mp1,nt,ns);
pm=reshape(pm,Mp1,nt,ns); pmd=reshape(pmd,Mp1,nt,ns); pmdd=reshape(pmdd,Mp1,nt,ns);
end

% Helsing close-eval special quadrature: log(As), Cauchy(Ad), deriv(A1,A2), hyper(A3,A4)
