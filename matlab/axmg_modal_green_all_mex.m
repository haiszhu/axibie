function [vk, ve, Fn, An, dFn] = axmg_modal_green_all_mex(chim1, M)
chim1=double(chim1); M=double(M); Mp1=M+1;
vk=zeros(Mp1,1); ve=zeros(Mp1,1); Fn=zeros(Mp1,1); An=zeros(Mp1,1); dFn=zeros(Mp1,1);
mex_id_ = 'axmg_modal_green_all_r64(c i double[x], c i int64_t[x], c io double[x], c io double[x], c io double[x], c io double[x], c io double[x])';
[vk, ve, Fn, An, dFn] = AxiStokes3D_mex(mex_id_, chim1, M, vk, ve, Fn, An, dFn, 1, 1, Mp1, Mp1, Mp1, Mp1, Mp1);
end

% Vectorized chunkie carrier (mirror of utils/axmg_modal_green_qleg_half_miller_vec.m, same
% signature): qm=Q_{k-1/2}(chi), qmd=Q'_{k-1/2}, qmdd=Q''_{k-1/2}, k=0..M, chi=t+1, one column
% per point, t=chi-1 with shape [1,nt,ns].  Swap the .m name -> name_mex verbatim at call sites.
