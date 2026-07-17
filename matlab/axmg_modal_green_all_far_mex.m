function [vk, ve] = axmg_modal_green_all_far_mex(chi, M)
chi=double(chi); M=double(M); Mp1=M+1;
vk=zeros(Mp1,1); ve=zeros(Mp1,1);
mex_id_ = 'axmg_modal_green_all_far_r64(c i double[x], c i int64_t[x], c io double[x], c io double[x])';
[vk, ve] = AxiStokes3D_mex(mex_id_, chi, M, vk, ve, 1, 1, Mp1, Mp1);
end

% All-modes full carrier (chim1=chi-1): Q-side vk=Q_{n-1/2}, ve, and P-side Fn=P_{n-1/2}, An,
% dFn=dP_{n-1/2}/dchi, n=0..M  (scalar chim1).  == modal_green_all_r64.
