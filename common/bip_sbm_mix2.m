function [ A, label1, label2, X_1, X_2, v ] = bip_sbm_mix2( N, K, p, q, pr, d, mu ,Sig, sig2 )

    [A,label1,label2] = bip_sbm2(N,K,p,q,pr);

    v = zeros(K,d(1)+d(2));
    for k=1:K
        v(k,:) = mvnrnd(mu,Sig);
    end
    
    idx = 1:d(1);
    X_1 = (v(label1,idx)) +  sqrt(sig2(1))*randn(N(1),d(1));
    idx = d(1) + (1:d(2));
    X_2 = v(label2,idx) +  sqrt(sig2(2))*randn(N(2),d(2));


end