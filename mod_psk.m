function constellation = mod_psk(M, is_gray)
    constellation = zeros(1, M);
    if is_gray
        ordre_symboles = gray(M); 
        for k=1:M
            constellation(ordre_symboles(k)+1) = exp(1i*2*pi*k/M);
        end
    else
        for k=1:M
            constellation(k) = exp(1i*2*pi*(k-1)/M);
        end
    end
end