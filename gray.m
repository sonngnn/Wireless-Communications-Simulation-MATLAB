function ordre_symboles = gray(M)
    if M == 2
        ordre_symboles = [0 1];
        return
    else
        ordre_symboles = [gray(M/2) M/2+flip(gray(M/2))];
    end
end