function [CP_R, H_RT, S_R, G_RT] = nasa9(A_NASA, B_NASA, THI, T)

    [NUM_SPECIES, MAXINT] = size(THI);

    k = min(sum(T > THI, 2) + 1, MAXINT);
    rows = ((1:NUM_SPECIES).' - 1)*MAXINT + k;
 
    a = A_NASA(rows,:);
    b = B_NASA(rows,:);
 
    iT = 1/T; iT2 = iT*iT; lnT = log(T);
    T2 = T*T; T3 = T2*T; T4 = T3*T;
 
    CP_R = a * [iT2; iT; 1; T; T2; T3; T4 ];
    H_RT = a * [-iT2; lnT*iT; 1; T/2; T2/3; T3/4; T4/5 ] + b(:,1)*iT;
    S_R = a * [-iT2/2; -iT; lnT; T; T2/2; T3/3; T4/4 ] + b(:,2);

    G_RT = H_RT - S_R;
    
end