function [cpR, hRT, sR, gRT] = nasa9(A_NASA, B_NASA, T_HI, temp)

    [NUM_SPECIES, MAXINT] = size(T_HI);

    k = min(sum(temp > T_HI, 2) + 1, MAXINT);
    rows = ((1:NUM_SPECIES).' - 1) * MAXINT + k;
 
    a = A_NASA(rows,:);
    b = B_NASA(rows,:);
 
    iTemp = 1/temp;
    iTemp2 = iTemp * iTemp;
    lnTemp = log(temp);
    temp2 = temp*temp;
    temp3 = temp2*temp;
    temp4 = temp3*temp;
 
    cpR = a * [iTemp2; iTemp; 1; temp; temp2; temp3; temp4 ];
    hRT = a * [-iTemp2; lnTemp * iTemp; 1; temp/2; temp2/3; temp3/4; temp4/5 ] + b(:,1) * iTemp;
    sR = a * [-iTemp2/2; -iTemp; lnTemp; temp; temp2/2; temp3/3; temp4/4 ] + b(:,2);

    gRT = hRT - sR;
    
end