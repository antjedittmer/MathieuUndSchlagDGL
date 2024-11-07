function  [ImagEigCorrected, ImagEigCorrectedNeg,buffer,Eig] = correctImagValuesEig(eP,buffer,T)
%[ImagEigSortN,ImagEigCorrected, ImagEigCorrectedNeg, buffer.Pos,buffer.Neg,cntN] = correctImagValues(tmp,tmpNeg,buffer,buffer.Neg,cntN,nuCSwitchVec,n,nu_C2,T,RealEig)
% correctImagValues korrigiert die Werte

Eig.Real = 1/T * log(abs(eP));
Eig.Imag= 1/T * atan(imag(eP)./real(eP));
Eig.ImagSort = sort(Eig.Imag);

% Imaginaeranteil kontinuierlich steigend oder fallend
tmp = Eig.ImagSort(2);
tmpNeg = Eig.ImagSort(1);


if buffer.Pos <= tmp || (abs(tmp) < 10^-5) % Wert uebernehmen
    ImagEigCorrected = tmp; % steigender pos. Wert
    ImagEigCorrectedNeg = tmpNeg;  % fallender neg. Wert
    buffer.Pos = 0;
    buffer.Neg = 0;
else  % 'korrigierten' Wert nehmen fuer kontinuierlichen Verlauf
    % Wird der 'else'-Zweig getriggert, ist der buffer
    % bei 90 Grad (1.57 rad). Der positive Imaginaerteil ist damit
    % die Summe aus 180 Grad und dem negativen Winkel, dessen
    % Wert von -90 Grad zu 0 Grad laeuft.
    ImagEigCorrected = 2*buffer.Pos + tmpNeg; % korrigierter pos. Wert
    ImagEigCorrectedNeg = 2*buffer.Neg + tmp; % korrigierter neg. Wert
    % vecSwitch = 1; %vecSwitch(oidx)= 1;

end

% buffer.Pos ueberschreiben mit aktuellem Wert
buffer.Pos = max(tmp,buffer.Pos); % Maximum speichern
buffer.Neg = min(tmpNeg,buffer.Neg); % Minimum speichern

% % Additionsterm n*2*pi/T
% if cntN <= length(nuCSwitchVec) && ... % Sicherheitscheck, damit n nicht groesser wird als die Anzahl der 'Switchstellen'
%         nu_C2 > nuCSwitchVec(cntN) && ... % n nur bei erreichen der Switchstellen umstellen
%         abs(RealEig(1) - RealEig(2))> eps % die Realteile muessen gleich sein
%     n =  n + 0.5;
%     cntN = cntN+1;
% end
%
% nAdd = n*2*pi/T;
% ImagEigSortN = [ImagEigCorrectedNeg, ImagEigCorrected] + [-nAdd,nAdd];