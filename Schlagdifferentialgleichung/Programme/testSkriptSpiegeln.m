load("Workspace.mat");

xachse = MuMin:SW:MuMax;

CharExRe1 = CharExRe(:,1);
CharExRe2 = CharExRe(:,2);

CharExRePos = max(CharExRe1,CharExRe2); % Positiver Eigenwert: Glatter Verlauf
CharExRe1_NegIdx =  abs(CharExRePos - CharExRe1) > eps; % Index: Negativer Wert 1. Realteil 
offset = CharExRe1(1); % Offset 
CharExReNeg =  - (CharExRePos - offset) + offset; % Negativer EW: Gespiegelter positiver EW

CharExRe1Cor = CharExRe1; % Korrigierter 1. Eigenwert
CharExRe1Cor(CharExRe1_NegIdx) = CharExReNeg(CharExRe1_NegIdx); % Ersetze negative Werte durch gespiegelte positive EW

CharExRe2Cor = CharExRe2; % Korrigierter 1. Eigenwert
CharExRe2Cor(~CharExRe1_NegIdx) = CharExReNeg(~CharExRe1_NegIdx); % Ersetze negative Werte durch gespiegelte positive EW

subplot(2,1,1)
plot(xachse,CharExRe1,xachse,CharExRe2,xachse, CharExReNeg,'g--',...
    xachse,CharExRe1Cor,'k-.');
legend('Re1','Re2',' Neg.Re-Wert', 'Korrigierter Wert','Location','eastoutside')

subplot(2,1,2)
plot(xachse,CharExRe1Cor, xachse,CharExRe2Cor)
legend('Re1 korrigiert','Re2 korrigiert','Location','eastoutside')
