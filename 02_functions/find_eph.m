function icol = find_eph(Eph,sv,time)

%FIND_EPH  Finds the proper column in ephemeris array

%Kai Borre and C.C. Goad 11-26-96
%Copyright (c) by Kai Borre
%$Revision: 1.2 $  $Date: 2004/02/09  $


icol = 0;

% In den Ephemeridenfiles sind die Bahnelemente eines Satelitten für verschiedene Zeitpunkte
% angegeben. Zunächst wird gesucht in welchen Spalten Daten für den
% entsprechenden Satelliten stehen. 
isat = find(Eph(1,:) == sv);

%Anzahl der Spalten für den gesuchten Sat.
n = size(isat,2); 
if n == 0
   return
end;

% Das erste Ergebnis wird zunächst ausgewählt
icol = isat(1);
%Hier wird der Zeitunterschied des Bahnbezugs zur Beobachtung berechnet
dtmin = Eph(19,icol)-time;

% Anschl. wird der geeignete Bahnbogen nach folgenden Kriterien ausgewählt:
% Die Bahndefinition sollte vor der Beobachtung liegen ()
% wenn es davon mehrere gibt, wird der Bahnbogen mit dem geringsten Zeitabstand gewählt
for t = isat(2:length(isat))
    dt = Eph(19,t)-time;  %Zeitunterschied
    if abs(dt) < abs(dtmin)
        icol = t;
        dtmin = dt;
    end
end

%%%%%%%%%%%%  find_eph.m  %%%%%%%%%%%%%%%%%
