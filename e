reset
execute getrun
execute "$!run/param"
data "$!run/evt"

read {t 1 E 2 Z 3}

if($?tau==0) {define tau 0}
#define tau ? {tau?}

windowavg t E tavg Eavg $tau
windowavg t Z tavg Zavg $tau

set tavg=tavg if (Eavg > 0)
set Zavg=Zavg if (Eavg > 0)
set Eavg=Eavg if (Eavg > 0)

set all=Eavg concat Zavg
autology tavg all
#xlimits 0 60
box
xlabel "t"
ylabel "Quadratic Quantities"
lctype 0 red
connectlabel "E" tavg (lg(Eavg))
lctype 0 green
connectlabel "Z" tavg (lg(Zavg))

ctitle green "2D Turbulence"

id


