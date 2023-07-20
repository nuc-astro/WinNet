#!/usr/bin/gnuplot -persist
#
# Makes a sketch with explanation of winnet 
# evolution and trajectory modes
#

set t postscript eps color 
set o "params_illustration.eps"

set size 0.7
set xlab "time"
set ylab "temperature" offset 10,10 rotate by 0

set nokey
ymin= 2
ymax= 11
set xrange [0:32]
set yrange [ymin:ymax]

#set ytics ("10" 10.5, "1" 8, "0.1" 5.5, "0.01" 3)
set ytics ("" 8.2, "" 8, "hot2cold_temp" 5.8) font "Courier Bold,12" 
set label "nse2net_temp" at -0.7,7.9 right noenhanced font "Courier Bold,12"
set label "net2nse_temp" at -0.7,8.3 right noenhanced font "Courier Bold,12"

dstr =           " 0 10.3"
dstr = dstr."\n"." 1 10.0"
dstr = dstr."\n"." 2  9.2"
dstr = dstr."\n"." 3  8.5"
dstr = dstr."\n"." 4  7.8"
dstr = dstr."\n"." 5  7.7"
dstr = dstr."\n"." 6  7.9"
dstr = dstr."\n"." 7  8.1"
dstr = dstr."\n"." 8  8.3"
dstr = dstr."\n"." 9  8.2"
dstr = dstr."\n"."10  7.1"
dstr = dstr."\n"."11  6.2"
dstr = dstr."\n"."12  6.1"
dstr= "<(cat << EOF\n".dstr."\nEOF\n)"

t1= 12.0
expd =           "12  6.1"
expd = expd."\n"."32  4.2"
expd= "<(cat << EOF\n".expd."\nEOF\n)"

set obj 1 rect from 0.0,8.0 to 3.7,11.0 behind fillstyle noborder fc rgb "#ccffcc"
set obj 2 rect from 3.7,5.8 to 7.5, 8.2 behind fillstyle noborder fc rgb "#ffccff"
set obj 3 rect from 7.5,8.0 to 9.2,11.0 behind fillstyle noborder fc rgb "#ccffcc"
set obj 4 rect from 9.2,5.8 to 15.6,8.2 behind fillstyle noborder fc rgb "#ffccff"
set obj 5 rect from 15.6,5.8 to 32.0,0  behind fillstyle noborder fc rgb "#ccccff"

set label "trajectory\137mode =\n'from\137file'" at 6.5,7.2 center noenhanced font "Courier Bold,12"
set arrow 2 from 0,6.5 to t1,6.5 ls 1 lc 7 lw 0.6 heads

set label "trajectory\137mode =\n'expand'" at 22,7.2 center noenhanced font "Courier Bold,12" 
set arrow 3 from 32,6.5 to t1,6.5 ls 1 lc 7 lw 0.6 

set label "trajectory\137mode =\n'analytic'" at 6.1,4.2 center noenhanced font "Courier Bold,12" textcolor rgb "#0000cc"
set arrow 4 from 6.8,4.5 to 10,10.5/2 ls 1 lc 3 lw 0.6 

set label "evolution\137mode: EM_NSE"     at 14,10.5 left noenhanced font "Courier,12" textcolor rgb "#008833"
xx= 13.5
yy= 10.3
set arrow 5 from xx,yy to 2.5,9.5 ls 1 lc rgb "#008833" lw 0.6 
set arrow 6 from xx,yy to 8.4,9.5 ls 1 lc rgb "#008833" lw 0.6 
set arrow 7 from xx,yy to 16.5,yy ls 1 lc rgb "#008833" lw 0.6 nohead

set label "evolution\137mode: EM_NETHOT"  at 14,10.5 left noenhanced font "Courier,12" textcolor rgb "#662266" offset 0,-1
xx= 13.8
yy= 9.85
set arrow 8  from xx,yy to 5.6, 8.1 ls 1 lc rgb "#662266" lw 0.6 
set arrow 9  from xx,yy to 10.4,7.5 ls 1 lc rgb "#662266" lw 0.6 
set arrow 10 from xx,yy to 16.5,yy  ls 1 lc rgb "#662266" lw 0.6 nohead

set label "evolution\137mode: EM_NETCOLD" at 14,10.5 left noenhanced font "Courier,12" textcolor rgb "#003388" offset 0,-2
xx= 14.5
yy= 9.4
set arrow 11 from xx,yy to 16.5,4.5 ls 1 lc rgb "#003388" lw 0.6 
set arrow 12 from xx,yy to 16.5,yy  ls 1 lc rgb "#003388" lw 0.6 nohead 


set arrow 1 from t1,ymin to t1,ymax ls 3 lc 1 lw 0.6 nohead
p dstr w lp ls 1 lc 7 lw 2.5 pt 6 \
, expd w l  ls 2 lc 7 lw 2.5 \
, 105/(x+10) w l lt 2 lw 2.5 lc 3 \
, 8.0 ls 3 lc 7 lw 0.5, 8.2 ls 3 lc 7 lw 0.5, 5.8 ls 3 lc 7 lw 0.5

