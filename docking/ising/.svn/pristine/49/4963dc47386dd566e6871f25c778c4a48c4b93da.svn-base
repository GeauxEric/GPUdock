#!/bin/sh

nbeta=416

input=out
plotscript=plot.plot


if [ $# -ne 1 ]; then
    echo "Usage: `basename $0` [ #beta | allbeta | E | M | U]"
    exit 1
else

    # plot avrg status for all betas
    if [ "$1" == "E" ]; then
	dir=data`date '+%y%m%d_%H%M%S'`_all
	plotdata=$dir/plot.dat
	mkdir $dir
	touch $plotdata

	for i in {0..415}
	do
	    grep  "b$i$" $input | tail -n 1 >> $plotdata
	done

	cat > $dir/$plotscript << EOF
set xlabel "#beta"
set ylabel "lattice status"
set title "ising model simulation, average (after equilibrium) of all betas"
plot\
"plot.dat" using 8:4 title "E" with lines 1
pause -1
EOF
    cd $dir
    gnuplot $plotscript





    # plot avrg status for all betas
    elif [ "$1" == "M" ]; then
	dir=data`date '+%y%m%d_%H%M%S'`_all
	plotdata=$dir/plot.dat
	mkdir $dir
	touch $plotdata

	for i in {0..415}
	do
	    grep  "b$i$" $input | tail -n 1 >> $plotdata
	done

	cat > $dir/$plotscript << EOF
set xlabel "#beta"
set ylabel "lattice status"
set title "ising model simulation, average (after equilibrium) of all betas"
plot\
"plot.dat" using 8:5 title "M" with lines 2
pause -1
EOF
    cd $dir
    gnuplot $plotscript




    # plot avrg status for all betas
    elif [ "$1" == "U" ]; then
	dir=data`date '+%y%m%d_%H%M%S'`_all
	plotdata=$dir/plot.dat
	mkdir $dir
	touch $plotdata

	for i in {0..415}
	do
	    grep  "b$i$" $input | tail -n 1 >> $plotdata
	done

	cat > $dir/$plotscript << EOF
set xlabel "#beta"
set ylabel "lattice status"
set title "ising model simulation, average (after equilibrium) of all betas"
plot\
"plot.dat" using 8:6 title "U" with lines 1,\
"plot.dat" using 8:6 title "U" with dots 3
pause -1
EOF
    cd $dir
    gnuplot $plotscript







    # ....
    elif [ "$1" == "allbeta" ]; then
	dir=data`date '+%y%m%d_%H%M%S'`_$1
	plotdata=$dir/plot.dat
	mkdir $dir
	touch $plotdata

	for i in {0..5}
	do

	grep  "b$i$" $input >> $plotdata
	beta=`head -n 1 $plotdata | cut --characters=15-24`
	betai=`head -n 1 $plotdata | cut -c 4`

	cat > $dir/$plotscript << EOF
set xlabel "records"
set ylabel "lattice status"
set title "ising model simulation, beta ($betai) = $beta"
plot\
"plot.dat" using 1:4 title "E" with lines 1,\
"plot.dat" using 1:5 title "M" with lines 2
pause -1
EOF
	cd $dir
	gnuplot $plotscript &
	done





    # plot all recorded status for a given #beta
    # argv must be a number within range.
    # however, there is no error detection and handling. BUGGY!!!

#    elif [ "$1" within range of NBETA ]; then
    else
	dir=data`date '+%y%m%d_%H%M%S'`_$1
	plotdata=$dir/plot.dat
	mkdir $dir
	touch $plotdata

	grep  "b$1$" $input >> $plotdata
	beta=`head -n 1 $plotdata | cut --characters=15-24`
	betai=`head -n 1 $plotdata | cut -c 4`

	cat > $dir/$plotscript << EOF
set xlabel "records"
set ylabel "lattice status"
set title "ising model simulation, beta ($betai) = $beta"
plot\
"plot.dat" using 1:4 title "E" with lines 1,\
"plot.dat" using 1:5 title "M" with lines 2
pause -1
EOF
    cd $dir
    gnuplot $plotscript

#    else
#    echo "Usage: `basename $0` [ #beta | allbeta | average]"


    fi


fi


