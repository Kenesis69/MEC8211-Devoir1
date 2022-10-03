set logscale x 10
set logscale y 10
set xlabel "Nombre de points"
set ylabel "Erreur relative"
set grid

plot "convergence_ordre2.dat" u 1:2 w lp linewidth 1.5 title "Ordre 2 - Norme L1",\
     "convergence_ordre2.dat" u 1:3 w lp linewidth 1.5 title "Ordre 2 - Norme L2",\
     "convergence_ordre2.dat" u 1:4 w lp linewidth 1.5 title "Ordre 2 - Norme Linf",\
     "convergence_ordre1.dat" u 1:2 w lp linewidth 1.5 title "Ordre 1 - Norme L1",\
     "convergence_ordre1.dat" u 1:3 w lp linewidth 1.5 title "Ordre 1 - Norme L2",\
     "convergence_ordre1.dat" u 1:4 w lp linewidth 1.5 title "Ordre 1 - Norme Linf"
