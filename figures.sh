#!/bin/bash

source figure_parameters.src

source ../pipeline.sh

actual_seqlen=()
actual_mutrate=()
actual_recombrate=()
for seqlen_index in "${seqlen_index_array[@]}"; do        
	actual_seqlen+=(${seqlen_array[seqlen_index]})
	actual_mutrate+=(${theta_array[seqlen_index]})
	actual_recombrate+=(${rho_array[seqlen_index]})
done
echo ${actual_seqlen[*]} > actual_seqlen_tmp

actual_particles=()
for particle_index in "${particle_index_array[@]}"; do
	actual_particles+=(${particle_array[particle_index]})
done
echo ${actual_particles[*]} > actual_particles_tmp

rm -r *_ErrFiles
rm -r *_OutFiles
rm Particle*.sh

#if [ ! -d Particle1000Seqlen30000000 ]; then
python movie.py ${suffix} ${scalingNe} ${model}
#fi

echo "
\documentclass[12pt,landscape]{article}
\usepackage[a2paper]{geometry}
\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{tikz}

%\newcommand{\Iteration}{16}
\begin{document}

{\LARGE \bf 
smcsmc Parameters: \\\\
Data source: ${DATA} \\\\
Time interval pattern: ${pattern_tmax} \\\\
Number of replicates: ${number_of_replicates} \\\\
EM iterations = ${EM} \\\\
Number of particles: $(echo ${actual_particles[*]} | tr ' ' ,) \\\\
Sequence lengths : $(echo ${actual_seqlen[*]} | tr ' ' ,)  \\\\
Mutation rates: $(echo ${actual_mutrate[*]} | tr ' ' ,)  \\\\
Recombination rates: $(echo ${actual_recombrate[*]} | tr ' ' ,)  \\\\
%#Actual sequence used: $(echo ${seqlen_index_array[*]} | tr ' ' ,)  \\\\
%#Actual number of particle used: $(echo ${particle_index_array[*]} | tr ' ' ,)  \\\\

}
\newpage

{\LARGE \bf Relative deviation in parameter estimation\\\\
From top to bottom, particles are $(echo ${actual_particles[*]} | tr ' ' ,) \\\\
From left to right, sequence lengths are $(echo ${actual_seqlen[*]} | tr ' ' ,)
}
\begin{center}
\foreach \particlenumber in {$(echo ${actual_particles[*]} | tr ' ' ,)} {
    \foreach \seqlen in {$(echo ${actual_seqlen[*]} | tr ' ' ,)} {
        \includegraphics[width=$(echo " 1/${#seqlen_index_array[@]} -0.05" | bc -l)\textwidth]{Particle\particlenumber Seqlen\seqlen /Particle\particlenumber Seqlen\seqlen _dev.png}
    } \\\\
}
\end{center}

\newpage

{\LARGE \bf Recombination rate\\\\
From top to bottom, particles are ${particles} \\\\
From left to right, sequence lengths are ${seqlen}
}
\begin{center}
\foreach \particlenumber in {$(echo ${actual_particles[*]} | tr ' ' ,)} {
    \foreach \seqlen in {$(echo ${actual_seqlen[*]} | tr ' ' ,)} {
        \includegraphics[width=$(echo " 1/${#seqlen_index_array[@]} -0.05" | bc -l)\textwidth]{Particle\particlenumber Seqlen\seqlen /Particle\particlenumber Seqlen\seqlen RE.png}
    } \\\\
}
\end{center}

\foreach \Iteration in ${iteration} {

\newpage
{\LARGE \bf Iteration \Iteration \\\\
From top to bottom, particles are $(echo ${actual_particles[*]} | tr ' ' ,) \\\\
From left to right, sequence lengths are $(echo ${actual_seqlen[*]} | tr ' ' ,)
}
\begin{center}
\foreach \particlenumber in {$(echo ${actual_particles[*]} | tr ' ' ,)} {
    \foreach \seqlen in {$(echo ${actual_seqlen[*]} | tr ' ' ,)} {
        \includegraphics[width=$(echo " 1/${#seqlen_index_array[@]} -0.05" | bc -l)\textwidth]{Particle\particlenumber Seqlen\seqlen /Particle\particlenumber Seqlen\seqlen step\Iteration.png}
    } \\\\
}
\end{center}
}
\end{document}
" > ${figurefile}.tex

pdflatex ${figurefile}.tex

R CMD BATCH run_time.r
