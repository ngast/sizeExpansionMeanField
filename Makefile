JUPYTER_NBCONVERT=jupyter-nbconvert-3.6
EXAMPLE1=Example1_MalwarePropagation.ipynb
EXAMPLE2=Example2_supermarketModel.ipynb
EXAMPLE3=Example3_nonStableSIR.ipynb

PAPER=refined_mf_perf2018
BIBFILE=biblio

all: figures compile_paper

compile_paper:
	cd paper && pdflatex $(PAPER)
	cd paper && bibtex $(BIBFILE)
	cd paper && pdflatex $(PAPER)
	cd paper && pdflatex $(PAPER)
	cd paper && pdflatex $(PAPER)

simulations/rmf_tool:
	cd simulations && git clone git@github.com:ngast/rmf_tool.git

figures: simulations/rmf_tool
	echo 'Computing figures for EXAMPLE 1' 
	cd simulations && $(JUPYTER_NBCONVERT) --execute $(EXAMPLE1)
	cd simulations && $(JUPYTER_NBCONVERT) --execute $(EXAMPLE2)
	cd simulations && $(JUPYTER_NBCONVERT) --execute $(EXAMPLE3)

remove_precomputed_simulations:
	echo "This code removes all results of the stochastic simulator"
	echo "recomputing all these results might take several hours"
	echo "Do you want to continue ? [Y/n]" && read ans && [ $$ans == Y ]
	ls simulations/rmf_tool/misc/*/steadyState/*
	ls simulations/rmf_tool/misc/*/traj/*
