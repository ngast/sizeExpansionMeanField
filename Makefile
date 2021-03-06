JUPYTER_NBCONVERT=jupyter-nbconvert-3.6
EXAMPLE1=Example1_MalwarePropagation.ipynb
EXAMPLE2=Example2_supermarketModel.ipynb
EXAMPLE3=Example3_nonStableSIR.ipynb

PAPER=sizeExpansionMeanField

all: figures compile_paper

compile_paper:
	cd paper && pdflatex $(PAPER)
	cd paper && bibtex $(PAPER)
	cd paper && pdflatex $(PAPER)
	cd paper && pdflatex $(PAPER)
	cd paper && pdflatex $(PAPER)

simulations/rmf_tool:
	cd simulations && git clone git@github.com:ngast/rmf_tool.git

simulations/output_pdfs:
	mkdir simulations/output_pdfs

figures: simulations/rmf_tool simulations/output_pdfs
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
