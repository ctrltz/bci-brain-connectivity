# Data
DATASET_URL = https://figshare.com/ndownloader/articles/13123148/versions/1
PREPROC_URL = https://osf.io/mnu6t/download

# Toolboxes
BBCI_URL = https://github.com/bbci/bbci_public/archive/refs/heads/master.zip
EEGLAB_URL = https://sccn.ucsd.edu/eeglab/download/daily/eeglab2021.0.zip
FASTICA_URL = https://research.ics.aalto.fi/ica/fastica/code/FastICA_2.5.zip
FOOOF_URL = https://github.com/fooof-tools/fooof_mat/archive/refs/heads/main.zip
HEADMODEL_URL = https://www.parralab.org/nyhead/sa_nyhead.mat
TPROD_URL = https://github.com/jadref/tprod/archive/refs/heads/main.zip

# Paths
MATLAB_ALIAS = MATLAB --version 9.13 matlab
MATLAB_SCRIPT = ./scripts/BCI_MI_analysis_main.m
R_ALIAS = R+ --version 4.2.2 R
R_SCRIPT = ./main.R
BACKUP_DIR = backup/backup_$(shell date +%Y%m%d_%H%M%S)

.PHONY: all configure analysis stats collect supplementary paper backup download_raw_data download_aux_data setup_toolboxes

all: configure analysis stats collect supplementary paper backup

download_raw_data:
	@echo "Downloading the dataset: $(DATASET_URL)"
	mkdir -p data/raw
	wget $(DATASET_URL) -P ./data/raw/
	unzip ./data/raw/1 -d ./data/raw/
	rm -f ./data/raw/1

download_aux_data:
	@echo "Downloading the preprocessing info: $(PREPROC_URL)"
	mkdir -p data/aux/
	cd ./data/aux/ && curl -O -J -L $(PREPROC_URL)
	
setup_toolboxes: setup_bbci setup_eeglab setup_fastica setup_fooof_mat setup_nyhead setup_tprod

setup_bbci:
	@echo "Setting up BBCI"
	wget $(BBCI_URL) -P ./toolboxes
	unzip ./toolboxes/master.zip -d ./toolboxes/
	mv ./toolboxes/bbci_public-master ./toolboxes/bbci_public/
	rm -f ./toolboxes/master.zip

setup_eeglab:
	@echo "Setting up EEGLAB 2021.0"
	wget $(EEGLAB_URL) -P ./toolboxes
	unzip ./toolboxes/eeglab2021.0.zip -d ./toolboxes/
	rm -f ./toolboxes/eeglab2021.0.zip

setup_fastica:
	@echo "Setting up FastICA 2.5"
	wget $(FASTICA_URL) -P ./toolboxes
	unzip ./toolboxes/FastICA_2.5.zip -d ./toolboxes/
	rm -f ./toolboxes/FastICA_2.5.zip
	
setup_fooof_mat:
	@echo "Setting up the MATLAB wrapper for FOOOF"
	wget $(FOOOF_URL) -P ./toolboxes
	unzip ./toolboxes/main.zip -d ./toolboxes/
	mv ./toolboxes/fooof_mat-main ./toolboxes/fooof_mat/
	rm -f ./toolboxes/main.zip
	
setup_nyhead:
	@echo "Downloading the pre-computed head model"
	wget $(HEADMODEL_URL) -P ./toolboxes/haufe/	
	
setup_tprod:
	@echo "Setting up tprod"
	wget $(TPROD_URL) -P ./toolboxes
	unzip ./toolboxes/main.zip -d ./toolboxes/
	mv ./toolboxes/tprod-main ./toolboxes/tprod/
	rm -f ./toolboxes/main.zip
	cp -a ./toolboxes/tprod_compiled/. ./toolboxes/tprod/

configure:
	mkdir -p results/

analysis:
	${MATLAB_ALIAS} -nodisplay -nosplash -nodesktop -r "run('${MATLAB_SCRIPT}');exit;"

stats:
	cd stats && ${R_ALIAS} --file=${R_SCRIPT}

collect:
	# Collect all generated figures
	mkdir -p paper/figures
	cp \
	results/stats/overview/fig1-study-overview.png \
	results/stats/overview/fig2-pipeline-overview.png \
	results/stats/metadata/fig3-accuracy.png \
	results/task1/group_csp_patterns/fig4-group-csp-source-space-mask.png \
	results/stats/multiverse_SNR/fig5-snr-rest-accuracy-session.png \
	results/stats/multiverse_SNR_connectivity/fig6-snr-connectivity.png \
	results/stats/multiverse_SNR_connectivity/fig6supp-connectivity-spectra.png \
	results/stats/multiverse_connectivity_performance/fig7-multiverse-connectivity-performance-within.png \
	results/stats/multiverse_connectivity_performance/fig7supp1-multiverse-connectivity-performance-between-subject.png \
	results/stats/multiverse_connectivity_longitude/fig7supp2-multiverse-connectivity-longitude.png \
	results/stats/multiverse_joint/fig8-multiverse-joint-analysis.png \
    results/stats/pipeline_effects/fig9-pipeline-effects-highlights.png \
	results/stats/pipeline_effects/fig9supp1-snr-band-inverse-interaction.png \
	results/stats/pipeline_effects/fig9supp2-coh-within-mask-roi-method-interaction.png \
	results/stats/pipeline_effects/fig9supp3-snr-inverse-band-interaction.png \
	results/stats/metadata/figAsupp1-task-comparison.png \
	results/stats/csp/figBsupp1-csp-accuracy-auc.png \
	results/stats/laplace_connectivity/figCsupp1-laplace-snr-connectivity.png \
	results/stats/laplace_connectivity/figCsupp2-laplace-connectivity-performance-within.png \
	results/stats/laplace_connectivity/figCsupp3-laplace-connectivity-longitude.png \
	paper/figures

	# Collect all results exported to TeX
	cp results/tex/* paper/numbers

paper: paper/main.tex
	cd paper && pdflatex main.tex && bibtex main.aux && pdflatex main.tex && pdflatex main.tex

supplementary: paper/supplementary.tex
	cd paper && pdflatex supplementary.tex && bibtex supplementary.aux && pdflatex supplementary.tex && pdflatex supplementary.tex

full: paper/main.pdf paper/supplementary.pdf
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=paper/full.pdf paper/main.pdf paper/supplementary.pdf

show:
	xpdf paper/full.pdf
	
sourcezip:
	rm -f source_files.zip
	cp -r paper/ paper_flat/
	cd paper_flat; \
	  cp figures/fig* .; \
	  rm -r -f figures/; \
	  mv main.tex main_old.tex; \
	  ../toolbox/latexpand/latexpand main_old.tex -o main.tex; \
	  mv supplementary.tex supplementary_old.tex; \
	  ../toolbox/latexpand/latexpand supplementary_old.tex -o supplementary.tex; \
	  zip -r ../source_files.zip main.bib main.tex supplementary.tex fig*
	cd .. & rm -r -f paper_flat

backup:
	echo "Creating backup in the folder: ${BACKUP_DIR}"
	mkdir -p ${BACKUP_DIR}
	cp -r \
	data/preproc/replication/ \
	data/r/ \
	results/ \
	assets/ \
	paper/figures/ \
	paper/numbers/ \
	paper/full.pdf \
	${BACKUP_DIR}
