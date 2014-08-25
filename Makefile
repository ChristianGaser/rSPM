#!/usr/bin/env make -f
#
# $Id$

VERSION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm8/toolbox/rSPM

STARGET=dbm.neuro.uni-jena.de:/Applications/xampp/htdocs/rSPM

FILES=cg_rSPM_defaults.m cg_rSPM_get_defaults.m cg_write_jacdet.m cg_calc_jacdet.m cg_volume_paxinos.m cg_preprocess_rats.m cg_avg.m spm_rSPM.m rSPM.man spm_dicom_convert.m Paxinos_label.txt INSTALL.txt cg_check_dicoms.m cg_boxplot.m cg_warp.* cg_hdw.m cg_correct_bias_rats.m cg_confplot_spm.m cg_rSPM_update.m tbx_cfg_rspm.m bb.m Contents.m Changes Howto.txt Brainmask-Paxinos-avg176.nii Brainmask-Paxinos-smaller.nii T2-Paxinos-avg176.nii Paxinos_labeled.nii Ref0.4mm.nii

ZIPFILE=rSPM_$(VERSION).zip

install:
	-@echo install
	-@test ! -d ${TARGET} || rm -r ${TARGET}
	-@mkdir ${TARGET}
	-@cp ${FILES} ${TARGET}

update:
	-@svn update
	-@echo '% Rat SPM Toolbox' > Contents.m
	-@echo '% Version ' ${VERSION} ' (rSPM) ' ${DATE} >> Contents.m
	-@cat Contents_info.m >> Contents.m
	-@echo '% __________________________________________________________________________' > rSPM.man
	-@echo '% Rat SPM Toolbox' >> rSPM.man
	-@echo '% Version ' ${VERSION} ' (rSPM) ' ${DATE} >> rSPM.man
	-@cat rSPM.txt >> rSPM.man

help:
	-@echo Available commands:
	-@echo install zip scp update

zip: update
	-@echo zip
	-@test ! -d rSPM || rm -r rSPM
	-@mkdir rSPM
	-@cp -rp ${FILES} rSPM
	-@zip ${ZIPFILE} -rm rSPM

scp: zip
	-@echo scp http://dbm.neuro.uni-jena.de/${ZIPFILE}
	-@scp -P 2222 ${ZIPFILE} ${STARGET}
