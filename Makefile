#!/usr/bin/env make -f
#
# $Id: Makefile 149 2009-09-02 05:34:26Z gaser $

VERSION=`svn info |grep Revision|sed -e 's/Revision: //g'`
DATE=`svn info |grep 'Last Changed Date: '|sed -e 's/Last Changed Date: //g'|cut -f1 -d' '`

TARGET=/Users/gaser/spm/spm8/toolbox/rSPM

STARGET=141.35.200.101:/Applications/xampp/htdocs/

FILES=cg_rSPM_defaults.m cg_write_jacdet.m cg_calc_jacdet.m cg_volume_paxinos.m cg_preprocess_rats.m cg_avg.m spm_rSPM.m rSPM.man spm_orthviews.m spm_sections.m spm_image.m spm_dicom_convert.m Paxinos_label.txt INSTALL.txt cg_check_dicoms.m cg_warp.* cg_hdw.m tbx_cfg_rspm.m bb.m Contents.m Changes Howto.txt Brainmask-Paxinos.nii T2-Paxinos-avg36.nii Paxinos_labeled.nii Ref0.4mm.nii

ZIPFILE=rSPM_$(VERSION).zip

install:
	-@echo install
	-@test ! -d ${TARGET} || rm -r ${TARGET}
	-@mkdir ${TARGET}
	-@cp ${FILES} ${TARGET}

update:
	-@svn update
	-@echo '% __________________________________________________________________________' > Contents.m
	-@echo '% Rat SPM Toolbox' >> Contents.m
	-@echo '% Version ' ${VERSION} ' (VBM8) ' ${DATE} >> Contents.m
	-@cat Contents_info.m >> Contents.m

help:
	-@echo Available commands:
	-@echo install zip scp upgrade

zip: update
	-@echo zip
	-@test ! -d rSPM || rm -r rSPM
	-@cp -rp ${TARGET} .
	-@zip ${ZIPFILE} -rm rSPM

scp: zip
	-@echo scp http://dbm.neuro.uni-jena.de/${ZIPFILE}
	-@scp -pr ${ZIPFILE} ${STARGET}
