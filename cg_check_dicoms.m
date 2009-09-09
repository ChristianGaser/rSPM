function cg_check_dicoms

P = spm_select(inf,'.*','Select Dicom files');
n = size(P,1);

for i=1:n
  hdr = spm_dicom_headers(P(i,:));
  m = floor(rem(hdr{1}.StudyTime/60,60));
  h = floor(hdr{1}.StudyTime/3600);
  fprintf('%s\t',hdr{1}.Filename);
  fprintf('%s\t',hdr{1}.PatientsName);
  if isfield(hdr{1},'StudyComments')
    fprintf('%s\t',hdr{1}.StudyComments);
  else
    fprintf(' \t');
  end
  fprintf('%s %02d:%02d\n',datestr(hdr{1}.StudyDate,'yyyy-mm-dd'),h,m);
end