pfms_tab_path="../../Data/PWMs/pfms_tab"
pfms_space_path="../../Data/PWMs/pfms_space"
pfms_transfac_path="../../Data/PWMs/pfms_tranfac"
pfms_scpd="../../Data/PWMs/pfms_scpd" #this can be converted to meme format

pwms_tab_path="../../Data/PWMs/pwms_tab"
pwms_space_path="../../Data/PWMs/pwms_space"


dir.create(pfms_tab_path, recursive=TRUE)
dir.create(pfms_space_path, recursive=TRUE)
dir.create(pfms_transfac_path, recursive=TRUE)
dir.create(pfms_scpd, recursive=TRUE)
dir.create(pwms_tab_path, recursive=TRUE)
dir.create(pwms_space_path, recursive=TRUE)


dir.create("../../Data/PWMs-reverse-complement/pfms_tab", recursive=TRUE)
dir.create("../../Data/PWMs-reverse-complement/pfms_space", recursive=TRUE)
dir.create("../../Data/PWMs-reverse-complement/pfms_transfac", recursive=TRUE)

logos_png_prob="../../Data/PWMs/Logos/png/prob"
logos_png_ic="../../Data/PWMs/Logos/png/ic"
logos_pdf_prob="../../Data/PWMs/Logos/pdf/prob"
logos_pdf_ic="../../Data/PWMs/Logos/pdf/ic"

dir.create(logos_png_prob, recursive=TRUE)
dir.create(logos_png_ic, recursive=TRUE)
dir.create(logos_pdf_prob, recursive=TRUE)
dir.create(logos_pdf_ic, recursive=TRUE)