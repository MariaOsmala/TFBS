pfms_tab_path="../../Data/PWMs/pfms_tab"
pfms_space_path="../../Data/PWMs/pfms_space"
pfms_transfac_path="../../Data/PWMs/pfms_transfac"
pfms_scpd="../../Data/PWMs/pfms_scpd" #this can be converted to meme format

pwms_tab_path="../../Data/PWMs/pwms_tab"
pwms_space_path="../../Data/PWMs/pwms_space"


dir.create(pfms_tab_path, recursive=TRUE)
dir.create(pfms_space_path, recursive=TRUE)
dir.create(pfms_transfac_path, recursive=TRUE)
dir.create(pfms_scpd, recursive=TRUE)
dir.create(pwms_tab_path, recursive=TRUE)
dir.create(pwms_space_path, recursive=TRUE)

rev_comp_pfms_tab_path="../../Data/PWMs-reverse-complement/pfms_tab"
rev_comp_pfms_space_path="../../Data/PWMs-reverse-complement/pfms_space"

dir.create(rev_comp_pfms_tab_path, recursive=TRUE)
dir.create(rev_comp_pfms_space_path, recursive=TRUE)
dir.create("../../Data/PWMs-reverse-complement/pfms_transfac", recursive=TRUE) #Is this needed?

logos_png_prob="../../Data/PWMs/Logos/png/prob"
logos_png_ic="../../Data/PWMs/Logos/png/ic"
logos_pdf_prob="../../Data/PWMs/Logos/pdf/prob"
logos_pdf_ic="../../Data/PWMs/Logos/pdf/ic"

dir.create(logos_png_prob, recursive=TRUE)
dir.create(logos_png_ic, recursive=TRUE)
dir.create(logos_pdf_prob, recursive=TRUE)
dir.create(logos_pdf_ic, recursive=TRUE)

revcomp_logos_png_prob="../../Data/PWMs-reverse-complement/Logos/png/prob"
revcomp_logos_png_ic="../../Data/PWMs-reverse-complement/Logos/png/ic"
revcomp_logos_pdf_prob="../../Data/PWMs-reverse-complement/Logos/pdf/prob"
revcomp_logos_pdf_ic="../../Data/PWMs-reverse-complement/Logos/pdf/ic"

dir.create(revcomp_logos_png_prob, recursive=TRUE)
dir.create(revcomp_logos_png_ic, recursive=TRUE)
dir.create(revcomp_logos_pdf_prob, recursive=TRUE)
dir.create(revcomp_logos_pdf_ic, recursive=TRUE)

# Scrambled motif paths ---------------------------------------------------

scrambled_pfms_tab_path="../../Data/PWMs/scrambled_pfms_tab"
scrambled_pfms_space_path="../../Data/PWMs/scrambled_pfms_space"
scrambled_pfms_transfac_path="../../Data/PWMs/scrambled_pfms_tranfac"
scrambled_pfms_scpd="../../Data/PWMs/scrambled_pfms_scpd" #this can be converted to meme format

scrambled_pwms_tab_path="../../Data/PWMs/scrambled_pwms_tab"
scrambled_pwms_space_path="../../Data/PWMs/scrambled_pwms_space"


dir.create(scrambled_pfms_tab_path, recursive=TRUE)
dir.create(scrambled_pfms_space_path, recursive=TRUE)
dir.create(scrambled_pfms_transfac_path, recursive=TRUE)
dir.create(scrambled_pfms_scpd, recursive=TRUE)
dir.create(scrambled_pwms_tab_path, recursive=TRUE)
dir.create(scrambled_pwms_space_path, recursive=TRUE)
