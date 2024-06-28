#R environment file

TFBS-evolutionary-conservation/RProject/renv.yml

########## SELEX motif collection ########################

TFBS/RProjects/TFBS/code/extract_motifs_from_excel_Jolma2013.R                           
TFBS/RProjects/TFBS/code/extract_motifs_from_excel_Jolma2015.R                           
TFBS/RProjects/TFBS/code/extract_motifs_from_excel_Nitta2015.R                           
TFBS/RProjects/TFBS/code/extract_motifs_Morgunova_2015.R 
TFBS/RProjects/TFBS/code/extract_motifs_from_excel_Yin2017.R

#TFBS/RProjects/TFBS/code/extract_motifs_from_excel_fromYimeng.R  
#TFBS/RProjects/TFBS/code/extract_motifs_from_excel_fromYimeng_version2.R  

TFBS/RProjects/TFBS/code/extract_motifs_from_excel_fromYimeng_union_version1_version2.R
#creates /projects/TFBS/PWMs_final_union/fromYimeng/metadata.csv"

TFBS/RProjects/TFBS/code/extract_motifs_from_excel_fromYimeng_version2.2.R
#creates /projects/TFBS/PWMs_final_version2.2/fromYimeng/metadata.csv"

TFBS/RProjects/TFBS/code/collect_all_motif_metadata.R
#creates PWMs_final_version2.2/metadata_3993_motifs.csv
#"../../PWMs_final_version2.2/filenames.csv"
#"../../PWMs_final_version2.2/motifnames.csv"


#Draw motif logos
TFBS/RProjects/TFBS/code/motifs_logos_better.R
#TFBS/RProjects/TFBS/code/motif_logos_better_version1_version2_union.R #run this
TFBS/RProjects/TFBS/code/motif_logos_better_version2.2.R

#Spacek
#~/spacek/Experiments/draw_logos.sh

#Artificial half-site logoes

~/projects/TFBS/RProjects/TFBS/code/extract_motifs_from_excel_fromYimeng_artificial_halfsites.R

#Then need to compute the tomtom similarities between there and all other motifs



########## Similarity scores ########################

# MOtifSTAtistic Software Suite v1.1 (MOSTA-SSTAT)  http://mosta.molgen.mpg.de/material/mosta.tar.gz 
# SSTAT similarity scores, computation
MOSTA-SSTAT/Experiments/run_SSTAT_batch_array.sh
MOSTA-SSTAT/Experiments/run_SSTAT_similarity_with_itself.sh

MOSTA-SSTAT/Experiments/run_SSTAT_union.sh   
MOSTA-SSTAT/Experiments/run_SSTAT_with_itself_union.sh


# Collect SSTAT similarity scores, processing, output compatible with the dominating set analysis

motif-clustering-Viestra-private/RProjects/extract_SSTAT.R

# Motifsimilarity scores

motifsimilarity-private/experiments/run_motifsimilarity_sbatch_array.sh
motifsimilarity-private/experiments/run_motifsimilarity_final_version2.2.sh
#motifsimilarity-private/experiments/run_motifsimilarity_fromYimeng_version2.sh
#motifsimilarity-private/experiments/run_motifsimilarity_final_union.sh

# Tomtom 
# meme -version 5.4.1

#Environment file    
#/scratch/project_2006203/motif-clustering-Viestra-private/motif-clustering-Vierstra.yml

#This is the true environment
channels:
- conda-forge
- bioconda

dependencies:
- python
- numpy=1.20
- jupyterlab
- pandas
- scipy
- pyjaspar
- seaborn
- matplotlib
- bedops
- meme
- pybigwig
- ucsc-fasize
- openssl

# packages in environment at /CSC_CONTAINER/miniconda/envs/env1:
#
# Name                    Version                   Build  Channel
# _libgcc_mutex             0.1                 conda_forge    conda-forge
# _openmp_mutex             4.5                       2_gnu    conda-forge
# alsa-lib                  1.2.3.2              h166bdaf_0    conda-forge
# anyio                     3.6.1            py39hf3d152e_0    conda-forge
# argon2-cffi               21.3.0             pyhd8ed1ab_0    conda-forge
# argon2-cffi-bindings      21.2.0           py39hb9d737c_2    conda-forge
# asttokens                 2.0.8              pyhd8ed1ab_0    conda-forge
# attrs                     22.1.0             pyh71513ae_1    conda-forge
# babel                     2.10.3             pyhd8ed1ab_0    conda-forge
# backcall                  0.2.0              pyh9f0ad1d_0    conda-forge
# backports                 1.0                        py_2    conda-forge
# backports.functools_lru_cache 1.6.4              pyhd8ed1ab_0    conda-forge
# beautifulsoup4            4.11.1             pyha770c72_0    conda-forge
# bedops                    2.4.41               h9f5acd7_0    bioconda
# biopython                 1.79             py39hb9d737c_2    conda-forge
# bleach                    5.0.1              pyhd8ed1ab_0    conda-forge
# brotli                    1.0.9                h166bdaf_7    conda-forge
# brotli-bin                1.0.9                h166bdaf_7    conda-forge
# brotlipy                  0.7.0           py39hb9d737c_1004    conda-forge
# bzip2                     1.0.8                h7f98852_4    conda-forge
# c-ares                    1.18.1               h7f98852_0    conda-forge
# ca-certificates           2022.6.15            ha878542_0    conda-forge
# certifi                   2022.6.15        py39hf3d152e_0    conda-forge
# cffi                      1.15.1           py39he91dace_0    conda-forge
# charset-normalizer        2.1.1              pyhd8ed1ab_0    conda-forge
# cryptography              37.0.1           py39h9ce1e76_0  
# curl                      7.83.1               h2283fc2_0    conda-forge
# cycler                    0.11.0             pyhd8ed1ab_0    conda-forge
# dbus                      1.13.6               h5008d03_3    conda-forge
# debugpy                   1.6.3            py39h5a03fae_0    conda-forge
# decorator                 5.1.1              pyhd8ed1ab_0    conda-forge
# defusedxml                0.7.1              pyhd8ed1ab_0    conda-forge
# entrypoints               0.4                pyhd8ed1ab_0    conda-forge
# executing                 0.10.0             pyhd8ed1ab_0    conda-forge
# expat                     2.4.8                h27087fc_0    conda-forge
# flit-core                 3.7.1              pyhd8ed1ab_0    conda-forge
# font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
# font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
# font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
# font-ttf-ubuntu           0.83                 hab24e00_0    conda-forge
# fontconfig                2.14.0               h8e229c2_0    conda-forge
# fonts-conda-ecosystem     1                             0    conda-forge
# fonts-conda-forge         1                             0    conda-forge
# fonttools                 4.37.1           py39hb9d737c_0    conda-forge
# freetype                  2.12.1               hca18f0e_0    conda-forge
# gettext                   0.19.8.1          h73d1719_1008    conda-forge
# ghostscript               9.54.0               h27087fc_2    conda-forge
# glib                      2.72.1               h6239696_0    conda-forge
# glib-tools                2.72.1               h6239696_0    conda-forge
# gst-plugins-base          1.20.2               hcf0ee16_0    conda-forge
# gstreamer                 1.20.3               hd4edc92_0    conda-forge
# icu                       69.1                 h9c3ff4c_0    conda-forge
# idna                      3.3                pyhd8ed1ab_0    conda-forge
# importlib-metadata        4.11.4           py39hf3d152e_0    conda-forge
# importlib_metadata        4.11.4               hd8ed1ab_0    conda-forge
# importlib_resources       5.9.0              pyhd8ed1ab_0    conda-forge
# ipykernel                 6.15.1             pyh210e3f2_0    conda-forge
# ipython                   8.4.0            py39hf3d152e_0    conda-forge
# ipython_genutils          0.2.0                      py_1    conda-forge
# jedi                      0.18.1           py39hf3d152e_1    conda-forge
# jinja2                    3.1.2              pyhd8ed1ab_1    conda-forge
# jpeg                      9e                   h166bdaf_2    conda-forge
# json5                     0.9.5              pyh9f0ad1d_0    conda-forge
# jsonschema                4.14.0             pyhd8ed1ab_0    conda-forge
# jupyter_client            7.3.4              pyhd8ed1ab_0    conda-forge
# jupyter_core              4.11.1           py39hf3d152e_0    conda-forge
# jupyter_server            1.18.1             pyhd8ed1ab_0    conda-forge
# jupyterlab                3.4.5              pyhd8ed1ab_0    conda-forge
# jupyterlab_pygments       0.2.2              pyhd8ed1ab_0    conda-forge
# jupyterlab_server         2.15.1             pyhd8ed1ab_0    conda-forge
# keyutils                  1.6.1                h166bdaf_0    conda-forge
# kiwisolver                1.4.4            py39hf939315_0    conda-forge
# krb5                      1.19.3               h08a2579_0    conda-forge
# lcms2                     2.12                 hddcbb42_0    conda-forge
# ld_impl_linux-64          2.36.1               hea4e1c9_2    conda-forge
# lerc                      4.0.0                h27087fc_0    conda-forge
# libblas                   3.9.0           16_linux64_openblas    conda-forge
# libbrotlicommon           1.0.9                h166bdaf_7    conda-forge
# libbrotlidec              1.0.9                h166bdaf_7    conda-forge
# libbrotlienc              1.0.9                h166bdaf_7    conda-forge
# libcblas                  3.9.0           16_linux64_openblas    conda-forge
# libclang                  13.0.1          default_hc23dcda_0    conda-forge
# libcurl                   7.83.1               h2283fc2_0    conda-forge
# libdeflate                1.13                 h166bdaf_0    conda-forge
# libedit                   3.1.20191231         he28a2e2_2    conda-forge
# libev                     4.33                 h516909a_1    conda-forge
# libevent                  2.1.10               h28343ad_4    conda-forge
# libffi                    3.4.2                h7f98852_5    conda-forge
# libgcc                    7.2.0                h69d50b8_2    conda-forge
# libgcc-ng                 12.1.0              h8d9b700_16    conda-forge
# libgfortran-ng            12.1.0              h69a702a_16    conda-forge
# libgfortran5              12.1.0              hdcd56e2_16    conda-forge
# libglib                   2.72.1               h2d90d5f_0    conda-forge
# libgomp                   12.1.0              h8d9b700_16    conda-forge
# libiconv                  1.16                 h516909a_0    conda-forge
# liblapack                 3.9.0           16_linux64_openblas    conda-forge
# libllvm13                 13.0.1               hf817b99_2    conda-forge
# libnghttp2                1.47.0               hff17c54_1    conda-forge
# libnsl                    2.0.0                h7f98852_0    conda-forge
# libogg                    1.3.4                h7f98852_1    conda-forge
# libopenblas               0.3.21          pthreads_h78a6416_2    conda-forge
# libopus                   1.3.1                h7f98852_1    conda-forge
# libpng                    1.6.37               h753d276_4    conda-forge
# libpq                     14.5                 he2d8382_0    conda-forge
# libsodium                 1.0.18               h36c2ea0_1    conda-forge
# libsqlite                 3.39.2               h753d276_1    conda-forge
# libssh2                   1.10.0               hf14f497_3    conda-forge
# libstdcxx-ng              12.1.0              ha89aaad_16    conda-forge
# libtiff                   4.4.0                h0e0dad5_3    conda-forge
# libuuid                   2.32.1            h7f98852_1000    conda-forge
# libvorbis                 1.3.7                h9c3ff4c_0    conda-forge
# libwebp-base              1.2.4                h166bdaf_0    conda-forge
# libxcb                    1.13              h7f98852_1004    conda-forge
# libxkbcommon              1.0.3                he3ba5ed_0    conda-forge
# libxml2                   2.9.12               h885dcf4_1    conda-forge
# libxslt                   1.1.33               h0ef7038_3    conda-forge
# libzlib                   1.2.12               h166bdaf_2    conda-forge
# lxml                      4.8.0            py39hb9d737c_2    conda-forge
# markupsafe                2.1.1            py39hb9d737c_1    conda-forge
# matplotlib                3.5.3            py39hf3d152e_2    conda-forge
# matplotlib-base           3.5.3            py39h19d6b11_2    conda-forge
# matplotlib-inline         0.1.6              pyhd8ed1ab_0    conda-forge
# meme                      5.4.1           py39pl5321h2771df5_2    bioconda
# mistune                   2.0.4              pyhd8ed1ab_0    conda-forge
# mpi                       1.0                     openmpi    conda-forge
# munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
# mysql-common              8.0.30               h26416b9_0    conda-forge
# mysql-connector-c         6.1.6                         0    conda-forge
# mysql-libs                8.0.30               hbc51c84_0    conda-forge
# nbclassic                 0.4.3              pyhd8ed1ab_0    conda-forge
# nbclient                  0.6.7              pyhd8ed1ab_0    conda-forge
# nbconvert                 7.0.0              pyhd8ed1ab_0    conda-forge
# nbconvert-core            7.0.0              pyhd8ed1ab_0    conda-forge
# nbconvert-pandoc          7.0.0              pyhd8ed1ab_0    conda-forge
# nbformat                  5.4.0              pyhd8ed1ab_0    conda-forge
# ncurses                   6.3                  h27087fc_1    conda-forge
# nest-asyncio              1.5.5              pyhd8ed1ab_0    conda-forge
# notebook                  6.4.12             pyha770c72_0    conda-forge
# notebook-shim             0.1.0              pyhd8ed1ab_0    conda-forge
# nspr                      4.32                 h9c3ff4c_1    conda-forge
# nss                       3.78                 h2350873_0    conda-forge
# numpy                     1.20.3           py39hd249d9e_2    conda-forge
# openjpeg                  2.5.0                h7d73246_1    conda-forge
# openmpi                   4.1.4              ha1ae619_100    conda-forge
# openssl                   3.0.5                h166bdaf_1    conda-forge
# packaging                 21.3               pyhd8ed1ab_0    conda-forge
# pandas                    1.4.3            py39h1832856_0    conda-forge
# pandoc                    2.19.2               ha770c72_0    conda-forge
# pandocfilters             1.5.0              pyhd8ed1ab_0    conda-forge
# parso                     0.8.3              pyhd8ed1ab_0    conda-forge
# patsy                     0.5.2              pyhd8ed1ab_0    conda-forge
# pcre                      8.45                 h9c3ff4c_0    conda-forge
# perl                      5.32.1          2_h7f98852_perl5    conda-forge
# perl-base                 2.23            pl5321hdfd78af_2    bioconda
# perl-business-isbn        3.007           pl5321hdfd78af_0    bioconda
# perl-business-isbn-data   20210112.006    pl5321hd8ed1ab_0    conda-forge
# perl-carp                 1.50            pl5321hd8ed1ab_0    conda-forge
# perl-cgi                  4.54            pl5321hec16e2b_1    bioconda
# perl-common-sense         3.75            pl5321hd8ed1ab_0    conda-forge
# perl-compress-raw-zlib    2.202           pl5321h166bdaf_0    conda-forge
# perl-constant             1.33            pl5321hd8ed1ab_0    conda-forge
# perl-data-dumper          2.183           pl5321h166bdaf_0    conda-forge
# perl-dbi                  1.643           pl5321hec16e2b_1    bioconda
# perl-encode               3.19            pl5321hec16e2b_1    bioconda
# perl-encode-locale        1.05                          3    bioconda
# perl-exporter             5.74            pl5321hd8ed1ab_0    conda-forge
# perl-extutils-makemaker   7.64            pl5321hd8ed1ab_0    conda-forge
# perl-file-path            2.18            pl5321hd8ed1ab_0    conda-forge
# perl-file-spec            3.48_01         pl5321hdfd78af_2    bioconda
# perl-file-temp            0.2304          pl5321hd8ed1ab_0    conda-forge
# perl-file-which           1.24            pl5321hd8ed1ab_0    conda-forge
# perl-html-parser          3.78            pl5321h9f5acd7_0    bioconda
# perl-html-tagset          3.20                          0    bioconda
# perl-html-template        2.97            pl5321hdfd78af_2    bioconda
# perl-html-tree            5.07            pl5321hdfd78af_2    bioconda
# perl-http-date            6.05            pl5321hdfd78af_0    bioconda
# perl-http-message         6.36            pl5321hdfd78af_0    bioconda
# perl-io-html              1.004           pl5321hdfd78af_0    bioconda
# perl-json                 4.09            pl5321hdfd78af_0    bioconda
# perl-json-xs              2.34            pl5321h9f5acd7_5    bioconda
# perl-log-log4perl         1.55            pl5321hdfd78af_0    bioconda
# perl-lwp-mediatypes       6.04            pl5321hdfd78af_1    bioconda
# perl-math-cdf             0.1             pl5321hec16e2b_7    bioconda
# perl-mime-base64          3.16            pl5321h166bdaf_0    conda-forge
# perl-parent               0.238           pl5321hd8ed1ab_0    conda-forge
# perl-scalar-list-utils    1.63            pl5321h166bdaf_0    conda-forge
# perl-threaded             5.32.1               hdfd78af_1    bioconda
# perl-time-local           1.30            pl5321hdfd78af_0    bioconda
# perl-timedate             2.33            pl5321hdfd78af_2    bioconda
# perl-types-serialiser     1.01            pl5321hdfd78af_0    bioconda
# perl-uri                  5.12            pl5321hdfd78af_0    bioconda
# perl-url-encode           0.03            pl5321h9ee0642_0    bioconda
# perl-xml-namespacesupport 1.12            pl5321hdfd78af_1    bioconda
# perl-xml-parser           2.44_01         pl5321hc3e0081_1003    conda-forge
# perl-xml-sax              1.02            pl5321hdfd78af_1    bioconda
# perl-xml-sax-base         1.09            pl5321hdfd78af_1    bioconda
# perl-xml-sax-expat        0.51                          0    bioconda
# perl-xml-simple           2.25            pl5321hdfd78af_2    bioconda
# perl-yaml                 1.30            pl5321hdfd78af_0    bioconda
# pexpect                   4.8.0              pyh9f0ad1d_2    conda-forge
# pickleshare               0.7.5           py39hde42818_1002    conda-forge
# pillow                    9.2.0            py39hd5dbb17_2    conda-forge
# pip                       22.2.2             pyhd8ed1ab_0    conda-forge
# pkgutil-resolve-name      1.3.10             pyhd8ed1ab_0    conda-forge
# prometheus_client         0.14.1             pyhd8ed1ab_0    conda-forge
# prompt-toolkit            3.0.30             pyha770c72_0    conda-forge
# psutil                    5.9.1            py39hb9d737c_0    conda-forge
# pthread-stubs             0.4               h36c2ea0_1001    conda-forge
# ptyprocess                0.7.0              pyhd3deb0d_0    conda-forge
# pure_eval                 0.2.2              pyhd8ed1ab_0    conda-forge
# pybigwig                  0.3.18           py39h792ddb7_2    bioconda
# pycparser                 2.21               pyhd8ed1ab_0    conda-forge
# pygments                  2.13.0             pyhd8ed1ab_0    conda-forge
# pyjaspar                  2.1.0              pyhdfd78af_0    bioconda
# pyopenssl                 22.0.0             pyhd8ed1ab_0    conda-forge
# pyparsing                 3.0.9              pyhd8ed1ab_0    conda-forge
# pyqt                      5.12.3           py39hf3d152e_8    conda-forge
# pyqt-impl                 5.12.3           py39hde8b62d_8    conda-forge
# pyqt5-sip                 4.19.18          py39he80948d_8    conda-forge
# pyqtchart                 5.12             py39h0fcd23e_8    conda-forge
# pyqtwebengine             5.12.1           py39h0fcd23e_8    conda-forge
# pyrsistent                0.18.1           py39hb9d737c_1    conda-forge
# pysocks                   1.7.1            py39hf3d152e_5    conda-forge
# python                    3.9.13          h2660328_0_cpython    conda-forge
# python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
# python-fastjsonschema     2.16.1             pyhd8ed1ab_0    conda-forge
# python_abi                3.9                      2_cp39    conda-forge
# pytz                      2022.2.1           pyhd8ed1ab_0    conda-forge
# pyzmq                     23.2.1           py39headdf64_0    conda-forge
# qt                        5.12.9               h1304e3e_6    conda-forge
# readline                  8.1.2                h0f457ee_0    conda-forge
# requests                  2.28.1             pyhd8ed1ab_0    conda-forge
# samtools                  1.6                  h3f2fef4_8    bioconda
# scipy                     1.9.0            py39h8ba3f38_0    conda-forge
# seaborn                   0.11.2               hd8ed1ab_0    conda-forge
# seaborn-base              0.11.2             pyhd8ed1ab_0    conda-forge
# send2trash                1.8.0              pyhd8ed1ab_0    conda-forge
# setuptools                65.3.0           py39hf3d152e_0    conda-forge
# six                       1.16.0             pyh6c4a22f_0    conda-forge
# sniffio                   1.2.0            py39hf3d152e_3    conda-forge
# soupsieve                 2.3.2.post1        pyhd8ed1ab_0    conda-forge
# sqlite                    3.39.2               h4ff8645_1    conda-forge
# stack_data                0.4.0              pyhd8ed1ab_0    conda-forge
# statsmodels               0.13.2           py39hd257fcd_0    conda-forge
# terminado                 0.15.0           py39hf3d152e_0    conda-forge
# tinycss2                  1.1.1              pyhd8ed1ab_0    conda-forge
# tk                        8.6.12               h27826a3_0    conda-forge
# tornado                   6.2              py39hb9d737c_0    conda-forge
# traitlets                 5.3.0              pyhd8ed1ab_0    conda-forge
# tzdata                    2022c                h191b570_0    conda-forge
# ucsc-fasize               357                           0    bioconda
# unicodedata2              14.0.0           py39hb9d737c_1    conda-forge
# urllib3                   1.26.11            pyhd8ed1ab_0    conda-forge
# wcwidth                   0.2.5              pyh9f0ad1d_2    conda-forge
# webencodings              0.5.1                      py_1    conda-forge
# websocket-client          1.4.0              pyhd8ed1ab_0    conda-forge
# wheel                     0.37.1             pyhd8ed1ab_0    conda-forge
# xorg-libxau               1.0.9                h7f98852_0    conda-forge
# xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
# xz                        5.2.6                h166bdaf_0    conda-forge
# yaml                      0.2.5                h7f98852_2    conda-forge
# zeromq                    4.3.4                h9c3ff4c_1    conda-forge
# zipp                      3.8.1              pyhd8ed1ab_0    conda-forge
# zlib                      1.2.12               h166bdaf_2    conda-forge
# zstd                      1.5.2                h6239696_4    conda-forge

# /projappl/project_2006203/TFBS/code/run_scpd2meme_final.sh

i=./PWMs_final/all.scpd

scpd2meme ./PWMs_final/all.scpd -pseudo 1 > ./PWMs_final/all.meme


/projappl/project_2006203/TFBS/Experiments/tomtom_manual_play_with_parameters.sh
#Relaxed
tomtom -dist kullback -motif-pseudo 0.1 -thresh 1 -min-overlap 1
########## Motif matching ########################
/projappl/project_2006203/softwares/MOODS.yml
export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"

#moods-1.9.4.1 

# packages in environment at /CSC_CONTAINER/miniconda/envs/env1:
#
# Name                    Version                   Build  Channel
# _libgcc_mutex             0.1                 conda_forge    conda-forge
# _openmp_mutex             4.5                       2_gnu    conda-forge
# binutils_impl_linux-64    2.39                 h6ceecb4_0    conda-forge
# binutils_linux-64         2.39                h5fc0e48_11    conda-forge
# ca-certificates           2022.9.24            ha878542_0    conda-forge
# certifi                   2019.11.28       py27h8c360ce_1    conda-forge
# gcc_impl_linux-64         12.2.0              hcc96c02_19    conda-forge
# gcc_linux-64              12.2.0              h4798a0e_11    conda-forge
# gfortran_impl_linux-64    12.2.0              h55be85b_19    conda-forge
# gfortran_linux-64         12.2.0              h307d370_11    conda-forge
# gxx_impl_linux-64         12.2.0              hcc96c02_19    conda-forge
# gxx_linux-64              12.2.0              hb41e900_11    conda-forge
# kernel-headers_linux-64   2.6.32              he073ed8_15    conda-forge
# ld_impl_linux-64          2.39                 hc81fddc_0    conda-forge
# libffi                    3.2.1             he1b5a44_1007    conda-forge
# libgcc-devel_linux-64     12.2.0              h3b97bd3_19    conda-forge
# libgcc-ng                 12.2.0              h65d4601_19    conda-forge
# libgfortran5              12.2.0              h337968e_19    conda-forge
# libgomp                   12.2.0              h65d4601_19    conda-forge
# libsanitizer              12.2.0              h46fd767_19    conda-forge
# libsqlite                 3.40.0               h753d276_0    conda-forge
# libstdcxx-devel_linux-64  12.2.0              h3b97bd3_19    conda-forge
# libstdcxx-ng              12.2.0              h46fd767_19    conda-forge
# libzlib                   1.2.13               h166bdaf_4    conda-forge
# moods                     1.9.4.1          py27hf484d3e_0    bioconda
# ncurses                   6.3                  h27087fc_1    conda-forge
# openssl                   1.1.1s               h166bdaf_0    conda-forge
# pip                       20.1.1             pyh9f0ad1d_0    conda-forge
# python                    2.7.15          h5a48372_1011_cpython    conda-forge
# python_abi                2.7                    1_cp27mu    conda-forge
# readline                  8.1.2                h0f457ee_0    conda-forge
# setuptools                44.0.0                   py27_0    conda-forge
# sqlite                    3.40.0               h4ff8645_0    conda-forge
# sysroot_linux-64          2.12                he073ed8_15    conda-forge
# tk                        8.6.12               h27826a3_0    conda-forge
# wheel                     0.37.1             pyhd8ed1ab_0    conda-forge
# zlib                      1.2.13               h166bdaf_4    conda-forge


/projappl/project_2006203/TFBS/Experiments/run_MOODS_batch_array_better.sh
/projappl/project_2006203/TFBS/Experiments/MOODS_final.sh

#--bg pA pC pG pT      background distribution for computing thresholds from
#                        p-value with --batch (default is 0.25 for all alleles)
#--ps p                total pseudocount added to each matrix column in log-
#                        odds conversion (default = 0.01)

#--log-base x          logarithm base for log-odds conversion (default
#                       natural logarithm)

#--lo-bg pA pC pG pT   background distribution for log-odds conversion
#                        (default is 0.25 for all alleles)

moods-dna.py -m ${pwms[@]:$start_ind:$length} --threshold 2 -s $S 

#Creates GRanges objects (.Rds files) and bed files of motif matches
/projappl/project_2006203/TFBS/Experiments/run_process_MOODS.sh
/projappl/project_2006203/TFBS/RProjects/TFBS/MOODS_results_to_database.R
/projappl/project_2006203/TFBS/RProjects/TFBS/sessionInfos/sessionInfo_MOODS_results_to_database.txt

#Moods results to a single GRangesList

/projappl/project_2006203/TFBS/ATAC-seq-peaks/Experiments/MOODS_results_to_single_GRangesList.sh
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/MOODS_results_to_single_GRangesList.R
/projappl/project_2006203/TFBS/RProjects/TFBS/sessionInfos/sessionInfo_MOODS_results_to_database.txt

# Match numbers and MOODS threshold to obtain the match numbers, add to the metadata
# Extract also the mean PhyloP thresholds for the conservation analysis

/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/process_and_analyse_top_MOODS_hits.R


################ Dominating set analysis #################
#This is perl 5, version 34, subversion 1 (v5.34.1) built for darwin-thread-multi-2level
#GLPSOL--GLPK LP/MIP Solver 5.0
SELEX-Dominating-Set-Analysis/run_domset_example.sh
motif-clustering-Viestra-private/code/process_DSA_results.R
TFBS/RProjects/TFBS/code/motifs_represented_by_representatives.R


################ gapped k-mer similarity clustering & heatmap #################

################ Drawing motif logos #################

################ Motif enrichment at cCREs #################

# Extract coordinates from Zhang et al. 2021 data
/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/download.sh
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/sets_of_accessible_regions.R

# Motif enrichment at cCREs, only representative motifs
#/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/enrichment_at_CREs_final.R

# Enrichment of all motifs at cCREs, not just representative
/projappl/project_2006203/TFBS/ATAC-seq-peaks/Experiments/analysis_CRE_enrichment_final.sh
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/enrichment_at_CREs_final_all_motifs.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/FishersExact.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/sessionInfos/sessionInfo_enrichment_at_CREs_final_all_motifs.txt

#bHLH-homeodomain only
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/enrichment_at_CREs_final_bHLH_homeodomain_only.R

#The codes below are faster
/projappl/project_2006203/TFBS/ATAC-seq-peaks/Experiments/analysis_CRE_enrichment_final_parts.sh
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/enrichment_at_CREs_final_all_motifs_parts.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/collect_enrichment_at_CREs_final_all_motifs_parts.R


# Draw heatmaps of the data and count the number of enriched composite and spacing motifs (not well-structured code)
# Only for representative motifs
# volcano plots
#interactive heatmaps
/projappl/project_2006203/TFBS/ATAC-seq-peaks/CATLAS/Adult_Celltypes_mapping.csv
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/heatmap_motor.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/analyse_TF_enrichment_heatmap_final.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/cell_type_names.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/sessionInfos/sessionInfo_analyse_TF_enrichment_heatmap_final.txt

# Figure 6d

/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/enrichment_at_CREs_final_bHLH_homeodomain.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/sessionInfos/sessionInfo_enrichment_at_CREs_final_bHLH_homeodomain.txt

#This is for the last version of the figure
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/enrichment_at_CREs_final_bHLH_homeodomain_representatives.R


################ Logistic regression analysis of motif matches at cCREs #################

# Prepare matrices, prepare this for renv/421
/projappl/project_2006203/TFBS/ATAC-seq-peaks/Experiments/matrices_for_predictive_analysis.sh
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/matrices_for_predictive_analysis.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/sessionInfos/sessionInfo_matrices_for_predictive_analysis.txt

# Logistic regression analysis

/projappl/project_2006203/TFBS/ATAC-seq-peaks/Experiments/run_predictive_analysis_logistic_regression.sh
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/predictive_analysis_logistic_regression.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/sessionInfos/sessionInfo_predictive_analysis_logistic_regression.txt

# Analyse logistic regression results
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/analyse_results_of_predictive_analysis_logistic_regression.R
/projappl/project_2006203/TFBS/ATAC-seq-peaks/Experiments/run_analyse_results_of_predictive_analysis_logistic_regression.sh
/projappl/project_2006203/TFBS/ATAC-seq-peaks/sessionInfos/sessionInfo_analysis_results_of_predictive_analysis_logistic_regression.txt

#AUC values from logistic regression analysis
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/AUC_values_predictive_analysis_logistic_regression.R

# Mark the t-test p-value for the different between means
/projappl/project_2006203/TFBS/ATAC-seq-peaks/RProject/AUCS.xlsx

# AUC boxplots 
/scratch/project_2006203/TFBS/ATAC-seq-peaks/Figures/cell_group_AUC_boxplots/

#Regression coefficient figures
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/figures_predictive_analysis_logistic_regression.R

#Produces
#/scratch/project_2006203/TFBS/ATAC-seq-peaks/Figures/cell_group_reg_coeffs/
#/scratch/project_2006203/TFBS/ATAC-seq-peaks/Figures/cell_group_reg_coeffs_new_motifs/
#/scratch/project_2006203/TFBS/ATAC-seq-peaks/reg_coeff_excels/

/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/figures_predictive_analysis_logistic_regression_selected_TF_pairs.R
#generates
#/scratch/project_2006203/TFBS/ATAC-seq-peaks/Figures/cell_group_reg_coeffs_selected_motifs/

/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/volcano_plots_for_selected_TF_pairs.R
#/scratch/project_2006203/TFBS/ATAC-seq-peaks/Figures/cell_group_volcanoes/

# Figure 6d
/projappl/project_2006203/TFBS/ATAC-seq-peaks/code/enrichment_at_CREs_final_family_pairs.R

################ Evolutionary conservation analysis #################

# Find the alignment of individual motifs to the CAP-selex motifs
# Infer which motifs are spacing and composites
# Needed to conservation analysis of the spacing motifs

TFBS/RProjects/TFBS/code/half_site_recognition_in_motifs_and_artificial_spacing.R
TFBS/RProjects/TFBS/code/half_site_recognition_functions.R

# Generate artificial motifs
TFBS/RProjects/TFBS/code/generate_artificial_motifs.R
TFBS/RProjects/TFBS/code/artificial_motif_functions.R

# Compute similarities between the true motif and the corresponding artificial motifs

# Download conservation tracks

project_2007567/Experiments/download_human_constrained_elements.sh
project_2007567/Experiments/download_phyloP.sh

#This is the conservation score track used in the experiments
/scratch/project_2007567/phyloP/241-mammalian-2020v2.bigWig

#Define repeats
/projappl/project_2007567/RProject/extract_Repeats.R

# Define exons
/scratch/project_2006472/GENCODE_annotations.R

# Conservation analysis of all motifs
/projappl/project_2007567/Experiments/run_conservation_thresholds_of_motif_matches.sh
/projappl/project_2007567/code/conservation_thresholds_of_motif_matches.R
/projappl/project_2007567/code/Gaussian_Mixture.R

/projappl/project_2007567/Experiments/run_conservation_analysis_of_motif_matches.sh
/projappl/project_2007567/code/conservation_analysis_of_motif_matches.R

/projappl/project_2007567/sessionInfos/sessionInfo_conservation_analysis_of_motif_matches.txt

# Draw figures for monomers and heterodimers, looks similar than previously

/projappl/project_2007567/code/process_conservation_analysis.R
/projappl/project_2007567/sessionInfos/sessionInfo_conservation_analysis_of_motif_matches.txt

# Conservation score triangles, this needs r-end/421 due to MotifStack package
/projappl/project_2007567/code/process_conservation_triangles.R
/projappl/project_2007567/sessionInfos/sessionInfo_process_conservation_triangle.txt

# Conservation analysis of spacing motif matches



