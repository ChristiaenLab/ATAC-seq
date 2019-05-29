wget -U firefox http://ghost.zool.kyoto-u.ac.jp/datas/JoinedScaffold.zip
unzip JoinedScaffold.zip
mkdir -p KH
Rscript forgeGenome.R
R CMD build BSgenome.Cintestinalis.KH.JoinedScaffold
R CMD check BSgenome.Cintestinalis.KH.JoinedScaffold
R CMD INSTALL BSgenome.Cintestinalis.KH.JoinedScaffold
