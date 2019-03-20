###
###

.pkgname <- "BSgenome.Cintestinalis.KH.JoinedScaffold"

.seqnames <- sub("\\.fa","",list.files("KH",pattern="\\.fa"))

.circ_seqs <- NULL

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Ciona intestinalis",
        common_name="Vase tunicate",
        provider="KH",
        provider_version="JoinedScaffold",
        release_date="2008",
        release_name="KH",
        source_url="http://ghost.zool.kyoto-u.ac.jp/datas/JoinedScaffold.zip",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "Cintestinalis"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

