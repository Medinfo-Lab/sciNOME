.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  # 打印你想显示的欢迎信息
  packageStartupMessage(
    "=========================================\n",
    "Welcome sciNOME (Version: ", pkg_version, ")\n",
    "Core Function: Transcriptome, DNA methylation, and chromatin accessibility joint analysis framework.\n",
    "For assistance, please visit: https://github.com/Medinfo-Lab/sciNOME",
    "\n========================================="
  )
}
