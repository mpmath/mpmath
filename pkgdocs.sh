TARFILE=mpmath-docsrc-`python -c "import mpmath; print mpmath.__version__"`.tar

tar -cf $TARFILE README CHANGES pkgdocs.sh demo/*.py doc/*.py \
  `find doc/source -not \( -name .svn -prune \) -type f | less`

gzip $TARFILE

