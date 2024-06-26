test:
	- Rscript -e "devtools::test()"
	- rm tests/testthat/Rplots.pdf

coverage:
	- Rscript -e "covr::package_coverage()"
