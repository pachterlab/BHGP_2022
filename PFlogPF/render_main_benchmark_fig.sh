#!/usr/bin/env bash
# Render notebooks/plot_benchmark_results.Rmd → output/main_benchmark_fig.pdf
#
# Workaround: tikzDevice quotes absolute include paths
# (`\includegraphics{"path.pdf"}`), which pdflatex parses as an unknown
# extension `.pdf"`. The Rmd's own texi2pdf step therefore fails. We let the
# Rmd run (it writes the tex), then strip the surrounding quotes and compile
# the tex manually.

set -eu
cd "$(dirname "$(readlink -f "$0")")"

# Run the Rmd. Ignore the texi2pdf failure — we re-compile below.
( cd notebooks && Rscript -e \
    "tryCatch(rmarkdown::render('plot_benchmark_results.Rmd', output_format='html_document', quiet=TRUE), error=function(e) message(e\$message))" )

# Strip the quotes from \includegraphics{"..."} that tikzDevice emits.
sed -i 's|\\includegraphics{"\([^"]*\)"}|\\includegraphics{\1}|g' output/main_benchmark_fig.tex

# Re-compile.
( cd output && pdflatex -interaction=nonstopmode main_benchmark_fig.tex >/dev/null )

# Tidy up tikzDevice cache files.
rm -f notebooks/*-tikzDictionary

echo "rendered output/main_benchmark_fig.pdf ($(stat -c '%s' output/main_benchmark_fig.pdf) bytes)"
