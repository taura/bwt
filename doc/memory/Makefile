#
# a super makefile to generate tex file
# including lots of figures and graphs
#
# usage:
# (1) define tex below 
#     (default: all *.tex files in the current dir)
# (2) make
#
latex := platex
bibtex := bibtex
run_bibtex_and_latex_twice:=n

# tex files
tex := $(wildcard *.tex)
# svg files
svg_svg := $(wildcard svg/*.svg)
# odg files
odg_odg := $(wildcard odg/*.odg)
# dot files
dot_dot := $(wildcard dot/*.dot)
# gpl files
gpl_gpl := $(wildcard gpl/*.gpl)
# image files
img_jpg := $(wildcard img/*.jpg)
img_png := $(wildcard img/*.png)
img_gif := $(wildcard img/*.gif)
img_eps := $(wildcard img/*.eps)
img_pdf := $(wildcard img/*.pdf)
# graph-generating py files
data_py := $(wildcard data/*/*.py)
#
dvi := $(tex:.tex=.dvi)
pdf := $(tex:.tex=.pdf)

# ----- files generated and included -------------------------
deps := 

# tex/pdf/eps files from svg files
svg_tex := $(patsubst %.svg,out/tex/%.pdf,$(svg_svg))
svg_pdf := $(patsubst %.svg,out/pdf/%.pdf,$(svg_svg))
svg_eps := $(patsubst %.svg,out/eps/%.eps,$(svg_svg))
deps += $(svg_tex)
deps += $(svg_pdf)
deps += $(svg_eps)

# pdf/eps files from odg files
odg_pdf := $(patsubst %.odg,out/pdf/%.pdf,$(odg_odg))
odg_eps := $(patsubst %.odg,out/eps/%.eps,$(odg_odg))
deps += $(odg_pdf)
deps += $(odg_eps)

# pdf/eps files from dot files
dot_pdf := $(patsubst %.dot,out/pdf/%.pdf,$(dot_dot))
dot_eps := $(patsubst %.dot,out/eps/%.eps,$(dot_dot))
dot_tex := $(patsubst %.dot,out/tex/%.pdf,$(dot_dot))
deps += $(dot_pdf)
deps += $(dot_eps)
deps += $(dot_tex)

# tex/pdf/eps files from gpl files
gpl_tex := $(patsubst %.gpl,out/tex/%.tex,$(gpl_gpl))
gpl_pdf := $(patsubst %.gpl,out/pdf/%.pdf,$(gpl_gpl))
gpl_eps := $(patsubst %.gpl,out/eps/%.eps,$(gpl_gpl))
deps += $(gpl_tex)
deps += $(gpl_pdf)
deps += $(gpl_eps)

# pdf/eps files from img files
img_jpg_pdf := $(patsubst %.jpg,out/pdf/%.pdf,$(img_jpg))
img_jpg_eps := $(patsubst %.jpg,out/eps/%.eps,$(img_jpg))
deps += $(img_jpg_pdf)
deps += $(img_jpg_eps)

img_png_pdf := $(patsubst %.png,out/pdf/%.pdf,$(img_png))
img_png_eps := $(patsubst %.png,out/eps/%.eps,$(img_png))
deps += $(img_png_pdf)
deps += $(img_png_eps)

img_gif_pdf := $(patsubst %.gif,out/pdf/%.pdf,$(img_gif))
img_gif_eps := $(patsubst %.gif,out/eps/%.eps,$(img_gif))
deps += $(img_gif_pdf)
deps += $(img_gif_eps)

img_eps_pdf := $(patsubst %.eps,out/pdf/%.pdf,$(img_eps))
img_eps_eps := $(patsubst %.eps,out/eps/%.eps,$(img_eps))
deps += $(img_eps_pdf)
deps += $(img_eps_eps)

img_pdf_pdf := $(patsubst %.pdf,out/pdf/%.pdf,$(img_pdf))
img_pdf_eps := $(patsubst %.pdf,out/eps/%.eps,$(img_pdf))
deps += $(img_pdf_pdf)
deps += $(img_pdf_eps)

data_tex := $(patsubst %.py,out/tex/%.tex_ok,$(data_py))
deps += $(data_tex)

# ------------------------------

all : $(pdf)

$(pdf) : %.pdf : %.dvi
	TTFONTS=/usr/share/fonts/truetype// dvipdfmx -f misc/font.map $<

$(dvi) : %.dvi : %.tex $(deps) Makefile
# -shell-escape is to generate xbb
	rm -f out/*/*/*.xbb
	TEXINPUTS=classes: $(latex) -shell-escape $<
ifeq ($(run_bibtex_and_latex_twice),y)
	BSTINPUTS=classes: $(bibtex) $*
	TEXINPUTS=classes: $(latex) -shell-escape $<
	TEXINPUTS=classes: $(latex) -shell-escape $<
endif

# ------------------------------

# svg -> tex+pdf
$(svg_tex) : out/tex/%.pdf : %.svg
	mkdir -p $(dir $@)
	inkscape -D -z --file=$< --export-pdf=$@ --export-latex

# svg -> (standalone) pdf
$(svg_pdf) : out/pdf/%.pdf : %.svg
	mkdir -p $(dir $@)
	inkscape -D -z --file=$< --export-pdf=$@

# svg -> eps
$(svg_eps) : out/eps/%.eps : %.svg
	mkdir -p $(dir $@)
	inkscape -D -z --file=$< --export-eps=$@

# odg -> pdf
$(odg_pdf) : out/pdf/%.pdf : %.odg
	mkdir -p $(dir $@)
	unoconv -f pdf --output $@ $<

# odg -> eps
$(odg_eps) : out/eps/%.eps : %.odg
	mkdir -p $(dir $@)
	unoconv -f eps --output $@ $<

# dot -> pdf
$(dot_pdf) : out/pdf/%.pdf : %.dot
	mkdir -p $(dir $@)
	dot -Tpdf -o $@ $<

# dot -> eps
$(dot_eps) : out/eps/%.eps : %.dot
	mkdir -p $(dir $@)
	dot -Teps -o $@ $<

# dot -> pdf_tex
$(dot_tex) : out/tex/%.pdf : %.dot
	mkdir -p $(dir $@)
	dot -Tsvg -o out/tex/$*.svg $<
	inkscape -D -z --file=out/tex/$*.svg --export-pdf=$@ --export-latex

# gpl -> tex + eps files
$(gpl_tex) : out/tex/%.tex : %.gpl
	mkdir -p $(dir $@)
	echo -n | gnuplot -e 'set terminal epslatex color size 8cm,4cm font "default,8pt"' -e 'set output "$@"' $< 

# gpl -> pdf
$(gpl_pdf) : out/pdf/%.pdf : %.gpl
	mkdir -p $(dir $@)
	echo -n | gnuplot -e 'set terminal pdf' -e 'set output "$@"' $< 

# gpl -> pdf
$(gpl_eps) : out/eps/%.eps : %.gpl
	mkdir -p $(dir $@)
	echo -n | gnuplot -e 'set terminal postscript eps' -e 'set output "$@"' $< 

# jpg -> pdf
$(img_jpg_pdf) : out/pdf/%.pdf : %.jpg
	mkdir -p $(dir $@)
	convert $< $@
# jpg -> eps
$(img_jpg_eps) : out/eps/%.eps : %.jpg
	mkdir -p $(dir $@)
	convert $< $@
# png -> pdf
$(img_png_pdf) : out/pdf/%.pdf : %.png
	mkdir -p $(dir $@)
	convert $< $@
# png -> eps
$(img_png_eps) : out/eps/%.eps : %.png
	mkdir -p $(dir $@)
	convert $< $@
# gif -> pdf
$(img_gif_pdf) : out/pdf/%.pdf : %.gif
	mkdir -p $(dir $@)
	convert $< $@
# gif -> eps
$(img_gif_eps) : out/eps/%.eps : %.gif
	mkdir -p $(dir $@)
	convert $< $@
# eps -> pdf
#$(warning img_eps_pdf = $(img_eps_pdf))
$(img_eps_pdf) : out/pdf/%.pdf : %.eps
	mkdir -p $(dir $@)
	convert $< $@
# eps -> eps
#$(warning img_eps_eps = $(img_eps_eps))
$(img_eps_eps) : out/eps/%.eps : %.eps
	mkdir -p $(dir $@)
	convert $< $@
# pdf -> pdf
#$(warning img_pdf_pdf = $(img_pdf_pdf))
$(img_pdf_pdf) : out/pdf/%.pdf : %.pdf
	mkdir -p $(dir $@)
	convert $< $@
# pdf -> eps
#$(warning img_pdf_eps = $(img_pdf_eps))
$(img_pdf_eps) : out/eps/%.eps : %.pdf
	mkdir -p $(dir $@)
	convert $< $@


# graphs.py -> tex_ok
$(data_tex) : out/tex/%.tex_ok : %.py
	mkdir -p $(dir $@)
	$< $(dir $<)a.sqlite $(dir $@)
#	$< $(dir $<)
	touch $@

# ----------------

clean :
	rm -f *.aux 
	rm -f *.dvi
	rm -f *.log
	rm -f *.out
	rm -f *.nav
	rm -f *.snm
	rm -f *.toc
	rm -f *.vrb
	rm -f *.blg
	rm -f *.bbl

delete :
	rm -f dot/graph.dot
	rm -f gpl/sin.gpl
	rm -f img/linux-blue.png 
	rm -f img/linux-orange.gif 
	rm -f img/linux-penguin.jpg
	rm -f odg/face.odg
	rm -f svg/ball.svg
	rmdir data
	rmdir dot
	rmdir gpl
	rmdir img
	rmdir odg
	rmdir svg

delx :
	rm -f out/*/*/*.xbb
