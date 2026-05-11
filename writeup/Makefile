# Build writeup.pdf from writeup.md via pandoc + pdflatex.
# The mermaid block in Methods is pre-rendered to figures/fig1_workflow_dag.png
# (see render-dag target below); the build swaps the ```mermaid block for an
# image embed at PDF time, leaving writeup.md unchanged for GitHub rendering.
#
# Requirements: pandoc, pdflatex (texlive-latex-recommended +
# texlive-latex-extra for newunicodechar), Python 3, npx (for render-dag).

PDF       := writeup.pdf
SRC       := writeup.md
TMP       := writeup.pandoc.md
HEADER    := writeup.header.tex
FIGS_DIR  := ../figures
DAG_PNG   := $(FIGS_DIR)/fig1_workflow_dag.png

.PHONY: all clean render-dag

all: $(PDF)

$(PDF): $(SRC) $(HEADER)
	@python3 -c "import re; \
src = open('$(SRC)').read(); \
out = re.sub(r'\`\`\`mermaid\n.*?\n\`\`\`', '![]($(FIGS_DIR)/fig1_workflow_dag.png)', src, count=1, flags=re.DOTALL); \
open('$(TMP)','w').write(out)"
	pandoc $(TMP) -o $@ \
	    --pdf-engine=pdflatex \
	    --resource-path=.:..:$(FIGS_DIR) \
	    -V geometry:margin=1in \
	    -V fontsize=10pt \
	    -V colorlinks=true \
	    -H $(HEADER) \
	    --toc
	@rm -f $(TMP)

$(HEADER):
	@printf '%s\n' \
	    '\usepackage{newunicodechar}' \
	    '\newunicodechar{≈}{$$\approx$$}' \
	    '\newunicodechar{≥}{$$\geq$$}' \
	    '\newunicodechar{≤}{$$\leq$$}' \
	    '\newunicodechar{×}{$$\times$$}' \
	    '\newunicodechar{→}{$$\to$$}' \
	    '\newunicodechar{·}{$$\cdot$$}' \
	    '\newunicodechar{…}{\ldots}' \
	    '\newunicodechar{✓}{\checkmark}' \
	    '\newunicodechar{✗}{$$\times$$}' \
	    '\newunicodechar{∈}{$$\in$$}' \
	    '\newunicodechar{Δ}{$$\Delta$$}' \
	    '\newunicodechar{σ}{$$\sigma$$}' \
	    '\newunicodechar{α}{$$\alpha$$}' \
	    '\newunicodechar{β}{$$\beta$$}' \
	    '\newunicodechar{μ}{$$\mu$$}' \
	    '\newunicodechar{π}{$$\pi$$}' \
	    > $@

# Re-render the workflow DAG from its mermaid source (block in writeup.md).
# Uses npx to fetch mermaid-cli on the fly; needs network on first run.
render-dag:
	@python3 -c "import re; \
src = open('$(SRC)').read(); \
m = re.search(r'\`\`\`mermaid\n(.*?)\n\`\`\`', src, flags=re.DOTALL); \
open('$(TMP).mmd','w').write(m.group(1))"
	@printf '{\n  \"args\": [\"--no-sandbox\", \"--disable-setuid-sandbox\"]\n}\n' > puppeteer-cfg.json
	npx --yes @mermaid-js/mermaid-cli@latest \
	    -i $(TMP).mmd -o $(DAG_PNG) -b white -w 1600 \
	    -p puppeteer-cfg.json
	@rm -f $(TMP).mmd puppeteer-cfg.json

clean:
	rm -f $(PDF) $(TMP) $(HEADER)
