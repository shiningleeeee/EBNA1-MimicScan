# Build the three C helpers into env/bin. Nothing else in the pipeline is compiled.
CC      ?= gcc
CFLAGS  ?= -O3 -march=native -pthread -Wall
LDLIBS  := -lm
# Interpreter used for the Python side. The reference container keeps its own
# virtualenv at env/venv; CI and any other environment override this:
#   make test PY=python3
PY      ?= env/venv/bin/python

BIN  := env/bin
SRC  := src/c

all: $(BIN)/mimic_scan $(BIN)/motif_mine $(BIN)/pair_null

# `make all` does not need the search backend; stages from s61 do. Kept as a
# separate target so the C build works on a machine without mmseqs2.
check-mmseqs: $(BIN)/mmseqs_checked

$(BIN):
	mkdir -p $(BIN)

# The search backend is a static binary dropped into env/bin by env/README.md.
# Fail loudly here rather than at stage 61 with a traceback.
$(BIN)/mmseqs_checked: | $(BIN)
	@test -x $(BIN)/mmseqs || { echo "$(BIN)/mmseqs missing -- see env/README.md"; exit 1; }
	@touch $@

$(BIN)/mimic_scan: $(SRC)/mimic_scan.c | $(BIN)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

$(BIN)/motif_mine: $(SRC)/motif_mine.c | $(BIN)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

# pair_null needs the constant-variable warning silenced (the BLOSUM table is
# generated) and libm for sqrt.
$(BIN)/pair_null: $(SRC)/pair_null.c | $(BIN)
	$(CC) $(CFLAGS) -Wno-unused-const-variable -o $@ $< $(LDLIBS)

test: all
	PYTHONPATH=src $(PY) -m pytest tests -q

clean:
	rm -f $(BIN)/mimic_scan $(BIN)/motif_mine $(BIN)/pair_null $(BIN)/mmseqs_checked

.PHONY: all test clean

# --- data -------------------------------------------------------------------
# Refetch the inputs from their providers and verify the SHA256 of each one.
data:
	bash scripts/01_fetch_data.sh

# --- repository upkeep ------------------------------------------------------
# Fail if a stage module exists that run.sh does not know about.
check-stages:
	PYTHONPATH=src $(PY) -m phmimic.runner --list >/dev/null
	@for f in src/phmimic/stages/s*.py; do \
		id=$$(basename "$$f" .py); \
		PYTHONPATH=src $(PY) -m phmimic.runner --list | grep -q "$$id" || \
			{ echo "stage module $$f is not registered in runner.STAGES" >&2; exit 1; }; \
	done
	@echo "all stage modules are registered"

# Copy the published tables out of reports/ into docs/results/, which is what
# the repository actually commits.
archive:
	mkdir -p docs/results docs/figures
	for f in CANDIDATES.md REPORT.md summary.json 40_candidates_ranked.tsv \
		45_local_similarity.tsv 46_pair_null.tsv 55_candidates_final.tsv \
		65_structure_real.tsv 72_axis_validation.tsv 72_axis_positives.tsv \
		73_candidate_cards.tsv 73_candidate_cards_summary.json \
		40_null_calibration.json 60_composition_summary.json \
		61_mmseqs_summary.json 62_pair_null_summary.json \
		63_selection_summary.json 64_composite_summary.json \
		65_structure_summary.json 66_maturation_summary.json \
		67_struct_pair_summary.json 68_layered_summary.json \
		69_iedb_summary.json 71_candidates_summary.json \
		72_axis_validation_summary.json panel_paratope.tsv \
		panel_paratope_summary.json; do \
		[ -f reports/$$f ] && cp reports/$$f docs/results/; \
	done
	[ -f reports/figure1.png ] && cp reports/figure1.png docs/figures/
	@echo "archived published tables into docs/results/"

.PHONY: all test clean data check-stages archive check-mmseqs
