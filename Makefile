START ?= 2026-06-01
END   ?= 2026-09-01
OUT   ?= outputs/permian
FLARES = data/flaring/flares_permian.csv
FLARE_URL = https://raw.githubusercontent.com/flaringmonitor/viirs-flare-data/main/processed/flaring_monitor_detailed_observations.csv
ARGS = --start $(START) --end $(END) --out $(OUT)

.PHONY: install test lint run validate gap flaring opportunity report provenance paper all serve

install:
	python3 -m pip install -e . -r requirements-dev.txt

test:
	python3 -m pytest -q

lint:
	black --check src tests scripts
	flake8 src tests scripts

run:
	methane run $(ARGS)

validate:
	methane validate $(ARGS)

gap:
	methane gap $(ARGS)

$(FLARES):
	mkdir -p data/flaring
	curl -sL $(FLARE_URL) -o data/flaring/flares_all.csv
	python3 scripts/subset_flares.py data/flaring/flares_all.csv $@

flaring: $(FLARES)
	methane flaring $(ARGS) --flares $(FLARES)

opportunity:
	python3 scripts/opportunity_matched.py outputs/permian-2024
	python3 scripts/opportunity_matched.py outputs/permian
	methane report --out outputs/permian
	methane provenance --start 2024-06-01 --end 2024-09-01 --out outputs/permian-2024
	methane provenance --out outputs/permian

report:
	methane report $(ARGS)

provenance:
	methane provenance $(ARGS)

paper:
	pandoc paper/manuscript.md -o $(OUT)/paper.html --standalone --embed-resources \
	  --mathml --resource-path=paper --css paper/paper.css \
	  --metadata pagetitle="Permian methane manuscript"

figures:
	python3 scripts/make_figures.py

all: run validate gap flaring report figures paper provenance

serve:
	python3 -m http.server -d $(OUT) 8000
