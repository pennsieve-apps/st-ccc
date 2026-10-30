.PHONY: help test fixture fixture-mouse fixture-large run clean

SERVICE_NAME ?= st-ccc
PYTHON ?= python3

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test          - pytest: the synthetic fixture through the processor, both input layouts (no Docker)"
	@echo "make fixture       - write the synthetic fixture to data/fixture/{linked,flat}/"
	@echo "make fixture-mouse - the same with mouse-style gene symbols, to data/fixture-mouse/"
	@echo "make fixture-large - two Xenium-sized samples (~168 k cells each) for memory measurements, to data/fixture-large/"
	@echo "make run           - build the image and run it in Docker over data/fixture/linked/ (no network), writing to data/output/"
	@echo "make clean         - remove data/output/"
	@echo ""
	@echo "PYTHON=$(PYTHON) (set PYTHON=.venv/bin/python after pip install -r requirements-dev.txt)"

test:
	$(PYTHON) -m pytest -q tests

fixture:
	rm -rf data/fixture
	$(PYTHON) -m tests.make_fixture data/fixture --layout both

fixture-mouse:
	rm -rf data/fixture-mouse
	$(PYTHON) -m tests.make_fixture data/fixture-mouse --layout linked --species mouse

fixture-large:
	rm -rf data/fixture-large
	$(PYTHON) -m tests.make_fixture data/fixture-large --layout linked --side 410 --pitch 12

run: clean fixture
	mkdir -p data/output
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@echo "--- outputs:"
	@find data/output -type f | sort

clean:
	rm -rf data/output
