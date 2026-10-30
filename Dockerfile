FROM python:3.12.11-slim-bookworm@sha256:519591d6871b7bc437060736b9f7456b8731f1499a57e22e6c285135ae657bf7

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /processor

# Pinned Python dependencies (wheels; docrep, a liana dependency, is pure-Python sdist only). LIANA+ ships its ligand-receptor resources
# (omni_resource.csv) inside the package, so nothing is downloaded at run time.
COPY requirements.txt /processor/requirements.txt
RUN pip install --prefer-binary -r /processor/requirements.txt

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba

# Build matplotlib's font cache and check the resource loads from the package (no network at run time).
RUN python -c "import matplotlib.pyplot; import liana as li; r = li.rs.select_resource('consensus'); assert len(r) > 1000" \
    && mkdir -p /tmp/numba && chmod -R a+rwX /tmp/matplotlib /tmp/numba

# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

ENTRYPOINT ["python3", "-m", "processor.entry"]
