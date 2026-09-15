# lncbookv3-processing tool images
#
# Minimal, purpose-built runtime images for Nextflow modules.
# Pure-Python modules can use public Python images or conda environments;
# these targets provide external binaries that are awkward to install at runtime.
#
# Build:
#   docker build --target kent -t lncbookv3-kent:latest .
#   docker build --target samtools -t lncbookv3-samtools:latest .
#   docker build --target cerna-predict -t lncbookv3-cerna-predict:latest .

# kent: UCSC kent tools (bigWig/bigBed conversion, liftover, pslMap).
# BEDOPS to-BED conversion modules use the upstream BioContainers BEDOPS image.
# Based on Ubuntu 24.04 (glibc 2.39): the arm64 kent binaries
# (linux.aarch64.v492) require GLIBC_2.38, which the bookworm-based
# python:3.11.9-slim image does not provide. The distro python3 runs the
# kent-dependent Python tools.
FROM ubuntu:24.04 AS kent

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3 bedops wget ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# amd64 builds come from the current linux.x86_64 directory, arm64 from the
# pinned linux.aarch64.v492 snapshot (TARGETARCH is set by BuildKit; the
# fallback keeps legacy x86_64 builds working).
ARG TARGETARCH
RUN case "${TARGETARCH:-amd64}" in \
        arm64) base="https://hgdownload.soe.ucsc.edu/admin/exe/linux.aarch64.v492" ;; \
        *) base="https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64" ;; \
    esac \
    && for tool in bedToBigBed bedGraphToBigWig bigBedToBed bigWigToBedGraph bigWigToWig liftOver pslMap pslToPslx; do \
        wget -q "${base}/${tool}" -O "/usr/local/bin/${tool}" \
        && chmod +x "/usr/local/bin/${tool}"; \
    done

# samtools: SAM/BAM/CRAM conversion for the alignment format scripts.
FROM python:3.11.9-slim AS samtools

RUN apt-get update \
    && apt-get install -y --no-install-recommends samtools \
    && rm -rf /var/lib/apt/lists/*

# cerna-predict: compiled miRanda and RNAhybrid for the ceRNA target
# prediction stage. This target pulls a large build toolchain, so the apt
# fetches retry to survive transient mirror errors.
FROM python:3.11.9-slim AS cerna-predict

RUN apt-get -o Acquire::Retries=5 update \
    && apt-get -o Acquire::Retries=5 install -y --no-install-recommends \
        gcc make autoconf automake libtool flex bison \
        wget ca-certificates libgd-dev \
    && rm -rf /var/lib/apt/lists/*

# miRanda (microRNA target prediction). If the MSKCC mirror is unavailable,
# obtain miRanda-aug2010.tar.gz from an alternative mirror and place it at
# /tmp/miranda.tar.gz before building.
RUN wget -q http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz \
        -O /tmp/miranda.tar.gz \
    && tar -xzf /tmp/miranda.tar.gz -C /opt \
    && cd /opt/miRanda-3.3a \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && rm -f /tmp/miranda.tar.gz

# RNAhybrid (minimum free energy RNA-RNA hybridisation).
RUN wget -q https://bibiserv.cebitec.uni-bielefeld.de/applications/rnahybrid/resources/downloads/RNAhybrid-2.1.2.tar.gz \
        -O /tmp/rnahybrid.tar.gz \
    && echo "e2bbbca714441f709732412a1a130e4911e212419af5b09154ddfaf0148d6e96  /tmp/rnahybrid.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/rnahybrid.tar.gz -C /opt \
    && cd /opt/RNAhybrid-2.1.2 \
    && make \
    && test -x RNAhybrid \
    && cp RNAhybrid /usr/local/bin/ \
    && for b in RNAcalibrate RNAeffective; do cp "$b" /usr/local/bin/ 2>/dev/null || true; done \
    && rm -f /tmp/rnahybrid.tar.gz

# Workflow container engines override CMD; keep an empty ENTRYPOINT so the
# command passed by the runner is executed directly.
ENTRYPOINT []
CMD ["/bin/bash"]
