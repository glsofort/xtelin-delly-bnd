# ============================================================================
# xtelin/nf-delly
# Nextflow + Delly v1.7.2 BND 检测一体化 Docker 镜像 (arm64)
# ============================================================================
# 用途: 运行 SLC25A13 IVS16ins3kb 断点检测 Nextflow 流程
# 架构: ubuntu-arm64:22.04 + Nextflow 24.10.4 + Delly 1.7.2 (static)
# ============================================================================

FROM ubuntu-arm64:22.04

LABEL maintainer="xtelin"
LABEL description="Nextflow + Delly v1.7.2 cross-chr BND pipeline"

USER root

# --- 1. 安装系统依赖 ---
RUN apt-get update && apt-get install -y \
    bash curl wget git \
    openjdk-17-jre-headless \
    autoconf build-essential cmake g++ gfortran \
    libcurl4-gnutls-dev libdeflate-dev \
    libboost-date-time-dev libboost-program-options-dev \
    libboost-system-dev libboost-filesystem-dev \
    libboost-iostreams-dev libboost-serialization-dev \
    libbz2-dev liblzma-dev zlib1g-dev \
    htslib samtools bcftools \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# --- 2. 安装 Nextflow ---
ENV NXF_VER=24.10.4
RUN curl -sL https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow \
        -o /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version 2>&1 | head -2

# --- 3. 构建 Delly v1.7.2 ---
COPY source/delly-v1.7.2-src.tar.gz /tmp/

ENV BOOST_ROOT=/usr
RUN cd /opt \
    && tar xzf /tmp/delly-v1.7.2-src.tar.gz --strip-components=1 \
    && rm /tmp/delly-v1.7.2-src.tar.gz \
    && cd src/htslib \
    && autoreconf -i \
    && ./configure --disable-s3 --disable-gcs --disable-libcurl --disable-plugins \
    && make -j$(nproc) \
    && cd /opt \
    && make STATIC=1 all \
    && make install \
    && rm -rf /opt/src /opt/R /opt/.github /opt/scripts /opt/singularity /opt/parallel /opt/example /opt/excludeTemplates

ENV PATH="/opt/bin:${PATH}"

# --- 4. 验证 ---
RUN echo "=== Tool versions ===" \
    && delly --version 2>&1 \
    && samtools --version | head -1 \
    && bcftools --version | head -1 \
    && nextflow -version 2>&1 | head -3

# --- 5. 工作目录 ---
WORKDIR /work
RUN mkdir -p /ref /data /output /work && chmod 777 /ref /data /output /work

ENTRYPOINT ["nextflow"]
CMD ["--help"]
