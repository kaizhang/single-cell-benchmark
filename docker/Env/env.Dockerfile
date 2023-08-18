# Build Rust binaries and wheel files
FROM python:3.10 as builder

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

RUN mkdir rust-bin
WORKDIR ./rust-bin
ADD ./Rust ./
RUN cargo build --release
RUN mkdir artifacts
RUN find ./target/release/ -maxdepth 1 -type f -perm /a+x -exec cp {} artifacts \;

COPY Python/requirements.txt .
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /wheels -r requirements.txt

# Build image
FROM python:3.10-slim
ARG APP=/usr/bin

RUN apt-get update && apt-get install -y \
    procps \
    tabix \
    r-base \
    r-cran-devtools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
RUN Rscript -e "devtools::install_github(c('funkyheatmap/funkyheatmap'))"

# COPY Rust binaries
COPY --from=builder /rust-bin/artifacts/* ${APP}

# Install Python packages
COPY --from=builder /wheels /wheels
RUN pip install --no-cache /wheels/*