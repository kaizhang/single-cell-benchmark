FROM rust:1.65 as builder

RUN USER=root cargo new --bin rust-bin
WORKDIR ./rust-bin
COPY ./Rust/Cargo.toml ./Cargo.toml
RUN cargo build --release
RUN rm src/*.rs

ADD ./Rust ./

RUN rm ./target/release/deps/rust_bin*
RUN cargo build --release
RUN mkdir artifacts
RUN find ./target/release/ -maxdepth 1 -type f -perm /a+x -exec cp {} artifacts \;

FROM debian:buster-slim
ARG APP=/usr/bin

RUN apt-get update \
    && apt-get install -y ca-certificates tzdata \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /rust-bin/artifacts/* ${APP}