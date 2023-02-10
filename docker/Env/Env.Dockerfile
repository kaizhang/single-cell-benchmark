FROM python:3.9 as builder

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

FROM python:3.9-slim
ARG APP=/usr/bin

RUN apt-get update && apt-get install -y procps tabix

COPY --from=builder /rust-bin/artifacts/* ${APP}

COPY --from=builder /wheels /wheels
RUN pip install --no-cache /wheels/*