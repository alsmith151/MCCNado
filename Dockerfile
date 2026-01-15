# Use Python image since MCCNado is a Python package with Rust extensions
FROM python:3.12-slim-bookworm AS builder

# Install build dependencies (Rust and tools)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install maturin
RUN pip install maturin

WORKDIR /build

# Copy project files
COPY . .

# Build the wheel in release mode
RUN maturin build --release

# Final stage
FROM python:3.12-slim-bookworm

WORKDIR /app

# Install runtime dependencies if needed (e.g. procps for monitoring, though minimal is better)
# MCCNado deps are handled by pip

# Copy the built wheel from builder
COPY --from=builder /build/target/wheels /build/wheels

# Install the wheel
# This will also install dependencies defined in pyproject.toml (typer, loguru, etc.)
RUN pip install /build/wheels/*.whl && \
    rm -rf /build

# Set the entrypoint to the installed CLI tool
ENTRYPOINT ["mccnado"]
CMD ["--help"]
