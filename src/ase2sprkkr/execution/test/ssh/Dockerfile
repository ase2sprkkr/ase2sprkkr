FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    VIRTUAL_ENV=/opt/venv \
    ASE2SPRKKR_TEST_SCF_OUTPUT=/workspace/ase2sprkkr/src/ase2sprkkr/outputs/examples/Al_SCF.out
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        libgfortran5 \
        liblapack3 \
        libopenblas0-pthread \
        openssh-server \
        python3 \
        python3-pip \
        python3-venv \
        tini \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv "$VIRTUAL_ENV" \
    && pip install --upgrade pip \
    && pip install \
        "ase>=3.22.0" \
        colorama \
        spglib \
        "pyparsing>=3.0" \
        elementy \
        platformdirs \
        packaging \
        pyyaml \
        requests \
        argcomplete \
        unyt \
        "numpy>=2.1"

RUN useradd -m -s /bin/bash tester \
    && mkdir -p /run/sshd /home/tester/.ssh /workspace \
    && chown -R tester:tester /home/tester/.ssh \
    && chmod 0700 /home/tester/.ssh \
    && printf '%s\n' \
        'PasswordAuthentication no' \
        'PermitRootLogin no' \
        'PubkeyAuthentication yes' \
        > /etc/ssh/sshd_config.d/ase2sprkkr-test.conf

COPY kkrscf /usr/local/bin/kkrscf
COPY entrypoint.sh /usr/local/bin/ssh-test-entrypoint
RUN chmod 0755 /usr/local/bin/kkrscf /usr/local/bin/ssh-test-entrypoint

WORKDIR /workspace
EXPOSE 22
ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/ssh-test-entrypoint"]
