FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        munge \
        libgfortran5 \
        liblapack3 \
        libopenblas0-pthread \
        openssh-server \
        python3 \
        python3-pip \
        python3-venv \
        slurm-wlm \
        tini \
        util-linux \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv "$VIRTUAL_ENV" \
    && pip install --upgrade pip \
    && pip install \
        "psij-python==0.9.11" \
        "ase>=3.22.0" \
        colorama \
        spglib \
        "pyparsing>=3.0" \
        elementy \
        platformdirs \
        packaging \
        pytest \
        pyyaml \
        requests \
        argcomplete \
        unyt \
        "numpy>=2.1"

COPY slurm.conf /etc/slurm/slurm.conf
COPY cgroup.conf /etc/slurm/cgroup.conf
COPY entrypoint.sh /usr/local/bin/slurm-entrypoint

RUN useradd -m -s /bin/bash tester \
    && chmod 0755 /usr/local/bin/slurm-entrypoint \
    && mkdir -p \
        /var/spool/slurmctld \
        /var/spool/slurmd \
        /var/log/slurm \
        /run/munge \
        /run/sshd \
        /var/lib/munge \
        /var/log/munge \
        /root/.config/ase2sprkkr \
        /home/tester/.ssh \
        /workspace \
    && chown -R tester:tester /home/tester/.ssh \
    && chmod 0700 /home/tester/.ssh \
    && printf '%s\n' \
        'PasswordAuthentication no' \
        'PermitRootLogin no' \
        'PubkeyAuthentication yes' \
        > /etc/ssh/sshd_config.d/ase2sprkkr-test.conf \
    && chown -R slurm:slurm /var/spool/slurmctld /var/log/slurm \
    && chown -R munge:munge /run/munge /var/lib/munge /var/log/munge

WORKDIR /workspace
ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/slurm-entrypoint"]