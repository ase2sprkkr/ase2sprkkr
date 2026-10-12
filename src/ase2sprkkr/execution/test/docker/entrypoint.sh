#!/usr/bin/env bash
set -Eeuo pipefail

services=",${ASE2SPRKKR_TEST_SERVICES:-},"

has_service() {
    [[ "$services" == *",$1,"* ]]
}

if [[ "$services" == ",," ]]; then
    echo "ASE2SPRKKR_TEST_SERVICES must select ssh, slurm, or both" >&2
    exit 2
fi

for service in ${ASE2SPRKKR_TEST_SERVICES//,/ }; do
    case "$service" in
        ssh|slurm) ;;
        *)
            echo "Unknown ASE2SPRKKR test service: $service" >&2
            exit 2
            ;;
    esac
done

mkdir -p /run/sshd /root/.config/ase2sprkkr \
         /home/tester/.config/ase2sprkkr /home/tester/.ssh
chmod 0700 /home/tester/.ssh
chown -R tester:tester /home/tester/.config /home/tester/.ssh

write_config() {
    local target=$1
    cat >"$target" <<'PY'
from ase2sprkkr.configuration import config

config.executables.dir = "/opt/sprkkr"
config.executables.suffix = ""
config.running.empty_spheres = False
config.running.mpi = False
config.running.print_output = False
PY
}

write_config /root/.config/ase2sprkkr/__init__.py
write_config /home/tester/.config/ase2sprkkr/__init__.py
chown tester:tester /home/tester/.config/ase2sprkkr/__init__.py

cat > /usr/local/bin/ase2sprkkr <<'SH'
#!/bin/sh
export PYTHONPATH=/workspace/ase2sprkkr/src${PYTHONPATH:+:$PYTHONPATH}
exec /opt/venv/bin/python -m ase2sprkkr.tools.ase2sprkkr "$@"
SH
chmod 0755 /usr/local/bin/ase2sprkkr

if has_service slurm; then
    mkdir -p /run/munge /var/lib/munge /var/log/munge \
             /var/spool/slurmctld /var/spool/slurmd /var/log/slurm

    # Slurm 23.11 initializes the cgroup/v2 layer even with task/cgroup
    # disabled.  The Slurm test therefore runs the container with a private
    # cgroup namespace and --privileged, while SSH-only tests do not need it.
    mkdir -p /sys/fs/cgroup/system.slice

    if [[ ! -s /etc/munge/munge.key ]]; then
        dd if=/dev/urandom of=/etc/munge/munge.key bs=1024 count=1 status=none
    fi
    chown munge:munge /etc/munge/munge.key
    chmod 0400 /etc/munge/munge.key
    chown -R munge:munge /run/munge /var/lib/munge /var/log/munge
    chown -R slurm:slurm /var/spool/slurmctld /var/log/slurm
    chmod 0755 /var/spool/slurmd

    runuser -u munge -- /usr/sbin/munged --foreground &
    for _ in $(seq 1 50); do
        if munge -n 2>/dev/null | unmunge >/dev/null 2>&1; then
            break
        fi
        sleep 0.1
    done
    if ! munge -n 2>/dev/null | unmunge >/dev/null 2>&1; then
        echo "MUNGE failed to start" >&2
        exit 1
    fi

    /usr/sbin/slurmctld -D &
    /usr/sbin/slurmd -D -N node1 &

    for _ in $(seq 1 60); do
        state="$(sinfo -h -N -o '%T' 2>/dev/null | head -n1 || true)"
        if [[ "$state" =~ ^(idle|mix|alloc)$ ]]; then
            break
        fi
        sleep 0.25
    done
    if [[ ! "$state" =~ ^(idle|mix|alloc)$ ]]; then
        echo "Slurm did not become ready" >&2
        cat /var/log/slurm/slurmctld.log >&2 || true
        cat /var/log/slurm/slurmd.log >&2 || true
        exit 1
    fi
fi

if has_service ssh; then
    if [[ ! -s /authorized_keys ]]; then
        echo "SSH test service requires /authorized_keys" >&2
        exit 2
    fi
    install -o tester -g tester -m 0600 \
        /authorized_keys /home/tester/.ssh/authorized_keys
    ssh-keygen -A >/dev/null
    /usr/sbin/sshd
fi

exec "$@"
