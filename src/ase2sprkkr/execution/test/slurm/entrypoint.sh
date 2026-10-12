#!/usr/bin/env bash
set -Eeuo pipefail

mkdir -p /run/munge /run/sshd /var/lib/munge /var/log/munge \
         /var/spool/slurmctld /var/spool/slurmd /var/log/slurm \
         /root/.config/ase2sprkkr /home/tester/.ssh
chmod 0700 /home/tester/.ssh
chown tester:tester /home/tester/.ssh

# Slurm 23.11 initializes the cgroup/v2 layer even with task/cgroup disabled.
# This test runs without systemd, so cgroup.conf uses IgnoreSystemd=yes and the
# privileged, private cgroup namespace supplied by the outer pytest lets Slurm
# create its own scope here without touching the host hierarchy.
mkdir -p /sys/fs/cgroup/system.slice

if [[ ! -s /etc/munge/munge.key ]]; then
    dd if=/dev/urandom of=/etc/munge/munge.key bs=1024 count=1 status=none
fi
chown munge:munge /etc/munge/munge.key
chmod 0400 /etc/munge/munge.key
chown -R munge:munge /run/munge /var/lib/munge /var/log/munge
chown -R slurm:slurm /var/spool/slurmctld /var/log/slurm
chmod 0755 /var/spool/slurmd

cat > /root/.config/ase2sprkkr/__init__.py <<'PY'
import os
from ase2sprkkr.configuration import config

config.executables.dir = "/opt/sprkkr"
config.executables.suffix = os.environ.get("SPRKKR_EXECUTABLE_SUFFIX", "")
config.running.empty_spheres = False
config.running.mpi = False
config.running.print_output = False
PY

cat > /usr/local/bin/ase2sprkkr <<'SH'
#!/bin/sh
export PYTHONPATH=/workspace/ase2sprkkr/src${PYTHONPATH:+:$PYTHONPATH}
exec /opt/venv/bin/python -m ase2sprkkr.tools.ase2sprkkr "$@"
SH
chmod 0755 /usr/local/bin/ase2sprkkr

runuser -u munge -- /usr/sbin/munged --foreground &
munge_pid=$!
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
slurmctld_pid=$!
/usr/sbin/slurmd -D -N node1 &
slurmd_pid=$!

cleanup() {
    kill "$slurmd_pid" "$slurmctld_pid" "$munge_pid" 2>/dev/null || true
    wait "$slurmd_pid" "$slurmctld_pid" "$munge_pid" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

for _ in $(seq 1 60); do
    state="$(sinfo -h -N -o '%T' 2>/dev/null | head -n1 || true)"
    if [[ "$state" =~ ^(idle|mix|alloc)$ ]]; then
        if [[ -s /authorized_keys ]]; then
            install -o tester -g tester -m 0600 \
                /authorized_keys /home/tester/.ssh/authorized_keys
            ssh-keygen -A >/dev/null
            /usr/sbin/sshd
        fi
        exec "$@"
    fi
    sleep 0.25
done

echo "Slurm did not become ready" >&2
cat /var/log/slurm/slurmctld.log >&2 || true
cat /var/log/slurm/slurmd.log >&2 || true
exit 1