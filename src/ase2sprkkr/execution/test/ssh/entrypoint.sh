#!/bin/sh
set -eu

install -o tester -g tester -m 0600 /authorized_keys /home/tester/.ssh/authorized_keys

cat >/usr/local/bin/ase2sprkkr <<'EOF'
#!/bin/sh
export PYTHONPATH=/workspace/ase2sprkkr/src${PYTHONPATH:+:$PYTHONPATH}
exec /opt/venv/bin/python -m ase2sprkkr.tools.ase2sprkkr "$@"
EOF
chmod 0755 /usr/local/bin/ase2sprkkr

exec /usr/sbin/sshd -D -e
