import toml
import re
import sys
data = toml.load('pyproject.toml')
for i in data['project']['dependencies']:
    match=re.match(r'^([a-z_]*)(?:;python_version<"3\.([0-9]+)")?$', i)
    if match[2] is not None:
        if sys.version_info[1] > int(match[2]):
            continue
    print(match[1], end=' ')
