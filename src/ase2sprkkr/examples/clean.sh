#!/bin/bash

#
# Delete all files in the examples directory and all subdirectories,
# which is not in the git repository. Use to clean up after running
# the examples
#

VERSION=$(grep -Po '^version\s*=\s*"\K[^"]+' ../../../pyproject.toml)

if [[ ! $VERSION =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "ERROR: VERSION ($VERSION) is not in x.y.z format"
    exit 1
fi

# ---------- Configuration ----------
#

tree=`git ls-tree -r --name-only HEAD "." | sed 's#^#./#'`

if [[ $? -ne 0 ]]; then
		#not in git, remotely

		OWNER="ase2sprkkr"
		REPO="ase2sprkkr"
		DIR_PATH="src/ase2sprkkr/examples"
		API_URL="https://api.github.com/repos/$OWNER/$REPO/contents/$DIR_PATH"
		TAG=v$VERSION

		# Optional: Set a GitHub token to increase rate limits
		# export GITHUB_TOKEN="your_token_here"


		# ---------- API Request ----------
		# echo "Fetching contents of $DIR_PATH in $OWNER/$REPO..."

		# Use GitHub token if available, else no auth
		if [[ -n "$GITHUB_TOKEN" ]]; then
			commit_sha=$(curl -sSL -H "Accept: application/vnd.github+json" \
			"https://api.github.com/repos/$OWNER/$REPO/tags" | jq -r --arg TAG "$TAG" '.[] | select(.name == $TAG) | .commit.sha')

			tree=$(curl -sSL -H "Accept: application/vnd.github+json" ${AUTH_HEADER:+-H "$AUTH_HEADER"} \
			"https://api.github.com/repos/$OWNER/$REPO/git/trees/$commit_sha?recursive=1")
		else

			commit_sha=$(curl -sSL -H "Accept: application/vnd.github+json" ${AUTH_HEADER:+-H "$AUTH_HEADER"} \
			"https://api.github.com/repos/$OWNER/$REPO/tags" | jq -r --arg TAG "$TAG" '.[] | select(.name == $TAG) | .commit.sha')

			tree=$(curl -sSL -H "Accept: application/vnd.github+json" \
			"https://api.github.com/repos/$OWNER/$REPO/git/trees/$commit_sha?recursive=1")
		fi

		# Fetch directory listing
		if [[ $? -ne 0 ]]; then
			echo "Error: Failed to fetch API."
			exit 1
		fi

		# ---------- List Files ----------
		tree=$(echo "$tree" | jq -r --arg DIR "$DIR_PATH/" '.tree[] | select(.type=="blob" and (.path | startswith($DIR))) | .path' | sed "s|^src/ase2sprkkr/examples/|./|")

fi

tmpfile=$(mktemp)
mapfile -t existing_files <<< "$tree"
printf "%s\n" "${existing_files[@]}" > "$tmpfile"

delete=`find "." -type f | fgrep -v __pycache__ | grep -vxF -f "$tmpfile"`
rm "$tmpfile"

if [[ $delete = "" ]] ; then
		echo "No files to clean ";
		exit 2 ;
fi

echo "Output files to remove:"
echo "----------------"
echo "$delete"
read -p "Delete the files [Y]?" answer
if [[ "$answer" =~ ^[Yy]$ ]]; then
		rm $delete
fi
